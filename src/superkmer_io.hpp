#pragma once

// On-disk superkmer format:
//   [uint8_t len_bases][uint8_t min_pos][uint64_t first_kmer][ceil(max(len-k,0)/4) packed suffix]
//
// Each superkmer stores the first k-mer explicitly, then packs only the suffix
// bases after that first k-mer (4 per byte, kache encoding A=0, C=1, G=2, T=3,
// MSB-first). This removes the phase-2 decode/rebuild of the first k bases.
//
//   len_bases — number of bases (max 255; superkmers are at most 2k−m bases).
//   min_pos   — 0-indexed start position of the minimizer m-mer within the
//               superkmer.  Stored so Phase 2 can compute ntHash(minimizer)
//               in O(m) without running MinimizerWindow::reset() in O(k).
//   first_kmer — the first k-mer of the superkmer encoded as Kmer<k>::as_int().
//
// Packed suffix encoding: suffix base i is at bits 7-2*(i%4) of byte i/4.
// Phase 2 initializes from `first_kmer` directly and only unpacks suffix bases
// while rolling through the superkmer.
//
// SuperkmerWriter — per-thread per-bucket buffered write.
//   Converts ASCII → packed on the fly; flushes to the shared ofstream under
//   its mutex when needs_flush() is true.
//
// SuperkmerReader — zero-copy sequential reader backed by mmap (Linux/POSIX).
//   Construction maps the entire file into the virtual address space.
//   next() reads the current record and advances the cursor; first_kmer(),
//   suffix_data(), size(), and min_pos() expose the current superkmer.

#include <fstream>
#include <mutex>
#include <string>
#include <cstdint>
#include <cstring>

#include "kache-hash/Kmer.hpp"

// Returns the path to the superkmer file for partition p under work_dir.
// Used in pipeline setup, counting, and cleanup — centralises the naming convention.
inline std::string partition_path(const std::string& work_dir, size_t p)
{
    return work_dir + "hash_" + std::to_string(p) + ".superkmers";
}

#include <sys/mman.h>
#include <sys/stat.h>
#include <fcntl.h>
#include <unistd.h>

#ifndef MAP_POPULATE
#  define MAP_POPULATE 0
#endif


// ─── Writer ───────────────────────────────────────────────────────────────────

struct SuperkmerWriter
{
    std::string buf;
    // flush_threshold is set per-writer based on n_parts so that total writer
    // memory across all partitions stays bounded to ~64 MB per thread:
    //   max(4 KB, 64 MB / n_parts)
    // Default (512 KB) is used when n_parts is small (≤128).
    const size_t flush_threshold;

    explicit SuperkmerWriter(size_t flush_thresh = 512u << 10)
        : flush_threshold(flush_thresh) {}

    static void pack_ascii_to_kache(const char* data, size_t len, uint8_t* packed)
    {
        std::memset(packed, 0, (len + 3u) / 4u);
        for (size_t i = 0; i < len; ++i) {
            const uint8_t b = ((uint8_t(data[i]) >> 2) ^ (uint8_t(data[i]) >> 1)) & 3u;
            packed[i >> 2] |= static_cast<uint8_t>(b << (6u - 2u * (i & 3u)));
        }
    }

    // Serialise one superkmer.
    // `data` is ASCII DNA (ACGT, any case); `len` is the number of bases;
    // `min_pos` is the 0-indexed start of the minimizer m-mer within the superkmer.
    template <uint16_t k>
    void append(const char* data, uint8_t len, uint8_t min_pos)
    {
        constexpr size_t FIRST_KMER_BYTES = sizeof(uint64_t);
        const size_t suffix_bases = len > k ? static_cast<size_t>(len - k)
                                            : static_cast<size_t>(len);
        const size_t suffix_bytes = (suffix_bases + 3u) / 4u;
        const size_t off = buf.size();
        buf.resize(off + 2u + FIRST_KMER_BYTES + suffix_bytes);

        // Two header bytes: length then minimizer position.
        buf[off]     = static_cast<char>(len);
        buf[off + 1] = static_cast<char>(min_pos);

        uint64_t first_kmer = 0;
        if (len >= k)
            first_kmer = kache_hash::Kmer<k>(data).as_int();
        std::memcpy(buf.data() + off + 2u, &first_kmer, FIRST_KMER_BYTES);

        uint8_t* suffix = reinterpret_cast<uint8_t*>(buf.data() + off + 2u + FIRST_KMER_BYTES);
        if (suffix_bytes > 0)
            pack_ascii_to_kache(data + (len > k ? k : 0), suffix_bases, suffix);
    }

    bool needs_flush() const { return buf.size() >= flush_threshold; }

    // Flush to the shared file under its mutex; no-op if buffer is empty.
    void flush_to(std::ofstream& file, std::mutex& mtx)
    {
        if (buf.empty()) return;
        std::lock_guard<std::mutex> g(mtx);
        file.write(buf.data(), static_cast<std::streamsize>(buf.size()));
        buf.clear();
    }

    // Flush to an in-memory string sink (streaming mode — avoids disk I/O).
    void flush_to_mem(std::string& dst, std::mutex& mtx)
    {
        if (buf.empty()) return;
        std::lock_guard<std::mutex> g(mtx);
        dst.append(buf);
        buf.clear();
    }
};


// ─── Reader ───────────────────────────────────────────────────────────────────

template <uint16_t k>
struct SuperkmerReader
{
    static constexpr size_t FIRST_KMER_BYTES = sizeof(uint64_t);

    explicit SuperkmerReader(const std::string& path)
    {
        fd_ = open(path.c_str(), O_RDONLY);
        if (fd_ < 0) return;

        struct stat sb;
        if (fstat(fd_, &sb) < 0 || sb.st_size == 0) return;
        size_ = static_cast<size_t>(sb.st_size);

        void* m = mmap(nullptr, size_, PROT_READ,
                       MAP_PRIVATE | MAP_POPULATE, fd_, 0);
        if (m == MAP_FAILED) return;
        map_ = static_cast<const char*>(m);

        madvise(const_cast<char*>(map_), size_,
                MADV_SEQUENTIAL | MADV_WILLNEED);

        cur_ = map_;
        end_ = map_ + size_;
    }

    ~SuperkmerReader()
    {
        if (map_) munmap(const_cast<char*>(map_), size_);
        if (fd_ >= 0) close(fd_);
    }

    SuperkmerReader(const SuperkmerReader&)            = delete;
    SuperkmerReader& operator=(const SuperkmerReader&) = delete;

    // Advance to the next superkmer.  Returns false at EOF.
    bool next()
    {
        if (cur_ + 2 > end_) return false;
        const uint8_t len8 = static_cast<uint8_t>(*cur_);
        if (len8 == 0) return false;
        min_pos_ = static_cast<uint8_t>(cur_[1]);
        cur_ += 2;

        if (cur_ + static_cast<ptrdiff_t>(FIRST_KMER_BYTES) > end_) return false;
        std::memcpy(&first_kmer_, cur_, FIRST_KMER_BYTES);
        cur_ += FIRST_KMER_BYTES;

        const size_t suffix_bases = len8 > k ? static_cast<size_t>(len8 - k)
                                             : static_cast<size_t>(len8);
        const size_t suffix_bytes = (suffix_bases + 3u) / 4u;
        if (cur_ + static_cast<ptrdiff_t>(suffix_bytes) > end_) return false;

        ptr_ = reinterpret_cast<const uint8_t*>(cur_);
        len_ = len8;
        cur_ += suffix_bytes;
        return true;
    }

    uint64_t first_kmer() const { return first_kmer_; }
    const uint8_t* suffix_data() const { return ptr_; }

    // Number of bases in the current superkmer.
    size_t size() const { return len_; }

    // 0-indexed start position of the minimizer m-mer within the current superkmer.
    uint8_t min_pos() const { return min_pos_; }

    // Total mapped bytes — available for diagnostics / capacity estimation.
    size_t file_size() const { return size_; }

    bool ok() const { return map_ != nullptr; }

private:
    int         fd_   = -1;
    size_t      size_ = 0;
    const char* map_  = nullptr;
    const char* cur_  = nullptr;
    const char* end_  = nullptr;
    const uint8_t* ptr_ = nullptr;
    size_t      len_     = 0;
    uint8_t     min_pos_ = 0;
    uint64_t    first_kmer_ = 0;
};


// ─── In-memory reader (streaming mode) ───────────────────────────────────────
//
// Same interface as SuperkmerReader but backed by an existing std::string
// rather than an mmap'd file.  Used by the streaming pipeline to avoid the
// disk write + mmap round-trip between Phase 1 and Phase 2.

template <uint16_t k>
struct MemoryReader
{
    static constexpr size_t FIRST_KMER_BYTES = sizeof(uint64_t);

    MemoryReader() = default;
    explicit MemoryReader(const std::string& data) noexcept
        : cur_(data.data()), end_(data.data() + data.size()) {}

    bool next() noexcept
    {
        if (cur_ + 2 > end_) return false;
        const uint8_t len8 = static_cast<uint8_t>(*cur_);
        if (len8 == 0) return false;
        min_pos_ = static_cast<uint8_t>(cur_[1]);
        cur_ += 2;
        if (cur_ + static_cast<ptrdiff_t>(FIRST_KMER_BYTES) > end_) return false;
        std::memcpy(&first_kmer_, cur_, FIRST_KMER_BYTES);
        cur_ += FIRST_KMER_BYTES;
        const size_t suffix_bases = len8 > k ? static_cast<size_t>(len8 - k)
                                             : static_cast<size_t>(len8);
        const size_t suffix_bytes = (suffix_bases + 3u) / 4u;
        if (cur_ + static_cast<ptrdiff_t>(suffix_bytes) > end_) return false;
        ptr_  = reinterpret_cast<const uint8_t*>(cur_);
        len_  = len8;
        cur_ += suffix_bytes;
        return true;
    }

    uint64_t       first_kmer() const noexcept { return first_kmer_; }
    const uint8_t* suffix_data() const noexcept { return ptr_; }
    size_t         size()        const noexcept { return len_; }
    uint8_t        min_pos()     const noexcept { return min_pos_; }
    bool           ok()          const noexcept { return cur_ != nullptr; }

private:
    const char*    cur_     = nullptr;
    const char*    end_     = nullptr;
    const uint8_t* ptr_     = nullptr;
    size_t         len_     = 0;
    uint8_t        min_pos_ = 0;
    uint64_t       first_kmer_ = 0;
};
