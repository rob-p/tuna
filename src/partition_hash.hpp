#pragma once

// Phase 1 — sequence input (FASTA/FASTQ, plain or gzip) → per-partition superkmer files.
//
// Partitioning: minimizer_hash % num_partitions.
//
// Parsing is delegated to SeqSource / GzInput + helicase SIMD parsers, which deliver
// ACTG-only chunks regardless of file format or compression.
//
// Internals:
//   extract_superkmers_from_actg<k, m, PartitionFn>  — pure ACTG sequence logic,
//                                                       no I/O, no threads.
//   partition_kmers_impl<k, m, PartitionFn>           — parallel harness.

#include "Config.hpp"
#include "superkmer_io.hpp"
#include "minimizer_window.hpp"
#include "seq_source.hpp"

#include <vector>
#include <deque>
#include <thread>
#include <mutex>
#include <condition_variable>
#include <atomic>
#include <exception>
#include <memory>
#include <string_view>
#include <utility>


// Per-writer flush threshold: max(4 KB, budget_per_thread / n_parts).
// Disk mode: budget is RAM-proportional — larger buffers yield fewer, bigger writes.
// In-memory mode: fixed 64 MB budget (writes go to RAM, no I/O benefit).
inline size_t writer_flush_threshold(size_t n_parts, size_t budget_per_thread)
{
    constexpr size_t MIN_FLUSH_BYTES = 4u << 10;  // 4 KB minimum
    return std::max(MIN_FLUSH_BYTES, budget_per_thread / n_parts);
}


// ─── Partition logic brick (ACTG-only) ────────────────────────────────────────
//
// Walk one ACTG-only sequence, detect superkmer boundaries on minimizer hash
// changes, and append each superkmer to the corresponding writer.
// Splits on hash change (not partition change) so every superkmer has one minimizer.
// min_pos comes from MinimizerWindow::min_lmer_pos() — no extra scan needed.
//
// FlushFn: void(std::vector<SuperkmerWriter<k,m>>&, size_t partition_id)
// Called after each append; O(1) per superkmer.
template <uint16_t k, uint16_t m>
static constexpr uint16_t phase1_signature_len_v =
#ifdef TUNA_PHASE1_SIG_LEN
    (TUNA_PHASE1_SIG_LEN < k ? TUNA_PHASE1_SIG_LEN : m)
#else
    m
#endif
;

template <uint16_t k, uint16_t m, uint16_t sig_m, typename PartitionFn, typename FlushFn>
void extract_superkmers_from_actg_sig(
    const char* const                    seq,
    const size_t                         seq_len,
    PartitionFn&&                        partition_fn,
    MinimizerWindow<k, sig_m>&           min_it,      // phase-1 signature iterator
    std::vector<SuperkmerWriter<k, m>>&  writers,
    uint64_t&                            kmer_count,
    uint64_t&                            sk_count,
    FlushFn&&                            flush_fn,
    std::vector<uint8_t>&                kache_buf)   // kache-encoded sequence buffer (A=0,C=1,G=2,T=3)
{
    using hdr_t = sk_hdr_t<k, m>;  // superkmer header type (local alias)
    static_assert(sig_m < k, "phase-1 signature length must be strictly less than k");
    static_assert(sig_m <= 32, "phase-1 signature length must be <= 32");
    // Max value of hdr_t: flush guard prevents sk_len from ever overflowing hdr_t.
    // Rare long same-minimizer runs can exceed 2k-m, so use the serialized
    // header limit rather than the natural two-window span.
    static constexpr size_t HDR_MAX = static_cast<size_t>(std::numeric_limits<hdr_t>::max());
    static constexpr hdr_t NO_MIN = sk_no_min<k, m>;

    if (seq_len < k) return;

    // Pre-encode ASCII→kache once; bases are read multiple times across superkmers.
    kache_buf.resize(seq_len);
    for (size_t i = 0; i < seq_len; ++i)
        kache_buf[i] = ((uint8_t(seq[i]) >> 2) ^ (uint8_t(seq[i]) >> 1)) & 3u;

    min_it.reset(seq);   // initialises from ASCII (called once per sequence — cheap)
    uint64_t prev_hash    = min_it.hash();
    uint64_t prev_min_pos = min_it.min_lmer_pos(); // absolute pos within seq
    size_t   pid          = partition_fn(prev_hash);  // partition id
    size_t   sk_start     = 0;                        // superkmer start position in current sequence

    for (size_t pos = k; pos < seq_len; ++pos) {
        min_it.advance_kache(kache_buf[pos]);
        const uint64_t new_hash = min_it.hash();
        if (__builtin_expect(new_hash != prev_hash || pos - sk_start >= HDR_MAX, 0)) {
            const auto sk_len  = static_cast<hdr_t>(pos - sk_start);
            const auto min_pos = [] (uint64_t p, size_t s) {
                if constexpr (sig_m == m) return static_cast<hdr_t>(p - s);
                else return NO_MIN;
            }(prev_min_pos, sk_start);
            writers[pid].append_kache(kache_buf.data() + sk_start, sk_len, min_pos);
            flush_fn(writers, pid);
            kmer_count += sk_len - k + 1;
            ++sk_count;
            prev_hash    = new_hash;
            prev_min_pos = min_it.min_lmer_pos();
            pid          = partition_fn(new_hash);
            sk_start     = pos - (k - 1);
        }
    }

    const auto sk_len  = static_cast<hdr_t>(seq_len - sk_start);
    const auto min_pos = [] (uint64_t p, size_t s) {
        if constexpr (sig_m == m) return static_cast<hdr_t>(p - s);
        else return NO_MIN;
    }(prev_min_pos, sk_start);
    writers[pid].append_kache(kache_buf.data() + sk_start, sk_len, min_pos);
    flush_fn(writers, pid);
    kmer_count += sk_len - k + 1;
    ++sk_count;
}

template <uint16_t k, uint16_t m, typename PartitionFn, typename FlushFn>
void extract_superkmers_from_actg(
    const char* const                    seq,
    const size_t                         seq_len,
    PartitionFn&&                        partition_fn,
    MinimizerWindow<k, m>&               min_it,
    std::vector<SuperkmerWriter<k, m>>&  writers,
    uint64_t&                            kmer_count,
    uint64_t&                            sk_count,
    FlushFn&&                            flush_fn,
    std::vector<uint8_t>&                kache_buf)
{
    static constexpr uint16_t sig_m = phase1_signature_len_v<k, m>;
    if constexpr (sig_m == m) {
        extract_superkmers_from_actg_sig<k, m, m>(
            seq, seq_len, std::forward<PartitionFn>(partition_fn), min_it, writers,
            kmer_count, sk_count, std::forward<FlushFn>(flush_fn), kache_buf);
    } else {
        MinimizerWindow<k, sig_m> sig_it;
        extract_superkmers_from_actg_sig<k, m, sig_m>(
            seq, seq_len, std::forward<PartitionFn>(partition_fn), sig_it, writers,
            kmer_count, sk_count, std::forward<FlushFn>(flush_fn), kache_buf);
    }
}

template <uint16_t k, uint16_t m, uint16_t sig_m, typename PackedSeq, typename PartitionFn, typename FlushFn>
void extract_superkmers_from_packed_nt_sig(
    const PackedSeq&                     seq,
    PartitionFn&&                        partition_fn,
    MinimizerWindow<k, sig_m>&           min_it,
    std::vector<SuperkmerWriter<k, m>>&  writers,
    uint64_t&                            kmer_count,
    uint64_t&                            sk_count,
    FlushFn&&                            flush_fn)
{
    using hdr_t = sk_hdr_t<k, m>;
    static_assert(sig_m < k, "phase-1 signature length must be strictly less than k");
    static_assert(sig_m <= 32, "phase-1 signature length must be <= 32");
    static constexpr size_t HDR_MAX = static_cast<size_t>(std::numeric_limits<hdr_t>::max());
    static constexpr hdr_t NO_MIN = sk_no_min<k, m>;

    const size_t seq_len = seq.len();
    if (seq_len < k) return;

    auto get_nt = [&](size_t i) noexcept -> uint8_t { return seq.get(i); };

    min_it.reset_nt(get_nt);
    uint64_t prev_hash    = min_it.hash();
    uint64_t prev_min_pos = min_it.min_lmer_pos();
    size_t   pid          = partition_fn(prev_hash);
    size_t   sk_start     = 0;

    for (size_t pos = k; pos < seq_len; ++pos) {
        min_it.advance_nt(get_nt(pos));
        const uint64_t new_hash = min_it.hash();
        if (__builtin_expect(new_hash != prev_hash || pos - sk_start >= HDR_MAX, 0)) {
            const auto sk_len  = static_cast<hdr_t>(pos - sk_start);
            const auto min_pos = [] (uint64_t p, size_t s) {
                if constexpr (sig_m == m) return static_cast<hdr_t>(p - s);
                else return NO_MIN;
            }(prev_min_pos, sk_start);
            writers[pid].append_nt(get_nt, sk_start, sk_len, min_pos);
            flush_fn(writers, pid);
            kmer_count += sk_len - k + 1;
            ++sk_count;
            prev_hash    = new_hash;
            prev_min_pos = min_it.min_lmer_pos();
            pid          = partition_fn(new_hash);
            sk_start     = pos - (k - 1);
        }
    }

    const auto sk_len  = static_cast<hdr_t>(seq_len - sk_start);
    const auto min_pos = [] (uint64_t p, size_t s) {
        if constexpr (sig_m == m) return static_cast<hdr_t>(p - s);
        else return NO_MIN;
    }(prev_min_pos, sk_start);
    writers[pid].append_nt(get_nt, sk_start, sk_len, min_pos);
    flush_fn(writers, pid);
    kmer_count += sk_len - k + 1;
    ++sk_count;
}

template <uint16_t k, uint16_t m, typename PackedSeq, typename PartitionFn, typename FlushFn>
void extract_superkmers_from_packed_nt(
    const PackedSeq&                     seq,
    PartitionFn&&                        partition_fn,
    MinimizerWindow<k, m>&               min_it,
    std::vector<SuperkmerWriter<k, m>>&  writers,
    uint64_t&                            kmer_count,
    uint64_t&                            sk_count,
    FlushFn&&                            flush_fn)
{
    static constexpr uint16_t sig_m = phase1_signature_len_v<k, m>;
    if constexpr (sig_m == m) {
        extract_superkmers_from_packed_nt_sig<k, m, m>(
            seq, std::forward<PartitionFn>(partition_fn), min_it, writers,
            kmer_count, sk_count, std::forward<FlushFn>(flush_fn));
    } else {
        MinimizerWindow<k, sig_m> sig_it;
        extract_superkmers_from_packed_nt_sig<k, m, sig_m>(
            seq, std::forward<PartitionFn>(partition_fn), sig_it, writers,
            kmer_count, sk_count, std::forward<FlushFn>(flush_fn));
    }
}


using PackedDNAWord = __uint128_t;

struct PackedReadRecord {
    size_t        word_offset = 0;
    size_t        word_count  = 0;
    PackedDNAWord tail        = 0;
    size_t        len         = 0;
};

struct PackedReadBatch {
    std::vector<PackedDNAWord> words;
    std::vector<PackedReadRecord> records;
    size_t bases = 0;

    bool empty() const noexcept { return records.empty(); }

    void reserve(size_t max_records, size_t target_bases)
    {
        records.reserve(max_records);
        words.reserve((target_bases + 63u) / 64u + max_records);
    }

    void append(const helicase::PackedDNA& dna)
    {
        auto [dna_words, tail] = dna.bits();
        PackedReadRecord rec;
        rec.word_offset = words.size();
        rec.word_count  = dna_words.size();
        rec.tail        = tail;
        rec.len         = dna.len();
        words.insert(words.end(), dna_words.begin(), dna_words.end());
        records.push_back(rec);
        bases += rec.len;
    }
};

struct PackedReadView {
    const PackedDNAWord* words = nullptr;
    size_t               word_count = 0;
    PackedDNAWord        tail = 0;
    size_t               len_bases = 0;

    size_t len() const noexcept { return len_bases; }

    uint8_t get(size_t i) const noexcept
    {
        const size_t word = i / 64u;
        const size_t bit  = i % 64u;
        const PackedDNAWord value = word < word_count ? words[word] : tail;
        return static_cast<uint8_t>((value >> (2u * bit)) & 0b11u);
    }
};

class AsyncPartitionWriters {
    struct Shard {
        std::mutex mtx;
        std::condition_variable cv;
        std::deque<SuperkmerWriteBlock> queue;
        size_t queued_bytes = 0;
        bool done = false;
    };

    std::vector<std::ofstream>& buckets_;
    std::vector<std::unique_ptr<Shard>> shards_;
    std::vector<std::thread> threads_;
    size_t max_queue_bytes_;
    std::exception_ptr error_ = nullptr;
    std::mutex error_mutex_;

    size_t shard_for(size_t partition) const noexcept
    {
        return partition % shards_.size();
    }

    void record_error()
    {
        std::lock_guard<std::mutex> lk(error_mutex_);
        if (!error_) error_ = std::current_exception();
    }

    void worker(size_t shard_id)
    {
        auto& shard = *shards_[shard_id];
        bool failed = false;
        while (true) {
            SuperkmerWriteBlock block;
            {
                std::unique_lock<std::mutex> lk(shard.mtx);
                shard.cv.wait(lk, [&]{ return shard.done || !shard.queue.empty(); });
                if (shard.queue.empty()) {
                    if (shard.done) break;
                    continue;
                }
                block = std::move(shard.queue.front());
                shard.queue.pop_front();
                shard.queued_bytes -= block.size;
            }
            shard.cv.notify_all();

            if (failed) continue;
            try {
                auto& out = buckets_[block.partition];
                out.write(block.data, static_cast<std::streamsize>(block.size));
                if (!out) throw std::runtime_error("tuna: failed while writing partition file");
            } catch (...) {
                failed = true;
                record_error();
            }
        }
    }

public:
    AsyncPartitionWriters(std::vector<std::ofstream>& buckets,
                          size_t n_shards,
                          size_t max_queue_bytes)
        : buckets_(buckets),
          max_queue_bytes_(max_queue_bytes)
    {
        n_shards = std::max<size_t>(1, std::min(n_shards, buckets_.size()));
        shards_.reserve(n_shards);
        for (size_t i = 0; i < n_shards; ++i)
            shards_.push_back(std::make_unique<Shard>());
        threads_.reserve(n_shards);
        for (size_t i = 0; i < n_shards; ++i)
            threads_.emplace_back([this, i]{ worker(i); });
    }

    AsyncPartitionWriters(const AsyncPartitionWriters&) = delete;
    AsyncPartitionWriters& operator=(const AsyncPartitionWriters&) = delete;

    void enqueue(SuperkmerWriteBlock&& block)
    {
        if (block.size == 0) return;
        auto& shard = *shards_[shard_for(block.partition)];
        {
            std::unique_lock<std::mutex> lk(shard.mtx);
            shard.cv.wait(lk, [&]{
                return shard.done ||
                       shard.queued_bytes + block.size <= max_queue_bytes_;
            });
            if (shard.done) return;
            shard.queued_bytes += block.size;
            shard.queue.emplace_back(std::move(block));
        }
        shard.cv.notify_one();
    }

    void finish()
    {
        for (auto& shard_ptr : shards_) {
            auto& shard = *shard_ptr;
            {
                std::lock_guard<std::mutex> lk(shard.mtx);
                shard.done = true;
            }
            shard.cv.notify_all();
        }
        for (auto& th : threads_)
            if (th.joinable()) th.join();
        if (error_) std::rethrow_exception(error_);
    }
};


// ─── Producer-consumer harness for a single gz file ──────────────────────────
//
// When there is only one input file and it is gzip-compressed, the standard
// work-stealing harness wastes all but one thread (min(n_threads,1)=1).
// This harness instead:
//   • 1 producer thread: gz decompression + split_actg → batches of ACTG chunks
//   • n_threads-1 consumer threads: extract_superkmers_from_actg per chunk
//
// The queue is bounded (MAX_QUEUE batches) for backpressure.  Consumers flush
// SuperkmerWriters to the shared bucket files under per-bucket mutexes.

template <uint16_t k, uint16_t m, typename PartitionFn>
PartitionStats partition_kmers_gz_pc(
    const Config&               cfg,
    const std::string&          gz_path,
    std::vector<std::ofstream>& buckets,
    PartitionFn                 partition_fn,
    size_t                      n_threads,          // ≥ 2 (1 producer + rest consumers)
    size_t                      write_budget_per_thread)
{
    using Batch = std::vector<std::string>;     // ACTG-only chunks pre-split by producer

    constexpr size_t MAX_QUEUE  = 32;           // max batches in flight
    constexpr size_t BATCH_SEQS = 512;          // sequences per batch

    const size_t n_parts     = cfg.num_partitions;
    const size_t n_consumers = n_threads - 1;

    std::deque<Batch>       queue;
    std::mutex              q_mutex;
    std::condition_variable q_cv;
    bool                    producer_done = false;
    std::exception_ptr      producer_error = nullptr;
    std::atomic<bool>       stop{false};
    std::exception_ptr      consumer_error = nullptr;
    std::mutex              consumer_error_mutex;

    const size_t writer_shards = std::min<size_t>(4, std::max<size_t>(1, n_threads / 4));
    AsyncPartitionWriters async_writers(
        buckets, writer_shards, std::max<size_t>(64u << 20, write_budget_per_thread));
    std::atomic<uint64_t>   total_seqs{0}, total_kmers{0}, total_superkmers{0};

    // Producer: decompress gz via GzInput, use helicase SIMD parser to deliver
    // ACTG-only chunks, accumulate into batches and push to the queue.
    auto producer_fn = [&]() {
        auto feed = [&](auto& parser) {
            while (true) {
                Batch batch;
                size_t chunk_count = 0;
                while (!stop.load(std::memory_order_relaxed)
                       && chunk_count < BATCH_SEQS
                       && parser.next()) {
                    auto [ptr, len] = parser.get_dna_raw();
                    batch.emplace_back(ptr, len);
                    ++chunk_count;
                }
                if (stop.load(std::memory_order_relaxed)) break;
                if (batch.empty()) break;
                {
                    std::unique_lock<std::mutex> lk(q_mutex);
                    q_cv.wait(lk, [&]{ return queue.size() < MAX_QUEUE; });
                    queue.push_back(std::move(batch));
                }
                q_cv.notify_one();
            }
            { std::lock_guard<std::mutex> lk(q_mutex); producer_done = true; }
            q_cv.notify_all();
        };
        try {
            GzInput inp(gz_path);
            if (inp.first_byte() == '@') {
                helicase::FastqParser<HELICASE_ACTG, GzInput> p(std::move(inp));
                feed(p);
            } else {
                helicase::FastaParser<HELICASE_ACTG, GzInput> p(std::move(inp));
                feed(p);
            }
        } catch (...) {
            {
                std::lock_guard<std::mutex> lk(q_mutex);
                producer_error = std::current_exception();
                producer_done = true;
            }
            stop.store(true, std::memory_order_relaxed);
            q_cv.notify_all();
        }
    };

    // Consumer: pull batches from the queue and extract superkmers.
    auto consumer_fn = [&]() {
        try {
            const size_t flush_thresh = writer_flush_threshold(n_parts, write_budget_per_thread);
            MinimizerWindow<k, m>               min_it;
            std::vector<SuperkmerWriter<k, m>>  writers(n_parts, SuperkmerWriter<k, m>(flush_thresh));
            std::vector<uint8_t>         kache_buf;
            uint64_t local_seqs = 0, local_kmers = 0, local_superkmers = 0;

            auto flush_fn = [&](std::vector<SuperkmerWriter<k, m>>& ws, size_t p) {
                if (ws[p].needs_flush()) async_writers.enqueue(ws[p].release_block(p));
            };

            while (true) {
                if (stop.load(std::memory_order_relaxed)) break;
                Batch batch;
                {
                    std::unique_lock<std::mutex> lk(q_mutex);
                    q_cv.wait(lk, [&]{ return !queue.empty() || producer_done; });
                    if (queue.empty()) break;       // producer_done && queue empty → done
                    batch = std::move(queue.front());
                    queue.pop_front();
                }
                q_cv.notify_one();                  // wake producer if it was waiting for space

                for (const auto& chunk : batch) {
                    if (stop.load(std::memory_order_relaxed)) break;
                    extract_superkmers_from_actg<k, m>(
                        chunk.data(), chunk.size(), partition_fn,
                        min_it, writers, local_kmers, local_superkmers, flush_fn, kache_buf);
                    ++local_seqs;
                }
            }

            for (size_t p = 0; p < n_parts; ++p) {
                if (!writers[p].empty())
                    async_writers.enqueue(writers[p].release_block(p));
            }

            total_seqs       .fetch_add(local_seqs,        std::memory_order_relaxed);
            total_kmers      .fetch_add(local_kmers,       std::memory_order_relaxed);
            total_superkmers .fetch_add(local_superkmers,  std::memory_order_relaxed);
        } catch (...) {
            {
                std::lock_guard<std::mutex> lk(consumer_error_mutex);
                if (!consumer_error) consumer_error = std::current_exception();
            }
            stop.store(true, std::memory_order_relaxed);
            q_cv.notify_all();
        }
    };

    std::vector<std::thread> threads;
    threads.reserve(n_threads);
    threads.emplace_back(producer_fn);
    for (size_t t = 0; t < n_consumers; ++t)
        threads.emplace_back(consumer_fn);
    for (auto& th : threads) th.join();
    async_writers.finish();
    if (producer_error) std::rethrow_exception(producer_error);
    if (consumer_error) std::rethrow_exception(consumer_error);

    return { total_seqs.load(), total_kmers.load(), total_superkmers.load() };
}


// Multiple gzipped inputs need finer-grained scheduling than file-level work
// stealing: each producer owns one gz stream at a time and consumers partition
// copied ACTG chunks from a shared queue.
template <uint16_t k, uint16_t m, typename PartitionFn>
PartitionStats partition_kmers_multi_gz_pc(
    const Config&               cfg,
    std::vector<std::ofstream>& buckets,
    PartitionFn                 partition_fn,
    size_t                      write_budget_per_thread)
{
    using Batch = PackedReadBatch;

    constexpr size_t MAX_QUEUE = 32;
    constexpr size_t MAX_BATCH_ITEMS = 524288;
    constexpr size_t TARGET_BATCH_BASES = 64u << 20;

    const size_t n_files = cfg.input_files.size();
    const size_t n_threads = static_cast<size_t>(cfg.num_threads);
    const size_t n_parts = cfg.num_partitions;
    const size_t n_producers = std::min(n_files, std::max<size_t>(1, n_threads / 2));
    const size_t n_consumers = std::max<size_t>(1, n_threads - n_producers);

    std::atomic<size_t> next_file{0};
    std::atomic<size_t> active_producers{n_producers};
    std::deque<Batch> queue;
    std::mutex q_mutex;
    std::condition_variable q_cv;
    std::atomic<bool> stop{false};
    std::exception_ptr producer_error = nullptr;
    std::exception_ptr consumer_error = nullptr;
    std::mutex error_mutex;

    const size_t writer_shards = std::min<size_t>(4, std::max<size_t>(1, n_threads / 4));
    AsyncPartitionWriters async_writers(
        buckets, writer_shards, std::max<size_t>(64u << 20, write_budget_per_thread));
    std::atomic<uint64_t> total_seqs{0}, total_kmers{0}, total_superkmers{0};

    auto push_batch = [&](Batch& batch, size_t& batch_bases) {
        if (batch.empty()) return;
        {
            std::unique_lock<std::mutex> lk(q_mutex);
            q_cv.wait(lk, [&]{
                return queue.size() < MAX_QUEUE || stop.load(std::memory_order_relaxed);
            });
            if (stop.load(std::memory_order_relaxed)) return;
            queue.push_back(std::move(batch));
        }
        q_cv.notify_one();
        batch = Batch{};
        batch.reserve(MAX_BATCH_ITEMS, TARGET_BATCH_BASES);
        batch_bases = 0;
    };

    auto producer_done = [&]() {
        active_producers.fetch_sub(1, std::memory_order_relaxed);
        q_cv.notify_all();
    };

    auto producer_fn = [&]() {
        try {
            while (!stop.load(std::memory_order_relaxed)) {
                const size_t fi = next_file.fetch_add(1, std::memory_order_relaxed);
                if (fi >= n_files) break;
                const auto& input_path = cfg.input_files[fi];

                Batch batch;
                batch.reserve(MAX_BATCH_ITEMS, TARGET_BATCH_BASES);
                size_t batch_bases = 0;
                GzInput inp(input_path);
                if (inp.first_byte() == '@') {
                    helicase::FastqParser<HELICASE_ACTG_PACKED, GzInput> p(std::move(inp));
                    while (!stop.load(std::memory_order_relaxed) && p.next()) {
                        const size_t len = p.get_dna_len();
                        if (len >= k) {
                            batch.append(p.get_dna_packed());
                            batch_bases += len;
                            if (batch_bases >= TARGET_BATCH_BASES ||
                                batch.records.size() >= MAX_BATCH_ITEMS)
                                push_batch(batch, batch_bases);
                        }
                    }
                } else {
                    helicase::FastaParser<HELICASE_ACTG_PACKED, GzInput> p(std::move(inp));
                    while (!stop.load(std::memory_order_relaxed) && p.next()) {
                        const size_t len = p.get_dna_len();
                        if (len >= k) {
                            batch.append(p.get_dna_packed());
                            batch_bases += len;
                            if (batch_bases >= TARGET_BATCH_BASES ||
                                batch.records.size() >= MAX_BATCH_ITEMS)
                                push_batch(batch, batch_bases);
                        }
                    }
                }
                if (!stop.load(std::memory_order_relaxed))
                    push_batch(batch, batch_bases);
            }
        } catch (...) {
            {
                std::lock_guard<std::mutex> lk(error_mutex);
                if (!producer_error) producer_error = std::current_exception();
            }
            stop.store(true, std::memory_order_relaxed);
            q_cv.notify_all();
        }
        producer_done();
    };

    auto consumer_fn = [&]() {
        try {
            const size_t flush_thresh = writer_flush_threshold(n_parts, write_budget_per_thread);
            MinimizerWindow<k, m> min_it;
            std::vector<SuperkmerWriter<k, m>> writers(n_parts, SuperkmerWriter<k, m>(flush_thresh));
            uint64_t local_seqs = 0, local_kmers = 0, local_superkmers = 0;

            auto flush_fn = [&](std::vector<SuperkmerWriter<k, m>>& ws, size_t p) {
                if (ws[p].needs_flush()) async_writers.enqueue(ws[p].release_block(p));
            };

            while (true) {
                if (stop.load(std::memory_order_relaxed)) break;
                Batch batch;
                {
                    std::unique_lock<std::mutex> lk(q_mutex);
                    q_cv.wait(lk, [&]{
                        return !queue.empty() ||
                               active_producers.load(std::memory_order_relaxed) == 0 ||
                               stop.load(std::memory_order_relaxed);
                    });
                    if (queue.empty()) break;
                    batch = std::move(queue.front());
                    queue.pop_front();
                }
                q_cv.notify_one();

                for (const auto& rec : batch.records) {
                    if (stop.load(std::memory_order_relaxed)) break;
                    PackedReadView chunk{
                        batch.words.data() + rec.word_offset,
                        rec.word_count,
                        rec.tail,
                        rec.len
                    };
                    extract_superkmers_from_packed_nt<k, m>(
                        chunk, partition_fn,
                        min_it, writers, local_kmers, local_superkmers, flush_fn);
                    ++local_seqs;
                }
            }

            for (size_t p = 0; p < n_parts; ++p) {
                if (!writers[p].empty())
                    async_writers.enqueue(writers[p].release_block(p));
            }

            total_seqs.fetch_add(local_seqs, std::memory_order_relaxed);
            total_kmers.fetch_add(local_kmers, std::memory_order_relaxed);
            total_superkmers.fetch_add(local_superkmers, std::memory_order_relaxed);
        } catch (...) {
            {
                std::lock_guard<std::mutex> lk(error_mutex);
                if (!consumer_error) consumer_error = std::current_exception();
            }
            stop.store(true, std::memory_order_relaxed);
            q_cv.notify_all();
        }
    };

    std::vector<std::thread> threads;
    threads.reserve(n_producers + n_consumers);
    for (size_t t = 0; t < n_producers; ++t)
        threads.emplace_back(producer_fn);
    for (size_t t = 0; t < n_consumers; ++t)
        threads.emplace_back(consumer_fn);
    for (auto& th : threads) th.join();
    async_writers.finish();
    if (producer_error) std::rethrow_exception(producer_error);
    if (consumer_error) std::rethrow_exception(consumer_error);

    return { total_seqs.load(), total_kmers.load(), total_superkmers.load() };
}


// One or a few large plain inputs need finer-grained scheduling than file-level
// work stealing.  A producer mmaps each file, publishes string_views into that
// mapping, and waits for consumers to drain them before advancing to the next
// file so the views never outlive their backing mmap.
template <uint16_t k, uint16_t m, typename PartitionFn>
PartitionStats partition_kmers_plain_pc(
    const Config&               cfg,
    std::vector<std::ofstream>& buckets,
    PartitionFn                 partition_fn,
    size_t                      write_budget_per_thread)
{
    using Batch = std::vector<std::string_view>;

    constexpr size_t MAX_QUEUE = 32;
    constexpr size_t MAX_BATCH_ITEMS = 8192;
    constexpr size_t TARGET_BATCH_BASES = 1u << 20;
    constexpr size_t TARGET_CHUNK_BASES = 1u << 20;

    const size_t n_threads = static_cast<size_t>(cfg.num_threads);
    const size_t n_parts = cfg.num_partitions;
    const size_t n_consumers = n_threads - 1;

    std::deque<Batch> queue;
    std::mutex q_mutex;
    std::condition_variable q_cv;
    bool producer_done = false;
    size_t active_batches = 0;
    std::atomic<bool> stop{false};
    std::exception_ptr producer_error = nullptr;
    std::exception_ptr consumer_error = nullptr;
    std::mutex error_mutex;

    const size_t writer_shards = std::min<size_t>(4, std::max<size_t>(1, n_threads / 4));
    AsyncPartitionWriters async_writers(
        buckets, writer_shards, std::max<size_t>(64u << 20, write_budget_per_thread));
    std::atomic<uint64_t> total_seqs{0}, total_kmers{0}, total_superkmers{0};

    auto push_batch = [&](Batch& batch, size_t& batch_bases) {
        if (batch.empty()) return;
        {
            std::unique_lock<std::mutex> lk(q_mutex);
            q_cv.wait(lk, [&]{
                return queue.size() < MAX_QUEUE || stop.load(std::memory_order_relaxed);
            });
            if (stop.load(std::memory_order_relaxed)) return;
            queue.push_back(std::move(batch));
        }
        q_cv.notify_one();
        batch = Batch{};
        batch_bases = 0;
    };

    auto wait_until_consumed = [&]() {
        std::unique_lock<std::mutex> lk(q_mutex);
        q_cv.wait(lk, [&]{
            return (queue.empty() && active_batches == 0) ||
                   stop.load(std::memory_order_relaxed);
        });
    };

    auto producer_fn = [&]() {
        auto feed = [&](auto& parser) {
            Batch batch;
            size_t batch_bases = 0;
            auto add_chunk = [&](const char* ptr, size_t len) {
                if (len < k) return;
                const size_t step = TARGET_CHUNK_BASES > k
                    ? TARGET_CHUNK_BASES - (k - 1)
                    : TARGET_CHUNK_BASES;
                for (size_t off = 0; off < len; ) {
                    const size_t sub_len = std::min(TARGET_CHUNK_BASES, len - off);
                    batch.emplace_back(ptr + off, sub_len);
                    batch_bases += sub_len;
                    if (batch_bases >= TARGET_BATCH_BASES ||
                        batch.size() >= MAX_BATCH_ITEMS)
                        push_batch(batch, batch_bases);
                    if (off + sub_len >= len) break;
                    off += step;
                }
            };

            while (!stop.load(std::memory_order_relaxed) && parser.next()) {
                auto [ptr, len] = parser.get_dna_raw();
                add_chunk(ptr, len);
            }
            if (!stop.load(std::memory_order_relaxed))
                push_batch(batch, batch_bases);
            wait_until_consumed();
        };

        try {
            for (const auto& input_path : cfg.input_files) {
                if (stop.load(std::memory_order_relaxed)) break;
                helicase::MmapInput inp(input_path);
                if (inp.first_byte() == '@') {
                    helicase::FastqParser<HELICASE_ACTG, helicase::MmapInput> p(std::move(inp));
                    feed(p);
                } else {
                    helicase::FastaParser<HELICASE_ACTG, helicase::MmapInput> p(std::move(inp));
                    feed(p);
                }
            }
            {
                std::lock_guard<std::mutex> lk(q_mutex);
                producer_done = true;
            }
            q_cv.notify_all();
        } catch (...) {
            {
                std::lock_guard<std::mutex> lk(q_mutex);
                producer_error = std::current_exception();
                producer_done = true;
            }
            stop.store(true, std::memory_order_relaxed);
            q_cv.notify_all();
        }
    };

    auto consumer_fn = [&]() {
        try {
            const size_t flush_thresh = writer_flush_threshold(n_parts, write_budget_per_thread);
            MinimizerWindow<k, m> min_it;
            std::vector<SuperkmerWriter<k, m>> writers(n_parts, SuperkmerWriter<k, m>(flush_thresh));
            std::vector<uint8_t> kache_buf;
            uint64_t local_seqs = 0, local_kmers = 0, local_superkmers = 0;

            auto flush_fn = [&](std::vector<SuperkmerWriter<k, m>>& ws, size_t p) {
                if (ws[p].needs_flush()) async_writers.enqueue(ws[p].release_block(p));
            };

            while (true) {
                if (stop.load(std::memory_order_relaxed)) break;
                Batch batch;
                {
                    std::unique_lock<std::mutex> lk(q_mutex);
                    q_cv.wait(lk, [&]{
                        return !queue.empty() || producer_done ||
                               stop.load(std::memory_order_relaxed);
                    });
                    if (queue.empty()) break;
                    batch = std::move(queue.front());
                    queue.pop_front();
                    ++active_batches;
                }
                q_cv.notify_all();

                for (const auto chunk : batch) {
                    if (stop.load(std::memory_order_relaxed)) break;
                    extract_superkmers_from_actg<k, m>(
                        chunk.data(), chunk.size(), partition_fn,
                        min_it, writers, local_kmers, local_superkmers, flush_fn, kache_buf);
                    ++local_seqs;
                }

                {
                    std::lock_guard<std::mutex> lk(q_mutex);
                    --active_batches;
                }
                q_cv.notify_all();
            }

            for (size_t p = 0; p < n_parts; ++p) {
                if (!writers[p].empty())
                    async_writers.enqueue(writers[p].release_block(p));
            }

            total_seqs.fetch_add(local_seqs, std::memory_order_relaxed);
            total_kmers.fetch_add(local_kmers, std::memory_order_relaxed);
            total_superkmers.fetch_add(local_superkmers, std::memory_order_relaxed);
        } catch (...) {
            {
                std::lock_guard<std::mutex> lk(error_mutex);
                if (!consumer_error) consumer_error = std::current_exception();
            }
            stop.store(true, std::memory_order_relaxed);
            q_cv.notify_all();
        }
    };

    std::vector<std::thread> threads;
    threads.reserve(n_threads);
    threads.emplace_back(producer_fn);
    for (size_t t = 0; t < n_consumers; ++t)
        threads.emplace_back(consumer_fn);
    for (auto& th : threads) th.join();
    async_writers.finish();
    if (producer_error) std::rethrow_exception(producer_error);
    if (consumer_error) std::rethrow_exception(consumer_error);

    return { total_seqs.load(), total_kmers.load(), total_superkmers.load() };
}


// ─── Parallel harness ─────────────────────────────────────────────────────────
//
// File-level work stealing: each worker atomically claims the next file and runs
// helicase + superkmer extraction on its own thread.
// Load balancing is coarse (file granularity) but negligible for equal-size files.

template <uint16_t k, uint16_t m, typename PartitionFn>
PartitionStats partition_kmers_impl(
    const Config&               cfg,
    std::vector<std::ofstream>& buckets,
    PartitionFn                 partition_fn,
    size_t                      write_budget_per_thread)
{
    const size_t n_files      = cfg.input_files.size();
    const size_t n_threads_req = static_cast<size_t>(cfg.num_threads);

    if (n_files >= 1 && n_threads_req > 1) {
        bool all_gz = true;
        bool all_plain = true;
        for (const auto& f : cfg.input_files) {
            if (!(f.size() > 3 && f.compare(f.size() - 3, 3, ".gz") == 0)) {
                all_gz = false;
            } else {
                all_plain = false;
            }
        }
        if (all_gz)
            return partition_kmers_multi_gz_pc<k, m>(
                cfg, buckets, partition_fn, write_budget_per_thread);
        if (all_plain && n_files <= n_threads_req)
            return partition_kmers_plain_pc<k, m>(
                cfg, buckets, partition_fn, write_budget_per_thread);
    }

    const size_t n_threads = std::min(n_threads_req, n_files);
    const size_t n_parts   = cfg.num_partitions;

    std::atomic<size_t>   next_file{0};
    const size_t writer_shards = std::min<size_t>(4, std::max<size_t>(1, n_threads_req / 4));
    AsyncPartitionWriters async_writers(
        buckets, writer_shards, std::max<size_t>(64u << 20, write_budget_per_thread));
    std::atomic<uint64_t>   total_seqs{0}, total_kmers{0}, total_superkmers{0};
    std::atomic<bool>       stop{false};
    std::exception_ptr      worker_error = nullptr;
    std::mutex              worker_error_mutex;

    auto worker = [&]() {
        try {
            const size_t flush_thresh = writer_flush_threshold(n_parts, write_budget_per_thread);
            SeqSource            source;
            MinimizerWindow<k, m>               min_it;
            std::vector<SuperkmerWriter<k, m>>  writers(n_parts, SuperkmerWriter<k, m>(flush_thresh));
            std::vector<uint8_t>         kache_buf;
            uint64_t local_seqs = 0, local_kmers = 0, local_superkmers = 0;

            auto flush_fn = [&](std::vector<SuperkmerWriter<k, m>>& ws, size_t p) {
                if (ws[p].needs_flush()) async_writers.enqueue(ws[p].release_block(p));
            };

            while (true) {
                if (stop.load(std::memory_order_relaxed)) break;
                const size_t fi = next_file.fetch_add(1, std::memory_order_relaxed);
                if (fi >= n_files) break;

                source.process(cfg.input_files[fi], [&](const char* chunk, size_t len) {
                    if (stop.load(std::memory_order_relaxed)) return;
                    extract_superkmers_from_actg<k, m>(
                        chunk, len, partition_fn, min_it, writers, local_kmers, local_superkmers, flush_fn, kache_buf);
                    ++local_seqs;
                });
            }

            for (size_t p = 0; p < n_parts; ++p) {
                if (!writers[p].empty())
                    async_writers.enqueue(writers[p].release_block(p));
            }

            total_seqs       .fetch_add(local_seqs,        std::memory_order_relaxed);
            total_kmers      .fetch_add(local_kmers,       std::memory_order_relaxed);
            total_superkmers .fetch_add(local_superkmers,  std::memory_order_relaxed);
        } catch (...) {
            {
                std::lock_guard<std::mutex> lk(worker_error_mutex);
                if (!worker_error) worker_error = std::current_exception();
            }
            stop.store(true, std::memory_order_relaxed);
        }
    };

    std::vector<std::thread> threads;
    threads.reserve(n_threads);
    for (size_t t = 0; t < n_threads; ++t)
        threads.emplace_back(worker);
    for (auto& th : threads) th.join();
    async_writers.finish();
    if (worker_error) std::rethrow_exception(worker_error);

    return { total_seqs.load(), total_kmers.load(), total_superkmers.load() };
}


// ─── Parallel harness (in-memory sinks) ───────────────────────────────────────
//
// Same logic as partition_kmers_impl but writes packed superkmers into
// per-partition std::string buffers instead of disk files.  Used by the
// streaming pipeline to avoid the Phase 1 disk write + Phase 2 mmap round-trip.

template <uint16_t k, uint16_t m, typename PartitionFn>
PartitionStats partition_kmers_mem_multi_gz_pc(
    const Config&             cfg,
    std::vector<std::string>& bufs,
    PartitionFn               partition_fn)
{
    using Batch = PackedReadBatch;

    constexpr size_t MAX_QUEUE = 32;
    constexpr size_t MAX_BATCH_ITEMS = 524288;
    constexpr size_t TARGET_BATCH_BASES = 64u << 20;

    const size_t n_files = cfg.input_files.size();
    const size_t n_threads = static_cast<size_t>(cfg.num_threads);
    const size_t n_parts = cfg.num_partitions;
    const size_t n_producers = std::min(n_files, std::max<size_t>(1, n_threads / 2));
    const size_t n_consumers = std::max<size_t>(1, n_threads - n_producers);

    std::atomic<size_t> next_file{0};
    std::atomic<size_t> active_producers{n_producers};
    std::deque<Batch> queue;
    std::mutex q_mutex;
    std::condition_variable q_cv;
    std::atomic<bool> stop{false};
    std::exception_ptr producer_error = nullptr;
    std::exception_ptr consumer_error = nullptr;
    std::mutex error_mutex;

    std::vector<std::mutex> buf_mutexes(n_parts);
    std::atomic<uint64_t> total_seqs{0}, total_kmers{0}, total_superkmers{0};

    auto push_batch = [&](Batch& batch, size_t& batch_bases) {
        if (batch.empty()) return;
        {
            std::unique_lock<std::mutex> lk(q_mutex);
            q_cv.wait(lk, [&]{
                return queue.size() < MAX_QUEUE || stop.load(std::memory_order_relaxed);
            });
            if (stop.load(std::memory_order_relaxed)) return;
            queue.push_back(std::move(batch));
        }
        q_cv.notify_one();
        batch = Batch{};
        batch.reserve(MAX_BATCH_ITEMS, TARGET_BATCH_BASES);
        batch_bases = 0;
    };

    auto producer_done = [&]() {
        active_producers.fetch_sub(1, std::memory_order_relaxed);
        q_cv.notify_all();
    };

    auto producer_fn = [&]() {
        try {
            while (!stop.load(std::memory_order_relaxed)) {
                const size_t fi = next_file.fetch_add(1, std::memory_order_relaxed);
                if (fi >= n_files) break;
                const auto& input_path = cfg.input_files[fi];

                Batch batch;
                batch.reserve(MAX_BATCH_ITEMS, TARGET_BATCH_BASES);
                size_t batch_bases = 0;
                GzInput inp(input_path);
                if (inp.first_byte() == '@') {
                    helicase::FastqParser<HELICASE_ACTG_PACKED, GzInput> p(std::move(inp));
                    while (!stop.load(std::memory_order_relaxed) && p.next()) {
                        const size_t len = p.get_dna_len();
                        if (len >= k) {
                            batch.append(p.get_dna_packed());
                            batch_bases += len;
                            if (batch_bases >= TARGET_BATCH_BASES ||
                                batch.records.size() >= MAX_BATCH_ITEMS)
                                push_batch(batch, batch_bases);
                        }
                    }
                } else {
                    helicase::FastaParser<HELICASE_ACTG_PACKED, GzInput> p(std::move(inp));
                    while (!stop.load(std::memory_order_relaxed) && p.next()) {
                        const size_t len = p.get_dna_len();
                        if (len >= k) {
                            batch.append(p.get_dna_packed());
                            batch_bases += len;
                            if (batch_bases >= TARGET_BATCH_BASES ||
                                batch.records.size() >= MAX_BATCH_ITEMS)
                                push_batch(batch, batch_bases);
                        }
                    }
                }
                if (!stop.load(std::memory_order_relaxed))
                    push_batch(batch, batch_bases);
            }
        } catch (...) {
            {
                std::lock_guard<std::mutex> lk(error_mutex);
                if (!producer_error) producer_error = std::current_exception();
            }
            stop.store(true, std::memory_order_relaxed);
            q_cv.notify_all();
        }
        producer_done();
    };

    auto consumer_fn = [&]() {
        try {
            const size_t flush_thresh = writer_flush_threshold(n_parts, 64u << 20);
            MinimizerWindow<k, m> min_it;
            std::vector<SuperkmerWriter<k, m>> writers(n_parts, SuperkmerWriter<k, m>(flush_thresh));
            uint64_t local_seqs = 0, local_kmers = 0, local_superkmers = 0;

            auto flush_fn = [&](std::vector<SuperkmerWriter<k, m>>& ws, size_t p) {
                if (ws[p].needs_flush()) ws[p].flush_to_mem(bufs[p], buf_mutexes[p]);
            };

            while (true) {
                if (stop.load(std::memory_order_relaxed)) break;
                Batch batch;
                {
                    std::unique_lock<std::mutex> lk(q_mutex);
                    q_cv.wait(lk, [&]{
                        return !queue.empty() ||
                               active_producers.load(std::memory_order_relaxed) == 0 ||
                               stop.load(std::memory_order_relaxed);
                    });
                    if (queue.empty()) break;
                    batch = std::move(queue.front());
                    queue.pop_front();
                }
                q_cv.notify_one();

                for (const auto& rec : batch.records) {
                    if (stop.load(std::memory_order_relaxed)) break;
                    PackedReadView chunk{
                        batch.words.data() + rec.word_offset,
                        rec.word_count,
                        rec.tail,
                        rec.len
                    };
                    extract_superkmers_from_packed_nt<k, m>(
                        chunk, partition_fn,
                        min_it, writers, local_kmers, local_superkmers, flush_fn);
                    ++local_seqs;
                }
            }

            for (size_t p = 0; p < n_parts; ++p)
                writers[p].flush_to_mem(bufs[p], buf_mutexes[p]);

            total_seqs.fetch_add(local_seqs, std::memory_order_relaxed);
            total_kmers.fetch_add(local_kmers, std::memory_order_relaxed);
            total_superkmers.fetch_add(local_superkmers, std::memory_order_relaxed);
        } catch (...) {
            {
                std::lock_guard<std::mutex> lk(error_mutex);
                if (!consumer_error) consumer_error = std::current_exception();
            }
            stop.store(true, std::memory_order_relaxed);
            q_cv.notify_all();
        }
    };

    std::vector<std::thread> threads;
    threads.reserve(n_producers + n_consumers);
    for (size_t t = 0; t < n_producers; ++t)
        threads.emplace_back(producer_fn);
    for (size_t t = 0; t < n_consumers; ++t)
        threads.emplace_back(consumer_fn);
    for (auto& th : threads) th.join();
    if (producer_error) std::rethrow_exception(producer_error);
    if (consumer_error) std::rethrow_exception(consumer_error);

    return { total_seqs.load(), total_kmers.load(), total_superkmers.load() };
}

template <uint16_t k, uint16_t m, typename PartitionFn>
PartitionStats partition_kmers_mem_impl(
    const Config&             cfg,
    std::vector<std::string>& bufs,
    PartitionFn               partition_fn)
{
    const size_t n_files      = cfg.input_files.size();
    const size_t n_threads_req = static_cast<size_t>(cfg.num_threads);
    const size_t n_parts      = cfg.num_partitions;

    // Single file + multiple threads → producer-consumer (reuse the same
    // consumers but flush to bufs instead of ofstreams).
    const size_t n_threads = (n_files == 1 && n_threads_req > 1)
        ? n_threads_req   // single-file producer-consumer: all threads participate
        : std::min(n_threads_req, n_files);

    if (n_files >= 1 && n_threads_req > 1) {
        bool all_gz = true;
        for (const auto& input_path : cfg.input_files) {
            if (!(input_path.size() > 3 &&
                  input_path.compare(input_path.size() - 3, 3, ".gz") == 0)) {
                all_gz = false;
                break;
            }
        }
        if (all_gz)
            return partition_kmers_mem_multi_gz_pc<k, m>(cfg, bufs, partition_fn);
    }

    // Single file: producer-consumer variant using in-memory sinks.
    if (n_files == 1 && n_threads_req > 1) {
        const auto& input_path = cfg.input_files[0];
        const bool is_gz = input_path.size() > 3 &&
                           input_path.compare(input_path.size() - 3, 3, ".gz") == 0;

        if (!is_gz) {
            using Batch = std::vector<std::string_view>;
            constexpr size_t MAX_QUEUE  = 32;
            constexpr size_t MAX_BATCH_ITEMS = 8192;
            constexpr size_t TARGET_BATCH_BASES = 1u << 20;
            constexpr size_t TARGET_CHUNK_BASES = 1u << 20;
            const size_t n_consumers = n_threads - 1;

            std::deque<Batch>       queue;
            std::mutex              q_mutex;
            std::condition_variable q_cv;
            bool                    producer_done = false;
            size_t                  active_batches = 0;
            std::exception_ptr      producer_error = nullptr;
            std::atomic<bool>       stop{false};
            std::exception_ptr      consumer_error = nullptr;
            std::mutex              consumer_error_mutex;
            std::vector<std::mutex> buf_mutexes(n_parts);
            std::atomic<uint64_t>   total_seqs{0}, total_kmers{0}, total_superkmers{0};

            auto producer_fn = [&]() {
                auto push_batch = [&](Batch& batch, size_t& batch_bases) {
                    if (batch.empty()) return;
                    {
                        std::unique_lock<std::mutex> lk(q_mutex);
                        q_cv.wait(lk, [&]{
                            return queue.size() < MAX_QUEUE ||
                                   stop.load(std::memory_order_relaxed);
                        });
                        if (stop.load(std::memory_order_relaxed)) return;
                        queue.push_back(std::move(batch));
                    }
                    q_cv.notify_one();
                    batch = Batch{};
                    batch_bases = 0;
                };

                auto wait_until_consumed = [&]() {
                    std::unique_lock<std::mutex> lk(q_mutex);
                    q_cv.wait(lk, [&]{
                        return (queue.empty() && active_batches == 0) ||
                               stop.load(std::memory_order_relaxed);
                    });
                };

                auto feed = [&](auto& parser) {
                    Batch batch;
                    size_t batch_bases = 0;
                    auto add_chunk = [&](const char* ptr, size_t len) {
                        if (len < k) return;
                        const size_t step = TARGET_CHUNK_BASES > k
                            ? TARGET_CHUNK_BASES - (k - 1)
                            : TARGET_CHUNK_BASES;
                        for (size_t off = 0; off < len; ) {
                            const size_t sub_len = std::min(TARGET_CHUNK_BASES, len - off);
                            batch.emplace_back(ptr + off, sub_len);
                            batch_bases += sub_len;
                            if (batch_bases >= TARGET_BATCH_BASES ||
                                batch.size() >= MAX_BATCH_ITEMS)
                                push_batch(batch, batch_bases);
                            if (off + sub_len >= len) break;
                            off += step;
                        }
                    };

                    while (!stop.load(std::memory_order_relaxed) && parser.next()) {
                        auto [ptr, len] = parser.get_dna_raw();
                        add_chunk(ptr, len);
                    }
                    if (!stop.load(std::memory_order_relaxed))
                        push_batch(batch, batch_bases);
                    {
                        std::lock_guard<std::mutex> lk(q_mutex);
                        producer_done = true;
                    }
                    q_cv.notify_all();
                    wait_until_consumed();
                };

                try {
                    helicase::MmapInput inp(input_path);
                    if (inp.first_byte() == '@') {
                        helicase::FastqParser<HELICASE_ACTG, helicase::MmapInput> p(std::move(inp));
                        feed(p);
                    } else {
                        helicase::FastaParser<HELICASE_ACTG, helicase::MmapInput> p(std::move(inp));
                        feed(p);
                    }
                } catch (...) {
                    {
                        std::lock_guard<std::mutex> lk(q_mutex);
                        producer_error = std::current_exception();
                        producer_done = true;
                    }
                    stop.store(true, std::memory_order_relaxed);
                    q_cv.notify_all();
                }
            };

            auto consumer_fn = [&]() {
                try {
                    const size_t flush_thresh = writer_flush_threshold(n_parts, 64u << 20);
                    MinimizerWindow<k, m>        min_it;
                    std::vector<SuperkmerWriter<k, m>>  writers(n_parts, SuperkmerWriter<k, m>(flush_thresh));
                    std::vector<uint8_t>         kache_buf;
                    uint64_t local_seqs = 0, local_kmers = 0, local_superkmers = 0;
                    auto flush_fn = [&](std::vector<SuperkmerWriter<k, m>>& ws, size_t p) {
                        if (ws[p].needs_flush()) ws[p].flush_to_mem(bufs[p], buf_mutexes[p]);
                    };
                    while (true) {
                        if (stop.load(std::memory_order_relaxed)) break;
                        Batch batch;
                        {
                            std::unique_lock<std::mutex> lk(q_mutex);
                            q_cv.wait(lk, [&]{ return !queue.empty() || producer_done; });
                            if (queue.empty()) break;
                            batch = std::move(queue.front());
                            queue.pop_front();
                            ++active_batches;
                        }
                        q_cv.notify_all();
                        for (const auto chunk : batch) {
                            if (stop.load(std::memory_order_relaxed)) break;
                            extract_superkmers_from_actg<k, m>(
                                chunk.data(), chunk.size(), partition_fn,
                                min_it, writers, local_kmers, local_superkmers, flush_fn, kache_buf);
                            ++local_seqs;
                        }
                        {
                            std::lock_guard<std::mutex> lk(q_mutex);
                            --active_batches;
                        }
                        q_cv.notify_all();
                    }
                    for (size_t p = 0; p < n_parts; ++p)
                        writers[p].flush_to_mem(bufs[p], buf_mutexes[p]);
                    total_seqs       .fetch_add(local_seqs,        std::memory_order_relaxed);
                    total_kmers      .fetch_add(local_kmers,       std::memory_order_relaxed);
                    total_superkmers .fetch_add(local_superkmers,  std::memory_order_relaxed);
                } catch (...) {
                    {
                        std::lock_guard<std::mutex> lk(consumer_error_mutex);
                        if (!consumer_error) consumer_error = std::current_exception();
                    }
                    stop.store(true, std::memory_order_relaxed);
                    q_cv.notify_all();
                }
            };

            std::vector<std::thread> threads;
            threads.reserve(n_threads);
            threads.emplace_back(producer_fn);
            for (size_t t = 0; t < n_consumers; ++t)
                threads.emplace_back(consumer_fn);
            for (auto& th : threads) th.join();
            if (producer_error) std::rethrow_exception(producer_error);
            if (consumer_error) std::rethrow_exception(consumer_error);
            return { total_seqs.load(), total_kmers.load(), total_superkmers.load() };
        }

            using Batch = std::vector<std::string>;
            constexpr size_t MAX_QUEUE  = 32;
            constexpr size_t MAX_BATCH_ITEMS = 8192;
            constexpr size_t TARGET_BATCH_BASES = 1u << 20;
            constexpr size_t TARGET_CHUNK_BASES = 1u << 20;
            const size_t n_consumers = n_threads - 1;

            std::deque<Batch>       queue;
            std::mutex              q_mutex;
            std::condition_variable q_cv;
            bool                    producer_done = false;
            std::exception_ptr      producer_error = nullptr;
            std::atomic<bool>       stop{false};
            std::exception_ptr      consumer_error = nullptr;
            std::mutex              consumer_error_mutex;
            std::vector<std::mutex> buf_mutexes(n_parts);
            std::atomic<uint64_t>   total_seqs{0}, total_kmers{0}, total_superkmers{0};

            auto producer_fn = [&]() {
                auto push_batch = [&](Batch& batch, size_t& batch_bases) {
                    if (batch.empty()) return;
                    {
                        std::unique_lock<std::mutex> lk(q_mutex);
                        q_cv.wait(lk, [&]{
                            return queue.size() < MAX_QUEUE ||
                                   stop.load(std::memory_order_relaxed);
                        });
                        if (stop.load(std::memory_order_relaxed)) return;
                        queue.push_back(std::move(batch));
                    }
                    q_cv.notify_one();
                    batch = Batch{};
                    batch_bases = 0;
                };
                try {
                    SeqSource source;
                    Batch batch;
                    size_t batch_bases = 0;
                    auto add_chunk = [&](const char* ptr, size_t len) {
                        if (len < k) return;
                        const size_t step = TARGET_CHUNK_BASES > k
                            ? TARGET_CHUNK_BASES - (k - 1)
                            : TARGET_CHUNK_BASES;
                        for (size_t off = 0; off < len; ) {
                            const size_t sub_len = std::min(TARGET_CHUNK_BASES, len - off);
                            batch.emplace_back(ptr + off, sub_len);
                            batch_bases += sub_len;
                            if (batch_bases >= TARGET_BATCH_BASES ||
                                batch.size() >= MAX_BATCH_ITEMS)
                                push_batch(batch, batch_bases);
                            if (off + sub_len >= len) break;
                            off += step;
                        }
                    };
                    source.process(input_path, [&](const char* ptr, size_t len) {
                        if (stop.load(std::memory_order_relaxed)) return;
                        add_chunk(ptr, len);
                    });
                    if (!stop.load(std::memory_order_relaxed))
                        push_batch(batch, batch_bases);
                    { std::lock_guard<std::mutex> lk(q_mutex); producer_done = true; }
                    q_cv.notify_all();
                } catch (...) {
                    {
                        std::lock_guard<std::mutex> lk(q_mutex);
                        producer_error = std::current_exception();
                        producer_done = true;
                    }
                    stop.store(true, std::memory_order_relaxed);
                    q_cv.notify_all();
                }
            };

            auto consumer_fn = [&]() {
                try {
                    const size_t flush_thresh = writer_flush_threshold(n_parts, 64u << 20);
                    MinimizerWindow<k, m>        min_it;
                    std::vector<SuperkmerWriter<k, m>>  writers(n_parts, SuperkmerWriter<k, m>(flush_thresh));
                    std::vector<uint8_t>         kache_buf;
                    uint64_t local_seqs = 0, local_kmers = 0, local_superkmers = 0;
                    auto flush_fn = [&](std::vector<SuperkmerWriter<k, m>>& ws, size_t p) {
                        if (ws[p].needs_flush()) ws[p].flush_to_mem(bufs[p], buf_mutexes[p]);
                    };
                    while (true) {
                        if (stop.load(std::memory_order_relaxed)) break;
                        Batch batch;
                        { std::unique_lock<std::mutex> lk(q_mutex);
                          q_cv.wait(lk, [&]{ return !queue.empty() || producer_done; });
                          if (queue.empty()) break;
                          batch = std::move(queue.front());
                          queue.pop_front(); }
                        q_cv.notify_one();
                        for (const auto& chunk : batch) {
                            if (stop.load(std::memory_order_relaxed)) break;
                            extract_superkmers_from_actg<k, m>(
                                chunk.data(), chunk.size(), partition_fn,
                                min_it, writers, local_kmers, local_superkmers, flush_fn, kache_buf);
                            ++local_seqs;
                        }
                    }
                    for (size_t p = 0; p < n_parts; ++p)
                        writers[p].flush_to_mem(bufs[p], buf_mutexes[p]);
                    total_seqs       .fetch_add(local_seqs,        std::memory_order_relaxed);
                    total_kmers      .fetch_add(local_kmers,       std::memory_order_relaxed);
                    total_superkmers .fetch_add(local_superkmers,  std::memory_order_relaxed);
                } catch (...) {
                    {
                        std::lock_guard<std::mutex> lk(consumer_error_mutex);
                        if (!consumer_error) consumer_error = std::current_exception();
                    }
                    stop.store(true, std::memory_order_relaxed);
                    q_cv.notify_all();
                }
            };

            std::vector<std::thread> threads;
            threads.reserve(n_threads);
            threads.emplace_back(producer_fn);
            for (size_t t = 0; t < n_consumers; ++t)
                threads.emplace_back(consumer_fn);
            for (auto& th : threads) th.join();
            if (producer_error) std::rethrow_exception(producer_error);
            if (consumer_error) std::rethrow_exception(consumer_error);
            return { total_seqs.load(), total_kmers.load(), total_superkmers.load() };
    }

    // Few large plain files: split parsed records/chunks across all worker
    // threads instead of assigning whole files to workers.  The producer keeps
    // each mmap-backed parser alive until consumers have drained all views for
    // that file.
    if (n_files > 1 && n_files <= n_threads_req && n_threads_req > 1) {
        bool all_plain = true;
        for (const auto& input_path : cfg.input_files) {
            if (input_path.size() > 3 &&
                input_path.compare(input_path.size() - 3, 3, ".gz") == 0) {
                all_plain = false;
                break;
            }
        }

        if (all_plain) {
            using Batch = std::vector<std::string_view>;
            constexpr size_t MAX_QUEUE  = 32;
            constexpr size_t MAX_BATCH_ITEMS = 8192;
            constexpr size_t TARGET_BATCH_BASES = 1u << 20;
            constexpr size_t TARGET_CHUNK_BASES = 1u << 20;
            const size_t n_consumers = n_threads_req - 1;

            std::deque<Batch>       queue;
            std::mutex              q_mutex;
            std::condition_variable q_cv;
            bool                    producer_done = false;
            size_t                  active_batches = 0;
            std::exception_ptr      producer_error = nullptr;
            std::atomic<bool>       stop{false};
            std::exception_ptr      consumer_error = nullptr;
            std::mutex              consumer_error_mutex;
            std::vector<std::mutex> buf_mutexes(n_parts);
            std::atomic<uint64_t>   total_seqs{0}, total_kmers{0}, total_superkmers{0};

            auto producer_fn = [&]() {
                auto push_batch = [&](Batch& batch, size_t& batch_bases) {
                    if (batch.empty()) return;
                    {
                        std::unique_lock<std::mutex> lk(q_mutex);
                        q_cv.wait(lk, [&]{
                            return queue.size() < MAX_QUEUE ||
                                   stop.load(std::memory_order_relaxed);
                        });
                        if (stop.load(std::memory_order_relaxed)) return;
                        queue.push_back(std::move(batch));
                    }
                    q_cv.notify_one();
                    batch = Batch{};
                    batch_bases = 0;
                };

                auto wait_until_consumed = [&]() {
                    std::unique_lock<std::mutex> lk(q_mutex);
                    q_cv.wait(lk, [&]{
                        return (queue.empty() && active_batches == 0) ||
                               stop.load(std::memory_order_relaxed);
                    });
                };

                auto feed = [&](auto& parser) {
                    Batch batch;
                    size_t batch_bases = 0;
                    auto add_chunk = [&](const char* ptr, size_t len) {
                        if (len < k) return;
                        const size_t step = TARGET_CHUNK_BASES > k
                            ? TARGET_CHUNK_BASES - (k - 1)
                            : TARGET_CHUNK_BASES;
                        for (size_t off = 0; off < len; ) {
                            const size_t sub_len = std::min(TARGET_CHUNK_BASES, len - off);
                            batch.emplace_back(ptr + off, sub_len);
                            batch_bases += sub_len;
                            if (batch_bases >= TARGET_BATCH_BASES ||
                                batch.size() >= MAX_BATCH_ITEMS)
                                push_batch(batch, batch_bases);
                            if (off + sub_len >= len) break;
                            off += step;
                        }
                    };

                    while (!stop.load(std::memory_order_relaxed) && parser.next()) {
                        auto [ptr, len] = parser.get_dna_raw();
                        add_chunk(ptr, len);
                    }
                    if (!stop.load(std::memory_order_relaxed))
                        push_batch(batch, batch_bases);
                    wait_until_consumed();
                };

                try {
                    for (const auto& input_path : cfg.input_files) {
                        if (stop.load(std::memory_order_relaxed)) break;
                        helicase::MmapInput inp(input_path);
                        if (inp.first_byte() == '@') {
                            helicase::FastqParser<HELICASE_ACTG, helicase::MmapInput> p(std::move(inp));
                            feed(p);
                        } else {
                            helicase::FastaParser<HELICASE_ACTG, helicase::MmapInput> p(std::move(inp));
                            feed(p);
                        }
                    }
                    {
                        std::lock_guard<std::mutex> lk(q_mutex);
                        producer_done = true;
                    }
                    q_cv.notify_all();
                } catch (...) {
                    {
                        std::lock_guard<std::mutex> lk(q_mutex);
                        producer_error = std::current_exception();
                        producer_done = true;
                    }
                    stop.store(true, std::memory_order_relaxed);
                    q_cv.notify_all();
                }
            };

            auto consumer_fn = [&]() {
                try {
                    const size_t flush_thresh = writer_flush_threshold(n_parts, 64u << 20);
                    MinimizerWindow<k, m>        min_it;
                    std::vector<SuperkmerWriter<k, m>>  writers(n_parts, SuperkmerWriter<k, m>(flush_thresh));
                    std::vector<uint8_t>         kache_buf;
                    uint64_t local_seqs = 0, local_kmers = 0, local_superkmers = 0;
                    auto flush_fn = [&](std::vector<SuperkmerWriter<k, m>>& ws, size_t p) {
                        if (ws[p].needs_flush()) ws[p].flush_to_mem(bufs[p], buf_mutexes[p]);
                    };
                    while (true) {
                        if (stop.load(std::memory_order_relaxed)) break;
                        Batch batch;
                        {
                            std::unique_lock<std::mutex> lk(q_mutex);
                            q_cv.wait(lk, [&]{ return !queue.empty() || producer_done; });
                            if (queue.empty()) break;
                            batch = std::move(queue.front());
                            queue.pop_front();
                            ++active_batches;
                        }
                        q_cv.notify_all();
                        for (const auto chunk : batch) {
                            if (stop.load(std::memory_order_relaxed)) break;
                            extract_superkmers_from_actg<k, m>(
                                chunk.data(), chunk.size(), partition_fn,
                                min_it, writers, local_kmers, local_superkmers, flush_fn, kache_buf);
                            ++local_seqs;
                        }
                        {
                            std::lock_guard<std::mutex> lk(q_mutex);
                            --active_batches;
                        }
                        q_cv.notify_all();
                    }
                    for (size_t p = 0; p < n_parts; ++p)
                        writers[p].flush_to_mem(bufs[p], buf_mutexes[p]);
                    total_seqs       .fetch_add(local_seqs,        std::memory_order_relaxed);
                    total_kmers      .fetch_add(local_kmers,       std::memory_order_relaxed);
                    total_superkmers .fetch_add(local_superkmers,  std::memory_order_relaxed);
                } catch (...) {
                    {
                        std::lock_guard<std::mutex> lk(consumer_error_mutex);
                        if (!consumer_error) consumer_error = std::current_exception();
                    }
                    stop.store(true, std::memory_order_relaxed);
                    q_cv.notify_all();
                }
            };

            std::vector<std::thread> threads;
            threads.reserve(n_threads_req);
            threads.emplace_back(producer_fn);
            for (size_t t = 0; t < n_consumers; ++t)
                threads.emplace_back(consumer_fn);
            for (auto& th : threads) th.join();
            if (producer_error) std::rethrow_exception(producer_error);
            if (consumer_error) std::rethrow_exception(consumer_error);
            return { total_seqs.load(), total_kmers.load(), total_superkmers.load() };
        }
    }

    // Multi-file: file-level work-stealing.
    std::atomic<size_t>     next_file{0};
    std::vector<std::mutex> buf_mutexes(n_parts);
    std::atomic<uint64_t>   total_seqs{0}, total_kmers{0}, total_superkmers{0};
    std::atomic<bool>       stop{false};
    std::exception_ptr      worker_error = nullptr;
    std::mutex              worker_error_mutex;

    auto worker = [&]() {
        try {
            const size_t flush_thresh = writer_flush_threshold(n_parts, 64u << 20);
            SeqSource            source;
            MinimizerWindow<k, m>               min_it;
            std::vector<SuperkmerWriter<k, m>>  writers(n_parts, SuperkmerWriter<k, m>(flush_thresh));
            std::vector<uint8_t>         kache_buf;
            uint64_t local_seqs = 0, local_kmers = 0, local_superkmers = 0;

            auto flush_fn = [&](std::vector<SuperkmerWriter<k, m>>& ws, size_t p) {
                if (ws[p].needs_flush()) ws[p].flush_to_mem(bufs[p], buf_mutexes[p]);
            };

            while (true) {
                if (stop.load(std::memory_order_relaxed)) break;
                const size_t fi = next_file.fetch_add(1, std::memory_order_relaxed);
                if (fi >= n_files) break;

                source.process(cfg.input_files[fi], [&](const char* chunk, size_t len) {
                    if (stop.load(std::memory_order_relaxed)) return;
                    extract_superkmers_from_actg<k, m>(
                        chunk, len, partition_fn, min_it, writers, local_kmers, local_superkmers, flush_fn, kache_buf);
                    ++local_seqs;
                });
            }

            for (size_t p = 0; p < n_parts; ++p)
                writers[p].flush_to_mem(bufs[p], buf_mutexes[p]);

            total_seqs       .fetch_add(local_seqs,        std::memory_order_relaxed);
            total_kmers      .fetch_add(local_kmers,       std::memory_order_relaxed);
            total_superkmers .fetch_add(local_superkmers,  std::memory_order_relaxed);
        } catch (...) {
            {
                std::lock_guard<std::mutex> lk(worker_error_mutex);
                if (!worker_error) worker_error = std::current_exception();
            }
            stop.store(true, std::memory_order_relaxed);
        }
    };

    std::vector<std::thread> threads;
    threads.reserve(n_threads);
    for (size_t t = 0; t < n_threads; ++t)
        threads.emplace_back(worker);
    for (auto& th : threads) th.join();
    if (worker_error) std::rethrow_exception(worker_error);
    return { total_seqs.load(), total_kmers.load(), total_superkmers.load() };
}


// ─── Public ───────────────────────────────────────────────────────────────────

template <uint16_t k, uint16_t m>
PartitionStats partition_kmers(
    const Config&               cfg,
    std::vector<std::ofstream>& buckets,
    size_t                      write_budget_per_thread)
{
    const size_t n    = cfg.num_partitions;
    const size_t mask = n - 1; // n is always a power of 2 (enforced in main.cpp)
    return partition_kmers_impl<k, m>(cfg, buckets,
        [mask](uint64_t h) -> size_t { return h & mask; },
        write_budget_per_thread);
}

template <uint16_t k, uint16_t m>
PartitionStats partition_kmers_mem(
    const Config&             cfg,
    std::vector<std::string>& bufs)
{
    const size_t n    = cfg.num_partitions;
    const size_t mask = n - 1; // n is always a power of 2 (enforced in main.cpp)
    return partition_kmers_mem_impl<k, m>(cfg, bufs,
        [mask](uint64_t h) -> size_t { return h & mask; });
}
