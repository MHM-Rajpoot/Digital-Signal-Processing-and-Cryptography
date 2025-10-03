// toy_hashes.cpp
// Educational toy implementations of Merkle-Damgard and Sponge-style hashes, plus birthday demo.
// NOT SECURE — for teaching only (Dr Markus Kuhn style).
//
// Compile: g++ -std=c++17 -O2 toy_hashes.cpp -o toy_hashes
// Run:     ./toy_hashes

#include <bits/stdc++.h>

using u8  = uint8_t;
using u32 = uint32_t;
using u64 = uint64_t;
using bytes = std::vector<u8>;

static std::string to_hex(const bytes &b) {
    std::ostringstream oss;
    oss << std::hex << std::setfill('0');
    for (u8 x : b) oss << std::setw(2) << (int)x;
    return oss.str();
}

// ------------------------- Merkle-Damgard-style toy hash -------------------------
// - Block size 8 bytes (64 bits) for demonstration
// - Chaining value 64-bit (8 bytes)
// - Compression function: simple mixing using xor, rotate, multiply (NOT cryptographically secure).
//
// Padding: MD-style pad: append 0x80, zeros, then 64-bit big-endian length in bits.

static inline u64 rotl64(u64 x, unsigned r) { return (x << r) | (x >> (64 - r)); }

u64 toy_compress(u64 h, u64 m_block) {
    // Very small toy compression: mix h and m_block
    // Not secure; only for didactic purpose.
    u64 x = h ^ 0x0123456789ABCDEFULL;
    x = x + m_block;
    x ^= rotl64(x, 23) ^ rotl64(m_block, 17);
    x = x * 0x9e3779b97f4a7c15ULL; // golden ratio multiplier
    x ^= (x >> 33);
    x += 0xA5A5A5A5A5A5A5A5ULL;
    return x;
}

bytes md_pad(const bytes &m, size_t block_bytes = 8) {
    // Merkle-Damgard pad: 0x80, zeros, 64-bit length in bits at end
    bytes p = m;
    uint64_t bitlen = (uint64_t)m.size() * 8ULL;
    // append 0x80
    p.push_back(0x80);
    // pad with zeros until final 8 bytes can contain length and total is multiple of block_bytes
    while ((p.size() + 8) % block_bytes != 0) p.push_back(0x00);
    // append 64-bit big-endian length in bits
    for (int i = 7; i >= 0; --i) p.push_back((u8)((bitlen >> (8*i)) & 0xFF));
    return p;
}

bytes merkle_damgard_hash(const bytes &message) {
    // IV: arbitrary constant
    u64 h = 0x0123456789ABCDEFULL;
    bytes padded = md_pad(message, 8);
    assert(padded.size() % 8 == 0);
    for (size_t off = 0; off < padded.size(); off += 8) {
        u64 m = 0;
        for (int i = 0; i < 8; ++i) m = (m << 8) | padded[off + i];
        h = toy_compress(h, m);
    }
    // produce 8-byte digest (big-endian)
    bytes out(8);
    for (int i = 7; i >= 0; --i) { out[7 - i] = (u8)((h >> (8*i)) & 0xFF); }
    return out;
}

// ------------------------- Sponge-style toy hash -------------------------
// - State: 5 * 8-byte words = 40 bytes (320 bits)
// - Rate: r bytes (we will use r = 16 bytes), capacity = b - r
// - Simple permutation: a few mix rounds (not cryptographically secure).
//
// Duplex-like behavior: absorb input blocks (rate-sized), apply permutation, then squeeze.

struct SpongeState {
    static const size_t W = 5; // number of 8-byte lanes
    u64 s[W];
    SpongeState() { for (size_t i=0;i<W;++i) s[i] = 0; }
};

void toy_permutation(SpongeState &st) {
    // A small mixing permutation across 5 lanes: rotate, xor, multiply
    for (int r = 0; r < 6; ++r) {
        // linear mixing
        u64 t0 = st.s[0] ^ (st.s[1] + 0x9e3779b97f4a7c15ULL);
        u64 t1 = st.s[1] ^ rotl64(st.s[2], 13);
        u64 t2 = st.s[2] ^ rotl64(st.s[3], 17);
        u64 t3 = st.s[3] ^ rotl64(st.s[4], 19);
        u64 t4 = st.s[4] ^ rotl64(st.s[0], 23);
        st.s[0] = (t0 * 0xC2B2AE3D27D4EB4FULL) ^ rotl64(t0, 7);
        st.s[1] = (t1 * 0x9E3779B185EBCA87ULL) ^ rotl64(t1, 11);
        st.s[2] = (t2 * 0xC3A5C85C97CB3127ULL) ^ rotl64(t2, 19);
        st.s[3] = (t3 * 0xB492B66FBE98F273ULL) ^ rotl64(t3, 3);
        st.s[4] = (t4 * 0xD6E8FEB86659FD93ULL) ^ rotl64(t4, 29);
        // small cross-xor
        st.s[0] ^= st.s[3];
        st.s[1] ^= st.s[4];
        st.s[2] ^= st.s[0];
        st.s[3] ^= st.s[1];
        st.s[4] ^= st.s[2];
    }
}

bytes sponge_hash(const bytes &message, size_t digest_bytes = 32) {
    const size_t rate = 16;      // bytes absorbed per round (r)
    SpongeState st;
    // pad the message with pad10*1 (simple: 0x01 then zeros to align to rate)
    bytes p = message;
    p.push_back(0x01);
    while (p.size() % rate != 0) p.push_back(0x00);
    // Absorb
    for (size_t off = 0; off < p.size(); off += rate) {
        // XOR rate bytes into the first lanes (little mapping)
        for (size_t i = 0; i < rate/8; ++i) {
            u64 lane = 0;
            for (int j = 0; j < 8; ++j) lane = (lane << 8) | p[off + i*8 + j];
            st.s[i] ^= lane;
        }
        toy_permutation(st);
    }
    // Squeeze
    bytes out;
    while (out.size() < digest_bytes) {
        // read rate bytes: first (rate/8) lanes
        for (size_t i = 0; i < rate/8; ++i) {
            u64 lane = st.s[i];
            for (int j = 7; j >= 0; --j) {
                if (out.size() >= digest_bytes) break;
                out.push_back((u8)((lane >> (8*j)) & 0xFF));
            }
        }
        if (out.size() < digest_bytes) toy_permutation(st);
    }
    return out;
}

// ------------------------- Birthday collision demo -------------------------
// Demonstrate that for d-bit digests, collisions are likely after ~2^{d/2} sample outputs.

void birthday_demo_truncated_md(size_t digest_bits_trunc, size_t samples = 1 << 16) {
    assert(digest_bits_trunc <= 64);
    std::mt19937_64 rng(123456789);
    std::unordered_map<uint64_t, size_t> seen;
    size_t collisions = 0;
    for (size_t i = 0; i < samples; ++i) {
        // sample random message (small)
        bytes msg(8);
        for (int j=0;j<8;++j) msg[j] = (u8)(rng() & 0xFF);
        bytes d = merkle_damgard_hash(msg);
        // interpret first digest_bits_trunc bits
        uint64_t val = 0;
        for (size_t k = 0; k < d.size() && k < 8; ++k) val = (val << 8) | d[k];
        if (digest_bits_trunc < 64) val >>= (64 - digest_bits_trunc);
        if (seen.find(val) != seen.end()) collisions++;
        seen[val] = i;
    }
    double expected = 1.0 - exp(- (double)samples * (samples - 1) / (2.0 * (1ULL << digest_bits_trunc)));
    std::cout << "[Birthday demo] " << digest_bits_trunc << "-bit truncated digest, samples=" << samples
              << ", collisions=" << collisions << ", approx_prob=" << expected << "\n";
}

// ------------------------- Entropy pool sketch -------------------------
// Very simple entropy pool: fold inputs into a small internal state using toy hash. Not cryptographically robust.

struct EntropyPool {
    bytes pool; // collect bytes
    EntropyPool(size_t capacity = 64) : pool() { pool.reserve(capacity); }
    void add(const bytes &b) {
        pool.insert(pool.end(), b.begin(), b.end());
        if (pool.size() > 1024) pool.erase(pool.begin(), pool.begin() + (pool.size() - 1024));
    }
    // extract: compress pool to seed using sponge
    bytes extract_seed(size_t nbytes) {
        bytes seed = sponge_hash(pool, nbytes);
        // clear pool (in a real design you'd keep some state and reseed)
        pool.clear();
        return seed;
    }
};

// ------------------------- Demonstration main -------------------------

int main() {
    std::cout << "Toy hash functions demo (Merkle-Damgard and Sponge)\n\n";

    std::string text = "The quick brown fox jumps over the lazy dog";
    bytes msg(text.begin(), text.end());

    // Merkle-Damgard toy hash
    bytes md_digest = merkle_damgard_hash(msg);
    std::cout << "Merkle-Damgard (toy) digest (8 bytes hex): " << to_hex(md_digest) << "\n";

    // Sponge toy hash (32-byte output)
    bytes sponge_digest = sponge_hash(msg, 32);
    std::cout << "Sponge (toy) digest (32 bytes hex): " << to_hex(sponge_digest) << "\n";

    // Compare with small change to show avalanche behavior (toy)
    bytes msg2 = msg;
    if (!msg2.empty()) msg2[0] ^= 1;
    bytes md2 = merkle_damgard_hash(msg2);
    bytes sp2 = sponge_hash(msg2, 32);
    std::cout << "Change one bit -> Merkle-Damgard digest: " << to_hex(md2) << "\n";
    std::cout << "Change one bit -> Sponge digest:         " << to_hex(sp2) << "\n";

    // Birthday demos: truncated digests to illustrate birthday bound
    std::cout << "\nBirthday collision demos (toy truncated digests):\n";
    birthday_demo_truncated_md(16, 1 << 12);  // 16-bit truncated digest, many samples
    birthday_demo_truncated_md(24, 1 << 16);  // 24-bit truncated digest

    // Entropy pool example
    EntropyPool pool;
    pool.add(msg);
    pool.add(sponge_digest);
    bytes seed = pool.extract_seed(16);
    std::cout << "\nEntropy pool extracted seed (16 bytes): " << to_hex(seed) << "\n";

    std::cout << "\nNotes:\n";
    std::cout << "- These toy functions illustrate structure & padding. They are NOT secure.\n";
    std::cout << "- Merkle-Damgard requires careful use: e.g. length-extension attacks exist on plain MD constructions.\n";
    std::cout << "- Sponge construction (Keccak/SHA-3) is widely used; it avoids length-extension and supports variable-length output.\n";
    std::cout << "- In practice prefer SHA-256 / SHA-3 / BLAKE2 / BLAKE3 and vetted implementations.\n";

    return 0;
}
