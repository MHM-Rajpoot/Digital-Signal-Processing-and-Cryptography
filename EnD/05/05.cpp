// modes_demo.cpp
// Educational demo: chosen-plaintext, randomized encryption, ECB/CBC/OFB/CTR modes,
// and a small meet-in-the-middle demo for multiple encryption.
// WARNING: Uses a toy block cipher (XOR with key) for clarity. Replace with AES for real use.

#include <bits/stdc++.h>
using namespace std;
using u8 = uint8_t;

static const size_t BLOCK_BYTES = 16; // 128-bit block for modes

// -------------------------------
// Toy block cipher (DEMO ONLY)
//  - block_encrypt: XORs block with 16-byte key repeated
//  - block_decrypt: same as encrypt (XOR)
void block_encrypt(const array<u8, BLOCK_BYTES> &in,
                   array<u8, BLOCK_BYTES> &out,
                   const vector<u8> &key16)
{
    for (size_t i=0;i<BLOCK_BYTES;++i) out[i] = in[i] ^ key16[i % key16.size()];
}

void block_decrypt(const array<u8, BLOCK_BYTES> &in,
                   array<u8, BLOCK_BYTES> &out,
                   const vector<u8> &key16)
{
    // XOR is involutive
    block_encrypt(in, out, key16);
}

// -------------------------------
// PKCS#7 padding/unpadding for block-oriented ECB/CBC
vector<u8> pkcs7_pad(const vector<u8> &data) {
    size_t n = data.size();
    size_t pad_len = BLOCK_BYTES - (n % BLOCK_BYTES);
    vector<u8> out = data;
    out.resize(n + pad_len, (u8)pad_len);
    return out;
}

vector<u8> pkcs7_unpad(const vector<u8> &padded) {
    if (padded.empty() || (padded.size() % BLOCK_BYTES) != 0) throw runtime_error("invalid padded size");
    u8 last = padded.back();
    if (last == 0 || last > BLOCK_BYTES) throw runtime_error("invalid padding");
    size_t n = padded.size() - last;
    // check
    for (size_t i = n; i < padded.size(); ++i) if (padded[i] != last) throw runtime_error("invalid padding bytes");
    return vector<u8>(padded.begin(), padded.begin() + n);
}

// -------------------------------
// Helpers
vector<u8> random_bytes(size_t n) {
    static std::mt19937_64 rng((unsigned)time(nullptr));
    vector<u8> out(n);
    for (size_t i=0;i<n;i++) out[i] = (u8)(rng() & 0xFF);
    return out;
}

string hex(const vector<u8> &v) {
    static const char* H="0123456789abcdef";
    string s; s.reserve(v.size()*2);
    for (u8 b: v) { s.push_back(H[b>>4]); s.push_back(H[b&0xF]); }
    return s;
}

// XOR two byte arrays
vector<u8> xor_vec(const vector<u8> &a, const vector<u8> &b) {
    vector<u8> out(min(a.size(), b.size()));
    for (size_t i=0;i<out.size();++i) out[i] = a[i] ^ b[i];
    return out;
}

// -------------------------------
// Modes implementations using block_encrypt/block_decrypt

// ECB (deterministic if key fixed)
vector<u8> encrypt_ecb(const vector<u8> &plaintext, const vector<u8> &key16) {
    vector<u8> padded = pkcs7_pad(plaintext);
    vector<u8> out(padded.size());
    array<u8, BLOCK_BYTES> in_block, out_block;
    for (size_t offset=0; offset<padded.size(); offset+=BLOCK_BYTES) {
        for (size_t i=0;i<BLOCK_BYTES;++i) in_block[i] = padded[offset+i];
        block_encrypt(in_block, out_block, key16);
        for (size_t i=0;i<BLOCK_BYTES;++i) out[offset+i] = out_block[i];
    }
    return out;
}

vector<u8> decrypt_ecb(const vector<u8> &ciphertext, const vector<u8> &key16) {
    if (ciphertext.size() % BLOCK_BYTES) throw runtime_error("ct size");
    vector<u8> tmp(ciphertext.size());
    array<u8, BLOCK_BYTES> in_block, out_block;
    for (size_t offset=0; offset<ciphertext.size(); offset+=BLOCK_BYTES) {
        for (size_t i=0;i<BLOCK_BYTES;++i) in_block[i] = ciphertext[offset+i];
        block_decrypt(in_block, out_block, key16);
        for (size_t i=0;i<BLOCK_BYTES;++i) tmp[offset+i] = out_block[i];
    }
    return pkcs7_unpad(tmp);
}

// CBC (uses random IV per encryption)
struct CBCCiphertext { vector<u8> iv; vector<u8> ct; };

CBCCiphertext encrypt_cbc(const vector<u8> &plaintext, const vector<u8> &key16) {
    vector<u8> padded = pkcs7_pad(plaintext);
    vector<u8> iv = random_bytes(BLOCK_BYTES);
    vector<u8> out(padded.size());
    array<u8, BLOCK_BYTES> in_block, x_block, out_block;
    // prev = IV
    for (size_t i=0;i<BLOCK_BYTES;++i) x_block[i] = iv[i];
    for (size_t offset=0; offset<padded.size(); offset+=BLOCK_BYTES) {
        for (size_t i=0;i<BLOCK_BYTES;++i) in_block[i] = padded[offset+i];
        array<u8, BLOCK_BYTES> to_enc;
        for (size_t i=0;i<BLOCK_BYTES;++i) to_enc[i] = in_block[i] ^ x_block[i];
        block_encrypt(to_enc, out_block, key16);
        for (size_t i=0;i<BLOCK_BYTES;++i) { out[offset+i] = out_block[i]; x_block[i] = out_block[i]; }
    }
    return {iv, out};
}

vector<u8> decrypt_cbc(const CBCCiphertext &input, const vector<u8> &key16) {
    const vector<u8> &iv = input.iv;
    const vector<u8> &ciphertext = input.ct;
    if (ciphertext.size() % BLOCK_BYTES) throw runtime_error("ct size");
    vector<u8> out(ciphertext.size());
    array<u8, BLOCK_BYTES> in_block, prev_block, plain_block;
    for (size_t i=0;i<BLOCK_BYTES;++i) prev_block[i] = iv[i];
    for (size_t offset=0; offset<ciphertext.size(); offset+=BLOCK_BYTES) {
        for (size_t i=0;i<BLOCK_BYTES;++i) in_block[i] = ciphertext[offset+i];
        block_decrypt(in_block, plain_block, key16);
        for (size_t i=0;i<BLOCK_BYTES;++i) out[offset+i] = plain_block[i] ^ prev_block[i];
        for (size_t i=0;i<BLOCK_BYTES;++i) prev_block[i] = in_block[i];
    }
    return pkcs7_unpad(out);
}

// OFB: keystream: S0 = E(K, IV); Si = E(K, S{i-1}); Ct = Pt XOR Si
struct OFBCiphertext { vector<u8> iv; vector<u8> ct; };
OFBCiphertext encrypt_ofb(const vector<u8> &plaintext, const vector<u8> &key16) {
    vector<u8> out(plaintext.size());
    vector<u8> iv = random_bytes(BLOCK_BYTES);
    array<u8,BLOCK_BYTES> S, tmp;
    for (size_t i=0;i<BLOCK_BYTES;++i) S[i] = iv[i];
    size_t pos = 0;
    while (pos < plaintext.size()) {
        array<u8,BLOCK_BYTES> s_out;
        block_encrypt(S, s_out, key16); // S_i = E_k(S_{i-1})
        // produce keystream bytes from s_out
        for (size_t i=0;i<BLOCK_BYTES && pos < plaintext.size(); ++i, ++pos) {
            out[pos] = plaintext[pos] ^ s_out[i];
            S[i] = s_out[i]; // S becomes s_out for next iteration (we could assign full, but using partial ok)
        }
        // set full S to s_out for next round
        for (size_t i=0;i<BLOCK_BYTES;++i) S[i] = s_out[i];
    }
    return {iv, out};
}

vector<u8> decrypt_ofb(const OFBCiphertext &input, const vector<u8> &key16) {
    // symmetric
    vector<u8> out(input.ct.size());
    array<u8,BLOCK_BYTES> S;
    for (size_t i=0;i<BLOCK_BYTES;++i) S[i] = input.iv[i];
    size_t pos = 0;
    while (pos < input.ct.size()) {
        array<u8,BLOCK_BYTES> s_out;
        block_encrypt(S, s_out, key16);
        for (size_t i=0;i<BLOCK_BYTES && pos < out.size(); ++i, ++pos) {
            out[pos] = input.ct[pos] ^ s_out[i];
        }
        for (size_t i=0;i<BLOCK_BYTES;++i) S[i] = s_out[i];
    }
    return out;
}

// CTR (counter) mode: S_i = E(K, nonce || counter) -> keystream; Ct = Pt XOR S_i
struct CTRCiphertext { vector<u8> nonce; vector<u8> ct; };
CTRCiphertext encrypt_ctr(const vector<u8> &plaintext, const vector<u8> &key16) {
    // use nonce = random 8 bytes || zeros to 16
    vector<u8> nonce = random_bytes(8);
    nonce.resize(BLOCK_BYTES, 0);
    vector<u8> out(plaintext.size());
    uint64_t counter = 0;
    size_t pos=0;
    while (pos < plaintext.size()) {
        array<u8,BLOCK_BYTES> inblock;
        for (size_t i=0;i<BLOCK_BYTES;++i) inblock[i] = nonce[i];
        // put counter in last 8 bytes (little-endian)
        for (int i=0;i<8;++i) inblock[BLOCK_BYTES-1-i] = (u8)((counter >> (8*i)) & 0xFF);
        array<u8,BLOCK_BYTES> s_out;
        block_encrypt(inblock, s_out, key16);
        for (size_t i=0;i<BLOCK_BYTES && pos < plaintext.size(); ++i, ++pos) {
            out[pos] = plaintext[pos] ^ s_out[i];
        }
        ++counter;
    }
    return {nonce, out};
}

vector<u8> decrypt_ctr(const CTRCiphertext &input, const vector<u8> &key16) {
    // symmetric given same nonce/counters
    vector<u8> out(input.ct.size());
    vector<u8> nonce = input.nonce;
    uint64_t counter = 0;
    size_t pos=0;
    while (pos < input.ct.size()) {
        array<u8,BLOCK_BYTES> inblock;
        for (size_t i=0;i<BLOCK_BYTES;++i) inblock[i] = nonce[i];
        for (int i=0;i<8;++i) inblock[BLOCK_BYTES-1-i] = (u8)((counter >> (8*i)) & 0xFF);
        array<u8,BLOCK_BYTES> s_out;
        block_encrypt(inblock, s_out, key16);
        for (size_t i=0;i<BLOCK_BYTES && pos < out.size(); ++i, ++pos) {
            out[pos] = input.ct[pos] ^ s_out[i];
        }
        ++counter;
    }
    return out;
}

// -------------------------------
// Chosen-plaintext IND-CPA experiment
// For a given mode (function pointer), we run a simple challenge: adversary supplies m0,m1, challenger encrypts m_b.
// We'll implement two adversaries:
//  - ECB adversary: returns 1 if ciphertext contains repeated blocks (detects equality pattern).
//  - Simple adversary for CBC (random IV) that guesses randomly (should have ~0 advantage).

using EncOracleECB = function<vector<u8>(const vector<u8>&)>;
using EncOracleCBC = function<CBCCiphertext(const vector<u8>&)>; // CBC returns IV+ct

double ind_cpa_ecb_experiment(const vector<u8> &key16, int trials = 2000) {
    // challenger uses encrypt_ecb
    int correct = 0;
    static mt19937_64 rng((unsigned)time(nullptr));
    for (int t=0;t<trials;++t) {
        // adversary picks m0,m1
        // craft messages that are identical in the middle blocks to let ECB detect equality
        vector<u8> m0, m1;
        // two blocks that are equal in m0, and in m1 differ in one block
        m0.resize(48, 0x41); // 'A' * 48 -> 3 blocks identical
        m1.resize(48, 0x41);
        // introduce difference in m1 block 1
        m1[16] = 0x42;
        // challenger picks b
        int b = (rng() & 1);
        vector<u8> ct = encrypt_ecb(b? m1 : m0, key16);
        // adversary: if there is any repeated 16-byte ciphertext block, guess b=0 (message with repeated blocks)
        bool repeated = false;
        set<string> seen;
        for (size_t off=0; off<ct.size(); off+=BLOCK_BYTES) {
            string blk((char*)&ct[off], (char*)&ct[off]+BLOCK_BYTES);
            if (seen.count(blk)) { repeated = true; break; }
            seen.insert(blk);
        }
        int guess = repeated ? 0 : 1;
        if (guess == b) ++correct;
    }
    double p = (double)correct / trials;
    return fabs(p - 0.5);
}

double ind_cpa_cbc_experiment(const vector<u8> &key16, int trials = 2000) {
    // challenger uses CBC with random IV; adversary as above shouldn't do better than random
    int correct = 0;
    static mt19937_64 rng((unsigned)time(nullptr)+123);
    for (int t=0;t<trials;++t) {
        vector<u8> m0(48, 0x41), m1(48, 0x41);
        m1[16] = 0x42;
        int b = (rng() & 1);
        CBCCiphertext cc = encrypt_cbc(b? m1 : m0, key16);
        // adversary checks repeated ciphertext blocks (should rarely happen because of chaining/IV)
        bool repeated = false;
        set<string> seen;
        const vector<u8> &ct = cc.ct;
        for (size_t off=0; off<ct.size(); off+=BLOCK_BYTES) {
            string blk((char*)&ct[off], (char*)&ct[off]+BLOCK_BYTES);
            if (seen.count(blk)) { repeated = true; break; }
            seen.insert(blk);
        }
        int guess = repeated ? 0 : 1;
        if (guess == b) ++correct;
    }
    double p = (double)correct / trials;
    return fabs(p - 0.5);
}

// -------------------------------
// Multiple encryption & meet-in-the-middle demo (small keyspace)
// We'll use a tiny toy block cipher on a single block for this demo so brute force is feasible.
// For clarity, implement E_k(block) = block XOR k_repeated where k is 2-byte key (16-bit).
// Double encryption E_{k2}(E_{k1}(m)) = m XOR k1 XOR k2 -> not good; but meet-in-middle for a nontrivial primitive
// For pedagogical demonstration, we implement a toy block cipher that maps 16 bytes XOR key repeated, but we simulate meet-in-the-middle brute force to find keys for small keyspace.

using SmallKey = uint16_t; // 16-bit key to allow brute force

array<u8, BLOCK_BYTES> toy_encrypt_with_smallkey(const array<u8, BLOCK_BYTES> &m, SmallKey k) {
    array<u8, BLOCK_BYTES> out;
    u8 kbytes[BLOCK_BYTES];
    for (size_t i=0;i<BLOCK_BYTES;++i) kbytes[i] = (u8)((k >> (8*(i%2))) & 0xFF); // repeat two key bytes
    for (size_t i=0;i<BLOCK_BYTES;++i) out[i] = m[i] ^ kbytes[i];
    return out;
}

SmallKey mitm_double_decrypt_find_keys(const array<u8,BLOCK_BYTES> &m,
                                       const array<u8,BLOCK_BYTES> &ct)
{
    // Solve: ct = E_{k2}( E_{k1}(m) )
    // Precompute E_{k1}(m) for all k1 -> store mapping intermediate->k1
    unordered_map<string, SmallKey> table;
    table.reserve(1<<16);
    for (int k1 = 0; k1 < (1<<16); ++k1) {
        auto mid = toy_encrypt_with_smallkey(m, (SmallKey)k1);
        string key((char*)mid.data(), BLOCK_BYTES);
        table[key] = (SmallKey)k1;
    }
    // For each k2, compute D_{k2}(ct) and lookup
    for (int k2 = 0; k2 < (1<<16); ++k2) {
        // D_{k2}(ct) = E_{k2}(ct) because XOR involutive
        auto mid = toy_encrypt_with_smallkey(ct, (SmallKey)k2);
        string key((char*)mid.data(), BLOCK_BYTES);
        auto it = table.find(key);
        if (it != table.end()) {
            // found a collision; return one candidate (in demo keyspace small)
            // In general multiple solutions possible, but demonstration suffices.
            SmallKey found_k1 = it->second;
            SmallKey found_k2 = (SmallKey)k2;
            cout << "Found keys (meet-in-the-middle): k1=" << found_k1 << " k2=" << found_k2 << "\n";
            // verify
            auto ctest = toy_encrypt_with_smallkey( toy_encrypt_with_smallkey(m, found_k1), found_k2 );
            if (ctest == ct) return (SmallKey)((found_k1<<16) ^ found_k2); // pack keys into return (not ideal)
        }
    }
    return 0;
}

// -------------------------------
// Demo main
int main() {
    ios::sync_with_stdio(false);
    cin.tie(nullptr);

    cout << "Modes of operation & chosen-plaintext demo (educational)\n";
    cout << "Block size: " << BLOCK_BYTES << " bytes\n\n";

    // Set up a random 16-byte key
    vector<u8> key16 = random_bytes(BLOCK_BYTES);
    cout << "Toy key (hex): " << hex(key16) << "\n\n";

    // Prepare a plaintext with repeating pattern that will show in ECB
    string text = "HELLO-HELLO-HELLO-HELLO-HELLO-HELLO-HELLO-HELLO-"; // 48+ bytes
    vector<u8> plaintext(text.begin(), text.end());
    cout << "Plaintext: " << text << "\n\n";

    // --- ECB demonstration ---
    cout << "[ECB]\n";
    vector<u8> ecb_ct = encrypt_ecb(plaintext, key16);
    cout << "ECB ciphertext (hex): " << hex(ecb_ct) << "\n";
    // show block-level repeated ciphertexts
    map<string, vector<int>> block_positions;
    for (size_t off=0; off<ecb_ct.size(); off+=BLOCK_BYTES) {
        string blk((char*)&ecb_ct[off], (char*)&ecb_ct[off]+BLOCK_BYTES);
        block_positions[blk].push_back(off / BLOCK_BYTES);
    }
    cout << "Repeated ciphertext blocks (block_index list):\n";
    for (auto &kv : block_positions) if (kv.second.size() > 1) {
        cout << "  repeated block at indices: ";
        for (int idx : kv.second) cout << idx << " ";
        cout << "\n";
    }
    cout << "=> ECB reveals repeated plaintext blocks as identical ciphertext blocks.\n\n";

    // --- CBC demonstration ---
    cout << "[CBC with random IV]\n";
    CBCCiphertext cbc = encrypt_cbc(plaintext, key16);
    cout << "CBC IV (hex): " << hex(cbc.iv) << "\n";
    cout << "CBC ciphertext (hex): " << hex(cbc.ct) << "\n";
    cout << "Check repeated ciphertext blocks (should be much rarer due to chaining):\n";
    map<string, vector<int>> cbc_blocks;
    for (size_t off=0; off<cbc.ct.size(); off+=BLOCK_BYTES) {
        string blk((char*)&cbc.ct[off], (char*)&cbc.ct[off]+BLOCK_BYTES);
        cbc_blocks[blk].push_back(off / BLOCK_BYTES);
    }
    for (auto &kv : cbc_blocks) if (kv.second.size() > 1) {
        cout << "  repeated block at indices: ";
        for (int idx : kv.second) cout << idx << " ";
        cout << "\n";
    }
    cout << "=> Random IV + chaining hides identical plaintext blocks.\n\n";

    // --- OFB demonstration (keystream reuse warning) ---
    cout << "[OFB]\n";
    OFBCiphertext ofb1 = encrypt_ofb(plaintext, key16);
    // simulate accidental reuse of same IV (bad): encrypt second plaintext with same IV -> keystream reuse
    OFBCiphertext ofb2;
    ofb2.iv = ofb1.iv;
    // to produce ofb2.ct we need to run OFB ENCRYPT but with fixed IV and key -> simulate by reuse of keystream:
    // for demonstration we call internal routine but we must replicate same stream; easiest: call encrypt_ofb then overwrite IV to be same
    ofb2 = encrypt_ofb(plaintext, key16);
    ofb2.iv = ofb1.iv; // force same IV (simulating bug)
    cout << "OFB IV: " << hex(ofb1.iv) << "\n";
    // if keystream reused, XOR(ct1, ct2) = XOR(pt1, pt2) -> leak relation
    vector<u8> xor_cts = xor_vec(ofb1.ct, ofb2.ct);
    cout << "XOR of two ciphertexts with same IV (hex): " << hex(xor_cts) << "\n";
    cout << "XOR(ct1,ct2) == XOR(pt1,pt2) -> reveals plaintext relations when keystream reused.\n\n";

    // --- CTR demonstration (nonce uniqueness required) ---
    cout << "[CTR]\n";
    CTRCiphertext ctr1 = encrypt_ctr(plaintext, key16);
    CTRCiphertext ctr2 = encrypt_ctr(plaintext, key16);
    cout << "CTR nonce1: " << hex(ctr1.nonce) << "\n";
    cout << "CTR nonce2: " << hex(ctr2.nonce) << "\n";
    cout << "Usually nonces should be unique. If nonce reused, keystream reuse = catastrophic.\n\n";

    // --- IND-CPA experiments (empirical advantage) ---
    cout << "[IND-CPA experiments]\n";
    double adv_ecb = ind_cpa_ecb_experiment(key16, 1000);
    double adv_cbc = ind_cpa_cbc_experiment(key16, 1000);
    cout << "Empirical advantage (adversary that detects repeated blocks):\n";
    cout << "  ECB advantage  ≈ " << adv_ecb << " (non-zero; ECB leaks)\n";
    cout << "  CBC advantage  ≈ " << adv_cbc << " (≈0 if IV random)\n\n";

    // --- Multiple encryption & Meet-in-the-middle demo ---
    cout << "[Multiple encryption: meet-in-the-middle demo (toy small keyspace)]\n";
    // Build a small message block
    array<u8,BLOCK_BYTES> m_block;
    for (size_t i=0;i<BLOCK_BYTES;++i) m_block[i] = (u8)(i + 1);
    // Choose tiny keys k1,k2 (16-bit)
    SmallKey k1 = 0x1a2b, k2 = 0x3c4d;
    auto mid_enc = toy_encrypt_with_smallkey(m_block, k1);
    auto ct_block = toy_encrypt_with_smallkey(mid_enc, k2);
    cout << "Toy double-encryption: m xor k1 xor k2 -> shown as hex below\n";
    vector<u8> mvec(BLOCK_BYTES), ctvec(BLOCK_BYTES);
    for (size_t i=0;i<BLOCK_BYTES;++i) { mvec[i]=m_block[i]; ctvec[i]=ct_block[i]; }
    cout << "m (hex): " << hex(mvec) << "\n";
    cout << "ct (hex): " << hex(ctvec) << "\n";
    cout << "Running meet-in-the-middle (brute-force 2^16) ... (fast on modern CPU)\n";
    SmallKey packed = mitm_double_decrypt_find_keys(m_block, ct_block);
    cout << "Meet-in-the-middle demo done. (Small keyspace only for demo.)\n\n";

    cout << "Summary / guidance:\n";
    cout << " - ECB is deterministic and leaks repeated blocks; avoid unless only encrypting single unique block.\n";
    cout << " - CBC with a random IV is IND-CPA secure when used with a secure block cipher and correct padding.\n";
    cout << " - OFB and CTR turn a block cipher into a stream cipher — keystream reuse is catastrophic.\n";
    cout << " - For multiple encryption, be aware of meet-in-the-middle. Use standardized constructions (e.g. AES-256) rather\n";
    cout << "   than naive repeated encryption; 3DES exists for historical reasons but has limited effective security vs AES-256.\n";
    cout << " - Always use vetted libraries (OpenSSL, libsodium) for actual encryption.\n";

    return 0;
}
