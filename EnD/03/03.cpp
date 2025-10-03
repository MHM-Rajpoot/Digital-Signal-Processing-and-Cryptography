// stream_prg_demo.cpp
// Educational demo: stream ciphers, PRGs (LCG, LFSR), attacks (LCG param recovery, Berlekamp-Massey),
// and an IND-CPA distinguishing experiment showing adversary advantage.
// Author: pedagogical code (in the style of concise university lecturing).
// NOT FOR PRODUCTION.

#include <bits/stdc++.h>
using namespace std;
using u64 = uint64_t;
using u32 = uint32_t;

// -------------------- Utilities --------------------
string bytes_xor(const string &a, const string &b) {
    string out;
    out.resize(min(a.size(), b.size()));
    for (size_t i = 0; i < out.size(); ++i) out[i] = a[i] ^ b[i];
    return out;
}

string str_from_random_bits(const vector<int> &bits) {
    // pack bits into ascii letters (a..z) for readability or return raw bytes.
    string out;
    for (size_t i = 0; i + 7 < bits.size(); i += 8) {
        unsigned char v = 0;
        for (int j = 0; j < 8; ++j) v = (v << 1) | (bits[i + j] & 1);
        out.push_back((char)v);
    }
    return out;
}

vector<int> bits_from_string(const string &s) {
    vector<int> bits;
    for (unsigned char c : s) {
        for (int i = 7; i >= 0; --i) bits.push_back((c >> i) & 1);
    }
    return bits;
}

// Pretty print hex for debugging
string to_hex(const string &s) {
    static const char *hex = "0123456789abcdef";
    string out;
    out.reserve(s.size() * 2);
    for (unsigned char c : s) {
        out.push_back(hex[c >> 4]);
        out.push_back(hex[c & 0xF]);
    }
    return out;
}

// -------------------- LCG (toy) --------------------
// X_{n+1} = (a * X_n + c) mod m
struct LCG {
    u64 a, c, m;
    u64 state;
    LCG(u64 a_, u64 c_, u64 m_, u64 seed) : a(a_), c(c_), m(m_), state(seed) {}
    u64 next() {
        state = (a * state + c) % m;
        return state;
    }
    // produce bytes by outputting lower bytes of next()
    string gen_bytes(size_t nbytes) {
        string out;
        out.reserve(nbytes);
        while (out.size() < nbytes) {
            u64 v = next();
            // take 8 bytes of v (little endian)
            for (int i = 0; i < 8 && out.size() < nbytes; ++i) out.push_back(char((v >> (8*i)) & 0xFF));
        }
        return out;
    }
};

// try to recover a,c mod m from three consecutive outputs X0, X1, X2
// returns pair<success, pair<a,c>>
bool lcg_recover_from3(u64 m, u64 x0, u64 x1, u64 x2, u64 &out_a, u64 &out_c) {
    // Solve: x1 = a*x0 + c (mod m)
    //        x2 = a*x1 + c (mod m)
    // => (x2 - x1) = a * (x1 - x0) mod m
    auto mod = [m](int64_t z)->u64 {
        int64_t r = z % (int64_t)m;
        if (r < 0) r += (int64_t)m;
        return (u64)r;
    };
    u64 dx1 = (x1 >= x0) ? (x1 - x0) : (m - (x0 - x1));
    u64 dx2 = (x2 >= x1) ? (x2 - x1) : (m - (x1 - x2));
    // need inverse of dx1 modulo m
    // compute gcd(dx1, m)
    u64 a = dx1, b = m;
    while (b) { u64 t = a % b; a = b; b = t; }
    u64 g = a;
    if (g != 1) {
        // not invertible; naive 3-sample method fails. For demo, we return false.
        return false;
    }
    // modular inverse via extended gcd
    auto egcd = [](long long A, long long B) {
        long long a=A, b=B;
        long long x0=1, y0=0, x1=0, y1=1;
        while (b) {
            long long q = a / b;
            long long t = a - q*b; a=b; b=t;
            t = x0 - q*x1; x0=x1; x1=t;
            t = y0 - q*y1; y0=y1; y1=t;
        }
        return tuple<long long,long long,long long>(a,x0,y0);
    };
    long long A = (long long)dx1, M = (long long)m;
    auto eg = egcd(A, M);
    long long inv = get<1>(eg);
    // ensure positive inverse
    inv %= (long long)m;
    if (inv < 0) inv += (long long)m;
    u64 a_candidate = ( ( (__int128)dx2 * (u64)inv ) % m );
    // then c = x1 - a*x0 mod m
    u64 c_candidate = (x1 + m - ( (__int128)a_candidate * x0 % m) ) % m;
    out_a = a_candidate; out_c = c_candidate;
    // quick sanity: check
    if ((( (__int128)out_a * x1 + out_c) % m) != x2) return false;
    return true;
}

// -------------------- LFSR (Fibonacci) --------------------
// Represent feedback polynomial by vector<int> taps: positions (0-based) where taps include bit_i.
// Sequence: s_{t+L} = XOR_{i in taps} s_{t + (L - 1 - i)}  (but we implement standard shift)
struct LFSR {
    vector<int> taps; // e.g., for polynomial x^5 + x^2 + 1 taps = {2, 0} or similar (positions from L-1..0)
    vector<int> state; // size L, state[0] is leftmost bit
    int L;
    LFSR() {}
    LFSR(const vector<int> &taps_, const vector<int> &seed_bits) : taps(taps_), state(seed_bits) {
        L = (int)seed_bits.size();
    }
    int step_bit() {
        // output is usually the leftmost bit
        int out = state[0];
        int fb = 0;
        // compute feedback as XOR of specific state bits (based on taps indices measured from left)
        for (int t : taps) {
            // tap index measured from left: t in [0..L-1], where 0 is leftmost
            fb ^= state[t];
        }
        // shift left
        for (int i = 0; i + 1 < L; ++i) state[i] = state[i+1];
        state[L-1] = fb;
        return out;
    }
    vector<int> gen_bits(size_t n) {
        vector<int> out;
        out.reserve(n);
        for (size_t i = 0; i < n; ++i) out.push_back(step_bit());
        return out;
    }
};

// -------------------- Berlekamp-Massey (over GF(2)) --------------------
// Input: bit sequence s (0/1). Output: minimal connection polynomial as vector<int> C (bits) where
// L = degree = connection length, and recurrence s_{n} = sum_{i=1..L} C[i]*s_{n-i} (mod 2)
vector<int> berlekamp_massey(const vector<int> &s) {
    vector<int> C(1,1), B(1,1);
    int L = 0, m = 1, b = 1;
    int N = s.size();
    for (int n = 0; n < N; ++n) {
        // compute discrepancy d
        int d = 0;
        for (int i = 0; i <= L; ++i) {
            if (C.size() > i) d ^= (C[i] & s[n - i]);
            else break;
        }
        if (d == 1) {
            vector<int> T = C;
            // C = C + x^m * B
            if (C.size() < B.size() + m) C.resize(B.size() + m);
            for (size_t i = 0; i < B.size(); ++i) C[i + m] ^= B[i];
            if (2*L <= n) {
                L = n + 1 - L;
                B = T;
                m = 1;
            } else {
                m += 1;
            }
        } else {
            m += 1;
        }
    }
    // Normalize: C[0] == 1. Return C (length L+1)
    C.resize(L+1);
    return C; // C[0] = 1, C[1..L] are connection coefficients
}

// Helper: derive LFSR taps from BM polynomial bits
vector<int> taps_from_BM(const vector<int> &C) {
    // C vector has length L+1, indexed 0..L, C[0]=1.
    int L = (int)C.size() - 1;
    vector<int> taps;
    for (int i = 1; i <= L; ++i) if (C[i]) taps.push_back(L - i); // position from left: L-i
    return taps;
}

// -------------------- Stream cipher using PRG --------------------
struct StreamCipherPRG {
    // PRG is abstracted as fun<size_t->string>
    function<string(size_t)> prg_bytes;
    StreamCipherPRG(function<string(size_t)> prg) : prg_bytes(prg) {}
    string encrypt(const string &plaintext) {
        string keystream = prg_bytes(plaintext.size());
        return bytes_xor(plaintext, keystream);
    }
};

// -------------------- IND-CPA game simulator --------------------
// Challenger holds key/PRG. Adversary gets encryption oracle that uses fresh PRG output each call.
struct INDCPA_Challenger {
    function<string(size_t)> prg; // maps nbytes -> keystream bytes
    INDCPA_Challenger(function<string(size_t)> prg_) : prg(prg_) {}
    // challenge: adversary supplies m0, m1 same length; challenger picks b and returns encryption of mb
    string challenge(const string &m0, const string &m1, int &out_b) {
        if (m0.size() != m1.size()) throw runtime_error("messages must be same length");
        out_b = (rand() & 1);
        string key = prg(m0.size());
        return bytes_xor( (out_b ? m1 : m0), key );
    }
};

// A simple adversary that tries to distinguish PRG-based keystream from random by an observable bias
// For demonstration: this adversary knows that real PRG (LCG with low entropy) tends to produce many zero bytes in low bits? (toy)
// We'll implement two adversaries: one naive random-guessing, one that computes parity-based distinguisher.
//
// Adversary interface: returns guessed b and we compute advantage empirically over multiple trials.
int adversary_trivial(const string &ct, const string &m0, const string &m1) {
    // trivial guess: random
    return (rand() & 1);
}

// Adversary that tries to guess by checking parity of ciphertext bytes (toy heuristic).
int adversary_parity(const string &ct, const string &m0, const string &m1) {
    // If keystream is from LFSR/LCG with biased low bits, parity of ciphertext may differ depending on plaintext.
    // Here we just compute parity of bytes and compare to parity of m0/m1 to guess.
    int parity_ct = 0;
    for (unsigned char c : ct) parity_ct ^= (c & 1);
    int parity_m0 = 0; for (unsigned char c : m0) parity_m0 ^= (c & 1);
    int parity_m1 = 0; for (unsigned char c : m1) parity_m1 ^= (c & 1);
    // guess whichever message parity matches ciphertext parity
    if (parity_ct == parity_m0) return 0;
    if (parity_ct == parity_m1) return 1;
    // else random
    return (rand() & 1);
}

// Run many trials and compute empirical advantage
double estimate_advantage(function<string(size_t)> prg_func, int trials = 2000) {
    INDCPA_Challenger chal(prg_func);
    int correct = 0;
    for (int t = 0; t < trials; ++t) {
        // random messages of 16 bytes
        string m0(16, 0), m1(16, 0);
        for (int i = 0; i < 16; ++i) {
            m0[i] = char(rand() & 0xFF);
            m1[i] = char(rand() & 0xFF);
        }
        int b;
        string ct = chal.challenge(m0, m1, b);
        int guess = adversary_parity(ct, m0, m1);
        if (guess == b) ++correct;
    }
    double p = (double)correct / trials;
    double adv = fabs(p - 0.5);
    return adv;
}

// -------------------- Demo / experiments --------------------
int main() {
    ios::sync_with_stdio(false);
    cin.tie(nullptr);
    srand((unsigned)time(nullptr));

    cout << "Private-key stream ciphers & PRG demo\n";
    cout << "-------------------------------------\n";

    // ---- LCG Demo ----
    cout << "\n[LCG PRG demo]\n";
    // choose m = 2^32 for demo
    u64 m = (1ULL << 32);
    // pick a odd to try to make differences invertible often
    u64 a = 1664525; // popular LCG multiplier (not secure)
    u64 c = 1013904223; // popular increment
    u64 seed = 123456789;
    LCG lcg(a,c,m,seed);

    // generate a few outputs
    u64 x0 = lcg.state;
    u64 x1 = lcg.next();
    u64 x2 = lcg.next();
    cout << "LCG outputs (3): " << x0 << " " << x1 << " " << x2 << "\n";

    // Try to recover a,c from these 3 outputs (demo assumes invertible difference)
    u64 rec_a, rec_c;
    bool ok = lcg_recover_from3(m, x0, x1, x2, rec_a, rec_c);
    if (ok) {
        cout << "Recovered LCG params: a=" << rec_a << " c=" << rec_c << "\n";
        cout << "Original params       : a=" << a << " c=" << c << "\n";
    } else {
        cout << "Could not recover LCG params from 3 outputs (difference not invertible modulo m). Try different seed.\n";
    }

    // Show vulnerability: if attacker recovers a,c and a single state sample, future keystream predictable.
    // Build a stream cipher using LCG as PRG
    LCG lcg_for_prg(a,c,m,seed);
    auto lcg_prg = [&](size_t nbytes){ return lcg_for_prg.gen_bytes(nbytes); };
    StreamCipherPRG sc_lcg(lcg_prg);
    string plaintext = "This is a short demo plaintext!";
    string ciphertext = sc_lcg.encrypt(plaintext);
    cout << "Plaintext:  " << plaintext << "\n";
    cout << "Ciphertext: " << to_hex(ciphertext) << "\n";
    // attacker who recovers (a,c, state) can reconstruct keystream and decrypt.

    // ---- LFSR Demo and attack via Berlekamp-Massey ----
    cout << "\n[LFSR PRG and Berlekamp-Massey attack]\n";
    // Example LFSR: degree L = 8, taps at positions {0,2,3,7} (just an example)
    int L = 8;
    vector<int> taps = {0,2,3,7}; // positions from left: 0..L-1
    // seed (non-zero)
    vector<int> seed_bits(L);
    for (int i = 0; i < L; ++i) seed_bits[i] = (rand() & 1);
    // ensure non-zero
    bool allzero = true; for (int b : seed_bits) if (b) allzero=false;
    if (allzero) seed_bits[0]=1;
    LFSR lfsr(taps, seed_bits);
    // generate keystream bits
    int needed_bits = 300;
    vector<int> ks_bits = lfsr.gen_bits(needed_bits);
    cout << "Generated " << ks_bits.size() << " bits from LFSR (degree " << L << ").\n";

    // Attacker receives keystream bits and runs Berlekamp-Massey to recover minimal polynomial
    vector<int> C = berlekamp_massey(ks_bits);
    int recovered_L = (int)C.size() - 1;
    cout << "Berlekamp-Massey recovered linear complexity L = " << recovered_L << "\n";
    vector<int> recovered_taps = taps_from_BM(C);
    cout << "Recovered taps (positions from left): ";
    for (int t : recovered_taps) cout << t << " ";
    cout << "\nOriginal taps: ";
    for (int t : taps) cout << t << " ";
    cout << "\n";
    // Once attacker has taps and L, they can solve linear system to find initial state from first L bits (or reconstruct future bits).

    // Demonstrate that recovered L equals original (should be <= L, often equals)
    cout << "If LFSR is linear, BM recovers degree ~ L given enough bits (~2L bits recommended).\n";

    // Build an LFSR-based byte-wise PRG (pack bits into bytes)
    LFSR lfsr_for_prg(taps, seed_bits);
    auto lfsr_prg = [&](size_t nbytes)->string {
        vector<int> bits = lfsr_for_prg.gen_bits(nbytes * 8);
        return str_from_random_bits(bits);
    };
    StreamCipherPRG sc_lfsr(lfsr_prg);
    string ct_lfsr = sc_lfsr.encrypt(plaintext);
    cout << "LFSR-based cipher ciphertext (hex): " << to_hex(ct_lfsr) << "\n";

    // ---- IND-CPA experiments: estimate advantage of parity adversary ----
    cout << "\n[IND-CPA experiments: estimating advantage of simple parity adversary]\n";
    // PRG1: LCG-based (as above)
    LCG lcg_prg_instance(a,c,m, (u64) (rand() & 0xFFFFFFFF));
    auto prg_lcg_func = [&](size_t nbytes){ return lcg_prg_instance.gen_bytes(nbytes); };
    // PRG2: true random (use rand())
    auto prg_rand_func = [&](size_t nbytes)->string {
        string out; out.resize(nbytes);
        for (size_t i = 0; i < nbytes; ++i) out[i] = char(rand() & 0xFF);
        return out;
    };

    int trials = 2000;
    double adv_lcg = estimate_advantage(prg_lcg_func, trials);
    double adv_rand = estimate_advantage(prg_rand_func, trials);
    cout << "Estimated advantage (parity adversary) vs LCG-PRG:  " << adv_lcg << " (empirical)\n";
    cout << "Estimated advantage (parity adversary) vs true-random: " << adv_rand << " (should be ~0)\n";
    cout << "Interpretation: larger advantage implies the adversary can distinguish PRG-based keystream from random.\n";

    // ---- Short discussion in comments printed to user ----
    cout << "\n[Short conceptual notes]\n";
    cout << " - Semantic security (IND-CPA): for any efficient adversary A making queries to encryption oracle, the\n";
    cout << "   advantage Adv_A = |Pr[A wins] - 1/2| should be negligible. In practice we estimate advantage empirically.\n";
    cout << " - Reduction idea (concrete security): If an adversary A distinguishes stream cipher built from PRG G with\n";
    cout << "   advantage eps, then we can build a distinguisher D that uses A to distinguish G's output from random with\n";
    cout << "   advantage roughly eps (plus small loss). Thus proving PRG indistinguishability implies IND-CPA for stream cipher.\n";
    cout << " - Computational vs information-theoretic: OTP with truly random key has perfect secrecy (information-theoretic).\n";
    cout << "   PRG-based stream ciphers only achieve computational security: secure against feasible adversaries, assuming\n";
    cout << "   the PRG is secure.\n";
    cout << " - Concrete-security proofs quantify advantage (e.g., eps <= q * Adv_G where q is number of queries/blocks).\n";
    cout << " - Attacks: LCGs are linear and often easy to recover (e.g., linear equations or lattice methods). LFSRs are linear\n";
    cout << "   over GF(2) and vulnerable to Berlekamp-Massey which recovers the minimal linear recurrence from 2L bits.\n";

    cout << "\nEnd of demo.\n";
    return 0;
}
