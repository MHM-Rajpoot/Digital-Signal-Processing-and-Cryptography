// rsa_tdp_demo.cpp
// -----------------------------------------------------------
// Educational demo: Trapdoor permutations, RSA as TDP, attacks on textbook RSA,
// common-factor attack, multiplicative malleability, broadcast (Hastad) low-e attack,
// and a toy randomized padding illustration.
//
// WARNING: This is pedagogical toy code ONLY. Uses tiny primes and simplistic padding.
// DO NOT use for real cryptography.
// -----------------------------------------------------------

#include <bits/stdc++.h>
using namespace std;
using u64 = unsigned long long;
using u128 = __uint128_t;

// -------------------- small number-theory helpers --------------------
u64 gcd_u64(u64 a, u64 b) {
    while (b) { u64 t = a % b; a = b; b = t; }
    return a;
}

long long egcd(long long a, long long b, long long &x, long long &y) {
    if (b==0) { x= (a>=0?1:-1); y=0; return llabs(a); }
    long long x1,y1;
    long long g = egcd(b, a%b, x1, y1);
    x = y1;
    y = x1 - (a/b) * y1;
    return g;
}

long long modinv_ll(long long a, long long m) {
    long long x,y;
    long long g = egcd(a,m,x,y);
    if (g != 1) return 0;
    long long res = x % m;
    if (res < 0) res += m;
    return res;
}

u64 modmul(u64 a, u64 b, u64 mod) {
    return (u64)((u128)a * b % mod);
}

u64 modpow(u64 a, u64 e, u64 mod) {
    u64 res = 1 % mod;
    u64 base = a % mod;
    while (e) {
        if (e & 1) res = modmul(res, base, mod);
        base = modmul(base, base, mod);
        e >>= 1;
    }
    return res;
}

// integer k-th root (floor) for unsigned __int128 - used for Hastad small-e attack
u64 integer_kth_root_u128(u128 x, unsigned k) {
    // binary search
    u128 lo = 0, hi = 1;
    while (pow((long double)hi, (long double)k) <= (long double)x) hi <<= 1;
    while (lo + 1 < hi) {
        u128 mid = (lo + hi) >> 1;
        // compute mid^k safely
        u128 acc = 1;
        for (unsigned i=0;i<k;i++) {
            acc *= mid;
            if (acc > x) break;
        }
        if (acc <= x) lo = mid;
        else hi = mid;
    }
    return (u64)lo;
}

// safer integer root using multiplication with overflow checks
u64 kth_root_small(u64 n, unsigned k) {
    if (n == 0) return 0;
    // for simplicity using 128-bit
    u128 N = n;
    u128 lo=0, hi = (u128)1 << 63;
    while (lo + 1 < hi) {
        u128 mid = (lo + hi) >> 1;
        // compute mid^k
        u128 acc=1;
        for (unsigned i=0;i<k;i++) {
            acc *= mid;
            if (acc > N) break;
        }
        if (acc <= N) lo = mid; else hi = mid;
    }
    return (u64)lo;
}

// compute integer cube root of u128 exactly if perfect cube, else floor
u64 integer_cuberoot_u128(u128 x) {
    u128 lo = 0, hi = ((u128)1)<<86; // plenty big
    while (lo + 1 < hi) {
        u128 mid = (lo + hi) >> 1;
        u128 acc = mid*mid*mid;
        if (acc == x) return (u64)mid;
        if (acc < x) lo = mid; else hi = mid;
    }
    return (u64)lo;
}

// -------------------- small RSA (toy) --------------------
// We'll use tiny primes from a small hard-coded list for deterministic runs.
// Note: primes are tiny to keep arithmetic in 64-bit for demonstration.

vector<u64> small_primes = {
    101, 103, 107, 109, 113, 127, 131, 137, 139, 149,
    151, 157, 163, 167, 173, 179, 181, 191, 193, 197,
    199, 211, 223, 227, 229, 233, 239, 241, 251, 257
};

// select two distinct primes from list (toy key generation)
void rsa_keygen_toy(u64 &n, u64 &e, u64 &d) {
    u64 p = small_primes[rand() % small_primes.size()];
    u64 q;
    do { q = small_primes[rand() % small_primes.size()]; } while (q == p);
    n = p * q;
    u64 phi = (p-1) * (q-1);
    e = 65537 % phi;
    if (gcd_u64(e, phi) != 1) e = 17;
    d = modinv_ll(e, phi);
    if (d == 0) {
        // fallback
        e = 3;
        d = modinv_ll(e, phi);
    }
    // store p and q not returned (in real keypair you'd keep them secret)
}

// RSA encrypt/decrypt (textbook)
u64 rsa_encrypt(u64 m, u64 e, u64 n) { return modpow(m, e, n); }
u64 rsa_decrypt(u64 c, u64 d, u64 n) { return modpow(c, d, n); }

// RSA private key from p,q and e -> compute d. (used if need explicit factors)
u64 rsa_compute_d_from_pq(u64 p, u64 q, u64 e) {
    u64 phi = (p-1)*(q-1);
    return modinv_ll(e, phi);
}

// -------------------- Trapdoor permutation & turning into PKE --------------------
/*
 A trapdoor permutation f_k(x) is easy to compute (public key) and hard to invert
 without trapdoor t (private key). RSA: f_{(n,e)}(x) = x^e mod n is a TDP if factoring n is hard.

 To turn a TDP into an encryption scheme:
  - deterministic: ciphertext = f_k(m)  (this is *not* IND-CPA secure).
  - randomized: apply randomized encoding/padding R to m to get M = R(m;r), then c = f_k(M).
  - for IND-CPA we need probabilistic encryption: e.g. RSA-OAEP (padding+hashing) or hybrid with ephemeral Diffie-Hellman.
*/

// toy randomized padding (NOT OAEP) - just append small random bytes to message integer
u64 toy_random_pad(u64 m) {
    // pack m (assume small) into 32-bit area and add 16-bit random pad
    u64 r = (u64)(rand() & 0xFFFF);
    return (m << 16) ^ r;
}
u64 toy_random_unpad(u64 M) {
    return M >> 16;
}

// -------------------- Attacks on textbook RSA --------------------

// 1) Deterministic ciphertext detection (IND-CPA failure demo):
//   if attacker can query encryption oracle, equality of ciphertexts reveals equality of messages.

// 2) Multiplicative malleability & CCA-style attack:
//   Given c = m^e mod n, attacker picks s, computes c' = c * s^e mod n.
//   If decryption oracle returns m' = m * s mod n, and attacker can invert s, m recovered.
//   We'll simulate an oracle with known d to show this.

void demo_malleability(u64 n, u64 e, u64 d, u64 message) {
    cout << ">> Textbook RSA malleability demo\n";
    cout << " original message m = " << message << "\n";
    u64 c = rsa_encrypt(message, e, n);
    cout << " ciphertext c = m^e mod n = " << c << "\n";

    // attacker chooses s
    u64 s = 5; // chosen small multiplier
    u64 s_e = rsa_encrypt(s, e, n);
    u64 c_prime = modmul(c, s_e, n);
    // decryption oracle (simulated: we have d)
    u64 mprime = rsa_decrypt(c_prime, d, n);
    cout << " attacker picked s = " << s << "\n";
    cout << " c' = c * s^e mod n = " << c_prime << "\n";
    cout << " oracle decrypts c' -> m' = " << mprime << "\n";
    // attacker computes m = m' * s^{-1} mod n
    long long invs = modinv_ll((long long)s, (long long)n);
    if (invs == 0) {
        cout << " s not invertible mod n (unexpected in demo)\n";
    } else {
        u64 recovered = modmul(mprime, (u64)invs, n);
        cout << " attacker computes recovered m = m' * s^{-1} mod n = " << recovered << "\n";
    }
    cout << "\n";
}

// 3) Common factor attack:
//    If two RSA moduli share a prime factor p (i.e., n1 and n2 not coprime), gcd(n1,n2) = p reveals p and hence factorization.
void demo_common_factor_attack(u64 n1, u64 n2) {
    cout << ">> Common-factor attack demo\n";
    cout << " n1=" << n1 << "\n";
    cout << " n2=" << n2 << "\n";
    u64 g = gcd_u64(n1, n2);
    if (g == 1) {
        cout << " gcd(n1,n2) = 1 -> no shared factor detected\n\n";
        return;
    }
    cout << " gcd(n1,n2) = " << g << " (non-trivial) -> shared prime factor found\n";
    u64 p = g;
    u64 q1 = n1 / p;
    u64 q2 = n2 / p;
    cout << " n1 = p * q1 -> p="<<p<<" q1="<<q1<<"\n";
    cout << " n2 = p * q2 -> p="<<p<<" q2="<<q2<<"\n";
    cout << " -> attacker can compute phi and recover private exponents for both keys\n\n";
}

// 4) Broadcast / Hastad low-exponent attack (e small, no padding):
//    If same message m is sent to k recipients with same small e and pairwise coprime moduli n_i,
//    and k >= e, attacker can CRT the ciphertexts to obtain m^e mod N (N = product n_i).
//    If m^e < N, taking exact e-th root yields m.

u64 crt_combine_pair(u64 a1, u64 n1, u64 a2, u64 n2) {
    // combine two congruences x = a1 (mod n1), x = a2 (mod n2); assume coprime for this demo
    long long s,t;
    long long g = egcd((long long)n1, (long long)n2, s, t);
    if (g != 1) {
        // not coprime; simple fallback (not expected in correct demo)
        return 0;
    }
    long long m1_inv = modinv_ll((long long)n1, (long long)n2);
    if (m1_inv == 0) return 0;
    long long k = ((long long)(a2 % n2) - (long long)(a1 % n2)) % (long long)n2;
    if (k < 0) k += n2;
    k = (k * m1_inv) % (long long)n2;
    // x = a1 + n1 * k
    __int128 x = ( __int128 ) a1 + ( __int128 ) n1 * k;
    return (u64)x;
}

u128 crt_combine_multiple(const vector<u64>& as, const vector<u64>& ns) {
    // combine all congruences; ns must be pairwise coprime ideally.
    u128 x = as[0];
    u128 N = ns[0];
    for (size_t i=1;i<as.size();++i) {
        // solve: x = current_x (mod N), x = as[i] (mod ns[i])
        // use extended gcd for N and ns[i]
        long long s,t;
        long long g = egcd((long long)N, (long long)ns[i], s, t);
        if ( ((long long)as[i] - (long long)(x % ns[i])) % g != 0 ) {
            // no solution (unexpected in demo)
            return 0;
        }
        // compute modular combine (safe since we use small numbers)
        // reduce to: x + N * ((as[i] - x) / g * s mod ns[i]/g)
        long long mod_i = (long long)(ns[i] / g);
        long long tmp = ( ( (long long)as[i] - (long long)(x % ns[i]) ) / (long long)g ) % mod_i;
        if (tmp < 0) tmp += mod_i;
        long long mul = ( (long long)s % mod_i + mod_i ) % mod_i;
        long long k = ( (__int128) tmp * mul ) % mod_i;
        x = x + (u128)N * (u128)k;
        N = N * (u128)(ns[i] / g);
        x %= N;
    }
    return x;
}

void demo_broadcast_hastad() {
    cout << ">> Hastad broadcast (low exponent) attack demo\n";
    // choose e = 3 and three small RSA moduli (pairwise coprime)
    u64 e = 3;
    // for demonstration choose small primes and compute n_i
    u64 p1 = 101, q1 = 103;
    u64 p2 = 107, q2 = 109;
    u64 p3 = 113, q3 = 127;
    u64 n1 = p1*q1;
    u64 n2 = p2*q2;
    u64 n3 = p3*q3;
    // Same plaintext message m (small so that m^3 < n1*n2*n3)
    u64 m = 42;
    u64 c1 = modpow(m, e, n1);
    u64 c2 = modpow(m, e, n2);
    u64 c3 = modpow(m, e, n3);
    cout << "moduli n1="<<n1<<" n2="<<n2<<" n3="<<n3<<"\n";
    cout << "ciphertexts c1="<<c1<<" c2="<<c2<<" c3="<<c3<<"\n";
    // combine with CRT to get M = m^e mod N, but since m^e < N, we get exact m^e
    vector<u64> as = {c1, c2, c3};
    vector<u64> ns = {n1, n2, n3};
    u128 M = crt_combine_multiple(as, ns); // M == m^e
    // compute integer e-th root
    u64 recovered = integer_cuberoot_u128(M);
    cout << "After CRT combined M = m^e (as integer) = " << (unsigned long long)M << "\n";
    cout << "Integer cube root recovered m = " << recovered << " (original " << m << ")\n\n";
}

// -------------------- Demonstrations & main driver --------------------

int main(){
    ios::sync_with_stdio(false);
    cin.tie(nullptr);
    srand((unsigned)time(nullptr) ^ (unsigned)getpid());

    cout << "Trapdoor permutations & RSA demo (educational toy)\n\n";

    // --- RSA key generation (toy) ---
    u64 n, e, d;
    rsa_keygen_toy(n,e,d);
    cout << "Toy RSA keypair generated (small primes):\n";
    cout << " n = " << n << "\n";
    cout << " e = " << e << "\n";
    cout << " d = " << d << "  (private)\n\n";

    // Show RSA as trapdoor permutation
    cout << "RSA as trapdoor permutation f_{(n,e)}(x) = x^e mod n\n";
    u64 m = 123; // small message
    u64 c = rsa_encrypt(m, e, n);
    cout << " f(m) = " << c << "\n";
    u64 m_rec = rsa_decrypt(c, d, n);
    cout << " f^{-1}(c) with trapdoor (d) recovers m = " << m_rec << "\n\n";

    // Demonstrate deterministic (textbook) CPA failure: same m -> same c
    cout << "IND-CPA failure demo (deterministic textbook RSA):\n";
    u64 m_same = 77;
    u64 c1 = rsa_encrypt(m_same, e, n);
    u64 c2 = rsa_encrypt(m_same, e, n);
    cout << " encrypt(m) -> c1="<<c1<<" c2="<<c2<<" (equal => attacker can detect equality)\n\n";

    // Show malleability attack
    demo_malleability(n,e,d,m);

    // Demonstrate common-factor attack (simulate two keys sharing factor)
    // Construct n1,n2 that share a prime p (we simulate by reusing a prime from small_primes)
    u64 p_shared = 127; // chosen small prime
    u64 qA = 131;
    u64 qB = 137;
    u64 n1 = p_shared * qA;
    u64 n2 = p_shared * qB;
    cout << "Now simulate two public moduli with a common prime factor (bad key generation):\n";
    demo_common_factor_attack(n1, n2);

    // Demonstrate Hastad broadcast attack (e=3)
    demo_broadcast_hastad();

    // Show toy randomized padding to mitigate deterministic problems (very simple)
    cout << "Toy randomized padding demo (NOT OAEP) — illustrates randomness prevents identical ciphertexts\n";
    u64 m_plain = 42;
    u64 M1 = toy_random_pad(m_plain);
    u64 M2 = toy_random_pad(m_plain);
    cout << " padded M1="<<M1<<" M2="<<M2<<"\n";
    u64 cM1 = rsa_encrypt(M1 % n, e, n);
    u64 cM2 = rsa_encrypt(M2 % n, e, n);
    cout << " ciphertexts cM1="<<cM1<<" cM2="<<cM2<<"\n";
    cout << " => padding caused ciphertexts to differ for same message\n\n";

    // Short discussion text (educational)
    cout << "Summary / guidance (short):\n";
    cout << "- A trapdoor permutation is easy to compute but hard to invert without a trapdoor. RSA (x^e mod n) is the classic example.\n";
    cout << "- Textbook (deterministic) RSA is NOT IND-CPA secure. Use randomized padding (OAEP) or hybrid KEM+DEM.\n";
    cout << "- RSA is multiplicatively malleable: raw RSA ciphertexts can be modified to produce predictable changes in plaintext.\n";
    cout << "- Common-factor (gcd) attacks occur if moduli share primes (bad key generation). Always ensure independent primes.\n";
    cout << "- Low-exponent broadcast attacks (Hastad) recover messages when same plaintext is sent under small e to multiple recipients without padding.\n";
    cout << "- In practice: use RSA-OAEP (encryption), PKCS#1 v1.5 is legacy, and use large key sizes and vetted libraries (OpenSSL, libsodium).\n";
    cout << "\nEnd of demo. Remember: this code is pedagogical only — do not use for real security.\n";

    return 0;
}
