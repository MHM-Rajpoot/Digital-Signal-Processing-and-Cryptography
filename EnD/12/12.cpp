// number_theory_and_dlog_demo.cpp
// -----------------------------------------------------------
// Educational demo for discrete-logarithm and number-theory concepts
// Inspired by Dr Markus Kuhn (Cambridge) lectures
//
// Features:
//  - Modular arithmetic utilities: gcd, extended gcd, modular inverse, modular exponentiation
//  - Factorization and Euler's totient function
//  - Baby-Step Giant-Step algorithm for discrete logarithms
//  - Diffie-Hellman key exchange
//  - ElGamal public-key encryption and toy hybrid encryption
//  - Schnorr subgroup generator example
//  - Elliptic-Curve (affine) point arithmetic and toy ECDH
//
// Purpose:
//  - Illustrates fundamental cryptographic building blocks
//  - Provides hands-on examples of number-theory applications in crypto
//  - Intended purely for educational use; not secure for real cryptography
//
// Usage:
//  g++ -std=c++17 -O2 number_theory_and_dlog_demo.cpp -o demo
//  ./demo
//
// -----------------------------------------------------------

#include <bits/stdc++.h>
using namespace std;
using u64 = unsigned long long;
using u128 = __uint128_t;
using i128 = __int128_t;

// -------------------- Basic utilities (from number_theory_demo_fixed) --------------------
u64 gcd_u64(u64 a, u64 b) {
    while (b) { u64 t = a % b; a = b; b = t; }
    return a;
}

tuple<u64, long long, long long> ext_gcd(long long a, long long b) {
    long long old_r = a, r = b;
    long long old_s = 1, s = 0;
    long long old_t = 0, t = 1;
    while (r != 0) {
        long long q = old_r / r;
        long long tmp = old_r - q * r; old_r = r; r = tmp;
        tmp = old_s - q * s; old_s = s; s = tmp;
        tmp = old_t - q * t; old_t = t; t = tmp;
    }
    return { (u64)old_r, old_s, old_t };
}

long long modinv_ll(long long a, long long m) {
    auto [g, x, y] = ext_gcd(a, m);
    if (g != 1) return 0;
    long long inv = x % m;
    if (inv < 0) inv += m;
    return inv;
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

// Miller-Rabin deterministic bases for 64-bit
bool is_prime_det_u64(u64 n) {
    if (n < 2) return false;
    static const u64 small_primes[] = {2,3,5,7,11,13,17,19,23,29,31,37};
    for (u64 p : small_primes) { if (n == p) return true; if (n % p == 0) return false; }
    u64 d = n - 1;
    int s = 0;
    while ((d & 1) == 0) { d >>= 1; ++s; }
    u64 bases[] = {2ULL, 325ULL, 9375ULL, 28178ULL, 450775ULL, 9780504ULL, 1795265022ULL};
    for (u64 a : bases) {
        if (a % n == 0) continue;
        u64 x = modpow(a, d, n);
        if (x == 1 || x == n-1) continue;
        bool composite = true;
        for (int r = 1; r < s; ++r) {
            x = modmul(x, x, n);
            if (x == n-1) { composite = false; break; }
        }
        if (composite) return false;
    }
    return true;
}

// Pollard-Rho factorization (used later)
u64 f_rho(u64 x, u64 c, u64 mod) {
    return (modmul(x, x, mod) + c) % mod;
}
u64 pollards_rho(u64 n) {
    if (n % 2 == 0) return 2;
    std::mt19937_64 gen((unsigned)chrono::high_resolution_clock::now().time_since_epoch().count());
    std::uniform_int_distribution<u64> dist(2, n-2);
    while (true) {
        u64 x = dist(gen);
        u64 y = x;
        u64 c = dist(gen);
        if (c >= n) c %= n;
        u64 d = 1;
        while (d == 1) {
            x = f_rho(x, c, n);
            y = f_rho(f_rho(y, c, n), c, n);
            u64 diff = x > y ? x - y : y - x;
            d = gcd_u64(diff, n);
            if (d == n) break;
        }
        if (d > 1 && d < n) return d;
    }
}
void factor_recursive(u64 n, map<u64,int> &out) {
    if (n == 1) return;
    if (is_prime_det_u64(n)) { out[n]++; return; }
    for (u64 p : {2ULL,3ULL,5ULL,7ULL,11ULL,13ULL,17ULL,19ULL,23ULL,29ULL,31ULL,37ULL}) {
        if (n % p == 0) {
            int cnt = 0;
            while (n % p == 0) { n /= p; ++cnt; }
            out[p] += cnt;
        }
    }
    if (n == 1) return;
    if (is_prime_det_u64(n)) { out[n]++; return; }
    u64 d = pollards_rho(n);
    factor_recursive(d, out);
    factor_recursive(n/d, out);
}
map<u64,int> factor(u64 n) { map<u64,int> res; if (n<=1) return res; factor_recursive(n,res); return res; }

u64 euler_phi(u64 n) {
    auto f = factor(n);
    u128 res = n;
    for (auto &kv : f) {
        u64 p = kv.first;
        res = res / p * (p - 1);
    }
    return (u64)res;
}

// -------------------- Baby-Step Giant-Step (BSGS) --------------------
// Solve g^x = h (mod p) where p is prime and g is generator (or element of order q).
// Returns x or -1 if not found.
long long baby_step_giant_step(u64 g, u64 h, u64 p) {
    u64 m = (u64)ceil(sqrt((double)(p - 1)));
    unordered_map<u64,u64> table;
    table.reserve(m*1.3);
    u64 baby = 1 % p;
    for (u64 j = 0; j < m; ++j) {
        if (!table.count(baby)) table[baby] = j;
        baby = modmul(baby, g, p);
    }
    // compute g^{-m}
    u64 g_inv_m = modpow(g, p-1 - (m % (p-1)), p); // g^{-m} == g^{p-1-m}
    u64 gamma = h % p;
    for (u64 i = 0; i <= m; ++i) {
        auto it = table.find(gamma);
        if (it != table.end()) {
            u64 j = it->second;
            long long x = (long long)i * (long long)m + (long long)j;
            return x;
        }
        gamma = modmul(gamma, g_inv_m, p);
    }
    return -1;
}

// -------------------- Diffie-Hellman (classic) --------------------
struct DH_Params { u64 p; u64 g; };

// generate small safe-ish parameters for demo (NOT SECURE)
DH_Params dh_demo_params() {
    // small prime for demo
    u64 p = 1019; // prime
    u64 g = 2;    // assume generator
    return {p,g};
}

// compute shared secret
u64 dh_shared(u64 p, u64 own_priv, u64 other_pub) {
    return modpow(other_pub, own_priv, p);
}

// -------------------- ElGamal encryption (over multiplicative group mod p) --------------------
struct ElGamal_Public {
    u64 p;
    u64 g;
    u64 y; // g^x mod p
};
struct ElGamal_Private {
    u64 x;
};

pair<u64,u64> elgamal_encrypt(const ElGamal_Public &pub, u64 message, u64 k) {
    // message must be < p
    u64 c1 = modpow(pub.g, k, pub.p);
    u64 s = modpow(pub.y, k, pub.p);
    u64 c2 = modmul(message % pub.p, s, pub.p);
    return {c1, c2};
}
u64 elgamal_decrypt(const ElGamal_Private &priv, const ElGamal_Public &pub, pair<u64,u64> ct) {
    u64 c1 = ct.first, c2 = ct.second;
    u64 s = modpow(c1, priv.x, pub.p);
    // compute inverse of s
    long long invs = modinv_ll((long long)s, (long long)pub.p);
    if (invs == 0) return 0;
    return modmul(c2, (u64)invs, pub.p);
}

// -------------------- Hybrid encryption (toy) --------------------
// We'll derive a symmetric key from DH/ElGamal shared secret and use XOR keystream (toy).
vector<uint8_t> simple_xor_encrypt(const vector<uint8_t>& data, u64 key56) {
    // key56 -> simple keystream by repeated bytes of key
    vector<uint8_t> out(data.size());
    uint8_t ks[8];
    for (int i=0;i<8;i++) ks[i] = (uint8_t)((key56 >> (8*i)) & 0xFF);
    for (size_t i=0;i<data.size();++i) out[i] = data[i] ^ ks[i % 8];
    return out;
}

// -------------------- Schnorr group demo (subgroup of order q in Z_p^*) --------------------
// Given p = k*q + 1, find generator g of subgroup of order q
u64 find_schnorr_generator(u64 p, u64 q) {
    // choose h in [2,p-2], compute g = h^{(p-1)/q} mod p, check g>1
    u64 exp = (p - 1) / q;
    for (u64 h = 2; h < p-1; ++h) {
        u64 g = modpow(h, exp, p);
        if (g != 1) return g;
    }
    return 0;
}

// -------------------- Elliptic Curve (affine) toy implementation --------------------
// Curve: y^2 = x^3 + ax + b over prime field p
struct EC_Point {
    long long x, y;
    bool inf;
    EC_Point(): x(0), y(0), inf(true) {}
    EC_Point(long long _x,long long _y): x(_x), y(_y), inf(false) {}
};

struct EC_Curve {
    long long p;
    long long a, b;
    EC_Point G; // base point
    long long n; // order of G (if known)
};

// modular inverse for positive modulus p (returns 0 on failure)
long long modinv_pos(long long a, long long p) {
    long long inv = modinv_ll((long long)a, (long long)p);
    return inv;
}

EC_Point ec_add(const EC_Curve &curve, const EC_Point &P, const EC_Point &Q) {
    if (P.inf) return Q;
    if (Q.inf) return P;
    long long p = curve.p;
    if (P.x == Q.x && ( (P.y + Q.y) % p == 0 )) return EC_Point(); // point at infinity
    long long lambda;
    if (P.x == Q.x && P.y == Q.y) {
        // lambda = (3 x^2 + a) / (2 y)
        long long num = ( ( (i128)3 * P.x % p) * P.x + curve.a ) % p;
        long long den = (2 * P.y) % p;
        if (den < 0) den += p;
        long long inv = modinv_pos(den, p);
        lambda = (i128)num * inv % p;
    } else {
        long long num = (Q.y - P.y) % p;
        if (num < 0) num += p;
        long long den = (Q.x - P.x) % p;
        if (den < 0) den += p;
        long long inv = modinv_pos(den, p);
        lambda = (i128)num * inv % p;
    }
    long long xr = ( (i128)lambda * lambda - P.x - Q.x ) % p;
    if (xr < 0) xr += p;
    long long yr = ( (i128)lambda * (P.x - xr) - P.y ) % p;
    if (yr < 0) yr += p;
    return EC_Point(xr, yr);
}

EC_Point ec_mul(const EC_Curve &curve, EC_Point P, long long k) {
    EC_Point R; // infinity
    long long K = k;
    while (K > 0) {
        if (K & 1) R = ec_add(curve, R, P);
        P = ec_add(curve, P, P);
        K >>= 1;
    }
    return R;
}

// -------------------- Demo driver --------------------
int main() {
    ios::sync_with_stdio(false);
    cin.tie(nullptr);
    cout << "Discrete log + DH + ElGamal + ECC demo (educational)\n\n";

    // 1) BSGS demo
    {
        cout << "*** Baby-step Giant-step demo ***\n";
        u64 p = 1019; // small prime for demo
        u64 g = 2;
        u64 secret = 123;
        u64 h = modpow(g, secret, p);
        cout << "Solve for x in " << g << "^x = " << h << " (mod " << p << ")\n";
        long long found = baby_step_giant_step(g, h, p);
        cout << "Found x = " << found << " (expected " << secret << ")\n\n";
    }

    // 2) Diffie-Hellman key exchange demo
    {
        cout << "*** Diffie-Hellman key exchange demo ***\n";
        auto params = dh_demo_params();
        u64 p = params.p, g = params.g;
        // small random private keys
        u64 a = 321, b = 789;
        u64 A = modpow(g, a, p);
        u64 B = modpow(g, b, p);
        u64 key1 = dh_shared(p, a, B);
        u64 key2 = dh_shared(p, b, A);
        cout << "p="<<p<<" g="<<g<<"\n";
        cout << "Alice priv a="<<a<<" pub A="<<A<<"\n";
        cout << "Bob   priv b="<<b<<" pub B="<<B<<"\n";
        cout << "Shared secrets (A^b, B^a) = "<<key1<<" , "<<key2<<"\n\n";
    }

    // 3) ElGamal encryption demo (and hybrid)
    {
        cout << "*** ElGamal encryption demo ***\n";
        u64 p = 2089; // prime
        // find generator g small
        u64 g = 2;
        while (modpow(g, (p-1)/2, p) == 1) ++g; // crude ensure not order 2
        // generate keypair
        u64 x = 123; // private
        u64 y = modpow(g, x, p);
        ElGamal_Public pub{p,g,y};
        ElGamal_Private priv{x};
        // message as integer mod p (toy)
        u64 message = 42;
        u64 k = 77; // ephemeral
        auto ct = elgamal_encrypt(pub, message, k);
        cout << "Public (p,g,y)=("<<p<<","<<g<<","<<y<<")\n";
        cout << "Plain message="<<message<<" -> ciphertext (c1,c2)=("<<ct.first<<","<<ct.second<<")\n";
        u64 rec = elgamal_decrypt(priv, pub, ct);
        cout << "Decrypted message = "<<rec<<"\n";

        // Hybrid: derive symmetric key = shared secret = y^k = g^{xk}
        u64 shared = modpow(pub.y, k, p); // equals modpow(ct.first, x, p)
        cout << "Derived symmetric key (toy) = "<<shared<<"\n";
        vector<uint8_t> data = { 'H','e','l','l','o','!','\n' };
        auto cipher = simple_xor_encrypt(data, shared);
        auto plain = simple_xor_encrypt(cipher, shared);
        string outplain(plain.begin(), plain.end());
        cout << "Hybrid (toy) symmetric roundtrip: " << outplain;
        cout << "\n\n";
    }

    // 4) Schnorr-group demo (find subgroup generator g of order q)
    {
        cout << "*** Schnorr group sketch/demo ***\n";
        // For demo pick small p=q*k+1 where q is prime.
        // Example: p=  467, q= 233, k=2  (233*2+1=467, 467 prime)
        u64 q = 233;
        u64 p = 2*q + 1;
        if (!is_prime_det_u64(p) || !is_prime_det_u64(q)) {
            cout << "Example p/q not prime; skipping Schnorr demo\n\n";
        } else {
            u64 g = find_schnorr_generator(p, q);
            cout << "p="<<p<<" q="<<q<<" found subgroup generator g="<<g<<"\n";
            // generator g has order q
            cout << "Check g^q mod p = "<<modpow(g,q,p)<<"\n\n";
        }
    }

    // 5) Elliptic-curve ECDH demo (toy small curve)
    {
        cout << "*** ECC ECDH demo (toy curve) ***\n";
        // Toy curve from classic examples:
        // y^2 = x^3 + 497*x + 1768 mod 9739
        EC_Curve curve;
        curve.p = 9739;
        curve.a = 497;
        curve.b = 1768;
        // base point G (from literature)
        curve.G = EC_Point(1804, 5368);
        curve.n = 0; // unknown (we won't use order)

        // Alice chooses a, Bob chooses b
        long long a_priv = 6534;
        long long b_priv = 7626;
        EC_Point A_pub = ec_mul(curve, curve.G, a_priv);
        EC_Point B_pub = ec_mul(curve, curve.G, b_priv);
        EC_Point S1 = ec_mul(curve, B_pub, a_priv);
        EC_Point S2 = ec_mul(curve, A_pub, b_priv);
        cout << "Alice pub: ("<<A_pub.x<<","<<A_pub.y<<")\n";
        cout << "Bob   pub: ("<<B_pub.x<<","<<B_pub.y<<")\n";
        if (!S1.inf && !S2.inf && S1.x == S2.x && S1.y == S2.y) {
            cout << "ECDH shared point: ("<<S1.x<<","<<S1.y<<")\n";
        } else {
            cout << "ECDH mismatch or infinity (unexpected in this demo)\n";
        }
        cout << "\n";
    }

    cout << "Notes (Discrete Log & DH / ECC Demo):\n";
    cout << "- Baby-Step Giant-Step (BSGS) solves discrete log problems in multiplicative groups; complexity O(sqrt(p)).\n";
    cout << "- Diffie-Hellman (DH) demonstrates shared-secret derivation via modular exponentiation.\n";
    cout << "- ElGamal encryption shows public-key encryption using DH principles; toy hybrid demo illustrates symmetric key derivation.\n";
    cout << "- Schnorr-group example demonstrates subgroup generation of prime order q inside Z_p*; used in signature schemes.\n";
    cout << "- Elliptic-Curve ECDH illustrates point multiplication and shared key derivation on a small prime-field curve.\n";
    cout << "- All parameters are small and educational; not secure for real cryptography.\n";
    cout << "- Modular arithmetic, inverse, and exponentiation routines are used extensively; pay attention to overflows for larger values.\n";
    cout << "- Understanding these building blocks is essential for cryptanalysis, protocol design, and hybrid cryptosystems.\n";

    return 0;
}
