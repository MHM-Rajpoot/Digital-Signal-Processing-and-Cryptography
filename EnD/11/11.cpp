// number_theory_demo_fixed.cpp
// Number theory demo (fixed variable-name collisions)
// Educational code: modular arithmetic, gcd/extgcd, factoring, CRT, Tonelli-Shanks, GF(2^n).
// NOT for production crypto use.

#include <bits/stdc++.h>
using namespace std;
using u64 = unsigned long long;
using u128 = __uint128_t;
using i128 = __int128_t;

// -------------------- Basic utilities --------------------
u64 gcd_u64(u64 a, u64 b) {
    while (b) { u64 t = a % b; a = b; b = t; }
    return a;
}

// extended gcd: returns (g,x,y) with ax + by = g = gcd(a,b)
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

long long modinv(long long a, long long m) {
    auto [g, x, y] = ext_gcd(a, m);
    if (g != 1) return 0; // inverse doesn't exist
    long long inv = x % m;
    if (inv < 0) inv += m;
    return inv;
}

// safe modular multiplication for 64-bit modulus using 128-bit intermediate
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

// -------------------- Miller-Rabin primality (deterministic for 64-bit) --------------------
bool is_prime_det_u64(u64 n) {
    if (n < 2) return false;
    static const u64 small_primes[] = {2,3,5,7,11,13,17,19,23,29,31,37};
    for (u64 p : small_primes) { if (n == p) return true; if (n % p == 0) return false; }
    // write n-1 = d * 2^s
    u64 d = n - 1;
    int s = 0;
    while ((d & 1) == 0) { d >>= 1; ++s; }
    // deterministic bases for 64-bit from research
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

// -------------------- Pollard Rho factorization --------------------
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

// recursively factor n into map {prime -> exponent}
void factor_recursive(u64 n, map<u64,int> &out) {
    if (n == 1) return;
    if (is_prime_det_u64(n)) {
        out[n]++;
        return;
    }
    // trial divide small primes first for speed
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

map<u64,int> factor(u64 n) {
    map<u64,int> res;
    if (n <= 1) return res;
    factor_recursive(n, res);
    return res;
}

// -------------------- Euler's totient --------------------
u64 euler_phi(u64 n) {
    auto f = factor(n);
    u128 res = n;
    for (auto &kv : f) {
        u64 p = kv.first;
        res = res / p * (p - 1);
    }
    return (u64)res;
}

// -------------------- Chinese Remainder Theorem --------------------
// CRT for pair (returns (x, lcm_mod) where x mod lcm_mod is solution), assumes moduli coprime if requested
bool crt_pair(u64 a1, u64 m1, u64 a2, u64 m2, u64 &out_a, u64 &out_m) {
    // Solve x = a1 (mod m1), x = a2 (mod m2)
    auto [g, s, t] = ext_gcd((long long)m1, (long long)m2);
    if ( (a2 - a1) % g != 0 ) return false;
    i128 l = (i128)m1 / g * m2;
    i128 tmp = (i128)(a2 - a1) / (long long)g;
    i128 mul = (i128)s * tmp % (m2 / g);
    i128 x = (i128)a1 + (i128)m1 * mul;
    i128 mod = l;
    i128 res = x % mod;
    if (res < 0) res += mod;
    out_a = (u64)res;
    out_m = (u64)mod;
    return true;
}

// general CRT for vector of congruences
// returns pair (x, M) or (0,0) if impossible
pair<u64,u64> crt(const vector<u64> &as, const vector<u64> &ms) {
    if (as.size() != ms.size() || as.empty()) return {0,0};
    u64 a = as[0], m = ms[0];
    for (size_t i=1;i<as.size();++i) {
        u64 a2 = as[i], m2 = ms[i];
        u64 na, nm;
        if (!crt_pair(a, m, a2, m2, na, nm)) return {0,0};
        a = na; m = nm;
    }
    return {a, m};
}

// -------------------- Multiplicative order and primitive root (mod prime) --------------------
u64 multiplicative_order(u64 g, u64 p) {
    // returns order of g modulo prime p (assumes 1 < g < p)
    if (gcd_u64(g, p) != 1) return 0;
    u64 phi = p - 1;
    auto f = factor(phi);
    u64 order = phi;
    for (auto &kv : f) {
        u64 q = kv.first;
        while (order % q == 0) {
            if (modpow(g, order / q, p) == 1) order /= q;
            else break;
        }
    }
    return order;
}

u64 find_primitive_root(u64 p) {
    if (p == 2) return 1;
    u64 phi = p - 1;
    auto f = factor(phi);
    vector<u64> factors;
    for (auto &kv : f) factors.push_back(kv.first);
    for (u64 g = 2; g < p; ++g) {
        bool ok = true;
        for (u64 q : factors) {
            if (modpow(g, phi / q, p) == 1) { ok = false; break; }
        }
        if (ok) return g;
    }
    return 0;
}

// -------------------- Legendre symbol and Tonelli-Shanks sqrt mod prime --------------------
int legendre_symbol(u64 a, u64 p) {
    if (a % p == 0) return 0;
    u64 ls = modpow(a, (p - 1) / 2, p);
    if (ls == 1) return 1;
    if (ls == p - 1) return -1;
    return 0; // shouldn't happen for prime p
}

// Tonelli-Shanks for square root modulo odd prime p
// returns (true, root) if root exists, else (false,0)
pair<bool,u64> tonelli_shanks(u64 n, u64 p) {
    if (n == 0) return {true, 0};
    if (p == 2) return {true, n % p};
    if (legendre_symbol(n, p) != 1) return {false, 0};
    u64 q = p - 1;
    int s = 0;
    while ((q & 1) == 0) { q >>= 1; ++s; }
    if (s == 1) {
        u64 r = modpow(n, (p + 1) / 4, p);
        return {true, r};
    }
    // find z which is a quadratic non-residue
    u64 z = 2;
    while (legendre_symbol(z, p) != -1) ++z;
    u64 c = modpow(z, q, p);
    u64 r = modpow(n, (q + 1) / 2, p);
    u64 t = modpow(n, q, p);
    int m = s;
    while (t != 1) {
        int i = 1;
        u64 tt = modmul(t, t, p);
        while (tt != 1) { tt = modmul(tt, tt, p); ++i; if (i == m) return {false,0}; }
        u64 b = modpow(c, 1ULL << (m - i - 1), p);
        r = modmul(r, b, p);
        c = modmul(b, b, p);
        t = modmul(t, c, p);
        m = i;
    }
    return {true, r};
}

// -------------------- GF(2^n) small-field arithmetic --------------------
// Represent elements as uint32_t (supports up to n<=32). irreducible polynomial given as bitmask (e.g., x^8 + x^4 + x^3 + x + 1 -> 0x11B for AES)
struct GF2N {
    int n;
    uint32_t mod_poly; // polynomial bits (degree n), e.g. for n=8 AES: 0x11B
    GF2N(int n_, uint32_t poly_) : n(n_), mod_poly(poly_) {}
    uint32_t add(uint32_t a, uint32_t b) const { return a ^ b; }
    uint32_t mul(uint32_t a, uint32_t b) const {
        uint64_t res = 0;
        uint64_t aa = a;
        for (int i = 0; i < n; ++i) {
            if ( (b >> i) & 1 ) res ^= (aa << i);
        }
        // reduce modulo mod_poly
        for (int i = 63; i >= n; --i) {
            if ( (res >> i) & 1 ) {
                res ^= ( (uint64_t)mod_poly << (i - n) );
            }
        }
        return (uint32_t)res;
    }
    uint32_t pow(uint32_t a, uint64_t e) const {
        uint32_t r = 1;
        while (e) {
            if (e & 1) r = mul(r, a);
            a = mul(a, a);
            e >>= 1;
        }
        return r;
    }
    // multiplicative inverse by exponentiation a^{2^n - 2}
    uint32_t inv(uint32_t a) const {
        // if a==0: no inverse
        if (a == 0) return 0;
        uint64_t exp = ((1ULL << n) - 1) - 1;
        return pow(a, exp);
    }
};

// -------------------- Examples in main --------------------
int main() {
    ios::sync_with_stdio(false);
    cin.tie(nullptr);
    cout << "Number theory demo (modular arithmetic, groups, CRT, Tonelli-Shanks, GF(2^n))\n\n";

    // gcd & inverse
    u64 a = 240, b = 46;
    cout << "gcd("<<a<<","<<b<<") = " << gcd_u64(a,b) << "\n";
    
    auto [eg_g, eg_x, eg_y] = ext_gcd((long long)a,(long long)b);
    cout << "extgcd: g="<<eg_g<<" x="<<eg_x<<" y="<<eg_y<<" -> " << a << "*" << eg_x << " + " << b << "*" << eg_y << " = " << eg_g << "\n";
    long long inv = modinv(17, 3120); // 17^{-1} mod 3120
    cout << "modinv(17,3120) = " << inv << "\n\n";

    // modpow & Euler phi
    cout << "modpow(3,100,1000000007) = " << modpow(3,100,1000000007ULL) << "\n";
    u64 n = 1001;
    cout << "factor("<<n<<") = ";
    auto facn = factor(n);
    for (auto &kv : facn) cout << kv.first << "^" << kv.second << " ";
    cout << "\nphi("<<n<<") = " << euler_phi(n) << "\n\n";

    // CRT example
    vector<u64> as = {2,3,2};
    vector<u64> ms = {3,5,7};
    auto [crt_x, crt_M] = crt(as, ms);
    if (crt_M) cout << "CRT solution x = " << crt_x << " (mod " << crt_M << ")\n";
    else cout << "CRT: no solution\n";
    cout << "\n";

    // primality test & factorization
    u64 testp = 104729; // prime
    cout << testp << " is prime? " << (is_prime_det_u64(testp) ? "yes" : "no") << "\n";
    u64 composite = 600851475143ULL;
    cout << "Factorizing " << composite << " ...\n";
    auto facc = factor(composite);
    for (auto &kv : facc) cout << kv.first << "^" << kv.second << " ";
    cout << "\n\n";

    // multiplicative order & primitive root (for prime)
    u64 p = 1019; // prime
    u64 prim_root = find_primitive_root(p);
    cout << "primitive root mod " << p << " is " << prim_root << "\n";
    cout << "order of g = " << multiplicative_order(prim_root, p) << " (should be p-1 = " << p-1 << ")\n\n";

    // Tonelli-Shanks sqrt mod prime
    u64 prime = 104729;
    u64 val = 12345;
    cout << "Legendre("<<val<<","<<prime<<") = " << legendre_symbol(val, prime) << "\n";
    auto [has_root, root] = tonelli_shanks(val % prime, prime);
    if (has_root) {
        cout << "sqrt("<<val<<") mod "<<prime<<" = "<<root<<" (check: "<<modmul(root,root,prime)<<")\n";
    } else cout << "No sqrt exists mod "<<prime<<"\n";
    cout << "\n";

    // GF(2^8) example (AES polynomial x^8 + x^4 + x^3 + x + 1 -> 0x11B)
    GF2N gf8(8, 0x11B);
    uint32_t a8 = 0x57, b8 = 0x83;
    uint32_t prod = gf8.mul(a8, b8);
    cout << "GF(2^8): 0x"<<hex<<a8<<" * 0x"<<b8<<" = 0x"<<prod<<dec<<"\n";
    uint32_t inv_a8 = gf8.inv(a8);
    cout << "inv(0x"<<hex<<a8<<dec<<") = 0x"<<hex<<inv_a8<<dec<<"\n\n";

    // multiplicative order example (small)
    u64 base = 2, modp = 13;
    cout << "order of "<<base<<" mod "<<modp<<" = "<<multiplicative_order(base, modp)<<"\n\n";

    // Show Euler theorem example
    u64 aa = 7, nn = 15;
    cout << aa << "^phi("<<nn<<") mod "<<nn<<" = "<< modpow(aa, euler_phi(nn), nn) << "\n\n";

    cout << "Notes:\n";
    cout << "- Pollard-Rho is probabilistic; Miller-Rabin bases chosen give deterministic correctness for 64-bit.\n";
    cout << "- Tonelli-Shanks finds modular square roots mod odd prime; for composite mod use different approach.\n";
    cout << "- GF(2^n) code is small-field polynomial arithmetic; for large n use specialized libraries.\n";
    cout << "- Hard problems used in crypto: factorization of large integers; discrete log in large prime fields or elliptic curves.\n";

    return 0;
}
