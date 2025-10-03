// signatures_demo.cpp
// -----------------------------------------------------------
// Educational demo: Digital signatures
// - Lamport one-time signatures
// - RSA signatures (textbook demo)
// - Schnorr identification -> Schnorr signatures
// - ElGamal signatures
// - DSA (toy parameters)
// - Simple certificate (self-signed) and verification
//
// WARNING: Toy parameters and toy hash. NOT secure.
// Use vetted crypto libs (OpenSSL, libsodium, BoringSSL) for real systems.
// -----------------------------------------------------------

#include <bits/stdc++.h>
using namespace std;
using u64 = unsigned long long;
using u128 = __uint128_t;
using i128 = __int128_t;

// ---------- small number-theory helpers ----------
u64 gcd_u64(u64 a, u64 b) { while (b) { u64 t = a % b; a = b; b = t; } return a; }

long long egcd(long long a, long long b, long long &x, long long &y) {
    if (b == 0) { x = (a >= 0 ? 1 : -1); y = 0; return llabs(a); }
    long long x1, y1;
    long long g = egcd(b, a % b, x1, y1);
    x = y1;
    y = x1 - (a / b) * y1;
    return g;
}

long long modinv_ll(long long a, long long m) {
    long long x, y;
    long long g = egcd(a, m, x, y);
    if (g != 1) return 0;
    long long res = x % m;
    if (res < 0) res += m;
    return res;
}

u64 modmul(u64 a, u64 b, u64 mod) { return (u64)((u128)a * b % mod); }

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

// simple (toy) hash: deterministic polynomial hash
u64 hash_to_u64(const string &m) {
    u64 h = 0;
    for (char c : m) {
        h = h * 31 + (unsigned char)c;
    }
    return h;
}

// ---------- Lamport one-time signatures ----------
struct LamportKeyPair {
    int kbits;
    vector<pair<u64,u64>> pub;
    vector<pair<u64,u64>> priv;
};

LamportKeyPair lamport_keygen(int kbits=16) {
    LamportKeyPair kp;
    kp.kbits = kbits;
    kp.priv.resize(kbits);
    kp.pub.resize(kbits);
    std::mt19937_64 rng(42ULL);
    for (int i = 0; i < kbits; ++i) {
        u64 s0 = rng(), s1 = rng();
        kp.priv[i] = {s0, s1};
        kp.pub[i].first = hash_to_u64(to_string(s0));
        kp.pub[i].second = hash_to_u64(to_string(s1));
    }
    return kp;
}

vector<u64> lamport_sign(const LamportKeyPair &kp, const string &message) {
    u64 h = hash_to_u64(message);
    vector<u64> sig;
    sig.reserve(kp.kbits);
    for (int i = 0; i < kp.kbits; ++i) {
        int bit = (h >> i) & 1;
        sig.push_back(bit ? kp.priv[i].second : kp.priv[i].first);
    }
    return sig;
}

bool lamport_verify(const LamportKeyPair &kp, const string &message, const vector<u64> &sig) {
    if ((int)sig.size() != kp.kbits) return false;
    u64 h = hash_to_u64(message);
    for (int i = 0; i < kp.kbits; ++i) {
        int bit = (h >> i) & 1;
        u64 expected = bit ? kp.pub[i].second : kp.pub[i].first;
        if (hash_to_u64(to_string(sig[i])) != expected) return false;
    }
    return true;
}

// ---------- RSA signatures (textbook/demo) ----------
struct RSAKeyPair {
    u64 n;
    u64 e;
    u64 d;
    u64 p, q;
};

vector<u64> small_primes = {
    101, 103, 107, 109, 113, 127, 131, 137, 139, 149,
    151, 157, 163, 167, 173, 179, 181, 191, 193, 197
};

void rsa_keygen_toy(RSAKeyPair &kp) {
    std::mt19937_64 rng(42ULL);
    u64 p = small_primes[rng() % small_primes.size()];
    u64 q;
    do { q = small_primes[rng() % small_primes.size()]; } while (q == p);
    kp.p = p; kp.q = q;
    kp.n = p * q;
    u64 phi = (p - 1) * (q - 1);
    kp.e = 65537 % phi;
    if (gcd_u64(kp.e, phi) != 1) kp.e = 17;
    kp.d = (u64)modinv_ll((long long)kp.e, (long long)phi);
    if (kp.d == 0) {
        kp.e = 3;
        kp.d = (u64)modinv_ll((long long)kp.e, (long long)phi);
    }
}

u64 rsa_sign_textbook(u64 m, const RSAKeyPair &kp) {
    return modpow(m % kp.n, kp.d, kp.n);
}

bool rsa_verify_textbook(u64 m, u64 sig, const RSAKeyPair &kp) {
    u64 recovered = modpow(sig, kp.e, kp.n);
    return recovered == (m % kp.n);
}

u64 rsa_sign_hashed(const string &message, const RSAKeyPair &kp) {
    u64 h = hash_to_u64(message) % kp.n;
    return rsa_sign_textbook(h, kp);
}

bool rsa_verify_hashed(const string &message, u64 sig, const RSAKeyPair &kp) {
    u64 h = hash_to_u64(message) % kp.n;
    return rsa_verify_textbook(h, sig, kp);
}

// ---------- Schnorr identification/signature ----------
struct SchnorrParams {
    u64 p; // large prime
    u64 q; // q divides p-1
    u64 g; // generator of subgroup order q
};

void schnorr_set_params(SchnorrParams &P) {
    // Deterministic parameters: p=607 (prime), q=101 (divides 606=6*101), g=64 (order 101)
    P.p = 607;
    P.q = 101;
    P.g = 64;
}

struct SchnorrKey {
    u64 x; // private
    u64 y; // public = g^x mod p
};

pair<u64,u64> schnorr_sign(const SchnorrParams &P, const SchnorrKey &K, const string &m) {
    u64 k = 7; // fixed for reproducibility
    u64 R = modpow(P.g, k, P.p);
    string comb = to_string(R) + "|" + m;
    u64 e = hash_to_u64(comb) % P.q;
    long long s = ((long long)k - (long long)K.x * (long long)e) % (long long)P.q;
    if (s < 0) s += P.q;
    return {e, (u64)s};
}

bool schnorr_verify(const SchnorrParams &P, const SchnorrKey &K, const string &m, pair<u64,u64> sig) {
    u64 e = sig.first, s = sig.second;
    u64 gs = modpow(P.g, s, P.p);
    u64 ye = modpow(K.y, e, P.p);
    u64 Rcalc = modmul(gs, ye, P.p);
    string comb = to_string(Rcalc) + "|" + m;
    u64 ecalc = hash_to_u64(comb) % P.q;
    return ecalc == e;
}

// ---------- ElGamal signatures ----------
struct ElGamalKey {
    u64 p, g;
    u64 x; // private
    u64 y; // public g^x mod p
};

pair<u64,u64> elgamal_sign(const ElGamalKey &K, const string &m) {
    u64 p = K.p;
    u64 q = p - 1;
    u64 k = 1; // fixed for reproducibility
    u64 r = modpow(K.g, k, p);
    u64 hm = hash_to_u64(m) % q;
    long long kinv = modinv_ll((long long)k, (long long)q);
    long long s = ((i128)kinv * ((i128)hm - (i128)K.x * (i128)r)) % (long long)q;
    if (s < 0) s += q;
    return {r, (u64)s};
}

bool elgamal_verify(const ElGamalKey &K, const string &m, pair<u64,u64> sig) {
    u64 r = sig.first, s = sig.second;
    if (r == 0 || r >= K.p) return false;
    u64 hm = hash_to_u64(m) % (K.p - 1);
    u64 lhs = modmul(modpow(K.y, r, K.p), modpow(r, s, K.p), K.p);
    u64 rhs = modpow(K.g, hm, K.p);
    return lhs == rhs;
}

// ---------- DSA (toy) ----------
pair<u64,u64> dsa_sign(const SchnorrParams &P, u64 x, const string &m) {
    u64 p = P.p, q = P.q, g = P.g;
    u64 k = 7; // fixed for reproducibility
    u64 r = modpow(g, k, p) % q;
    if (r == 0) return {0,0};
    u64 hm = hash_to_u64(m) % q;
    long long kinv = modinv_ll((long long)k, (long long)q);
    long long s = ((i128)kinv * ((i128)hm + (i128)x * (i128)r)) % (long long)q;
    if (s < 0) s += q;
    return {r, (u64)s};
}

bool dsa_verify(const SchnorrParams &P, u64 y, const string &m, pair<u64,u64> sig) {
    u64 p = P.p, q = P.q, g = P.g;
    u64 r = sig.first, s = sig.second;
    if (r == 0 || r >= q || s == 0 || s >= q) return false;
    u64 hm = hash_to_u64(m) % q;
    long long w = modinv_ll((long long)s, (long long)q);
    if (w == 0) return false;
    u64 u1 = modmul(hm, (u64)w, q);
    u64 u2 = modmul(r, (u64)w, q);
    u64 v = (modmul(modpow(g, u1, p), modpow(y, u2, p), p)) % q;
    return v == r;
}

// ---------- Simple certificate / PKI demo ----------
struct Certificate {
    string subject;
    string pubkey_desc;
    u64 signature;
    string issuer;
};

Certificate self_signed_certificate_rsa(const string &name, const RSAKeyPair &issuer_kp) {
    Certificate cert;
    cert.subject = name;
    cert.pubkey_desc = "rsa_pub_n=" + to_string(issuer_kp.n) + "_e=" + to_string(issuer_kp.e);
    string payload = cert.subject + "|" + cert.pubkey_desc;
    cert.signature = rsa_sign_hashed(payload, issuer_kp);
    cert.issuer = "Self";
    return cert;
}

bool verify_cert_rsa(const Certificate &cert, const RSAKeyPair &issuer_kp) {
    string payload = cert.subject + "|" + cert.pubkey_desc;
    return rsa_verify_hashed(payload, cert.signature, issuer_kp);
}

// ---------- Demo driver ----------
int main() {
    ios::sync_with_stdio(false);
    cin.tie(nullptr);
    cout << "Digital Signatures demo (toy implementations)\n\n";

    // ---- Lamport one-time signature demo ----
    {
        cout << "--- Lamport one-time signature ---\n";
        auto kp = lamport_keygen(16);
        string msg = "Hello Lamport";
        auto sig = lamport_sign(kp, msg);
        cout << "Verify lamport sig: " << (lamport_verify(kp, msg, sig) ? "OK" : "FAIL") << "\n\n";
    }

    // ---- RSA signatures demo ----
    {
        cout << "--- RSA signature (textbook, toy) ---\n";
        RSAKeyPair rsa_kp;
        rsa_keygen_toy(rsa_kp);
        cout << "Generated RSA n=" << rsa_kp.n << " e=" << rsa_kp.e << " (toy)\n";
        string message = "Message to sign";
        u64 sig = rsa_sign_hashed(message, rsa_kp);
        cout << "Signature (textbook sign of H(m)) = " << sig << "\n";
        cout << "Verify RSA signature: " << (rsa_verify_hashed(message, sig, rsa_kp) ? "OK" : "FAIL") << "\n\n";
    }

    // ---- Schnorr identification -> signature demo ----
    {
        cout << "--- Schnorr signature demo ---\n";
        SchnorrParams P;
        schnorr_set_params(P);
        SchnorrKey K;
        K.x = 23; // fixed for reproducibility
        K.y = modpow(P.g, K.x, P.p);
        string msg = "Schnorr message";
        auto sig = schnorr_sign(P, K, msg);
        cout << "r = " << sig.first << ", s = " << sig.second << "\n";
        cout << "Schnorr verify: " << (schnorr_verify(P, K, msg, sig) ? "PASS" : "FAIL") << "\n\n";
    }

    // ---- ElGamal signature demo ----
    {
        cout << "--- ElGamal signature demo ---\n";
        u64 p = 563;
        u64 g = 2;
        ElGamalKey K; K.p = p; K.g = g;
        K.x = 123;
        K.y = modpow(g, K.x, p);
        string msg = "ElGamal message";
        auto sig = elgamal_sign(K, msg);
        cout << "r = " << sig.first << ", s = " << sig.second << "\n";
        cout << "ElGamal verify: " << (elgamal_verify(K, msg, sig) ? "PASS" : "FAIL") << "\n\n";
    }

    // ---- DSA demo (toy) ----
    {
        cout << "--- DSA (toy) demo ---\n";
        SchnorrParams P;
        schnorr_set_params(P);
        u64 x = 23;
        u64 y = modpow(P.g, x, P.p);
        string msg = "DSA message";
        auto sig = dsa_sign(P, x, msg);
        cout << "r = " << sig.first << ", s = " << sig.second << "\n";
        cout << "DSA verify: " << (dsa_verify(P, y, msg, sig) ? "PASS" : "FAIL") << "\n";
        cout << "Note: reusing k in DSA leaks x (private key).\n\n";
    }

    // ---- Simple certificate / PKI demo ----
    {
        cout << "--- Simple certificate (self-signed, toy) ---\n";
        RSAKeyPair ca; rsa_keygen_toy(ca);
        Certificate cert = self_signed_certificate_rsa("example.com", ca);
        cout << "Cert subject=" << cert.subject << " issuer=" << cert.issuer << "\n";
        cout << "Certificate verification: " << (verify_cert_rsa(cert, ca) ? "OK" : "FAIL") << "\n";
        cout << "In real PKI, a CA signs other people's public keys, forming trust chains.\n\n";
    }

    // ---- Closing notes ----
    cout << "Closing notes (educational):\n";
    cout << "- One-time signatures (Lamport) are information-theoretically secure but single-use.\n";
    cout << "- Textbook RSA signing (m^d mod n) is deterministic and has practical pitfalls; use RSA-PSS.\n";
    cout << "- Schnorr is compact and forms the basis for many modern signature schemes.\n";
    cout << "- ElGamal/DSA require fresh randomness 'k' per signature; nonce reuse leaks private keys.\n";
    cout << "- Certificates bind identities to public keys; PKI must handle issuance, revocation, and trust anchors.\n";
    cout << "- For production, always use vetted implementations and standards (RSA-PSS, ECDSA/EdDSA, X.509/PKIX rules).\n";
    cout << "\nDemo finished.\n";
    return 0;
}