#include <iostream>
#include <string>
#include <vector>
#include <cmath>

// ---------------- Symmetric Cipher (XOR) ----------------
std::string xor_encrypt(const std::string &msg, const std::string &key) {
    std::string out = msg;
    for (size_t i = 0; i < msg.size(); i++) {
        out[i] = msg[i] ^ key[i % key.size()];
    }
    return out;
}

// ---------------- Simple MAC (checksum) ----------------
// Not secure! Just illustrates the idea of a tag with a shared secret
int simple_mac(const std::string &msg, const std::string &key) {
    int mac = 0;
    for (char c : msg) mac += c;
    for (char k : key) mac += k * 3;
    return mac;
}

// ---------------- Toy RSA-like signature ----------------
// (Insecure! Just demonstrates private/public key idea)

long long mod_pow(long long base, long long exp, long long mod) {
    long long result = 1;
    base = base % mod;
    while (exp > 0) {
        if (exp % 2 == 1) result = (result * base) % mod;
        exp = exp >> 1;
        base = (base * base) % mod;
    }
    return result;
}

struct RSAKeys {
    long long n, e, d;
};

// super simple keygen (not real RSA security)
RSAKeys generate_keys() {
    long long p = 61, q = 53;   // small primes
    long long n = p * q;        // modulus
    long long phi = (p - 1) * (q - 1);
    long long e = 17;           // public exponent
    long long d = 2753;         // private exponent (precomputed for demo)
    return {n, e, d};
}

// sign: sig = hash^d mod n
long long sign_message(const std::string &msg, const RSAKeys &keys) {
    long long h = 0;
    for (char c : msg) h += c;
    return mod_pow(h, keys.d, keys.n);
}

// verify: check hash == sig^e mod n
bool verify_signature(const std::string &msg, long long sig, const RSAKeys &keys) {
    long long h = 0;
    for (char c : msg) h += c;
    long long check = mod_pow(sig, keys.e, keys.n);
    return h == check;
}

// ---------------- Main demo ----------------
int main() {
    std::string message = "Attack at dawn";
    std::string key = "secret";

    // --- Symmetric encryption ---
    std::string ciphertext = xor_encrypt(message, key);
    std::string decrypted  = xor_encrypt(ciphertext, key);
    std::cout << "Original:  " << message << "\n";
    std::cout << "Encrypted: ";
    for (char c : ciphertext) std::cout << (int)(unsigned char)c << " ";
    std::cout << "\nDecrypted: " << decrypted << "\n\n";

    // --- Message Authentication Code ---
    int tag = simple_mac(message, key);
    std::cout << "Message: " << message << "\nMAC tag: " << tag << "\n";
    int verify = simple_mac(message, key);
    std::cout << "MAC verification: " << (verify == tag ? "OK" : "FAIL") << "\n\n";

    // --- Digital Signature ---
    RSAKeys keys = generate_keys();
    long long signature = sign_message(message, keys);
    std::cout << "RSA signature: " << signature << "\n";
    bool valid = verify_signature(message, signature, keys);
    std::cout << "Signature verification: " << (valid ? "OK" : "FAIL") << "\n";

    return 0;
}
