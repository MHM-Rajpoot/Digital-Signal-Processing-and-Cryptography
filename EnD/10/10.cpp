// key_dist_demo.cpp
// Educational demo: Key distribution, Needham-Schroeder (toy), Kerberos (toy), HSM, RSA-based PKE vulnerabilities.
// NOT SECURE — toy code for lecture/demonstration only.

#include <bits/stdc++.h>
using namespace std;
using u64 = unsigned long long;
using u128 = __uint128_t;

// ----------------------- small/insecure RSA helpers (toy) -----------------------
u64 gcd_u64(u64 a, u64 b) { while (b) { u64 t = a % b; a = b; b = t; } return a; }

u64 modinv(u64 a, u64 m) {
    // extended gcd
    long long t = 0, newt = 1;
    long long r = m, newr = a;
    while (newr != 0) {
        long long q = r / newr;
        long long tmp = newt; newt = t - q * newt; t = tmp;
        tmp = newr; newr = r - q * newr; r = tmp;
    }
    if (r > 1) return 0; // not invertible
    if (t < 0) t += m;
    return (u64)t;
}

u64 modexp(u64 a, u64 e, u64 m) {
    u128 res = 1;
    u128 base = a % m;
    while (e) {
        if (e & 1) res = (res * base) % m;
        base = (base * base) % m;
        e >>= 1;
    }
    return (u64)res;
}

struct RSAKeyPair {
    u64 n, e, d;
    // For demonstration, small primes used
};

RSAKeyPair rsa_generate_small() {
    // VERY small primes for demonstration only
    u64 p = 10007; // small prime
    u64 q = 10009; // small prime
    u64 n = p * q;
    u64 phi = (p-1)*(q-1);
    u64 e = 65537 % phi;
    if (gcd_u64(e, phi) != 1) e = 17;
    u64 d = modinv(e, phi);
    if (d == 0) { e = 3; d = modinv(e, phi); }
    return {n,e,d};
}

u64 rsa_encrypt(u64 m, const RSAKeyPair &pub) {
    // assume m < n
    return modexp(m, pub.e, pub.n);
}
u64 rsa_decrypt(u64 c, const RSAKeyPair &priv) {
    return modexp(c, priv.d, priv.n);
}

// ----------------------- HSM simulation -----------------------
struct HSM {
    RSAKeyPair keypair;
    HSM(const RSAKeyPair &kp) : keypair(kp) {}
    // HSM exposes decrypt/sign but never private key directly
    u64 decrypt(u64 c) const {
        // simulating secure private decryption inside HSM
        return rsa_decrypt(c, keypair);
    }
    u64 sign(u64 m) const {
        // sign = m^d mod n
        return rsa_decrypt(m, keypair);
    }
    // No method to extract d or raw key
};

// ----------------------- IND-CPA experiment: deterministic RSA fails -----------------------
void demo_rsa_cpa_failure(const RSAKeyPair &pub) {
    cout << "\n--- IND-CPA experiment with deterministic (textbook) RSA ---\n";
    // Adversary chooses m0,m1 (distinct small messages)
    u64 m0 = 42, m1 = 99;
    // Challenger picks random b and returns ciphertext = Enc(m_b)
    int b = rand() & 1;
    u64 challenge = rsa_encrypt(b ? m1 : m0, pub);
    // Adversary: query encryption oracle on m0 (since deterministic)
    u64 c0 = rsa_encrypt(m0, pub);
    int guess = (challenge == c0) ? 0 : 1;
    cout << "m0="<<m0<<" m1="<<m1<<" challenge_cipher="<<challenge<<" oracle_c0="<<c0<<"\n";
    cout << "Adversary guess="<<guess<<" actual b="<<b;
    cout << " => advantage = " << (guess==b ? "wins" : "loses") << "\n";
    cout << "Conclusion: deterministic encryption leaks equality -> not IND-CPA secure.\n";
}

// ----------------------- RSA malleability / CCA demonstration -----------------------
void demo_rsa_malleability_and_cca(const RSAKeyPair &pub, const RSAKeyPair &priv) {
    cout << "\n--- RSA multiplicative malleability & CCA-style oracle demonstration ---\n";
    // Suppose we have ciphertext c = m^e mod n. Attacker wants m but only has decryption oracle for chosen ciphertexts (except challenge).
    // With multiplicative property attacker can transform c' = c * s^e mod n. Decrypt(c') = m * s mod n.
    // If attacker can choose s and query decryption oracle for c' (and get back plaintext), recover m.
    u64 m = 1234 % pub.n;
    u64 c = rsa_encrypt(m, pub);
    cout << "Original message m="<<m<<" ciphertext c="<<c<<"\n";
    // Attacker picks s
    u64 s = 7;
    u64 s_e = modexp(s, pub.e, pub.n);
    u64 c_prime = ( (u128)c * s_e ) % pub.n;
    // Attacker sends c_prime to decryption oracle and receives m*s mod n
    u64 decrypted_prime = rsa_decrypt(c_prime, priv);
    // compute original m = decrypted_prime * inv(s) mod n
    u64 invs = modinv(s, pub.n);
    u64 recovered = ( (u128)decrypted_prime * invs ) % pub.n;
    cout << "Attacker created c'="<<c_prime<<" oracle returned "<<decrypted_prime<<"\n";
    cout << "Recovered m = "<<recovered<<" (original m="<<m<<")\n";
    cout << "Conclusion: raw RSA is malleable; with decryption oracle attacker can recover plaintext -> not IND-CCA secure.\n";
}

// ----------------------- Needham-Schroeder public-key protocol (toy) -----------------------
struct Principal {
    string name;
    RSAKeyPair kp; // public/private
    HSM *hsm;      // optional HSM to hold private key
    Principal(const string &n, const RSAKeyPair &k): name(n), kp(k), hsm(nullptr) {}
    void attach_hsm(HSM *h) { hsm = h; }
};

struct Message {
    string from;
    string to;
    u64 ciphertext; // toy: encrypt integer-coded payload
};

// Utility: encode pair (nonce, id) as integer for toy
u64 encode_pair(u64 nonce, const string &id) {
    // naive encoding: (nonce << 16) ^ (hash(id) & 0xFFFF)
    hash<string> hh;
    u64 h = hh(id) & 0xFFFFULL;
    return (nonce << 16) ^ h;
}
pair<u64,string> decode_pair(u64 v) {
    u64 nonce = v >> 16;
    u64 hid = v & 0xFFFFULL;
    // we can't invert id from hid in toy; we'll carry identity externally in simulation
    return {nonce, to_string(hid)};
}

// Simulate the Needham-Schroeder public-key protocol (classic)
// A -> B : {Na, A}_Kb
// B -> A : {Na, Nb}_Ka
// A -> B : {Nb}_Kb
// We'll show how an active intruder M can impersonate A to B by relaying messages to A and B (toy simulation).
void demo_needham_schroeder_attack(Principal &A, Principal &B) {
    cout << "\n--- Needham-Schroeder public-key protocol (toy) and simple active attack ---\n";
    // Clean run without intruder
    cout << "Normal run (no intruder):\n";
    u64 Na = 1001;
    u64 msg1_plain = encode_pair(Na, A.name);
    u64 msg1 = rsa_encrypt(msg1_plain, B.kp); // A->B
    // B decrypts:
    u64 dec1 = rsa_decrypt(msg1, B.kp);
    auto p = decode_pair(dec1); (void)p;
    u64 Nb = 2002;
    u64 msg2_plain = (Na << 32) ^ Nb; // combine in toy
    u64 msg2 = rsa_encrypt(msg2_plain, A.kp);
    u64 dec2 = rsa_decrypt(msg2, A.kp);
    u64 recomputed = (dec2 >> 32);
    cout << "A sees Na="<<Na<<" B responded with Nb, processed normally.\n";

    // Now active intruder M tries to impersonate A to B by relaying
    cout << "\nActive intruder M relaying attack (toy):\n";
    Principal M("Mallory", rsa_generate_small());
    // Step 1: A->B (A sends to B), M intercepts and stores ciphertext C1
    Na = 1111;
    u64 C1 = rsa_encrypt(encode_pair(Na, A.name), B.kp);
    cout << "A -> B (intercepted by M): ciphertext C1="<<C1<<"\n";
    // M now opens a session with A pretending to be B (or uses C1 to trick A) - classic Lowe attack
    // For simplicity we simulate the relay: M forwards C1 to B but under its own identity, then relays B->A responses back to B
    // In this toy simulation, we show that without binding identities inside messages this relaying can let M impersonate A.
    // (Rigorous Lowe attack requires more detailed identity handling.)
    cout << "Result: Without identity binding and proper freshness verification, active relay may allow impersonation.\n";
    cout << "Fix: include explicit identities inside encrypted messages and/or authenticate messages (signatures/HSM) and use nonce checks.\n";
}

// ----------------------- Kerberos-style simplified simulation -----------------------
// We'll simulate: Client C authenticates to AS (trusted) and gets TGT (encrypted to TGS) and session key Kc_tgs.
// Then C requests a service ticket for S from TGS by presenting TGT; TGS returns ticket encrypted to S and a session key Kc_s.
// Then C presents ticket to S to use service.

using Bytes = vector<unsigned char>;
Bytes xor_encrypt(const Bytes &plain, const Bytes &key) {
    Bytes out(plain.size());
    for (size_t i=0;i<plain.size();++i) out[i] = plain[i] ^ key[i % key.size()];
    return out;
}
string bytes_to_hex(const Bytes &b) {
    string s; char hexch[]="0123456789abcdef";
    for (auto x : b) { s.push_back(hexch[x>>4]); s.push_back(hexch[x&0xF]); }
    return s;
}

struct ServiceTicket {
    Bytes encrypted_ticket_for_server; // contains client id and session key, encrypted with server key
    Bytes encrypted_for_client; // bomb: client's copy containing session key, encrypted with client key
};

void demo_kerberos() {
    cout << "\n--- Simplified Kerberos-style exchange (toy) ---\n";
    // Long-term symmetric keys:
    Bytes Kc = {1,2,3,4}; // client's key shared with AS
    Bytes Ktgs = {9,9,9,9}; // TGS key
    Bytes Ks = {7,7,7,7}; // service S key

    // Client C asks AS for TGT
    string Cid = "clientA";
    // AS generates session key Kc_tgs
    Bytes Kc_tgs = {5,5,5,5};
    // AS creates TGT: encrypt {Cid, Kc_tgs} under Ktgs
    string ticket_plain = Cid + ":" + bytes_to_hex(Kc_tgs);
    Bytes ticket_plain_b(ticket_plain.begin(), ticket_plain.end());
    Bytes TGT = xor_encrypt(ticket_plain_b, Ktgs); // encrypted to TGS
    // AS also sends Kc_tgs encrypted to client under Kc
    Bytes Kc_tgs_b = Kc_tgs;
    Bytes enc_for_client = xor_encrypt(Kc_tgs_b, Kc);
    cout << "AS->C: TGT(encrypted for TGS)="<<bytes_to_hex(TGT)<<"  enc_for_client="<<bytes_to_hex(enc_for_client)<<"\n";

    // Client presents TGT to TGS and requests service ticket for S
    // TGS decrypts TGT using Ktgs, gets Cid and Kc_tgs, then issues service ticket
    // TGS creates Ks_session = {random}; encrypts {Cid, Kc_s} under Ks (server key) to make ticket for server
    Bytes Kc_s = {8,8,8,8};
    string t_for_s = Cid + ":" + bytes_to_hex(Kc_s);
    Bytes t_for_s_b(t_for_s.begin(), t_for_s.end());
    Bytes ticket_for_server = xor_encrypt(t_for_s_b, Ks); // encrypted to server
    // TGS also sends Kc_s to client encrypted under Kc_tgs
    Bytes enc_kc_s = xor_encrypt(Kc_s, Kc_tgs);
    cout << "TGS->C: ticket_for_server="<<bytes_to_hex(ticket_for_server)<<" enc_kc_s="<<bytes_to_hex(enc_kc_s)<<"\n";

    // Client presents ticket_for_server to S; S decrypts using Ks and checks client id and issues service
    Bytes recovered = xor_encrypt(ticket_for_server, Ks);
    string recstr(recovered.begin(), recovered.end());
    cout << "Server decrypted ticket: '"<<recstr<<"' -> session key for client-server established.\n";
    cout << "Conclusion: Kerberos centralizes auth: AS/TGS issue tickets encrypted under long-term keys; clients never learn server's long-term key.\n";
}

// ----------------------- Main demo -----------------------
int main() {
    ios::sync_with_stdio(false);
    cin.tie(nullptr);
    srand((unsigned)time(nullptr));

    cout << "Key distribution, Needham-Schroeder, Kerberos, HSM, PKE security demo\n";

    // Generate small RSA pairs for A and B
    RSAKeyPair A_kp = rsa_generate_small();
    RSAKeyPair B_kp = rsa_generate_small();
    cout << "\nGenerated toy RSA keypairs (small insecure numbers) for demonstration.\n";
    cout << "A.n="<<A_kp.n<<" A.e="<<A_kp.e<<"\n";
    cout << "B.n="<<B_kp.n<<" B.e="<<B_kp.e<<"\n";

    // HSM demonstration: B stores private key in HSM
    HSM b_hsm(B_kp);
    Principal A("Alice", A_kp), B("Bob", B_kp);
    B.attach_hsm(&b_hsm);

    // Demo IND-CPA vulnerability of deterministic RSA
    demo_rsa_cpa_failure(A_kp);

    // Demo malleability & CCA-style attack
    demo_rsa_malleability_and_cca(A_kp, A_kp); // using A's keypair as both pub/priv for demo (small)

    // Needham-Schroeder protocol demo (toy)
    demo_needham_schroeder_attack(A, B);

    // Kerberos simplified demo
    demo_kerberos();

    // Closing notes
    cout << "\nSummary notes:\n";
    cout << "- Key distribution: symmetric keys require secure channel; PKI or trusted directories help authenticate public keys.\n";
    cout << "- Needham-Schroeder public-key protocol needs identity binding to prevent replay/impersonation attacks (Lowe fix).\n";
    cout << "- Kerberos uses a trusted AS/TGS with tickets encrypted under long-term keys; careful nonce/timestamp handling prevents replays.\n";
    cout << "- HSMs protect private keys by exposing only operations (decrypt/sign) without key export.\n";
    cout << "- Use randomized, non-malleable encryption (RSA-OAEP, ElGamal+KDF+AE) for IND-CPA/IND-CCA security.\n";

    return 0;
}
