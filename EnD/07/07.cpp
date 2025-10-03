// mac_demo.cpp
// Educational demo: MACs, CBC-MAC, CMAC sketch, Carter-Wegman one-time MAC, birthday attack.
// Author: Educational toy implementation, inspired by Dr Markus Kuhn’s Cambridge cryptography lectures.
// *** DO NOT USE for real cryptography ***

#include <bits/stdc++.h>
using namespace std;

// ---------------- Toy block cipher (Feistel from earlier) ----------------
uint8_t Ffunc(uint8_t half, uint8_t subkey) {
    static uint8_t Sbox[16] = {0xE,4,0xD,1,2,0xF,0xB,8,3,0xA,6,0xC,5,9,0,7};
    uint8_t x = half ^ subkey;
    uint8_t lo = Sbox[x & 0xF];
    uint8_t hi = Sbox[(x>>4) & 0xF];
    return (hi<<4)|lo;
}

uint16_t feistel_encrypt(uint16_t block, const vector<uint8_t> &roundkeys) {
    uint8_t L = block >> 8;
    uint8_t R = block & 0xFF;
    for (uint8_t rk : roundkeys) {
        uint8_t newL = R;
        uint8_t newR = L ^ Ffunc(R, rk);
        L = newL; R = newR;
    }
    return (uint16_t(L) << 8) | R;
}

// ---------------- CBC-MAC ----------------
// Toy block cipher with 16-bit blocks
uint16_t cbc_mac(const vector<uint16_t> &messageBlocks, const vector<uint8_t> &key) {
    uint16_t X = 0;
    for (uint16_t m : messageBlocks) {
        X = feistel_encrypt(X ^ m, key);
    }
    return X; // final tag
}

// ---------------- Carter-Wegman one-time MAC ----------------
// Universal hash: h(m) = (a*m + b) mod p
uint64_t carter_wegman_mac(const vector<uint64_t>& message, uint64_t a, uint64_t b, uint64_t p, uint64_t key) {
    uint64_t h = 0;
    for (uint64_t m : message) {
        h = (a * m + h + b) % p;
    }
    return (h + key) % p; // add secret pad
}

// ---------------- Birthday attack demo ----------------
void birthday_attack_demo(int nbits, int trials) {
    int N = 1<<nbits;
    unordered_map<int,int> seen;
    int collisions=0;
    mt19937 rng(42);
    uniform_int_distribution<int> dist(0,N-1);
    for (int t=0;t<trials;t++) {
        int tag=dist(rng);
        if (seen.count(tag)) { collisions++; }
        seen[tag]=t;
    }
    cout << "Birthday attack demo: "<<nbits<<"-bit tags, "<<trials<<" trials\n";
    cout << " Unique tags: "<<seen.size()<<", collisions: "<<collisions<<"\n";
    double approx = 1 - exp(-(double)trials*(trials-1)/(2*N));
    cout << " Approx collision probability ≈ "<<approx<<"\n";
}

// ---------------- Demo main ----------------
int main(){
    cout<<"MAC demo (CBC-MAC, Carter-Wegman, Birthday attack)\n";

    // Example CBC-MAC
    vector<uint8_t> key = {0x3A,0xC5,0x9F,0x12}; // toy subkeys
    vector<uint16_t> msg = {0x1111, 0x2222, 0x3333};
    uint16_t tag = cbc_mac(msg,key);
    cout<<"CBC-MAC tag for message [1111,2222,3333] = "<<hex<<tag<<dec<<"\n";

    // Example Carter-Wegman
    vector<uint64_t> msg2 = {123,456,789};
    uint64_t a=17, b=23, p=1000003; // small prime modulus
    uint64_t secretPad = 99991;
    uint64_t cwtag = carter_wegman_mac(msg2,a,b,p,secretPad);
    cout<<"Carter-Wegman one-time MAC tag = "<<cwtag<<"\n";

    // Birthday attack demonstration
    birthday_attack_demo(16, 300); // collisions in 16-bit tags
    birthday_attack_demo(32, 50000); // collisions in 32-bit tags

    // Sketch of CMAC/ECBC idea:
    cout<<"[Sketch] CMAC would derive subkeys (K1,K2) from main key and tweak the last block.\n";
    cout<<"[Sketch] ECBC uses independent keys for different message lengths.\n";
}
