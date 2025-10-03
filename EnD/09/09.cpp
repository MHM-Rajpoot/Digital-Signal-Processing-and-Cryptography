// hash_apps.cpp
// Educational demo of applications of hash functions:
// HMAC, Merkle Tree, Commitments, Blockchain mining (toy).
// Author: Toy cryptography code, inspired by Dr Markus Kuhn (Cambridge).
// *** DO NOT USE for real security ***

#include <bits/stdc++.h>
using namespace std;

// ---------------- Toy hash function ----------------
// Very simple (not secure!): based on std::hash
string toy_hash(const string &input) {
    std::hash<string> h;
    size_t v = h(input);
    stringstream ss;
    ss << hex << v;
    return ss.str();
}

// ---------------- HMAC ----------------
string hmac(const string &key, const string &message) {
    string k = key;
    // Pad or trim key
    if (k.size() > 64) k = toy_hash(k);
    if (k.size() < 64) k.append(64 - k.size(), '0');

    string inner = toy_hash((k + string(64, 0x36)) + message);
    string outer = toy_hash((k + string(64, 0x5c)) + inner);
    return outer;
}

// ---------------- Commitment scheme ----------------
struct Commitment {
    string commit;
    string message;
    string nonce;
};

Commitment commit_message(const string &msg) {
    string nonce = to_string(rand()%1000000);
    Commitment c;
    c.message = msg;
    c.nonce = nonce;
    c.commit = toy_hash(msg + "|" + nonce);
    return c;
}

bool open_commitment(const Commitment &c, const string &msg, const string &nonce) {
    return c.commit == toy_hash(msg + "|" + nonce);
}

// ---------------- Merkle tree ----------------
string merkle_tree_root(const vector<string> &leaves) {
    vector<string> level;
    for (auto &m : leaves) level.push_back(toy_hash(m));

    while (level.size() > 1) {
        vector<string> next;
        for (size_t i=0;i<level.size();i+=2) {
            if (i+1 < level.size())
                next.push_back(toy_hash(level[i] + level[i+1]));
            else
                next.push_back(toy_hash(level[i] + level[i])); // duplicate last
        }
        level = next;
    }
    return level[0];
}

// ---------------- Blockchain toy ----------------
struct Block {
    int index;
    string prevHash;
    string data;
    string hash;
    int nonce;
};

string compute_block_hash(const Block &b) {
    return toy_hash(to_string(b.index) + b.prevHash + b.data + to_string(b.nonce));
}

Block mine_block(int index, const string &prevHash, const string &data, int difficulty) {
    Block b;
    b.index = index;
    b.prevHash = prevHash;
    b.data = data;
    b.nonce = 0;
    string target(difficulty, '0'); // require hash prefix
    do {
    b.hash = compute_block_hash(b);
    b.nonce++;
    } while (stoull(b.hash, nullptr, 16) % (1 << difficulty*4) != 0);
    return b;
}

// ---------------- Main ----------------
int main() {
    srand(time(NULL));
    cout<<"Hash Applications Demo\n";

    // HMAC demo
    string tag = hmac("secretkey","Hello World");
    cout<<"HMAC tag = "<<tag<<"\n";

    // Commitment demo
    Commitment c = commit_message("vote=Alice");
    cout<<"Commitment = "<<c.commit<<"\n";
    cout<<"Revealing... "<<(open_commitment(c,"vote=Alice",c.nonce)?"Valid":"Invalid")<<"\n";

    // Merkle tree demo
    vector<string> leaves = {"tx1","tx2","tx3","tx4"};
    cout<<"Merkle root = "<<merkle_tree_root(leaves)<<"\n";

    // Blockchain demo
    Block genesis = mine_block(0,"0","Genesis Block",1);
    Block b1 = mine_block(1,genesis.hash,"Alice pays Bob",1);
    Block b2 = mine_block(2,b1.hash,"Bob pays Carol",1);
    cout<<"Blockchain:\n";
    cout<<" Block0 hash="<<genesis.hash<<"\n";
    cout<<" Block1 hash="<<b1.hash<<" prev="<<b1.prevHash<<"\n";
    cout<<" Block2 hash="<<b2.hash<<" prev="<<b2.prevHash<<"\n";
}
