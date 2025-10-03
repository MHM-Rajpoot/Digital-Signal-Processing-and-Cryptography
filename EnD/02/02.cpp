// classic_crypto_demo_fixed.cpp
// Educational implementations of classic ciphers and simple attacks.
// NOT secure. For teaching only.

#include <bits/stdc++.h>
using namespace std;

// ---------- Utilities ----------
string normalize(const string &s) {
    string out;
    for (char c : s) if (isalpha((unsigned char)c)) out.push_back(tolower(c));
    return out;
}

void print_div() { cout << "----------------------------------------\n"; }

// ---------- Caesar / Substitution ----------
string caesar_encrypt(const string &pt, int shift) {
    string out = pt;
    for (size_t i = 0; i < pt.size(); ++i) {
        char c = tolower(pt[i]);
        if (isalpha((unsigned char)c)) {
            out[i] = char('a' + (c - 'a' + shift + 26) % 26);
        } else out[i] = pt[i];
    }
    return out;
}

// monoalphabetic substitution: provide 26-char key: mapping plaintext 'a'..'z' -> key[0]..key[25]
string subst_encrypt(const string &pt, const string &map26) {
    string out;
    for (char ch : pt) {
        if (isalpha((unsigned char)ch)) {
            char c = tolower(ch);
            out.push_back(map26[c - 'a']);
        } else out.push_back(ch);
    }
    return out;
}

string subst_decrypt(const string &ct, const string &map26) {
    // build inverse map
    string inv(26, '?');
    for (int i = 0; i < 26; ++i) inv[map26[i]-'a'] = char('a'+i);
    string out;
    for (char ch : ct) {
        if (isalpha((unsigned char)ch)) out.push_back(inv[tolower(ch)-'a']);
        else out.push_back(ch);
    }
    return out;
}

// ---------- Frequency analysis helper ----------
vector<int> freq(const string &s) {
    vector<int> f(26,0);
    for (char c : s) if (isalpha((unsigned char)c)) ++f[tolower(c)-'a'];
    return f;
}

void print_freq(const vector<int> &f) {
    int total = 0;
    for (int v : f) total += v;
    cout << "letter : count (percent)\n";
    for (int i=0;i<26;++i) {
        double p = total? (100.0 * f[i] / total) : 0.0;
        cout << char('a'+i) << " : " << f[i] << " (" << fixed << setprecision(1) << p << "%)\n";
    }
}

// Simple substitution attack: map most frequent cipher letters to most frequent English letters.
// It's a heuristic; may partially recover plaintext.
string attack_subst_freq(const string &ct) {
    string english_freq = "etaoinshrdlcumwfgypbvkjxqz"; // order of frequency
    vector<int> f = freq(ct);
    vector<pair<int,char>> v;
    for (int i=0;i<26;++i) v.push_back({f[i], char('a'+i)});
    sort(v.begin(), v.end(), greater<>());
    // build naive mapping
    string mapping(26,'?');
    for (int i=0;i<26;++i) mapping[v[i].second - 'a'] = english_freq[i];
    // decrypt using mapping (only letters mapped)
    string out;
    for (char ch : ct) {
        if (isalpha((unsigned char)ch)) {
            char lower = tolower(ch);
            char dec = mapping[lower - 'a'];
            out.push_back(dec);
        } else out.push_back(ch);
    }
    return out;
}

// ---------- Columnar Transposition ----------
string columnar_encrypt(const string &pt, int cols) {
    string n = normalize(pt);
    int r = (n.size() + cols -1)/cols;
    vector<string> grid(r, string(cols,'x')); // pad with 'x'
    for (size_t i=0;i<n.size();++i) {
        int rr = i/cols, cc = i%cols;
        grid[rr][cc] = n[i];
    }
    string out;
    for (int c=0;c<cols;++c)
        for (int rr=0; rr<r; ++rr) out.push_back(grid[rr][c]);
    return out;
}

string columnar_decrypt(const string &ct, int cols) {
    int r = (ct.size()+cols-1)/cols;
    vector<string> grid(r, string(cols,'?'));
    int idx=0;
    for (int c=0;c<cols;++c) for (int rr=0; rr<r; ++rr) {
        if (idx < (int)ct.size()) grid[rr][c] = ct[idx++];
    }
    string out;
    for (int rr=0; rr<r; ++rr) for (int c=0;c<cols;++c) out.push_back(grid[rr][c]);
    return out;
}

// Score text by english bigram/simple word count heuristic
double english_score(const string &s) {
    // count common digrams and common words
    vector<string> common = {"th","he","in","er","an","re","ed","on","es","st","en","at","to","it","is","or","as","te","be","of"};
    double score=0;
    string lower;
    for (char c : s) if (isalpha((unsigned char)c)) lower.push_back(tolower(c)); else lower.push_back(' ');
    for (auto &d : common) {
        size_t pos = lower.find(d,0);
        while (pos != string::npos) { score += 1.0; pos = lower.find(d, pos+1); }
    }
    // common words
    vector<string> words = {"the","and","that","have","for","not","with","you","this"};
    for (auto &w: words) {
        size_t pos = lower.find(w,0);
        while (pos != string::npos) { score += 2.0; pos = lower.find(w, pos+1); }
    }
    return score;
}

// Attempt to break transposition by trying column counts and scoring
void attack_transposition(const string &ct) {
    cout << "Attempting column count attack (try cols 2..20)\n";
    for (int cols=2; cols<=20; ++cols) {
        string cand = columnar_decrypt(ct, cols);
        double sc = english_score(cand);
        if (sc > 5.0) {
            cout << "cols=" << cols << " score=" << sc << " -> " << cand.substr(0,80) << "\n";
        }
    }
}

// ---------- Vigenere ----------
string vigenere_encrypt(const string &pt, const string &key) {
    string n = normalize(pt);
    string out;
    for (size_t i=0;i<n.size();++i) {
        int p = n[i]-'a';
        int k = key[i % key.size()]-'a';
        out.push_back(char('a' + (p + k) % 26));
    }
    return out;
}

string vigenere_decrypt(const string &ct, const string &key) {
    string out;
    for (size_t i=0;i<ct.size();++i) {
        int c = ct[i]-'a';
        int k = key[i % key.size()]-'a';
        out.push_back(char('a' + (c - k + 26) % 26));
    }
    return out;
}

// Index of Coincidence (IoC)
double index_of_coincidence(const string &s) {
    vector<int> f = freq(s);
    int N = 0; for (int x : f) N += x;
    if (N <= 1) return 0.0;
    double sum = 0;
    for (int x : f) sum += x * (x - 1);
    return sum / (double)(N * (N - 1));
}

// Fixed Kasiski-like repeated substring distances collector
vector<int> kasiski_distances(const string &ct, int minlen=3) {
    int n = (int)ct.size();
    vector<int> dists;
    for (int L = minlen; L <= 5; ++L) {
        unordered_map<string, vector<int>> pos;
        for (int i = 0; i + L <= n; ++i) {
            pos[ct.substr(i, L)].push_back(i);
        }
        for (auto &kv : pos) {
            const vector<int> &positions = kv.second;
            if (positions.size() >= 2) {
                for (size_t i = 1; i < positions.size(); ++i) {
                    int dist = positions[i] - positions[i-1];
                    if (dist > 0) dists.push_back(dist);
                }
            }
        }
    }
    return dists;
}

// gcd helper
int gcd_int(int a,int b){ return b==0? a: gcd_int(b,a%b); }

// Estimate key length by Kasiski distances gcd & IoC fallback
int estimate_vigenere_keylen(const string &ct) {
    auto dists = kasiski_distances(ct, 3);
    if (!dists.empty()) {
        int g = dists[0];
        for (size_t i=1;i<dists.size();++i) g = gcd_int(g, dists[i]);
        if (g > 1 && g <= 20) return g;
    }
    // fallback: try IoC per keylen guess
    double english_ioc = 0.066;
    int best_k = 1;
    double best_diff = 1e9;
    for (int k=1;k<=12;++k) {
        double avg = 0;
        for (int i=0;i<k;++i) {
            string col;
            for (int j=i;j<(int)ct.size(); j+=k) col.push_back(ct[j]);
            avg += index_of_coincidence(col);
        }
        avg /= k;
        double diff = fabs(avg - english_ioc);
        if (diff < best_diff) { best_diff = diff; best_k = k; }
    }
    return best_k;
}

// Break Vigenere given key length by frequency-shifting each column to align 'e' as highest freq heuristic
string break_vigenere_simple(const string &ct, int keylen) {
    string key;
    for (int i=0;i<keylen;++i) {
        string col;
        for (int j=i;j<(int)ct.size(); j+=keylen) col.push_back(ct[j]);
        vector<int> f = freq(col);
        int maxidx = max_element(f.begin(), f.end()) - f.begin(); // assume maps to 'e'
        int k = (maxidx - ('e'-'a') + 26) % 26;
        key.push_back(char('a' + k));
    }
    return key;
}

// ---------- One-Time Pad ----------
string otp_encrypt(const string &pt, const string &key) {
    string out;
    for (size_t i=0;i<pt.size();++i) {
        char p = tolower(pt[i]);
        if (isalpha((unsigned char)p)) {
            int c = (p - 'a') ^ (key[i] - 'a'); // XOR on 0..25 domain (toy)
            out.push_back(char('a' + (c % 26)));
        } else out.push_back(pt[i]);
    }
    return out;
}
string otp_decrypt(const string &ct, const string &key) {
    // symmetric for this toy XOR
    return otp_encrypt(ct, key);
}

// ---------- Demo main ----------
int main() {
    cout << "Classic ciphers demo & simple attacks\n";
    print_div();

    string plaintext = "Attack at dawn! The quick brown fox jumps over the lazy dog.";
    cout << "Plaintext: " << plaintext << "\n\n";

    // --- Caesar ---
    cout << "Caesar (shift=3)\n";
    string ca = caesar_encrypt(plaintext, 3);
    cout << ca << "\n";
    cout << "Decrypt: " << caesar_encrypt(ca, -3) << "\n";
    print_div();

    // --- Substitution (monoalphabetic) ---
    string map26 = "phqgiumeaylnofdxjkrcvstzwb"; // example permutation
    string sub_ct = subst_encrypt(plaintext, map26);
    cout << "Substitution ciphertext:\n" << sub_ct << "\n\n";
    cout << "Frequency analysis (ciphertext):\n";
    print_freq(freq(normalize(sub_ct)));
    cout << "\nNaive freq attack output (mapping top-freq -> etaoin...):\n";
    string naive = attack_subst_freq(normalize(sub_ct));
    cout << naive << "\n";
    cout << "Note: frequency mapping gives partial recovery but not perfect.\n";
    print_div();

    // --- Columnar transposition ---
    cout << "Columnar transposition (cols=6)\n";
    string trans_ct = columnar_encrypt(plaintext, 6);
    cout << "Ciphertext: " << trans_ct << "\n";
    cout << "Attempt attack by trying column counts (show candidate decryptions with score):\n";
    attack_transposition(trans_ct);
    print_div();

    // --- Vigenere ---
    string vkey = "lemon";
    cout << "Vigenere key: " << vkey << "\n";
    string vig = vigenere_encrypt(plaintext, vkey);
    cout << "Ciphertext (vigenere):\n" << vig << "\n";
    cout << "Index of Coincidence (full cipher): " << fixed << setprecision(4) << index_of_coincidence(vig) << "\n";
    auto dists = kasiski_distances(vig);
    cout << "Kasiski distances (some): ";
    for (size_t i=0;i<dists.size() && i<10; ++i) cout << dists[i] << " ";
    cout << "\nEstimating key length...\n";
    int klen = estimate_vigenere_keylen(vig);
    cout << "Estimated key length: " << klen << "\n";
    string guessed_key = break_vigenere_simple(vig, klen);
    cout << "Guessed key (heuristic): " << guessed_key << "\n";
    cout << "Decrypted with guessed key: " << vigenere_decrypt(vig, guessed_key) << "\n";
    cout << "Actual decrypt: " << vigenere_decrypt(vig, vkey) << "\n";
    print_div();

    // --- One-time pad (perfect secrecy if key truly random and used once) ---
    cout << "One-Time Pad (toy XOR on letters)\n";
    string otp_key;
    string norm = normalize(plaintext);
    // very simple pseudorandom key (not secure)
    mt19937 rng(12345);
    for (size_t i=0;i<norm.size();++i) otp_key.push_back(char('a' + (rng() % 26)));
    string otp_ct = otp_encrypt(plaintext, otp_key);
    cout << "OTP ciphertext: " << otp_ct << "\n";
    cout << "OTP decrypt (same key): " << otp_decrypt(otp_ct, otp_key) << "\n";
    cout << "Important: If OTP key reused for two plaintexts, XOR of ciphertexts cancels key and leaks plaintext relations.\n";
    print_div();

    cout << "End of demo.\n";
    return 0;
}
