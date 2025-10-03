# Modern Cryptography Course Implementation in C++

## Course Overview
---------------

This C++ implementation is based on the Modern Cryptography course taught by Principal Lecturers, covering essential concepts in modern cryptography. The course introduces students to both theoretical and practical aspects of cryptographic standards, algorithms, and secure system design.

[Course Webpage](https://www.cl.cam.ac.uk/teaching/2223/Crypto/)

## Course Details
----------------

* **Number of Lectures:** 14+
* **Prerequisite Courses:** Discrete Mathematics, Probability, Algorithms, and Basic Programming

## Aims
-----

This course aims to provide students with a solid understanding of modern cryptographic techniques, including private-key and public-key encryption, authentication, secure hash functions, and digital signatures. Students will learn how to implement cryptographic algorithms securely and understand potential vulnerabilities in naive schemes.

## Lectures
----------

The course covers the following topics:

1. **Cryptography Overview**: Private vs. public-key ciphers, MACs vs. signatures, certificates, adversary capabilities, Kerckhoffs’ principle.  
2. **Classic Ciphers**: Substitution and transposition ciphers, Vigenère, perfect secrecy, and one-time pads.  
3. **Private-Key Encryption**: Stream ciphers, pseudo-random generators, attacking linear-congruential RNGs and LFSRs, semantic security, oracle queries, computational security, concrete-security proofs.  
4. **Block Ciphers**: Pseudo-random functions and permutations, birthday problem, Feistel/Luby-Rackoff structure, DES, TDES, AES.  
5. **Chosen-Plaintext Attack Security**: Security with multiple encryptions, randomized encryption, modes of operation: ECB, CBC, OFB, CTR.  
6. **Message Authenticity**: Malleability, MACs, existential unforgeability, CBC-MAC, ECBC-MAC, CMAC, birthday attacks, Carter-Wegman one-time MAC.  
7. **Authenticated Encryption**: Chosen-ciphertext security, ciphertext integrity, encrypt-and-authenticate, authenticate-then-encrypt, encrypt-then-authenticate, padding oracle, GCM.  
8. **Secure Hash Functions**: One-way functions, collision resistance, padding, Merkle-Damgård construction, sponge function, duplex construct, entropy pool, SHA standards.  
9. **Applications of Secure Hash Functions**: HMAC, stream authentication, Merkle trees, commitment protocols, blockchains, Bitcoin.  
10. **Key Distribution Problem**: Needham-Schroeder protocol, Kerberos, hardware-security modules, public-key encryption schemes, CPA and CCA security.  
11. **Number Theory, Finite Groups, and Fields**: Modular arithmetic, Euclid’s algorithm, inversion, groups, rings, fields, GF(2ⁿ), subgroup order, cyclic groups, Euler’s theorem, Chinese remainder theorem, modular roots, quadratic residues, modular exponentiation.  
12. **Discrete Logarithm Problem**: Baby-step-giant-step algorithm, computational and decision Diffie-Hellman problems, DH key exchange, ElGamal encryption, hybrid cryptography, Schnorr groups, elliptic-curve systems, key sizes.  
13. **Trapdoor Permutations**: Security definitions, turning trapdoor permutations into public-key encryption, RSA, attacks on textbook RSA, optimal asymmetric encryption padding, common factor attacks.  
14. **Digital Signatures**: One-time signatures, RSA signatures, Schnorr identification scheme, ElGamal signatures, DSA, certificates, PKI, case study: PS3 hack.

## Objectives
------------

By the end of the course, students should be able to:

* Understand commonly used standardized cryptographic building blocks  
* Match application requirements with concrete security definitions and identify insecure schemes  
* Understand adversarial capabilities, attack algorithms, and implications for key sizes  
* Compare finite groups used in discrete-logarithm schemes  
* Understand number theory underlying common public-key schemes and efficient implementations  

## Recommended Reading
---------------------

* Katz, J., Lindell, Y. (2015). _Introduction to Modern Cryptography_ (2nd ed.). Chapman and Hall/CRC  

## Implementation
---------------

This C++ implementation provides a practical approach to cryptography, covering:

* Private-key and public-key encryption  
* Stream and block cipher implementations  
* MACs and digital signatures  
* Secure hash functions  
* Authenticated encryption and key distribution  

## Usage
-----

Clone the repository and compile the code using a C++ compiler. The implementation provides functions and classes to perform cryptographic operations securely.

## Contributing
------------

Contributions are welcome. Submit a pull request or issue for improvements or suggestions.

## License
-------

This implementation is licensed under the MIT License. See the LICENSE file for details.
