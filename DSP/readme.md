# Signal Processing Course Implementation in C++
=============================================

## Course Overview
---------------

This C++ implementation is based on the Signal Processing course taught by Principal Lecturer Dr. Markus Kuhn, taken by Part II students. The course covers fundamental signal-processing principles, digital communications examples, and practical experience with MATLAB-based programming assignments.

## Course Details
----------------

* **Number of Lectures:** 12
* **Prerequisite Courses:** Probability, Mathematical Methods for Computer Science, Unix Tools (MATLAB introduction), and Floating-Point Computation

## Aims
-----

This course aims to teach students the basic signal-processing principles necessary to understand modern high-tech systems, with a focus on digital communications examples. Students will gain practical experience through numerical experiments in MATLAB-based programming assignments.

## Lectures
----------

The course covers the following topics:

1. **Signals and Systems**: Discrete sequences and systems, their types and properties, linear time-invariant systems, convolution.
2. **Phasors**: Eigen functions of linear time-invariant systems, review of complex arithmetic, examples from electronics, optics, and acoustics.
3. **Fourier Transform**: Phasors as orthogonal base functions, forms of the Fourier transform, convolution theorem, Dirac's delta function, impulse combs in time and frequency domains.
4. **Discrete Sequences and Spectra**: Periodic sampling of continuous signals, periodic signals, aliasing, sampling and reconstruction of low-pass and band-pass signals, spectral inversion.
5. **Discrete Fourier Transform**: Continuous vs. discrete Fourier transform, symmetry, linearity, review of FFT, real-valued FFT.
6. **Spectral Estimation**: Leakage and scalloping phenomena, windowing, zero padding.
7. **Finite and Infinite Impulse-Response Filters**: Properties of filters, implementation forms, window-based FIR design, frequency-inversion, modulation, FFT-based convolution, polynomial representation, z-transform, zeros and poles.
8. **Digital Modulation**: IQ representation of band-pass signals, AM, FM, MSK, QAM, and OFDM signals, clock recovery, symbol detection, matched filter, software-defined radio.
9. **Random Sequences and Noise**: Random variables, stationary processes, autocorrelation, crosscorrelation, filtered random sequences, white noise, exponential averaging.
10. **Correlation Coding**: Random vectors, dependence vs. correlation, covariance, decorrelation, matrix diagonalization, eigen decomposition, Karhunen-Loève transform, principal component analysis.
11. **Lossy vs. Lossless Compression**: Perceptual scales, masking, spatial resolution, color coordinates, demonstration experiments.
12. **Quantization, Image Coding Standards**: A/mu-law coding, delta coding, JPEG.

## Objectives
------------

By the end of the course, students should be able to:

* Apply basic properties of time-invariant linear systems
* Understand sampling, aliasing, convolution, filtering, and spectral estimation
* Explain concepts in time and frequency domain representations
* Use filter-design software
* Visualize and discuss digital filters in the z-domain
* Use FFT for convolution, deconvolution, filtering
* Implement and evaluate simple DSP applications in MATLAB
* Apply transforms to reduce correlation between signal sources
* Understand basic principles of modulation and image coding techniques

## Recommended Reading
---------------------

* Lyons, R.G. (2010). Understanding digital signal processing. Prentice Hall (3rd ed.).
* Oppenheim, A.V. & Schafer, R.W. (2007). Discrete-time digital signal processing. Prentice Hall (3rd ed.).
* Stein, J. (2000). Digital signal processing - a computer science perspective. Wiley.
* Salomon, D. (2002). A guide to data compression methods. Springer.

## Implementation
---------------

This C++ implementation provides a practical approach to signal processing, covering topics such as:

* Signal generation and manipulation
* Filtering and convolution
* Fourier transform and spectral analysis
* Modulation and demodulation
* Image processing and compression

## Usage
-----

To use this implementation, simply clone the repository and compile the code using a C++ compiler. The implementation provides a set of functions and classes that can be used to perform various signal processing tasks.

## Contributing
------------

Contributions to this implementation are welcome. If you have any suggestions or improvements, please submit a pull request or issue on the repository.

## License
-------

This implementation is licensed under the MIT License. See the LICENSE file for more information.