# SnFFT

The original SnFFT package is Julia package designed to facilitate harmonic analysis on the symmetric group of order n, denoted Sn .
SnFFT can be downloaded [here](https://github.com/GDPlumb/SnFFT.jl) . This version is based on representations of Sn with values in a finite field.
The main differences are the following :
  * Young seminormal representations instead of Young orthogonal representations
  * Working with BigInt modulo M (M prime) instead of Float64

Out of the box, this SnFFT version implements:

  *  Group operations and factorizations for Sn
  *  Functionality to set up functions over Sn
  *  The fast Fourier transform
  *  The inverse fast Fourier transform
  *  The fast convolution product of functions over Sn
  *  The fast exponentiation for the convolution product over Sn
