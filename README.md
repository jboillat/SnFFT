# &#x1D54A;<sub>n</sub>FFT&#x1D53D;<sub>p</sub>

The original &#x1D54A;<sub>n</sub>FFT package is a Julia package designed to facilitate harmonic analysis on the symmetric group of order n, denoted &#x1D54A;<sub>n</sub>. &#x1D54A;<sub>n</sub>FFT can be downloaded [here](https://github.com/GDPlumb/SnFFT.jl). The version &#x1D54A;<sub>n</sub>FFT&#x1D53D;<sub>p</sub> is based on representations of &#x1D54A;<sub>n</sub> with values in a finite field &#x1D53D;<sub>p</sub>.

The two main differences between the packages are the following:

  * Young seminormal representations instead of Young orthogonal representations
  * Working with `BigInt` modulo p (p prime) instead of `Float64`

Out of the box, this &#x1D54A;<sub>n</sub>FFT&#x1D53D;<sub>p</sub> version implements:

  *  Group operations and factorizations for &#x1D54A;<sub>n</sub>
  *  Functionality to set up &#x1D53D;<sub>p</sub> valued functions over &#x1D54A;<sub>n</sub>
  *  The fast Fourier transform
  *  The inverse fast Fourier transform
  *  The fast convolution product of functions over &#x1D54A;<sub>n</sub>
  *  The fast exponentiation for the convolution product over &#x1D54A;<sub>n</sub>
