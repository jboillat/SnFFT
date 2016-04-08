# SnFFT

The SnFFT package is Julia package designed to facilitate harmonic analysis on the symmetric group of order n, denoted Sn . Out of the box, SnFFT implements:

  *  Group operations and factorizations for Sn
  *  Functionality to set up functions over Sn
  *  The fast Fourier transform with additional options if the function is sparse or bandlimited
  *  The inverse fast Fourier transform with additional options if the function is bandlimited or the user is only interested in the result from the top few components
  *  The convolution and correlation of two Fourier transforms
