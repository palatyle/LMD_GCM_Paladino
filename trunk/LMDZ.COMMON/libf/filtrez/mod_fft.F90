MODULE mod_fft

#ifdef FFT_MATHKEISAN
  USE mod_fft_mathkeisan
#else
#ifdef FFT_FFTW
  USE mod_fft_fftw
#else
#ifdef FFT_MKL
  USE mod_fft_mkl
#else
  USE mod_fft_wrapper
#endif
#endif
#endif

END MODULE mod_fft
