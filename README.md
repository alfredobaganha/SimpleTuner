# SimpleTuner
Simple Guitar Tuner

This is the early stage of the SimpleTunner Project that will eventually become a VST3 plugin

At the current stage we are modeling the overall algorithm in matlab.

The core principle here (at this moment) is using Hamonic Product Spectrum,
where multiple filtered/downsampled copies of a given input are coherently 
multiplied to obtain a peak at the fundamental pitch frequency when observing the
power density function in the frequency domain.

Notes:
- Resampling requires windowing to avoid high frequency noise due to discontinuities
- Need to research a better window function with better properties for now Hanning will work
- The FFTs are computed over a real vector, thus hermitian symmetry is kept and for our
  purposes we only need the upper half of the spectrum, which corresponds to the first 
  half of the output vector
- The objectiv at this point is to explore the algorithm and adjust its complexity and details
  before c++ implementation

Improvements:
- Automatically determine the best (tradeoff) prime for the size of FFT based on the sampling rate.
  The tradeoff here is the frequency resolution and FFT complexity.
- Print out a log file with the detected pitch, expected pitch, error and note name
- Error is not compensated for quantisation error yet
- Reference pitch is now an input parameter