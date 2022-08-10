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
- The object at this point is to explore the algorithm and adjust its complexity and details
  before c++ implementation

