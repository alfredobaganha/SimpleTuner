tuner_filt is the current working implementation.
- added new test samples (freesound.org)
- implemented new sets of FIR filters with different assumptions
- settled with simple fir filtering without decimation (very good results)
- expanded the current small prime test set to include base 10 and 11
- fixed minor logic errors from tuner_mul
- current tests have been consistent

Observed a some logic error on tuner_mul. Plus it is also very labor intensive
with the resampling.

[DEPRECATED] tuner_mul is the actual working implementation

tuner was an early approach a and used a poor interpretation of the properties 
desired. just dont use it. keeping it for the sake of tracking progress