# solensim
Solenoid simulation project

### TODO:
 - Plan
 - learn biblatex
 - ~MATLAB/Python integration~ done.

### Notes:
 - ~MATLAB engine for Python is not yet updated for 3.8.~ Did a dirty trick to "make" it also "compatible" with Python 3.8.
 - Oct2Py seems to work fine for simpler scripts. It is also much faster, surprisingly, than the MATLAB interface.
 - there is no apparent slowdown from calling Octave through Python, as long as same session objects are used (no Python object recreation/redeletion); also, passing a lot of data through the interface causes slowdowns - better use file buffers for large amounts of num. data?
