# solensim
Solenoid simulation project

### TODO:
 - Plan
 - learn biblatex

### Notes:
 - MATLAB engine for Python is not yet updated for 3.8. Oct2Py seems to work fine for simpler scripts.
 - there is no apparent slowdown from calling Octave through Python, as long as same session objects are used (no Python object recreation/redeletion); also, passing a lot of data through the interface causes slowdowns - better use file buffers for large amounts of num. data
