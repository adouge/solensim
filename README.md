# solensim
Solenoid simulation project - DEV branch

### TODO:
 - Plan
 - learn biblatex
 - ~MATLAB/Python integration~ done.

### Notes:
 - ~MATLAB engine for Python is not yet updated for 3.8.~ Did a dirty trick to "make" it also "compatible" with Python 3.8.
 - Oct2Py seems to work fine for simpler scripts. It is also much faster, surprisingly, than the MATLAB interface.
 - there is no apparent slowdown from calling Octave through Python, as long as same session objects are used (no Python object recreation/redeletion); also, passing a lot of data through the interface causes slowdowns - better use file buffers for large amounts of num. data?
 - might have to do a simple interactive looping CLI myself.

### Dependencies:
 - oct2py package

### Structure:
 - User <-CLI, data in/out, plots out-> [Frontend,Wrappers] <--> Backend (mcode)
 - Corresponding files:
   - solensim.py - interface and main script
   - frontend.py - user-side functionality handles, link to mcode
   - mdoce - backend
   - wrapper.py - Oc2Py/MATLAB interfaces, I/O wrappers, plotter?
