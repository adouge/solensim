# solensim
Solenoid simulation project

### TODO:
- General:
   - **Разделение работы?**
   - *Антон прочитать теорию* -> Problem, model, design Zusammenfassung 
   - Прописать концепт работы, от указанных параметров до результата
 
- Präsentation:
   - Теор. выкладка
   - Выжимка см. выше
   - *TBA*

- Software:
   - Концепт работы --> общий алгоритм *Антон*
   - mcode *(Андрій)*
   - pycode *(Антон)*
     - ~backend base~
     - review parameter model 
     - Constraints definition
     - Optimization
     - Output
   - Wrapper: *(Aнтон)*
     - flexible interface for iPython integration --> don't overcomplicate; provide both wrapped API for future CLI, as well as a more pythonic, naked user-facing methods for iPython


### Notes:
### Problem:
 - Parameters:
   - Geometry: mean radius r, radial width a, axial width b; 
   - I*N scaling factor
 - Constraints: 
   - FWHM L ~50 mm
   - focal distance f > 50 cm
   - Peak field on axis ~100 mT
   - Energy 3.5 MeV, beam RMS radius 1 mm
   - *2pi transversal rotation phase???*
 - Optimization:
   - Minimal spheric aberration
   - limit power, material usage

### Design:
 - Backends:
   - Oct2Py seems to work fine for simpler scripts. It is also much faster, surprisingly, than the MATLAB interface.
   - there is no apparent slowdown from calling Octave through Python, as long as same session objects are used (no Python object recreation/redeletion); also, passing a lot of data through the interface causes slowdowns - might be better to use file buffers for large amounts of num. data?
   - *А не проще ли всё запитонить?* **Проще.**
 - Interface:
   - iPython command API
   - rudimentary looping CLI


### Dependencies:
 - *oct2py package, MATLAB Python interface (optional, doubted)*
 - numpy, scipy

### Structure:
 - User <-CLI, data in/out, plots out-> [Frontend,Wrappers] <--> Backend (mcode, python?)
 - Corresponding files:
   - solensim.py - interface and main script
   - frontend.py - user-side functionality handles, link to mcode
   - mdoce/pycode - backend
   - wrapper.py - Oc2Py/MATLAB interfaces, I/O wrappers, plotter?

### Branches: *(fucked here up)*
 - master - main branch - only chunk in "done" fragments, no more merging master-->dev
 - dev - development branch
