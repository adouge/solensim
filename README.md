# solensim
Solenoid simulation project

### TODO:
- General:
  - Problem, model, design Zusammenfassung
  - Прописать концепт работы, от указанных параметров до результата:
    - цільові характеристики, відносні ваги обмежень --> оптимізація --> аналіз (чутливість, характеристики солющн), презентація характ, чутливості
  - ~replace NonlinearConstraint with LinearConstraint where possible~
  - finish backend
  - bridge backend to wrapper
  - iPython API, frontend (output)


- Präsentation: **Формули у model.tex, потім переведу у слайди**
   1. Общая теория e в соленоїдах, aберрації
   2. Цель, target Parameters
   3. Phys. model (формули, інтеграли)
   4. present software:
      - Example - calc with AREAL params
      - Implementation
      - Optimization
      - Results

- Software:
   - mcode *(Андрій)* - *а надо?*
   - pycode *(Антон)*
     - ~backend base~
     - ~get rid of bloat in backend~
     - review parameter model?
     - Flexible constraints definition - they influence CTR output strongly:
        - **ALLOW UNCONSTRAINED**
        - ~starting parameter values from bounds~
     - Optimization:
        - **Optimize in focal spot size?**
        - ctr
        - COBYLA, SLSQP?
        - try Rsq, c instead of r, a, b?
     - Output
     - *Documentation!*
   - Wrapper: *(Aнтон)*
     - **iPython API**
     - flexible interface for iPython integration --> don't overcomplicate; provide both wrapped API for future CLI, as well as a more pythonic, naked user-facing methods for iPython


### Notes:
Optimization:
 - https://docs.scipy.org/doc/scipy/reference/generated/scipy.optimize.minimize.html
 - https://docs.scipy.org/doc/scipy/reference/optimize.minimize-cobyla.html#optimize-minimize-cobyla
 - https://docs.scipy.org/doc/scipy/reference/optimize.minimize-slsqp.html#optimize-minimize-slsqp
 - https://docs.scipy.org/doc/scipy/reference/generated/scipy.optimize.OptimizeResult.html#scipy.optimize.OptimizeResult
 - https://docs.scipy.org/doc/scipy/reference/generated/scipy.optimize.NonlinearConstraint.html#scipy.optimize.NonlinearConstraint

### Problem:
 - Parameters:
   - Geometry: mean radius r, radial width a, axial width b;
   - I*N scaling factor
   - Alt. geometry: *optimize in Rsq, c?*
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
 - numpy, scipy (found in e.g. Conda)

### Structure:
 - User <-CLI, data in/out, plots out-> [Frontend,Wrappers] <--> Backend (mcode, python?)
 - Corresponding files:
   - solensim.py - interface and main script
   - frontend.py - user-side functionality handles, link to mcode
   - mdoce/pycode - backend
   - wrapper.py - Oc2Py/MATLAB interfaces, pycode handle, I/O wrappers, plotter?

### Branches: *Неактуально*
 - master - main branch - only chunk in "done" fragments, no more merging master-->dev
 - dev - development branch
