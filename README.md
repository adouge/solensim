# solensim
Solenoid simulation project

### TODO:
- General:
  - Problem, model, design Zusammenfassung
  - Прописать концепт работы, от указанных параметров до результата:
    - цільові характеристики, відносні ваги обмежень --> оптимізація --> аналіз (чутливість, характеристики солющн), презентація характ, чутливості
  - API, plots
  - *chromatic aberrations?*

- Präsentation:
   1. Общая теория e в соленоїдах, aберрації
   2. Цель, target Parameters
   3. Phys. model (формули, інтеграли)
   4. present software:
      - Example - calc with AREAL params
      - Implementation
      - Optimization
      - Results discussion

- Software:
   - backend: *chromatic aberrations?*
   - Output: plots
   - *Documentation!*
   - Finish API

### Problem formulation:
 - Parameters:
   - Geometry: mean radius r, radial width a, axial width b / Rsq, c;
   - I*N scaling factor / N Windungen @ 8 A
 - Provided constraints:
   - FWHM L ~50 mm
   - focal distance f > 50 cm
   - Peak field on axis ~100 mT
   - Energy 3.5 MeV, beam RMS radius 1 mm
 - Optimization:
   - Minimal spheric aberration
   - limit power, material usage

### Dependencies:
 - numpy, scipy (found in e.g. Conda)
 - ~oct2py package, MATLAB Python interface~ - phased out

### Structure:
 - User <-commands in, results, plots out-> [Frontend,Wrapper] <--> Backend (pycode, mcode?)

### Branches: *Неактуально*
 - master - main branch - only chunk in "done" fragments, no more merging master-->dev
 - dev - development branch
