# solensim
Solenoid simulation project

### TODO:
- General:
   - Разделение работы
   - **Антон прочитать теорию**
   - Прописать концепт работы, от указанных параметров до результата *Антон*
   - ~Форма катушки~ AREAL (see notes)
   - *To Be Added* Нічого не забув такого, що ми говорили по плануванню? *Aндрій*
 
- Präsentation:
   -  Plan

- Software:
   - Концепт работы --> общий алгоритм *Антон*
   - **mcode:** *(Андрій)*
     - *TBA* *(Антон)*
   - *TBA* *(Антон)*


### Notes:
 - Parameters:
   - Geometry: М слоёв по N витков, abstand d, первый слой на радиусе R
 - Backends:
   - Oct2Py seems to work fine for simpler scripts. It is also much faster, surprisingly, than the MATLAB interface.
   - there is no apparent slowdown from calling Octave through Python, as long as same session objects are used (no Python object recreation/redeletion); also, passing a lot of data through the interface causes slowdowns - might be better to use file buffers for large amounts of num. data?
   - *А не проще ли всё запитонить?*
 - Interface:
   - rudimentary looping CLI


### Dependencies:
 - oct2py package, MATLAB Python interface
 - numpy, pandas

### Structure:
 - User <-CLI, data in/out, plots out-> [Frontend,Wrappers] <--> Backend (mcode, python?)
 - Corresponding files:
   - solensim.py - interface and main script
   - frontend.py - user-side functionality handles, link to mcode
   - mdoce - backend (possibly some python)
   - wrapper.py - Oc2Py/MATLAB interfaces, I/O wrappers, plotter?

### Branches:
 - master - main branch - only chunk in "done" fragments, no more merging master-->dev
 - dev - development branch
