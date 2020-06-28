# solensim
Solenoid simulation project

### TODO:

### Post-submission:

**Presentation:**

0. Project formulation, motivation _(TODO: rewrite)_
1. Necessary theory beyond basics: _(TODO: condense, remove redundant)_
2. Project methodology: _(TODO: Rewrite)_
    1. Model recap **_(TODO: avoid unneccessary redundancy & repetitions!!!)_**
    2. Algorithm:
        1. General schematic
        2. Optimization, briefly
3. Software presentation: _(TODO: **Rewrite completely**)_
    1. General structure
    2. Interface
    3. Example output - compare to calc
    4. Example design study - tradeoffs **_(TODO: more different examples)_**
4. Conclusions:
    1. **Further development** - pitch ideas
    2. Summary **_(Rewrite!)_**

### Further improvements:
1. Expand model - yoke, FEM etc.
2. Better optimization study
3. Improve interface
4. normalized parameter space
-------

### TODO Documentation:
- This Readme File
- Software manual

### Dependencies:
 - numpy, scipy (found in e.g. Conda)
 - _oct2py package, MATLAB Python interface_ - not used currently

### Structure:
 - User <-commands in, results, plots out-> [Frontend,Wrapper] <--> Backend (pycode, mcode?)
 - File tree:
    - src - software code
    - tex - docus, präsi

### Branches:
 - master - main branch - only chunk in "done" fragments, no more merging master-->dev
 - dev_D - development branch Anton
 - dev_Y - development branch Andrii
