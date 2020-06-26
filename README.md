# solensim
Solenoid simulation project

### TODO:
 - Documentation:
  - This Readme File
  - Software manual

 - Improvements:
  1. Model:
    - chromatic aberrations?
    - yoke modeling?
    - better field model?
  2. Software:
    - finish API
    - comment code
    - package/port mcode properly

 - Presentation structure, _TODO_:
  0. Project formulation, motivation
  1. Necessary theory beyond basics: _(TODO: condense, remove redundant)_
    1. Електронная оптика, общая теория, почему именно она?:
       1. Классическая оптика -> Электронная оптика (_1 слайд_)
    2. Фокусна відстань, Общая теория e в соленоїдах:
       1. Формулировка требований к линзе, вытекающих из поставленной задачи (_2 слайд_)
       2. Представление соленоида, как подходящего кандидата (что такое соленоид, два графика: теор. схема и пример чертежа реального соленоида в ускорителе) (_3 слайд_)
       3. Рассчёт траектории полёта одного електрона через линзу:
         1. Исходные уравнения и следующая из них система диф. уравнений (_4 слайд_)
         2. Решение системы с применением свойств поля и закона и сох. энерг. (_5 слайд_)
       4. Выведение формулы фокусного расстояния (_6 слайд_)
       5. Применение обобщений и презентация конечной формулы Z компоненты магн. индукции для рассчёта фок. рас. и дефектов (_7 слайд_)
    3. Дефекты линзы, сфер. аберрация, эмиссия, размер пучка в фокусе:
       1. Откуда берутся дефекты, математическое объяснение (_8 слайд_)
       2. Типы дефетков (кратко, и сразу переход на эмиссию) (_9 слайд_)
       3. Эмиссия (кратко, изначальная и конечные формулы на слайде) (_10 слайд_)
       4. Сфер аберрация как самый важный деффект в нашем случае (объяснить важность) (_11 слайд_)
    4. Выводы:
       1. Полное теоретическое описание проблемы поиска параметров с минимальными дефектами, конечные формулы с которыми придётся работать (_12 слайд_)
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


### Dependencies:
 - numpy, scipy (found in e.g. Conda)
 - ~oct2py package, MATLAB Python interface~ - phased out

### Structure:
 - User <-commands in, results, plots out-> [Frontend,Wrapper] <--> Backend (pycode, mcode?)

### Branches:
 - master - main branch - only chunk in "done" fragments, no more merging master-->dev
 - dev_d - development branch Anton
 - dev_y - development branch Andrii
