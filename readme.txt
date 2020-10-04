solensim v0.3.4
Comprehensive solenoid lens design tool.
-------------------------------------------------------------------
Aims:
 - flexibly and capably calculate solenoid fields with yoke, and analyze their effect on various electron
   beams via tracking and analytic methods;
 - optimize magnet design in accordance with given search regions and target functionality demands.

Currently implemented:
 - ASTRA interface
 - Basic field calculations and characterization
 - Field characterization via track module

Currently 100% usable:
 - ASTRA interface

Dependencies:
 - numpy, pandas, scipy, f90nml
 - ASTRA
 - Has not been tested/configured to work with Windows, i.e. linux only

To get ASTRA executables, sh get_astra.sh inside the (main) solensim directory

To use, please run solensim.py from within this directory in an iPython environment of choice.
-------------------------------------------------------------------

License notice:

    Copyright 2020 Anton Douginets, Andrii Yanovets

    solensim is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    solensim is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with solensim.  If not, see <https://www.gnu.org/licenses/>.
