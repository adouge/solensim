### solensim
Comprehensive solenoid lens design tool.
-------------------------------------------------------------------
**Aims**:
 - flexibly and capably calculate solenoid fields with yoke, and analyze their effect on various electron
   beams via tracking and analytic methods;
 - optimize magnet design in accordance with given search regions and target functionality demands.

**Currently implemented**:
 - ASTRA interface
 - Basic field calculations and characterization

**Currently 100% usable**:
 - ASTRA interface

**Dependencies**:
 - numpy, pandas, scipy
 - ASTRA
 - Has not been tested/configured to work with Windows, i.e. linux only

To get ASTRA executables, sh get_astra.sh inside the (main) solensim directory

To use, please run solensim.py from within this directory in an iPython environment of choice.
