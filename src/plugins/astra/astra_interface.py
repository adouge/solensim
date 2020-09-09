#########################################################################
#    Copyright 2020 Anton Douginets
#    This file is part of solensim.
#
#    solensim is free software: you can redistribute it and/or modify
#    it under the terms of the GNU General Public License as published by
#    the Free Software Foundation, either version 3 of the License, or
#    (at your option) any later version.
#
#    solensim is distributed in the hope that it will be useful,
#    but WITHOUT ANY WARRANTY; without even the implied warranty of
#    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#    GNU General Public License for more details.
#
#    You should have received a copy of the GNU General Public License
#    along with solensim.  If not, see <https://www.gnu.org/licenses/>.
#########################################################################

import numpy as np
import pandas as pd
import f90nml
import os.path
import subprocess

import sscode.wrapper as wrapper
from sscode.units import *


class Core():
    """
    TODO
    """
    
    exedir = os.path.abspath("./plugins/astra/")

    def __init__(self, exename="test"):
        self.exename = exename

    def get_exename(self):
        return self._exename
    def set_exename(self, exename):
        self._exename = exename
        self.exepath = os.path.join(self.exedir, exename)
    exename = property(get_exename, set_exename)

    def run(self, namelist):
        """
            Runs Astra/generator with namelist provided in argument
            Returns stdout, to be printed
        """
        arg = [self.exepath, namelist]
        Exe = subprocess.Popen(arg,
            stdout=subprocess.PIPE,
            stderr=subprocess.STDOUT)
        stdout, stderr = Exe.communicate()
        lines = str(stdout).split("\\n")
        lines[0] = lines[0][3:]
        output = "\n".join(lines)
        return output

    def read_output(self, todo):
        pass

class Newrun():
    """
    Instances of Newrun will hold information on a particular run
    and manage it's files
    """
    def __init__(self):
        pass

    def run(self):
        pass

    def clean(self):
        pass
