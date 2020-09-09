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

import sscode.wrapper as wrapper
from sscode.units import *


class Core():
    """
    Shared methods?
    """
    exedir = os.path.abspath("./plugins/astra/")

    def __init__(self):
        self.exename = "test"
        self.exepath = os.path.join(self.exedir, self.exename)

    def run(self, namelist):
        pass

    def read_output(self, todo):
        pass


class Generator(Core):
    """
    Generator-related API
    """
    def __init__(self):
        Core.__init__(self)
        self.exename = "generator"

class Astra(Core):
    """
    ASTRA-related API
    """
    def __init__(self):
        Core.__init__(self)
        self.exename = "Astra"

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
