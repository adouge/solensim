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
from os import listdir
import subprocess

import sscode.wrapper as wrapper
from sscode.units import *


class Core():
    """
    TODO
    """
# Astra/generator path setup
    plugindir = os.path.abspath("./plugins/astra/")
    workdir  = os.path.abspath("./plugins/astra/workspace")
    ssdir = os.path.abspath(".")
    def get_exename(self):
        return self._exename
    def set_exename(self, exename):
        self._exename = exename
        self.exepath = os.path.join(self.workdir, exename)
        self.presetsdir = os.path.join(self.plugindir, "presets", exename)
    exename = property(get_exename, set_exename)

    def __init__(self, exename="test"):
        self.exename = exename
        self.preset = "example"

    def workspace(self):
        "Returns workspace contents"
        ls = subprocess.Popen(["ls", self.workdir],
            stdout=subprocess.PIPE)
        files, stderr = ls.communicate()
        files = str(files).split("\\n")[0:-1]
        files[0] = files[0][2:]
        return files

# Presets
    def presets(self):
        return listdir(self.presetsdir)

    def load_preset(self, preset):
        self.clean()
        toload = os.path.join(self.presetsdir, preset)
        cp = "cp %s/* %s"%(toload, self.workdir)
        os.system(cp)
        self._preset = preset

    def get_preset(self):
        return self._preset
    preset = property(get_preset, load_preset)

    def clean(self):
        def mop(filename):
            file = os.path.join(self.workdir, filename)
            os.system("rm %s"%(file))
        if self.exename == "generator":
            mop("generate.in")
            mop("beam.ini")
        else:
            mop("beam.ini")
            mop("run.*")
            mop("solenoid.dat")

# Current namelist access
## TODO

# Run
    def run(self, namelist):
        """
            Runs Astra/generator with namelist provided in argument
            Returns stdout, to be printed
        """
        cmd = "cd %s; %s %s"%(self.workdir, self.exepath, namelist)
        Exe = subprocess.Popen(cmd,
            shell=True,
            stdout=subprocess.PIPE,
            stderr=subprocess.STDOUT)
        stdout, stderr = Exe.communicate()
        lines = str(stdout).split("\\n")
        lines[0] = lines[0][2:]
        lines[-1] = lines[-1][0:-2]
        output = "\n".join(lines)
        return output

# Output
# Column headers: TODO!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    _beam_labels = ["x", "y", "z", "px", "py", "pz", "t", "q", "type", "flag"]
# Beam:
    def get_beam(self):
        path = os.path.join(self.workdir, "beam.ini")
        beam = pd.read_table(path, names=self._beam_labels, skipinitialspace=True, sep=" +", engine="python")
        return beam

# Astra output:
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
