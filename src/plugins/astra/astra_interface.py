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
    beamsdir = os.path.abspath("./plugins/astra/presets/beams")

    def get_exename(self):
        return self._exename
    def set_exename(self, exename):
        self._exename = exename
        self.exepath = os.path.join(self.workdir, exename)
        self.presetsdir = os.path.join(self.plugindir, "presets", exename)
    exename = property(get_exename, set_exename)

    def __init__(self, exename="Astra"):
        self.exename = exename
        #self.preset = "example"
        #self.beam_preset = "example"

    def workspace(self):
        files =  listdir(self.workdir)
        files.remove("Astra")
        files.remove("generator")
        files.remove("NORRAN")
        return files

# Presets, beams
    def beams(self):
        return listdir(self.beamsdir)

    def read_beam(self):
        self.beam = self.get_beam()

    def load_beam_preset(self, beam):
        toload = os.path.join(self.beamsdir, beam)
        cp = "cp %s/* %s"%(toload, self.workdir)
        os.system(cp)
        self._beam_preset = beam
        self.beam = self.get_beam()

    def loaded_beam_preset(self):
        return self._beam_preset

    beam_preset = property(loaded_beam_preset, load_beam_preset)


    def presets(self):
        return listdir(self.presetsdir)

    def load_preset(self, preset):
        toload = os.path.join(self.presetsdir, preset)
        cp = "cp %s/* %s"%(toload, self.workdir)
        os.system(cp)
        self._preset = preset
        if preset in self.beams() and self.exename=="Astra": self.beam_preset = preset
        self.read_runfile()

    def get_preset(self):
        return self._preset

    preset = property(get_preset, load_preset)




    def clean(self):
        def mop(filename):
            file = os.path.join(self.workdir, filename)
            os.system("rm %s"%(file))
        mop("beam.ini")
        mop("run.*")
        mop("solenoid.dat")

# Current namelist access
    def read_nml(self, file):
        return f90nml.read(os.path.join(self.workdir, file))

    def write_nml(self, nml, file):
        outpath = os.path.join(self.workdir, file)
        nml.write(outpath, force=True)

    def read_runfile(self):
        self.runfile = self.read_nml("run.in")

    def update_runfile(self):
        self.runfile.write(os.path.join(self.workdir,"run.in"), force=True)


# Run
    def run(self, namelist="run.in"):
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
