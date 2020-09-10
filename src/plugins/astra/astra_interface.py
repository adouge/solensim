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
    _plugindir = os.path.abspath("./plugins/astra/")
    _workdir  = os.path.abspath("./plugins/astra/workspace")
    _ssdir = os.path.abspath(".")
    _beamsdir = os.path.abspath("./plugins/astra/presets/beams")
    _tpresetsdir = os.path.abspath("./plugins/astra/presets/track")
    _genpresetsdir = os.path.abspath("./plugins/astra/presets/generator")


#    def get_exename(self):
#        return self._exename
#    def set_exename(self, exename):
#        self._exename = exename
#        self.exepath = os.path.join(self._workdir, exename)
#        self.presetsdir = os.path.join(self._plugindir, "presets", exename)
#    exename = property(get_exename, set_exename)

    def __init__(self):
        pass

# Presets, beams
    def beams(self):
        return listdir(self._beamsdir)

    def track_presets(self):
        return listdir(self._tpresetsdir)

    def gen_presets(self):
        return listdir(self._genpresetsdir)

    def load_beam_preset(self, beam):
        toload = os.path.join(self._beamsdir, beam)
        cp = "cp %s/* %s"%(toload, self._workdir)
        os.system(cp)
        self._beam_preset = beam
        self.beam = self.get_beam()

    def loaded_beam_preset(self):
        return self._beam_preset

    beam_preset = property(loaded_beam_preset, load_beam_preset)

    def load_track_preset(self, preset):
        toload = os.path.join(self._tpresetsdir, preset)
        cp = "cp %s/* %s"%(toload, self._workdir)
        os.system(cp)
        self._track_preset = preset
        if preset in self.beams(): self.beam_preset = preset
        self.read_runfile()

    def get_track_preset(self):
        return self._track_preset

    track_preset = property(get_track_preset, load_track_preset)

    def load_gen_preset(self, preset):
        toload = os.path.join(self._genpresetsdir, preset)
        cp = "cp %s/* %s"%(toload, self._workdir)
        os.system(cp)
        self._gen_preset = preset
        self.read_genfile()

    def get_gen_preset(self):
        return self._gen_preset

    gen_preset = property(get_gen_preset, load_gen_preset)




# Current workspace access
    def write_field(self, z, Bz):
        field = {"z":z, "Bz":Bz}
        fielddf = pd.DataFrame.from_dict(field)
        fielddf.to_csv(os.path.join(self._workdir, "solenoid.dat"), sep="\t", index=False, header=False)

    def get_field(self):
        fielddf = pd.read_table(os.path.join(self._workdir,"solenoid.dat"), names=["z", "Bz"])
        z = fielddf["z"].values
        Bz = fielddf["Bz"].values
        return z, Bz

    def read_nml(self, file):
        return f90nml.read(os.path.join(self._workdir, file))

    def write_nml(self, nml, file):
        outpath = os.path.join(self._workdir, file)
        nml.write(outpath, force=True)

    def read_runfile(self):
        self.runfile = self.read_nml("run.in")

    def read_genfile(self):
        self.genfile = self.read_nml("generator.in")

    def update_runfile(self):
        self.runfile.write(os.path.join(self._workdir,"run.in"), force=True)

    def update_genfile(self):
        self.genfile.write(os.path.join(self._workdir,"generator.in"), force=True)

    def read_beam(self):
        self.beam = self.get_beam()

    def mop(self, filename):
        file = os.path.join(self._workdir, filename)
        os.system("rm %s"%(file))

    def clean(self):
        self.mop("beam.ini")
        self.mop("run.*")
        self.mop("generator.in")
        self.mop("solenoid.dat")

    def workspace(self):
        files =  listdir(self._workdir)
        files.remove("Astra")
        files.remove("generator")
        files.remove("NORRAN")
        return files

# Run
    def run(self, namelist="run.in", exe="Astra"):
        """
            Runs Astra/generator with namelist provided in argument
            Returns stdout, to be printed
        """
        exepath = os.path.join(self._workdir, exe)
        cmd = "cd %s; %s %s"%(self._workdir, exepath, namelist)
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

    def generate(self, namelist="generator.in"):
        return self.run(namelist=namelist, exe="generator")

# Output
# Column headers: TODO!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    _beam_labels = ["x", "y", "z", "px", "py", "pz", "t", "q", "type", "flag"]

# Beam:
    def get_beam(self):
        path = os.path.join(self._workdir, "beam.ini")
        beam = pd.read_table(path, names=self._beam_labels, skipinitialspace=True, sep=" +", engine="python")
        return beam

# Astra output:
    def read_output(self, todo):
        pass
