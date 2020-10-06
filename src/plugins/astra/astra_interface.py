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
import re

import solensim.wrapper as wrapper
from solensim.units import *


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

    def __init__(self):
        pass

# Presets, beams
    def beam_presets(self):
        return listdir(self._beamsdir)

    def track_presets(self):
        return listdir(self._tpresetsdir)

    def gen_presets(self):
        return listdir(self._genpresetsdir)

    def load_beam_preset(self, beam):
        toload = os.path.join(self._beamsdir, beam)
        cp = "cp -f %s/* %s"%(toload, self._workdir)
        os.system(cp)
        self._beam_preset = beam
        self.read_beam()

    def loaded_beam_preset(self):
        return self._beam_preset

    beam_preset = property(loaded_beam_preset, load_beam_preset)

    def load_track_preset(self, preset):
        toload = os.path.join(self._tpresetsdir, preset)
        cp = "cp -f %s/* %s"%(toload, self._workdir)
        os.system(cp)
        self._track_preset = preset
        self.read_runfile()

    def get_track_preset(self):
        return self._track_preset

    track_preset = property(get_track_preset, load_track_preset)

    def load_gen_preset(self, preset):
        toload = os.path.join(self._genpresetsdir, preset)
        cp = "cp -f %s/* %s"%(toload, self._workdir)
        os.system(cp)
        self._gen_preset = preset
        self.read_genfile()

    def get_gen_preset(self):
        return self._gen_preset

    gen_preset = property(get_gen_preset, load_gen_preset)

    def save_preset(self, preset, type="track"):
        "Save a preset/beam (currently loaded), modify existing one"
        new = True
        flag = False

        if type=="track":
            target = self._tpresetsdir
            source = os.path.join(self._workdir, "run.in")
        elif type=="gen":
            target = self._genpresetsdir
            source = os.path.join(self._workdir, "generator.in")
        elif type=="beam":
            target = self._beamsdir
            source = os.path.join(self._workdir, "beam.ini")
        else:
            flag = True
            return flag, new  # return error status

        path = os.path.join(target, preset)
        if preset in listdir(target):
            new = False
        else:
            mkdir = "mkdir %s"%path
            os.system(mkdir)

        cp = "cp -f %s %s"%(source, path)
        os.system(cp)

        return flag, new

    def delete_preset(self, preset, type):
        flag = False
        existed = True

        if type=="track":
            target = self._tpresetsdir
        elif type=="gen":
            target = self._genpresetsdir
        elif type=="beam":
            target = self._beamsdir
        else:
            flag = True
            return flag, existed  # return error status

        if preset not in listdir(target):
            existed = False
            flag = True
            return flag, existed

        path = os.path.join(target, preset)
        rm = "rm -r %s"%(path)
        os.system(rm)

        return flag, existed


# Current workspace access
    def write_field(self, z, Bz, file="solenoid.dat"):
        field = {"z":z, "Bz":Bz}
        fielddf = pd.DataFrame.from_dict(field)
        self.field = fielddf
        fielddf.to_csv(os.path.join(self._workdir, file), sep="\t", index=False, header=False)

    def read_field(self, file="solenoid.dat", sep="\t"):
        fielddf = pd.read_table(os.path.join(self._workdir, file), names=["z", "Bz"], engine="python", sep=sep)
        self.field = fielddf
        z = fielddf["z"].values
        Bz = fielddf["Bz"].values
        return z, Bz

    def read_nml(self, file):
        return f90nml.read(os.path.join(self._workdir, file))

    def write_nml(self, nml, file):
        outpath = os.path.join(self._workdir, file)
        nml.write(outpath, force=True, sort=True)

    def read_runfile(self):
        self.runfile = self.read_nml("run.in")

    def read_genfile(self):
        self.genfile = self.read_nml("generator.in")

    def update_runfile(self):
        self.runfile.write(os.path.join(self._workdir,"run.in"), force=True, sort=True)

    def update_genfile(self):
        self.genfile.write(os.path.join(self._workdir,"generator.in"), force=True, sort=True)

    def read_beam(self):
        self.beam = self.get_beam()

    def write_beam(self, beam):
        beam.to_csv(os.path.join(self._workdir, "beam.ini"), sep=" ", index=False, header=False)

    def mop(self, filename):
        file = os.path.join(self._workdir, filename)
        os.system("rm %s"%(file))

    def clean(self):
        self.mop("beam.ini")
        self.mop("run.*")
        self.mop("generator.in")
        self.track_preset = self.track_preset
        self.gen_preset = self.gen_preset
        self.beam_preset = self.beam_preset

    def workspace(self):
        files =  listdir(self._workdir)
        files.remove("Astra")
        files.remove("generator")
        files.remove("NORRAN")
        files = sorted(files)
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
        out = self.run(namelist=namelist, exe="generator")
        self.read_beam()
        return out

# Output
# Column headers:
    _beam_labels = ["x", "y", "z", "px", "py", "pz", "t", "q", "type", "flag"]
    _zemit_labels = ["z", "t", "E", "zrms", "dErms", "epszrms", "z*dEbar/dz"]
    _tracking_labels = ["n", "flag", "z", "x", "y", "Fz", "Fx", "Fy"]
# Beam:
    def get_beam(self):
        path = os.path.join(self._workdir, "beam.ini")
        beam = pd.read_table(path, names=self._beam_labels, skipinitialspace=True, sep=" +", engine="python")
        return beam

# Astra output:
    def read_states(self):
        """
        Reads beam states from all run.pos.001 files;
        outputs a massive dataframe, indexed by state z positions and particle indices
        """
        pattern = re.compile("run.[0-9]+.001")
        files = list(filter(pattern.match, self.workspace()))
        states = []
        zpos = []
        zstop = self.runfile["output"]["zstop"]
        zstart = self.runfile["output"]["zstart"]

        divisor = float(files[-1].split(".")[1])/(zstop-zstart)
        for file in files:
            ident = file.split(".")[1]
            z = float(ident)/divisor + zstart
            path = os.path.join(self._workdir, file)
            state = pd.read_table(path, names=self._beam_labels, skipinitialspace=True, sep=" +", engine="python")
            states.append(state)
            zpos.append(z)
        states.append(self.beam.copy())
        zpos.append(0.0)
        states = pd.concat(states, keys=zpos, names=["zpos", "particle"]).sort_index()
        return states

#    def read_last(self):
#        """
#        Reads beam state last zpos defined in runfile
#        """
#        zstop = self.runfile["output"]["zstop"]
#        ident = str(zstop/cm)[0:-2]
#        if len(ident) == 2: ident = "00"+ident
#        elif len(ident) == 3: ident = "0"+ident
#        ident = "run."+ident+".001"
#        path = os.path.join(self._workdir, ident)
#        screen = pd.read_table(path, names=self._beam_labels, skipinitialspace=True, sep=" +", engine="python")
#        screens = pd.concat([self.beam.copy(), screen], keys = [0.0, zstop], names=["zpos", "particle"])
#        return screens



    def read_zemit(self):
        path = os.path.join(self._workdir, "run.Zemit.001")
        zemit = pd.read_table(path, names=self._zemit_labels, skipinitialspace=True, sep=" +", engine="python")
        return zemit

    def read_trajectories(self):
        path = os.path.join(self._workdir, "run.track.001")
        traj = pd.read_table(path, names=self._tracking_labels, skipinitialspace=True, sep=" +", engine="python")
        return traj
