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
import matplotlib.pyplot as plt
from os import listdir
import os.path

from solensim.units import *
import solensim.wrapper as wrapper
import plugins.astra.astra_interface as astra_interface


class Core(wrapper.CoreHandle):
    """
        Main Interface
    """
    def __init__(self):
        wrapper.CoreHandle.__init__(self)

    _helptext = """
        This is a helptext.
    """

    def help(self):
        print(self._helptext)


class Tracker(wrapper.TrackHandle):
    """
        Dedicated tracking functionality interface
    """
    def __init__(self, astra):
        wrapper.TrackHandle.__init__(self, astra)

    _helptext = """
        This is a helptext.
    """

    def help(self):
        print(self._helptext)


class Astra_Interface(astra_interface.Core):
    """
        Direct access to the astra interface
    """

    def __init__(self):
        astra_interface.Core.__init__(self);
        self.track_preset = "default"
        self.gen_preset = "default"
        self.beam_preset = "default"
        self.verbose = True

    _helptext = """
        General options, presets:
            .verbose - True default, set to False to supress ASTRA stdout piping.
            .presets() - list available and loaded presets
            .track_preset = "preset" to set ASTRA runfile preset,
            .beam_preset = "beam" to choose beam,
            .gen_preset = "preset" to load generator preset.

        Main commands:
            .run(namelist) - run a particular ASTRA input deck, defaults to run.in
            .generate(namelist) - run generator on a distro specification, defaults to generator.in

        Data handles:
            .runfile - loaded ASTRA runfile in dict form, see current setup editing
            .genfile - analogously
            .beam - currently loaded beam file as pandas dataframe
            .field - currently loaded field as pandas dataframe

        Current run setup editing:
            .read_runfile(),
            .read_genfile() - load input namelists into handle object;
                editable as dictionary-type objects.

            .update_genfile(),
            .update_runfile() to save changes into loaded setup

            .get_beam() - returns starting beam distribution as pd dataframe
            .read_beam() - loads beam into handle.beam
            .beam - direct access to loaded starting distribution (beam)
            .write_beam(beam_dataframe) - wrie beam_dataframe to current beam (for, say, tracking same beam multiple times)

        Preset saving/deletion:
            .save_preset(preset_name, preset_type) to create new or overwrite existing preset from loaded setup.
            .delete_preset(preset_name, preset_type) to delete;

        Workspace control:
            .workspace() to view ASTRA workdir
            .clean() to clean everything in workspace and reload presets
            .mop(filename) to delete particular file; works with asterisk patterns (e.g. *.001)

            .read_nml(file) - return contents of file (in workspace) as namelist object
            .write_nml(nml, file) - write nml namelist to file

            .read_field() - return z, Bz from solenoid.dat file, update stored .field
            .write_field(z, Bz) - write solenoid to solenoid.dat, update stored .field

        Output readin:
            .read_screens() - returns dictionary of screen output based on screens specified in &OUTPUT namelist;
                keyed according to screen positions; essentially a collection of .beam dataframes
            .read_trajectories() - returns contents of the trajectory tracking output
            .read_zemit() - returns contents of the Zemit file output
    """

    def help(self):
        print(self._helptext)

    def run(self, namelist="run.in", exe="Astra"):
        out = astra_interface.Core.run(self, namelist, exe)
        if self.verbose: print(out)

    def generate(self, namelist="generator.in"):
        out = astra_interface.Core.run(self, namelist=namelist, exe="generator")
        self.read_beam()
        if self.verbose: print(out)

    def presets(self):
        print("Available presets:")
        print("beam:", self.beam_presets())
        print("track:", self.track_presets())
        print("gen:", self.gen_presets())
        print("\nSet presets:")
        print("beam:", self.beam_preset)
        print("track:", self.track_preset)
        print("gen:", self.gen_preset)

    def save_preset(self, preset, type="track"):
        if type=="track": pot_targets = self.track_presets()
        elif type=="gen": pot_targets = self.gen_presets()
        elif type=="beam": pot_targets = self.beam_presets()
        else: pot_targets = []  # let later error core handle it

        if preset in pot_targets:
            yn = input("Are you sure you want to overwrite %s preset %s? [yes/no]\n"%(type, preset))
        else: yn = "yes"

        if yn == "yes":
            flag, new = astra_interface.Core.save_preset(self, preset, type)
            newstr = "new" if new else "existing"
            if flag:
                print("Failed to save %s %s preset \"%s\""%(newstr, type, preset))
            else:
                print("Successfully saved %s %s preset \"%s\""%(newstr, type, preset))
        else: pass

    def delete_preset(self, preset, type):
        yn = input("Are you sure you want to delete %s preset %s? [yes/no]\n"%(type, preset))
        if yn != "yes": pass
        else:
            flag, existed = astra_interface.Core.delete_preset(self, preset, type)
            if flag:
                exstr = "" if existed else "- preset does not exist!"
                print("Failed to delete %s preset \"%s\" %s"%(type, preset, exstr))
            else:
                print("Successfully deleted %s preset \"%s\""%(type, preset))
