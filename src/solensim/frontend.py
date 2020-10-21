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

from solensim.aux import *
import solensim.wrapper as wrapper
import plugins.astra.astra_interface as astra_interface


class Core(wrapper.CoreHandle):
    """
        Main Interface
    """
    def __init__(self, trace="core"):
        self.verbose = True  # verbose by default
        self.trace = trace
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
    def __init__(self, astra, trace="track"):
        self.verbose = True  # verbose by default
        self.trace = trace
        self.E = 3.5  # default energy
        wrapper.TrackHandle.__init__(self, astra)


    _helptext = """
        This is a helptext.
    """

    def help(self):
        print(self._helptext)

    # Various plots
    def check_ray_fitting(self, label=None):
        lbl = self.resolve_label(label)
        if lbl not in self.runs.index:
            raise(ValueError("Unknown label: %s"%lbl))
        self.msg("Illustrating trajectory fitting at focal region at label %s"%lbl)
        plt.figure(figsize=(9,9))
        p = self.data[lbl]["s_f"].swaplevel()
        z_solenoid = self.runs.loc[lbl, "z_solenoid"]
        p = p.query("z>@z_solenoid")
        fits = self.data[lbl]["fits"]
        parts = self.data[lbl]["parts"]
        z = np.linspace(p.get("z").min(), p.get("z").max(), num=1000)
        if "drdz" in fits.columns:
            for part in parts:
                plt.plot(z, self.model.axial_trajectory(z, fits.loc[part, "z_f"], fits.loc[part, "drdz"])/mm, "--g")
            plt.plot(z, self.model.axial_trajectory(z, fits.loc[part, "z_f"], fits.loc[part, "drdz"])/mm, "--g", label="Axial approx.")
        else:
            for part in parts:
                plt.plot(z, self.model.off_axis_trajectory(z, fits.loc[part, "z_f"], fits.loc[part, "r_min"], fits.loc[part, "dxdz"], fits.loc[part, "dydz"], fits.loc[part, "cos0"])/mm, "--g")
            plt.plot(z, self.model.off_axis_trajectory(z, fits.loc[part, "z_f"], fits.loc[part, "r_min"], fits.loc[part, "dxdz"], fits.loc[part, "dydz"], fits.loc[part, "cos0"])/mm, "--g", label="Approx. with offset")
        for part in parts:
            plt.plot(p.loc[part, "z"].values, p.loc[part, "r"].values/mm, "+r")
        plt.plot(p.loc[part, "z"].values, p.loc[part, "r"].values/mm, "+r", label="data")
        plt.xlabel("Axial position [m]", fontsize=16)
        plt.ylabel("Radial position [mm]", fontsize=16)
        plt.axis([self.runs.loc[label, "z_focal_left"], self.runs.loc[label, "z_focal_right"], 0, p.get("r").max()*1.05/mm])
        plt.legend(loc="upper center", fontsize=16)

    def check_felddurchgang(self, label=None):
        lbl = self.resolve_label(label)
        if lbl not in self.runs.index:
            raise(ValueError("Unknown label: %s"%lbl))
        self.msg("Illustrating beam state vs. position at label %s"%lbl)
        s = self.data[lbl]["s"]
        z_solenoid = self.runs.loc[lbl, "z_solenoid"]
        zpos = self.data[lbl]["zpos"]
        if "heads" in self.runs.columns and self.runs.loc[lbl, "heads"] == True:
            heads = self.data[lbl]["heads"]
        else:
            heads = self.make_heads(s, zpos)
            self.data[lbl]["heads"] = heads
            self.runs.loc[lbl, "heads"] = True

        plt.figure(figsize=(9,9))
        plt.plot(heads.get("z").values, heads.get("r").values*1/heads["r"].max(), "-k", label="Avg. beam radius [max(r)]")
        plt.plot(heads.get("z").values, heads.get("pr").values*1/heads["pr"].abs().max(), "-r", label="Avg. radial momentum [max(pr)]")
        plt.plot(heads.get("z").values, heads.get("pphi").values*1/heads["pr"].abs().max(), "-b", label="Avg. rot. momentum [max(pr)]")
        plt.plot(heads.get("z").values, heads.get("turn").values, "-g", label="Avg. cum. turn [rad/pi]")
        z, Bz = self.astra.read_field()
        plt.plot(z+z_solenoid, Bz/np.max(Bz), "--k", label="Axial field component [max(Bz)]")
        plt.xlabel("Axial position [m]", fontsize=16)
        plt.ylabel("Arbitrary units", fontsize=16)
        plt.axis([heads["z"].min(), heads["z"].max(), -1, np.max((2, heads["pphi"].abs().max()/heads["pr"].abs().max()))])
        plt.legend(loc="best", fontsize=14)
        plt.grid()
        plt.show()

class Astra_Interface(astra_interface.Core):
    """
        Direct access to the astra interface
    """

    def __init__(self):
        astra_interface.Core.__init__(self);
        self.track_preset = "default"
        self.gen_preset = "uniform"
        self.verbose = False  # not verbose by default
        self.clean()

    _helptext = """
        General options, presets:
            .verbose - False default, set to True to show ASTRA stdout piping.
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
            .field - currently loaded field as pandas dataframe, if loaded

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
            .clean() to clean everything in workspace and reload presets (Attention: overwrites current runfiles!)
                - does not delete user-made files
            .mop(filename) to delete particular file; works with asterisk patterns (e.g. *.001)

            .read_nml(file) - return contents of file (in workspace) as namelist object
            .write_nml(nml, file) - write nml namelist to file

            .read_field(file) - return z, Bz from file (defaults to solenoid.dat), update stored .field
            .write_field(z, Bz, file) - write solenoid to file (defaults to solenoid.dat), update stored .field

        Output readin:
            .read_states() - returns dataframe with screen output based on screens & zphase/start/stop specified in &OUTPUT namelist,
                as well as the initial beam state; keyed according to beam positions
            .read_trajectories() - returns contents of the trajectory tracking output
            .read_zemit() - returns contents of the Zemit file output (WIP)
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
