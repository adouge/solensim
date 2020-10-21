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

import solensim.backend.core as core
import solensim.backend.optim as optim
import solensim.backend.track as track
from solensim.aux import *

import time
import pandas as pd
import numpy as np

def test_load_mcode_plugin():
    import plugins.mcode.wrapper as mwrapper
    o = mwrapper.OWrapper()
    return o

def load_ini():
    wip()

class TrackHandle(track.TrackModule):
    """
    Interlayer to tracking functionality
    """
    def __init__(self, astra):
        track.TrackModule.__init__(self, astra)
        self.linked_core = None

    # Result storage:
        self.runs = pd.DataFrame()  # run info container
        self.data = {}
        self._run_ticker = 0

    def resolve_label(self, label):
        if type(label)==type(None):
            lbl = self._run_ticker
        else:
            lbl = label
        return lbl

    def init_run(self, rel_decrement, label=None):
        if label not in self.runs.index:
            self._run_ticker += 1
        lbl = self.resolve_label(label)
        track.TrackModule.init_run(self, label=lbl, rel_decrement=rel_decrement)


    # Interaction with core:
    def bind_to_core(self, core):
        self._linked_core = core
        if core != None:
            core.register_track_module(self)
    def get_link_to_bound_core(self):
        return self._linked_core
    linked_core = property(get_link_to_bound_core, bind_to_core)

    # Logging:
    def msg(self, msg):
        dt = pd.Timestamp.fromtimestamp(time.time())
        if self.verbose:  print("%s %s :  %s"%(dt.time(), self.trace, msg))


    def cut_field_edges(self, Bz):
        if Bz[0] != 0 and Bz[-1] != 0:
            background = np.min(Bz)
            Bz_new = Bz-background
            rel_decrement = background/np.max(Bz)
            self.msg("Non-zero field edge, making decrement: %.2e relative to max(B_z)"%rel_decrement)
            return Bz_new, rel_decrement
        else:
            self.msg("Field seems to cut off on itself, no decrement made.")
            return Bz, 0

    def use_field(self, z, Bz, label=None):
        """
        Requirements:
            symmetrical field (Bmax at z=0), ~0 at bounds;
            [z] - m, [Bz] - T
        """
        Bz2, rel_decrement = self.cut_field_edges(Bz)
        self.msg("Updating currently used field.")
        self.field_z = z
        self.field_Bz = Bz2
        self.field_width = (z[-1] - z[0])
        self.astra.write_field(z, Bz2)
        self.init_run(rel_decrement, label=label)

    def use_dat(self, file, sep="\t", label=None):
        """
        Read from file in current dir (temporary)
        Requirements:
            symmetrical field (Bmax at z=0), ~0 at bounds;
            [z] - m, [Bz] - T
        """
        self.msg("Using field from %s."%file)
        fielddf = pd.read_table(file, names=["z", "Bz"], engine="python", sep=sep)
        z = fielddf["z"].values
        Bz = fielddf["Bz"].values
        self.use_field(z, Bz, label=label)

class CoreHandle(core.Core):
    """
    Pre-optim core interlayer
    """
    def __init__(self):
        core.Core.__init__(self)
        self.track = None

    def register_track_module(self, track_module):
        self.track = track_module

    # Logging:
    def msg(self, msg):
        dt = pd.Timestamp.fromtimestamp(time.time())
        if self.verbose:  print("%s %s : %s"%(dt.time(), self.trace, msg))
