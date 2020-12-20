"""solensim wrapper module."""
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
    """Test."""
    import plugins.mcode.wrapper as mwrapper
    o = mwrapper.OWrapper()
    return o


def load_ini():
    """Load program configuration."""  # TODO: config loading
    wip()


class TrackHandle(track.TrackModule):
    """Interlayer to tracking functionality."""

    def __init__(self, astra):
        """Init interlayer."""
        track.TrackModule.__init__(self, astra)
        self.linked_core = None

    # Result storage:
        self.runs = pd.DataFrame()  # run info container
        self.data = {}
        self._run_ticker = 0

    def drop_labels(self, *args):
        """Delete sotred results at labels specified."""
        L = len(args)
        if L == 1:
            labels = args[0]
            self.msg("Dropping label \"%s\" from result DF." % labels)
            self.data.pop(labels)
        else:
            labels = [*args]
            self.msg("Dropping %d labels from result DF." % L)
            for label in labels:
                self.data.pop(label)
        self.runs = self.runs.drop(index=labels)

    def resolve_label(self, label):
        """Check if label exists and assign a run counter if needed."""
        if type(label) == type(None):
            lbl = self._run_ticker
        else:
            lbl = label
        return lbl

    def init_run(self, rel_decrement, label=None):
        """Initialize new result entry."""
        if label not in self.runs.index:
            self._run_ticker += 1
        lbl = self.resolve_label(label)
        track.TrackModule.init_run(self, label=lbl, rel_decrement=rel_decrement)

# Interaction with core:
    def bind_to_core(self, core):
        """Connect with an instance of calc module."""
        self._linked_core = core
        if core is not None:
            core.register_track_module(self)

    def get_link_to_bound_core(self):
        """Get connected calc module instance."""
        return self._linked_core

    linked_core = property(get_link_to_bound_core, bind_to_core)

    # Logging:
    def msg(self, msg):
        """Log msg with caller signature and timestamp."""
        dt = pd.Timestamp.fromtimestamp(time.time())
        if self.verbose:
            print("%s %s :  %s" % (dt.time(), self.trace, msg))

    def cut_field_edges(self, Bz):
        """Cut field edges."""  # TODO: docstring
        if Bz[0] != 0 and Bz[-1] != 0:
            background = np.min(Bz)
            Bz_new = Bz-background
            rel_decrement = background/np.max(Bz)
            self.msg("Non-zero field edge, making decrement: %.2e relative to max(B_z)" % rel_decrement)
            return Bz_new, rel_decrement
        else:
            self.msg("Field seems to cut off on itself, no decrement made.")
            return Bz, 0

    def use_field(self, z, Bz, label=None, normalize=False):
        """
        Select field for furhter handling.

        Requirements:
            symmetrical field (Bmax at z=0), ~0 at bounds;
            [z] - m, [Bz] - T
        """
        self.msg(">>> New input field!")
        if normalize:
            self.msg("Normalizing input field...")
            self.linked_core.sample_field(z, Bz)
            k = self.linked_core.get_scale_factor(self.baseline_f, self.E)
            Bz2 = Bz*k
        else:
            Bz2 = Bz
        Bz2, rel_decrement = self.cut_field_edges(Bz2)
        self.field_z = z
        self.field_Bz = Bz2
        self.field_width = (z[-1] - z[0])
        self.astra.write_field(z, Bz2)
        self.init_run(rel_decrement, label=label)

    def use_dat(self, file, sep="\t", label=None, normalize=False):
        """
        Read from file in current dir (temporary).

        Requirements:
            symmetrical field (Bmax at z=0), ~0 at bounds;
            [z] - m, [Bz] - T
        """
        self.msg("Using field from %s." % file)
        fielddf = pd.read_table(file, names=["z", "Bz"], engine="python", sep=sep)
        z = fielddf["z"].values
        Bz = fielddf["Bz"].values
        self.use_field(z, Bz, label=label, normalize=normalize)


class CoreHandle(core.Core):
    """Pre-optim core interlayer."""

    def __init__(self):
        """Init the interlayer."""
        core.Core.__init__(self)
        self.track = None

    def register_track_module(self, track_module):
        """Connect to an instance of track module."""
        self.track = track_module

    # Logging:
    def msg(self, msg):
        """Log msg with caller signature and timestamp."""
        dt = pd.Timestamp.fromtimestamp(time.time())
        if self.verbose:
            print("%s %s : %s" % (dt.time(), self.trace, msg))
