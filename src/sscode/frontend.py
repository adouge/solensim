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

from sscode.units import *
import sscode.wrapper as wrapper
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
    def __init__(self):
        wrapper.TrackHandle.__init__(self)

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
        print("Loaded default track & generator presets")

    _helptext = """
        This is a helptext.
    """

    def help(self):
        print(self._helptext)

    def run(self, namelist="run.in", exe="Astra"):
        out = astra_interface.Core.run(self, namelist, exe)
        print(out)

    def generate(self, namelist="generator.in"):
        out = astra_interface.Core.run(self, namelist, "generator")
        print(out)
