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


class Core(wrapper.CoreHandle):
    """
        Main Interface
    """
    def __init__(self):
        wrapper.CoreHandle.__init__(self)

class Tracker(wrapper.TrackHandle):
    """
        Dedicated tracking functionality interface
    """
    def __init__(self):
        wrapper.TrackHandle.__init__(self)
