import moshpit_utils.units.cgs as cgs
import numpy as np
from moshpit_utils import get_header, get_parameter

def number_density(h5file):
    abar = get_header(h5file,"abar")
    return h5file["density"][:]/(abar * cgs.mH)

