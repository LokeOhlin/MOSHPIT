import h5py as hp
import numpy as np


def get_parameter(h5file, param):
    return h5file["/Parameters"].attrs.get(param)

def get_header(h5file, head):
    return h5file["/Headers"].attrs.get(head)
def current_time(h5file):
    return h5file["/Headers"].attrs.get("time")
