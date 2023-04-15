import os
import pandas as pd
import numpy as np
import time

class Figure(object):
    """
    Superclass of figures . Each figure is a concrete subclass.
    """
    name = None
    def __init__(self, outpath = os.getcwd() +"/output"):
        self.outpath = outpath
        datafile_name = self.outpath + "/{}.h5".format(self.name)
        self.data = pd.read_hdf(datafile_name, mode="r")

    def save(self, figure_name=None, bbox_inches="tight"):
        print("AAAAAAA")

    def load_true_values(self,filename = "true_values.npy"):
        data_filename = self.outpath + "/{}".format(filename)
        return np.load(data_filename)

class plot_summary(Figure):
    name = "summary"
    def plot(self,  true_values= False):
        print("NNNNNNN")

class Trace(Figure):
    name = "summary"
    def arginfer_trace(self,  true_values= False):
        print("HHHHHHHHH")
