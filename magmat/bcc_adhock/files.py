import numpy as np
import sys

class DispForceFile:
    def __init__(self):
        self.filename = ""
        self.ndata = 0
        self.nstart = 0
        self.nend = 0
        self.skip_s = 0
        self.skip_e = 0

class Files:
    def __init__(self):
        self.job_title = ""
        self.datfile_train = DispForceFile()
        self.datfile_validation = DispForceFile()

    # Files::init() in alm
    def initialization(self, prefix_in):
        self.job_title = prefix_in

    def set_prefix(self, prefix_in):
        self.job_title = prefix_in

    def get_prefix(self):
        return self.job_title

    def set_datfile_train(self, dat_in):
        self.datfile_train = dat_in

    def set_datfile_validation(self, dat_in):
        self.datfile_validation = dat_in

    def get_datfile_train(self):
        return self.datfile_train

    def get_datfile_validation(self):
        return self.datfile_validation