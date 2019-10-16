
import matplotlib.pyplot as plt

from . import fileio
from . import plot
from . import config

class Float(object):
    """Class to read and plot float data"""
    def __init__(self, float_no="lovbio014b"):
        """Create instance of float class"""
        self.float_no = float_no
        self.mld_time, self.mld = fileio.read_mld(float_no=float_no)

    def load(self, varname="chl_smooth"):
        """Load variable from float dat file"""
        self.time,self.depth,self.vals = fileio.read_variable(
            float_no=self.float_no, varname=varname)
        
    def plot(self, varname="chl_smooth"):
        """Plot variable as time-depth panel"""
        self.load(varname=varname)
        ckws = config.get_colordict(varname, vals=self.vals)
        plot.time_depth(self.time, self.depth, self.vals, **ckws)
        plot.mld(self.mld_time, self.mld)
        plot.colorbar(vmin=ckws["vmin"], vmax=ckws["vmax"], cmap=ckws["cmap"])

def process_float(float_no="lovbio014b", varname="chl_smooth"):
    fl = Float(float_no=float_no)
    fl.plot(varname=varname)
    plt.savefig(f"{float_no}_{varname}.png")
    
                
