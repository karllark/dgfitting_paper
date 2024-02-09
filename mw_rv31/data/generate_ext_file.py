import numpy as np
import astropy.units as u
from astropy.table import QTable
from dust_extinction.parameter_averages import G23


if __name__ == "__main__":

    waves = np.logspace(np.log10(0.092), np.log10(32.0), num=200) * u.micron
    emod = G23(Rv=3.1)
    ext = emod(waves)
    ext_uncs = 0.01 * ext

    otab = QTable()
    otab["wave"] = waves
    otab["A(l)/A(V)"] = ext
    otab["unc"] = ext_uncs
    otab.write("MW_diffuse_Gordon23_ext.dat", format="ascii.commented_header", overwrite=True)