import argparse
import numpy as np
import matplotlib.pyplot as plt

from astropy.wcs import WCS
from astropy.io import fits
from astropy.table import QTable
from astropy.coordinates import SkyCoord
import astropy.units as u

if __name__ == "__main__":

    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--gal", help="get the LMC ave SED", choices=["lmc", "smc"], default="smc"
    )
    parser.add_argument("--png", help="save figure as a png file", action="store_true")
    parser.add_argument("--pdf", help="save figure as a pdf file", action="store_true")
    args = parser.parse_args()

    if args.gal == "smc":

        filenames = [
            "SAGE_SMC_IRAC3.6_2_resid.fits",
            "SAGE_SMC_IRAC4.5_2_resid.fits",
            "SAGE_SMC_IRAC5.8_2_resid.fits",
            "SAGE_SMC_IRAC8.0_2_resid.fits",
            "SAGE_SMC_MIPS24_E012_resid.fits",
            "SAGE_SMC_MIPS70_E012_resid.fits",
            "SMC_Herschel-PACS_100_Feathered.fits",
            "SAGE_SMC_MIPS160_E012_resid.fits",
            "SMC_Herschel-PACS_160_Feathered.fits",
            "SMC_Herschel-SPIRE_250_Feathered.fits",
            "SMC_Herschel-SPIRE_350_Feathered.fits",
            "SMC_Herschel-SPIRE_500_Feathered.fits",
            "SMC_askap_parkes_nhi.fits",
        ]
        datapath = "SMC_Images"

        # fmt: off
        waves = [3.6, 4.5, 5.8, 8.0, 24.0, 70.0, 100.0, 155.0, 165.0, 250., 350., 500.]
        # fmt: on
        posfile = f"~/Python/hst_smc_ext/data/smc_positions_nobump.dat"
    else:
        filenames = [
            "SAGE_LMC_IRAC3.6_2_resid.fits",
            "SAGE_LMC_IRAC4.5_2_resid.fits",
            "SAGE_LMC_IRAC5.8_2_resid.fits",
            "SAGE_LMC_IRAC8.0_2_resid.fits",
            "SAGE_LMC_MIPS24_E12.fits",
            "SAGE_LMC_MIPS70_E12.fits",
            "LMC_Herschel-PACS_100_Feathered.fits",
            "SAGE_LMC_MIPS160_E12.fits",
            "LMC_Herschel-PACS_160_Feathered.fits",
            "LMC_Herschel-SPIRE_250_Feathered.fits",
            "LMC_Herschel-SPIRE_350_Feathered.fits",
            "LMC_Herschel-SPIRE_500_Feathered.fits",
            "lmc_hi_160proj.fits",
        ]
        datapath = "LMC_Images"

        # fmt: off
        waves = [3.6, 4.5, 5.8, 8.0, 24.0, 70.0, 100.0, 155.0, 165.0, 250., 350., 500.]
        # fmt: on
        posfile = f"~/Python/hst_smc_ext/data/lmc_positions_avg.dat"

    ptab = QTable.read(posfile, format="ascii.commented_header")

    n_files = len(filenames)
    n_stars = len(ptab)
    n_waves = len(waves)
    fluxes = np.zeros((n_stars, n_files))

    for kk, filename in enumerate(filenames):

        hdu = fits.open(f"{datapath}/{filename}")[0]
        data = hdu.data
        if "CTYPE3" in list(hdu.header.keys()):
            del hdu.header["CTYPE3"]
            del hdu.header["CDELT3"]
            del hdu.header["CRPIX3"]
            del hdu.header["CRVAL3"]
        wcs = WCS(hdu.header)

        for k in range(len(ptab)):
            coord = SkyCoord(
                ptab["ra"][k],
                ptab["dec"][k],
                unit=(u.hourangle, u.deg),
            )

            pix = wcs.world_to_pixel(coord)
            if (
                (pix[1] >= 0)
                & (pix[1] < data.shape[0] - 1)
                & (pix[0] >= 0)
                & (pix[0] < data.shape[1] - 1)
            ):
                fluxes[k, kk] = data[int(pix[1]), int(pix[0])]

    # normalize to the HI column density
    for k in range(n_files - 1):
        fluxes[:, k] = fluxes[:, k] / fluxes[:, -1]

    # plot
    fontsize = 16

    font = {"size": fontsize}

    plt.rc("font", **font)

    plt.rc("lines", linewidth=1)
    plt.rc("axes", linewidth=2)
    plt.rc("xtick.major", width=2)
    plt.rc("xtick.minor", width=2)
    plt.rc("ytick.major", width=2)
    plt.rc("ytick.minor", width=2)

    figsize = (8, 6)
    fig, ax = plt.subplots(figsize=figsize)

    # normalize to the HI column density
    for k in range(n_files - 1):
        fluxes[:, k] = fluxes[:, k] / fluxes[:, -1]

    for k in range(n_stars):
        ax.plot(waves, fluxes[k, 0:-1], "ko")

    # determine the average
    ave_fluxes = np.average(fluxes[:, 0:n_waves], axis=0)
    unc_fluxes = np.std(fluxes[:, 0:n_waves], axis=0) / np.sqrt(n_stars)
    ax.errorbar(waves, ave_fluxes, fmt="rs", yerr=unc_fluxes)

    ax.set_xlabel("$\lambda$ [$\mu$m]")
    ax.set_ylabel("SB/N(HI) 10$^{-20}$ MJy sr$^{-1}$ (HI atom)$^{-1}$")

    ax.set_xscale("log")
    ax.set_yscale("log")

    fig.tight_layout()

    save_str = "smc_iremission"
    if args.png:
        fig.savefig(f"{save_str}.png")
    elif args.pdf:
        fig.savefig(f"{save_str}.pdf")
    else:
        plt.show()
