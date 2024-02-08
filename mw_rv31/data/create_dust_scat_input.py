import numpy as np
from astropy.table import QTable

if __name__ == "__main__":
    scat_path = "Scat_Data/"

    files_dgl = [
        "mathis73",
        "morgan76",
        "lillie76",
        "toller81",
        "murthy93",
        "murthy95",
        "petersohn97",
        "witt97",
        "schiminovich01",
        "shalima04",
        "sujatha05",
        "sujatha07",
        "sujatha10",
    ]

    scat_waves = []
    scat_albedo = []
    scat_albedo_unc = []
    scat_g = []
    scat_g_unc = []
    scat_ref = []
    for sfile in files_dgl:
        f = open(scat_path + sfile + ".dat", "r")
        ref = f.readline().rstrip()
        f.close()

        t = QTable.read(
            scat_path + sfile + ".dat", format="ascii", header_start=1
        )
        for k in range(len(t)):
            scat_waves.append(t["wave,"][k])
            scat_albedo.append(t["albedo,"][k])
            scat_albedo_unc.append(t["delta,"][k])
            scat_g.append(t["g,"][k])
            scat_g_unc.append(t["delta"][k])
            scat_ref.append(ref)

    # remove all the measurements with zero uncertainty
    (gindxs,) = np.where(np.array(scat_albedo_unc) > 0.0)
    scat_a_waves = np.array(scat_waves)[gindxs] * 1e-4
    scat_albedo = np.array(scat_albedo)[gindxs]
    scat_albedo_unc = np.array(scat_albedo_unc)[gindxs]
    scat_a_ref = np.array(scat_ref)[gindxs]

    # sort
    sindxs = np.argsort(scat_a_waves)
    scat_a_waves = scat_a_waves[sindxs]
    scat_albedo = scat_albedo[sindxs]
    scat_albedo_unc = scat_albedo_unc[sindxs]
    scat_a_ref = scat_a_ref[sindxs]

    atab = QTable()
    atab["wave"] = scat_a_waves
    atab["albedo"] = scat_albedo
    atab["unc"] = scat_albedo_unc
    atab["ref"] = scat_a_ref

    atab.write("MW_diffuse_DGL_Gordon04_albedo.dat", overwrite=True, format="ascii.commented_header")

    # remove all the measurements with zero uncertainty
    (gindxs,) = np.where(np.array(scat_g_unc) > 0.0)
    scat_g_waves = np.array(scat_waves)[gindxs] * 1e-4
    scat_g = np.array(scat_g)[gindxs]
    scat_g_unc = np.array(scat_g_unc)[gindxs]
    scat_g_ref = np.array(scat_ref)[gindxs]

    # sort
    sindxs = np.argsort(scat_g_waves)
    scat_g_waves = scat_g_waves[sindxs]
    scat_g = scat_g[sindxs]
    scat_g_unc = scat_g_unc[sindxs]
    scat_g_ref = scat_g_ref[sindxs]

    gtab = QTable()
    gtab["wave"] = scat_g_waves
    gtab["g"] = scat_g
    gtab["unc"] = scat_g_unc
    gtab["ref"] = scat_g_ref

    gtab.write("MW_diffuse_DGL_Gordon04_g.dat", overwrite=True, format="ascii.commented_header")
