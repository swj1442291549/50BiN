import pickle
from pathlib import Path
from operator import itemgetter
from matplotlib import pyplot as plt
import numpy as np
import astropy.units as u
from astropy.time import Time
from astropy.coordinates import SkyCoord, EarthLocation, AltAz

if __name__ == "__main__":
    file_name = "N2301ALL_B.0agcat_cal.pkl"
    mergecat_dict = pickle.load(open(file_name, "rb"))
    (
        nframe,
        medframe_index,
        nstar,
        ndate,
        frame_info,
        nomatch,
        coord,
        psfmagmatch,
        apmagmatch,
        magtype,
        magx,
        ommag,
        ommag_err,
    ) = itemgetter(
        "nframe",
        "medframe_index",
        "medframe_nstar",
        "ndate",
        "frame_info",
        "nomatch",
        "coord",
        "psfmagmatch",
        "apmagmatch",
        "magtype",
        "magx",
        "ommag",
        "ommag_err",
    )(
        mergecat_dict
    )

    noc = 20
    with open("mdev.dat", "r") as f:
        lines = f.readlines()
        i1 = [int(list(filter(None, line.split(" ")))[0]) for line in lines]
        i2 = [int(list(filter(None, line.split(" ")))[1]) for line in lines]
        ncs = [int(list(filter(None, line.split(" ")))[2]) for line in lines]
    if len(ncs) < noc:
        print("too few std stars selected")
        exit()
    else:
        ncs = ncs[:noc]
        print("Std tot: {0:d}".format(noc))
    magtype = "a"
    if magtype == "a":
        magmatch = apmagmatch
    else:
        magmatch = psfmagmatch

    a = np.subtract(magmatch[ncs, :, 0].T, magmatch[ncs, medframe_index, 0])
    b = np.subtract(magx[ncs, :, 0].T, magx[ncs, medframe_index, 0])
    amjd = frame_info.amjd
    ut = frame_info.mid_time
    cc = frame_info[frame_info.mjd == 16677].mid_time
    airmass = frame_info[frame_info.mjd == 16677].airmass
    aa = a[frame_info.mjd == 16677]
    bb = b[frame_info.mjd == 16677]

    cc = list()
    fig, ax = plt.subplots(1, 1)
    for i in range(len(ncs)):
        mm = magx[ncs[i], :, 0]
        ax.scatter(ut, mm, s=3)
        c = np.std(mm[mm > 1])
        cc.append(c)
    ax.set_xlabel("x")
    ax.set_ylabel("y")
    plt.show()
