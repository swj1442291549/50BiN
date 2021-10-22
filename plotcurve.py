import pickle
from sys import exit
from pathlib import Path
from operator import itemgetter
from matplotlib import pyplot as plt
import click
import numpy as np

@click.command()
@click.argument("file_name", type=str)
def cli(target):
    """Input catalog file from catmerge (.pkl)"""
    pass





if __name__ == "__main__":
    mergecat_file_name = "N2301ALL_B.gcat0.pkl"
    magtype = "a"
    noc = 20
    plot_flag = True

    if not Path(mergecat_file_name).is_file():
        print("File not found!")
        exit()

    mergecat_dict = pickle.load(open(mergecat_file_name, "rb"))

    nframe, medframe_index, nstar, ndate, df_info, nomatch, coord, psfmagmatch, apmagmatch = itemgetter("nframe", "medframe_index", "medframe_nstar", "ndate", "df_info", "nomatch", "coord", "psfmagmatch", "apmagmatch")(mergecat_dict)

    df_info = df_info.assign(amjd=df_info.mjd + df_info.mid_time/24 - 0.5) # convert observing time to modified julian day AMJD(1-nf)

    print("total number of stars: {0:d}".format(nstar))

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
        for i in range(noc):
            print("{0:3d} Std star ID: {1:4d}".format(i, ncs[i]))

    ut1 = df_info.mid_time.iloc[0] - (df_info.mid_time.iloc[nframe - 1] - df_info.mid_time.iloc[0]) / 30
    ut2 = df_info.mid_time.iloc[nframe - 1] + (df_info.mid_time.iloc[nframe - 1] - df_info.mid_time.iloc[0]) / 30

    if magtype == "a":
        magmatch = apmagmatch
    else:
        magmatch = psfmagmatch


    mag2 = np.zeros(nframe)
    err = np.zeros(nframe)
    magx = np.zeros_like(magmatch)
    ommag = np.zeros(nstar)
    ut = df_info.mid_time.values
    magmatch_medframe = magmatch[ncs, medframe_index, 0]
    for ipg in range(nstar):
        js1 = ipg
        for k in range(nframe):
            if magmatch[js1, k, 0] > 30:
                magmatch[js1, k, 0] = 0
                magmatch[js1, k, 1] = 0
            mag2[k] = magmatch[js1, k, 0]
            err[k] = magmatch[js1, k, 1]
            magx[js1, k, 0] = magmatch[js1, k, 0]
            magx[js1, k, 1] = magmatch[js1, k, 1]

            if k != medframe_index:
                y = magmatch[ncs, k, 0]
                dm0 = (y - magmatch_medframe)
                dm0 = dm0 * (y > 1)
                if magmatch[js1, k, 0] > 1:
                    mag2[k] = magmatch[js1, k, 0] - np.sum(dm0) / len(dm0)
                    magx[js1, k, 0] = mag2[k]
                    err[k] = magmatch[js1, k, 1]
                    magx[js1, k, 1] = magmatch[js1, k, 1]
            else:
                mag2[k] = magmatch[js1, medframe_index, 0]
                err[k] = magmatch[js1, medframe_index, 1]

        mag0 = magmatch[js1, :, 0]
        err = magmatch[js1, :, 1]

        znone = mag2[(mag2 > 1) & (mag2 < 25)]
        zmean = np.mean(znone)
        sig = np.sqrt(np.mean((znone - zmean) ** 2))

        magb = mag2[np.abs(mag2-zmean) < 3 * sig]
        utb = ut[(mag2 > 1) & (mag2 < 25) & (np.abs(mag2-zmean) < 3 * sig)]
        ommag[ipg] = np.mean(magb)

        ra = coord[js1, 0]
        rah = int(ra / 15)
        ram = int((ra * 60 / 15 - rah * 60))
        ras = (ra * 240 - rah * 3600 - ram * 60)

