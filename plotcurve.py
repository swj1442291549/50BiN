import pickle
from sys import exit
from pathlib import Path
from operator import itemgetter
from matplotlib import pyplot as plt
import click
import numpy as np


@click.command()
@click.argument("input_file_name", type=str)
@click.option(
    "--magtype", type=str, default="a", help="magnitude type. a: aperture; p: psf"
)
@click.option("--noc", type=int, default=20, help="Number of selected standard stars")
@click.option("--plot_flag", type=bool, default=True, help="Enter interactive plot")
def cli(input_file_name, magtype, noc, plot_flag):
    """Input catalog file from catmerge (.pkl)"""
    if not Path(input_file_name).is_file():
        print("File not found!")
        exit()


    mergecat_dict = pickle.load(open(input_file_name, "rb"))

    if "magx" in mergecat_dict.keys():
        print("This catalog has alreadly contained corrected photometry")
        if plot_flag:
            plot_lc(input_file_name)
            return
        else:
            exit()

    (
        nframe,
        medframe_index,
        nstar,
        ndate,
        df_info,
        nomatch,
        coord,
        psfmagmatch,
        apmagmatch,
    ) = itemgetter(
        "nframe",
        "medframe_index",
        "medframe_nstar",
        "ndate",
        "df_info",
        "nomatch",
        "coord",
        "psfmagmatch",
        "apmagmatch",
    )(
        mergecat_dict
    )

    df_info = df_info.assign(
        amjd=df_info.mjd + df_info.mid_time / 24 - 0.5
    )  # convert observing time to modified julian day AMJD(1-nf)

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

    ut1 = (
        df_info.mid_time.iloc[0]
        - (df_info.mid_time.iloc[nframe - 1] - df_info.mid_time.iloc[0]) / 30
    )
    ut2 = (
        df_info.mid_time.iloc[nframe - 1]
        + (df_info.mid_time.iloc[nframe - 1] - df_info.mid_time.iloc[0]) / 30
    )
    amjd = df_info.mjd + df_info.mid_time / 24 - 0.5

    if magtype == "a":
        magmatch = apmagmatch
    else:
        magmatch = psfmagmatch

    magx = np.zeros_like(magmatch)
    ommag = np.zeros(nstar)
    ommag_err = np.zeros(nstar)
    ut = df_info.mid_time.values
    magmatch_medframe = magmatch[ncs, medframe_index, 0]
    for ipg in range(nstar):
        magmatch[ipg, magmatch[ipg, :, 0] > 30, 1] = 0
        magmatch[ipg, magmatch[ipg, :, 0] > 30, 0] = 0
    for ipg in range(nstar):
        mag2 = np.zeros(nframe)
        err = np.zeros(nframe)
        js1 = ipg
        for k in range(nframe):
            mag2[k] = magmatch[js1, k, 0]
            err[k] = magmatch[js1, k, 1]
            magx[js1, k, 0] = magmatch[js1, k, 0]
            magx[js1, k, 1] = magmatch[js1, k, 1]

            if k != medframe_index:
                y = magmatch[ncs, k, 0]
                dm0 = y - magmatch_medframe
                dm0 = dm0[y > 1]
                if magmatch[js1, k, 0] > 1:
                    mag2[k] = magmatch[js1, k, 0] - np.sum(dm0) / len(dm0)
                    magx[js1, k, 0] = mag2[k]
                    err[k] = magmatch[js1, k, 1]
                    magx[js1, k, 1] = magmatch[js1, k, 1]

            else:
                mag2[k] = magmatch[js1, medframe_index, 0]
                err[k] = magmatch[js1, medframe_index, 1]

        znone = mag2[(mag2 > 1) & (mag2 < 25)]
        if len(znone) != 0:
            zmean = np.mean(znone)
            sig = np.sqrt(np.mean((znone - zmean) ** 2))

            magb = mag2[np.abs(mag2 - zmean) < 3 * sig]
            errb = magmatch[js1, :, 1][np.abs(mag2 - zmean) < 3 * sig]
            utb = ut[(mag2 > 1) & (mag2 < 25) & (np.abs(mag2 - zmean) < 3 * sig)]
            if len(magb) != 0:
                ommag[ipg] = np.mean(magb)
                ommag_err[ipg] = np.mean(err)

    # Save average magnitude for each star
    mmag_catfile_name = "{0}.{1}{2}gcat_mmag".format(
        input_file_name.split(".")[0], input_file_name.split(".")[1][0], magtype
    )
    with open(mmag_catfile_name, "w") as f:
        for i in range(nstar):
            f.write(
                "{0:5d} {1:15.8f} {2:15.8f} {3:10.5f} {4:10.5f}\n".format(
                    i, coord[i, 0], coord[i, 1], ommag[i], ommag_err[i]
                )
            )

    # Save final catalog
    final_catfile_name = "{0}.{1}{2}gcat_cal.pkl".format(
        input_file_name.split(".")[0], input_file_name.split(".")[1][0], magtype
    )

    mergecat_dict = {
        "nframe": nframe,
        "medframe_index": medframe_index,
        "medframe_nstar": nstar,
        "ndate": ndate,
        "df_info": df_info,
        "nomatch": nomatch,
        "coord": coord,
        "psfmagmatch": psfmagmatch,
        "apmagmatch": apmagmatch,
        "magtype": magtype,
        "magx": magx,
        "ommag": ommag,
        "ommag_err": ommag_err,
    }
    pickle.dump(mergecat_dict, open(final_catfile_name, "wb"))

    if plot_flag:
        plot_lc(final_catfile_name)


def coord_to_str(ra, dec):
    rah = int(ra / 15)
    ram = int((ra * 60 / 15 - rah * 60))
    ras = ra * 240 - rah * 3600 - ram * 60
    decd = int(dec)
    decm = int(dec * 60 - decd * 60)
    decs = dec * 3600 - decd * 3600 - decm * 60
    ra_str = "{0:2d}:{1:0>2d}:{2:0>5.2f}".format(rah, ram, ras)
    dec_str = "{0:3d}:{1:0>2d}:{2:0>5.2f}".format(decd, decm, decs)
    return ra_str, dec_str


def save_single_phot(input_file_name, magtype, star_index, amjd, coord, magmatch, magx):
    ra_str, dec_str = coord_to_str(coord[star_index][0], coord[star_index][1])
    orig_file_name = "{0}.{1}{3}S{2:0>4d}.orig".format(
        input_file_name.split(".")[0],
        input_file_name.split(".")[1][0],
        star_index,
        magtype,
    )
    mag_file_name = "{0}.{1}{3}S{2:0>4d}.dat".format(
        input_file_name.split(".")[0],
        input_file_name.split(".")[1][0],
        star_index,
        magtype,
    )
    print(
        "Save photometry data of #{0:3d} to {1}, {2}".format(
            star_index, orig_file_name, mag_file_name
        )
    )
    with open(orig_file_name, "w") as f:
        f.write(
            "{0} {1} {2} (J2000) {3}\n".format(star_index, ra_str, dec_str, magtype)
        )
        for i in range(magmatch.shape[0]):
            f.write(
                "{0:8d} {1:7d} {2:10.5f} {3:10.5f} {4:10.5f}\n".format(
                    i,
                    int(amjd[star_index]),
                    (amjd[star_index] - int(amjd[star_index]) + 0.5) * 24,
                    magmatch[i, star_index, 0],
                    magmatch[i, star_index, 1],
                )
            )
    with open(mag_file_name, "w") as f:
        f.write(
            "{0} {1} {2} (J2000) {3}\n".format(star_index, ra_str, dec_str, magtype)
        )
        for i in range(magmatch.shape[0]):
            f.write(
                "{0:8d} {1:7d} {2:10.5f} {3:10.5f} {4:10.5f}\n".format(
                    i,
                    int(amjd[star_index]),
                    (amjd[star_index] - int(amjd[star_index]) + 0.5) * 24,
                    magx[i, star_index, 0],
                    magx[i, star_index, 1],
                )
            )


def prepare_lc_df(star_index, df_info, magmatch, magx):
    lc = df_info.copy()
    lc = lc.assign(mag=magmatch[star_index, :, 0])
    lc = lc.assign(mag_err=magmatch[star_index, :, 1])
    lc = lc.assign(magx=magx[star_index, :, 0])
    lc = lc.assign(magx_err=magx[star_index, :, 1])
    lc = lc[~np.isnan(lc.magx) & (lc.magx > 1)]
    return lc


def plot_lc(file_name):
    if "s" in plt.rcParams["keymap.save"]:
        plt.rcParams["keymap.save"].remove("s")
    # if ''
    mergecat_dict = pickle.load(open(file_name, "rb"))
    (
        nframe,
        medframe_index,
        nstar,
        ndate,
        df_info,
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
        "df_info",
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
    if magtype == "a":
        magmatch = apmagmatch
    else:
        magmatch = psfmagmatch

    mjd_list = list(set(df_info.mjd))

    global star_index
    global mday_flag
    global mjd_index
    star_index = 0
    mjd_index = -1
    mday_flag = False
    lc = prepare_lc_df(star_index, df_info, magmatch, magx)

    def on_key(event):
        global star_index
        global mday_flag
        global mjd_index
        if event.key == "n":
            if mday_flag:
                mjd_index += 1
                mjd_index = mjd_index % len(mjd_list)
            else:
                star_index += 1
        elif event.key == "N":
            if mday_flag:
                mjd_index -= 1
                mjd_index = mjd_index % len(mjd_list)
            else:
                star_index -= 1
        elif event.key == "s":
            save_single_phot(
                file_name, magtype, star_index, df_info["amjd"], coord, magmatch, magx
            )
            return
        elif event.key == "d":
            mday_flag = not mday_flag
            if not mday_flag:
                mjd_index = -1

        plt.clf()
        lc = prepare_lc_df(star_index, df_info, magmatch, magx)
        if mjd_index != -1:
            lc = lc[lc.mjd == mjd_list[mjd_index]]
        if mday_flag:
            x = lc.amjd - int(min(lc.amjd))
            y1 = lc.mag
            y2 = lc.magx
            xlabel = "MJD"
            if mjd_index == -1:
                title = "#: {0:4d}  M = {1:.2f}  RA: {3}  DEC: {4}  MJD+{2:d}".format(
                    star_index, ommag[star_index], int(min(lc.amjd)), *coord_to_str(coord[star_index][0], coord[star_index][1])
                )
            else:
                title = "#: {0:4d}  M = {1:.2f}  RA: {5}  DEC: {6}    MJD+{2:d} ({3}/{4})".format(
                    star_index,
                    ommag[star_index],
                    int(min(lc.amjd)),
                    mjd_index + 1,
                    len(mjd_list), *coord_to_str(coord[star_index][0], coord[star_index][1])
                )
            plt.scatter(x, y1, s=4, c="C0")
            plt.scatter(x, y2, s=4, c="C1")
            plt.hlines(
                ommag[star_index], min(x), max(x), linestyles="dashed", colors="red"
            )
            plt.title(title)
            plt.draw()
        else:
            x = lc.mid_time
            y1 = lc.mag
            y2 = lc.magx
            xlabel = "UT"
            title = "#: {0:4d}  M = {1:.2f}  RA: {2}  DEC: {3}".format(star_index, ommag[star_index], *coord_to_str(coord[star_index][0], coord[star_index][1]))
        ylabel = "Magnitude"
        plt.scatter(x, y1, s=4, c="C0")
        plt.scatter(x, y2, s=4, c="C1")
        plt.hlines(ommag[star_index], min(x), max(x), linestyles="dashed", colors="red")
        plt.ylabel(ylabel)
        plt.xlabel(xlabel)
        plt.title(title)
        plt.draw()

    fig, ax = plt.subplots(1, 1)
    x = lc.mid_time
    y1 = lc.mag
    y2 = lc.magx
    xlabel = "UT"
    ylabel = "Magnitude (mag)"
    title = "#: {0:4d}  M = {1:.2f}  RA: {2}  DEC: {3}".format(star_index, ommag[star_index], *coord_to_str(coord[star_index][0], coord[star_index][1]))
    ax.scatter(x, y1, s=4, c="C0")
    ax.scatter(x, y2, s=4, c="C1")
    ax.hlines(ommag[star_index], min(x), max(x), linestyles="dashed", colors="red")
    ax.set_xlabel(xlabel)
    ax.set_ylabel("Magnitude")
    ax.set_title(title)
    fig.canvas.mpl_connect("key_press_event", on_key)
    plt.show()


if __name__ == "__main__":
    cli()
