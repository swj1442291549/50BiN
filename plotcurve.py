import pickle
from pathlib import Path
from operator import itemgetter
from matplotlib import pyplot as plt
import click
import numpy as np


@click.command()
@click.argument("input_file_name", type=str)
@click.option(
    "-i", "--init_star_index", type=int, default=0, help="Index of star to plot"
)
def cli(input_file_name, init_star_index):
    """Input catalog file from catmerge (.pkl)"""
    if not Path(input_file_name).is_file():
        print("File not found!")
        return

    mergecat_dict = pickle.load(open(input_file_name, "rb"))

    if "magx" in mergecat_dict.keys():
        print("This catalog has alreadly contained corrected photometry")
        plot_lc(input_file_name, init_star_index)
    else:
        print(
            "This catalog does not have corrected photometry! Please run command `correctphot` in advance!"
        )
        return


def coord_to_str(ra, dec):
    """Transfer coord into str

    Args:
        ra (float): RA in deg
        dec (float): DEC in deg

    Returns:
        ra_str (str): RA str
        dec_str (str): DEC str
    """
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
    """Save single star's photometry into file

    Args:
        input_file_name (str): input catalog pkl file name
        magtype (str): magtype input
        star_index (int): index of the star
        amjd (array): AMJD
        coord (array): 2D coords array
        magmatch (array): raw photometry array
        magx (array): corrected photometry array
    """
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


def prepare_lc_df(star_index, frame_info, magmatch, magx):
    """Prepare cleaned light curve data

    Add mag, mag_err, magx, and magx_err to info
    Remove nan values or too bright values in magx

    Args:
        star_index (int): index of the star
        frame_info (DataFrame): info data
        magmatch (array): raw photometry array
        magx (array): corrected photometry array

    Returns:
        lc (array): light curve data
    """
    lc = frame_info.copy()
    lc = lc.assign(mag=magmatch[star_index, :, 0])
    lc = lc.assign(mag_err=magmatch[star_index, :, 1])
    lc = lc.assign(magx=magx[star_index, :, 0])
    lc = lc.assign(magx_err=magx[star_index, :, 1])
    lc = lc[~np.isnan(lc.magx) & (lc.magx > 1)]
    return lc


def plot_lc(file_name, init_star_index):
    """Plot light curve

    Keyboard interaction:
    n: next star; next day (day mode)
    N: previous star; previous day (day mode)
    d: switch day mode
    s: save the data of this star 

    Args:
        file_name (str): pkl file name with magx
    """
    print("Plotting ...")
    if "s" in plt.rcParams["keymap.save"]:
        plt.rcParams["keymap.save"].remove("s")
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
        nframe_date_list,
        mjd_date_list,
        ncs,
    ) = itemgetter(
        "nframe",
        "medframe_index",
        "nstar",
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
        "nframe_date_list",
        "mjd_date_list",
        "ncs",
    )(
        mergecat_dict
    )
    if magtype == "a":
        magmatch = apmagmatch
    else:
        magmatch = psfmagmatch

    mjd_list = list(set(frame_info.mjd))

    global star_index
    global mday_flag
    global mjd_index
    star_index = init_star_index
    mjd_index = -1
    mday_flag = False
    lc = prepare_lc_df(star_index, frame_info, magmatch, magx)

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
                file_name,
                magtype,
                star_index,
                frame_info["amjd"],
                coord,
                magmatch,
                magx,
            )
            return
        elif event.key == "d":
            mday_flag = not mday_flag
            if not mday_flag:
                mjd_index = -1

        plt.clf()
        lc = prepare_lc_df(star_index, frame_info, magmatch, magx)
        if mjd_index != -1:
            lc = lc[lc.mjd == mjd_list[mjd_index]]
        if mday_flag:
            x = lc.amjd - int(min(lc.amjd))
            y1 = lc.mag
            y2 = lc.magx
            xlabel = "MJD"
            if mjd_index == -1:
                title = "#: {0:4d}  M = {1:.2f}  RA: {3}  DEC: {4}  MJD+{2:d}".format(
                    star_index,
                    ommag[star_index],
                    int(min(lc.amjd)),
                    *coord_to_str(coord[star_index][0], coord[star_index][1])
                )
            else:
                title = "#: {0:4d}  M = {1:.2f}  RA: {5}  DEC: {6}    MJD+{2:d} ({3}/{4})".format(
                    star_index,
                    ommag[star_index],
                    int(min(lc.amjd)),
                    mjd_index + 1,
                    len(mjd_list),
                    *coord_to_str(coord[star_index][0], coord[star_index][1])
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
            title = "#: {0:4d}  M = {1:.2f}  RA: {2}  DEC: {3}".format(
                star_index,
                ommag[star_index],
                *coord_to_str(coord[star_index][0], coord[star_index][1])
            )
        ylabel = "Magnitude"
        plt.scatter(x, y1, s=4, c="C0")
        plt.scatter(x, y2, s=4, c="C1")
        plt.hlines(ommag[star_index], min(x), max(x), linestyles="dashed", colors="red")
        plt.ylabel(ylabel)
        plt.xlabel(xlabel)
        plt.title(title)
        print(
            "#: {0:4d}  M = {1:.2f}  RA: {2}  DEC: {3}   ##: {4:d}".format(
                star_index,
                ommag[star_index],
                *coord_to_str(coord[star_index][0], coord[star_index][1]),
                len(lc)
            )
        )
        plt.ylim(min(y2) - 0.05, max(y2) + 0.05)
        plt.gca().invert_yaxis()
        plt.draw()

    fig, ax = plt.subplots(1, 1)
    x = lc.mid_time
    y1 = lc.mag
    y2 = lc.magx
    xlabel = "UT"
    ylabel = "Magnitude (mag)"
    title = "#: {0:4d}  M = {1:.2f}  RA: {2}  DEC: {3}".format(
        star_index,
        ommag[star_index],
        *coord_to_str(coord[star_index][0], coord[star_index][1])
    )
    print(
        "#: {0:4d}  M = {1:.2f}  RA: {2}  DEC: {3}   ##: {4:d}".format(
            star_index,
            ommag[star_index],
            *coord_to_str(coord[star_index][0], coord[star_index][1]),
            len(lc)
        )
    )
    ax.scatter(x, y1, s=4, c="C0")
    ax.scatter(x, y2, s=4, c="C1")
    ax.hlines(ommag[star_index], min(x), max(x), linestyles="dashed", colors="red")
    ax.set_xlabel(xlabel)
    ax.set_ylabel("Magnitude")
    ax.set_title(title)
    ax.set_ylim(min(y2) - 0.05, max(y2) + 0.05)
    ax.invert_yaxis()
    fig.canvas.mpl_connect("key_press_event", on_key)
    plt.show()


if __name__ == "__main__":
    cli()
