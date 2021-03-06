import pickle
from pathlib import Path
from operator import itemgetter
from matplotlib import pyplot as plt
from glob import glob
from matplotlib import cm
from mergecat import convert_coord_str_float
import click
import numpy as np


@click.command()
@click.option("-f", "--file_name", type=str, help="Input pkl file name")
@click.option(
    "-i", "--init_star_index", type=int, default=0, help="Index of star to plot"
)
@click.option(
    "--coord",
    type=str,
    help="Coordinate of the query star. Example: '0:0:0.0 0:0:0.0' or '0 0 0.0 0 0 0.0' or '0.0 0.0'",
)
def cli(file_name, init_star_index, coord):
    """Plot light curves"""

    dmatch = 1
    if file_name is None:
        candidate_file_list = glob("*gcat_cal.pkl")
        if len(candidate_file_list) == 1:
            file_name = candidate_file_list[0]
            print("Find {0}".format(file_name))
        elif len(candidate_file_list) == 0:
            print(
                "No *gcat_cal.pkl file is found! Please run command `correctphot` in advance!"
            )
            return
        else:
            print(
                "More than one *gcat_cal.pkl is found! Please specify which file to use by `plotcurve -f FILE_NAME`"
            )
            return
    if not Path(file_name).is_file():
        print("WARNING: File not found!")
        return

    mergecat_dict = pickle.load(open(file_name, "rb"))

    if coord is not None:
        print("Searching coordinate ...")
        c_split = coord.split()
        try:
            if len(c_split) == 2:
                if ":" not in c_split[0] and ":" not in c_split[1]:
                    ra_q = float(c_split[0].strip())
                    dec_q = float(c_split[1].strip())
                else:
                    ra_q, dec_q = convert_coord_str_float(
                        c_split[0].strip().split(":"), c_split[1].strip().split(":")
                    )
            else:
                ra_q, dec_q = convert_coord_str_float(c_split[0:3], c_split[3:6])
        except Exception:
            print("WARNING: Unable to parse the coordinate!")
            return
        else:
            coord = mergecat_dict["coord"]
            sep = np.sqrt((ra_q - coord[:, 0]) ** 2 + (dec_q - coord[:, 1]) ** 2)
            if np.min(sep) < dmatch / 3600:
                init_star_index = sep.argmin()
                print("Find target star with index: {0}".format(init_star_index))
            else:
                print("WARNING: No matching stars within {} arcsec".format(dmatch))
                return

    if "magx" in mergecat_dict.keys():
        plot_lc(file_name, init_star_index)
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


def save_single_phot(
    file_name, magtype, star_index, amjd, coord, frame_info, magmatch, magx
):
    """Save single star's photometry into file

    Args:
        file_name (str): input catalog pkl file name
        magtype (str): magtype input
        star_index (int): index of the star
        amjd (array): AMJD
        coord (array): 2D coords array
        magmatch (array): raw photometry array
        magx (array): corrected photometry array
    """
    ra_str, dec_str = coord_to_str(coord[star_index][0], coord[star_index][1])
    orig_file_name = "{0}.{1}{3}S{2:0>4d}.orig".format(
        file_name.split(".")[0], file_name.split(".")[1][0], star_index, magtype,
    )
    mag_file_name = "{0}.{1}{3}S{2:0>4d}.dat".format(
        file_name.split(".")[0], file_name.split(".")[1][0], star_index, magtype,
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
        for i in range(magmatch.shape[1]):
            f.write(
                "{0:5d} {1:7d} {2:10.5f} {3:10.5f} {4:10.5f} {5}\n".format(
                    i,
                    int(amjd[i]),
                    (amjd[i] - int(amjd[i]) + 0.5) * 24,
                    magmatch[star_index, i, 0],
                    magmatch[star_index, i, 1],
                    frame_info.iloc[i]["band"],
                )
            )
    with open(mag_file_name, "w") as f:
        f.write(
            "{0} {1} {2} (J2000) {3}\n".format(star_index, ra_str, dec_str, magtype)
        )
        for i in range(magmatch.shape[1]):
            f.write(
                "{0:5d} {1:7d} {2:10.5f} {3:10.5f} {4:10.5f} {5}\n".format(
                    i,
                    int(amjd[i]),
                    (amjd[i] - int(amjd[i]) + 0.5) * 24,
                    magx[star_index, i, 0],
                    magx[star_index, i, 1],
                    frame_info.iloc[i]["band"],
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


def handle_margin(mag_surronding, mag_zoom_index):
    if mag_zoom_index >= -2:
        return mag_surronding * 1.4 ** mag_zoom_index
    else:
        return (mag_zoom_index + 2) * mag_surronding / 2


def plot_lc(file_name, init_star_index):
    """Plot light curve

    Keyboard interaction:
    n: next star; next day (day mode)
    N: previous star; previous day (day mode)
    d: switch day mode
    a: switch airmass / ligt curve
    z: zoom out
    Z: zoom in
    r: reset scale
    s: save the data of this star 
    q: quit

    Args:
        file_name (str): pkl file name with magx
    """
    print("Plotting ...")
    # Remove "s", "a", "f", "l" from keymap
    for key in plt.rcParams:
        if key[:6] == "keymap":
            for key_to_remove in ["s", "a", "f", "l"]:
                if key_to_remove in plt.rcParams[key]:
                    plt.rcParams[key].remove(key_to_remove)

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
        nband,
        band_list,
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
        "nband",
        "band_list",
    )(
        mergecat_dict
    )
    if magtype == "a":
        magmatch = apmagmatch
    else:
        magmatch = psfmagmatch

    global star_index
    global mday_flag
    global mjd_index
    global plot_flag
    global plot_flag_old
    global mag_zoom_index
    global fwhm_zoom_index
    star_index = init_star_index
    mjd_index = -1
    mday_flag = False
    plot_flag = "l"
    plot_flag_old = "l"
    mag_surronding = 0.02
    fwhm_surronding = 10
    mag_zoom_index = 0
    fwhm_zoom_index = 0

    plot_label_dict = {"l": "Light curve", "a": "Airmass", "f": "FWHM"}

    marker_list = ["o", "v", "^", "D", "s"]
    color_list = ["C0", "C1", "C2", "C3", "C4"]

    def on_key(event):
        global star_index
        global mday_flag
        global mjd_index
        global plot_flag
        global plot_flag_old
        global mag_zoom_index
        global fwhm_zoom_index
        plt.clf()
        text_x = 0.05
        text_y = 0.02

        if event.key == "n":
            fwhm_zoom_index = 0
            if mday_flag:
                mjd_index += 1
                mjd_index = mjd_index % ndate
                plt.text(
                    text_x,
                    text_y,
                    "Day: ({0}/{1})".format(mjd_index + 1, ndate),
                    transform=plt.gca().transAxes,
                )
            else:
                star_index += 1
                print(
                    "#: {0:4d}  M = {1:5.2f}  RA: {2}  DEC: {3}".format(
                        star_index,
                        ommag[star_index],
                        *coord_to_str(coord[star_index][0], coord[star_index][1]),
                    )
                )
        elif event.key == "N":
            fwhm_zoom_index = 0
            if mday_flag:
                mjd_index -= 1
                mjd_index = mjd_index % ndate
                plt.text(
                    text_x,
                    text_y,
                    "Day: ({0}/{1})".format(mjd_index + 1, ndate),
                    transform=plt.gca().transAxes,
                )
            else:
                star_index -= 1
                print(
                    "#: {0:4d}  M = {1:5.2f}  RA: {2}  DEC: {3}".format(
                        star_index,
                        ommag[star_index],
                        *coord_to_str(coord[star_index][0], coord[star_index][1]),
                    )
                )
        elif event.key == "s":
            save_single_phot(
                file_name,
                magtype,
                star_index,
                frame_info["amjd"],
                coord,
                frame_info,
                magmatch,
                magx,
            )
            return
        elif event.key == "d":
            mday_flag = not mday_flag
            fwhm_zoom_index = 0
            if not mday_flag:
                mjd_index = -1
                plt.text(text_x, text_y, "NORMAL mode", transform=plt.gca().transAxes)
            else:
                plt.text(text_x, text_y, "DAY mode", transform=plt.gca().transAxes)
        elif event.key == "z":
            if plot_flag == "f":
                fwhm_zoom_index += 1
            else:
                mag_zoom_index += 1
            plt.text(text_x, text_y, "Zoom Out", transform=plt.gca().transAxes)
        elif event.key == "Z":
            if plot_flag == "f":
                fwhm_zoom_index -= 1
            else:
                mag_zoom_index -= 1
            plt.text(text_x, text_y, "Zoom In", transform=plt.gca().transAxes)
        elif event.key == "a":
            if plot_flag != "a":
                plot_flag_old = plot_flag
                plot_flag = "a"
            else:
                plot_flag_old, plot_flag = plot_flag, plot_flag_old
            plt.text(
                text_x,
                text_y,
                plot_label_dict[plot_flag],
                transform=plt.gca().transAxes,
            )
        elif event.key == "f":
            if plot_flag != "f":
                plot_flag_old = plot_flag
                plot_flag = "f"
            else:
                plot_flag_old, plot_flag = plot_flag, plot_flag_old
            plt.text(
                text_x,
                text_y,
                plot_label_dict[plot_flag],
                transform=plt.gca().transAxes,
            )
        elif event.key == "l":
            if plot_flag != "l":
                plot_flag_old = plot_flag
                plot_flag = "l"
            else:
                plot_flag_old, plot_flag = plot_flag, plot_flag_old
            plt.text(
                text_x,
                text_y,
                plot_label_dict[plot_flag],
                transform=plt.gca().transAxes,
            )
        elif event.key == "r":
            mag_zoom_index = 0
            fwhm_zoom_index = 0
            plt.text(text_x, text_y, "Reset", transform=plt.gca().transAxes)
        else:
            return

        inner_plot(
            frame_info,
            magmatch,
            magx,
            star_index,
            mday_flag,
            mjd_index,
            plot_flag,
            plot_flag_old,
            mag_zoom_index,
            fwhm_zoom_index,
        )
        plt.draw()

    def inner_plot(
        frame_info,
        magmatch,
        magx,
        star_index,
        mday_flag,
        mjd_index,
        plot_flag,
        plot_flag_old,
        mag_zoom_index,
        fwhm_zoom_index,
    ):
        lc = prepare_lc_df(star_index, frame_info, magmatch, magx)
        if plot_flag == "l":
            if mjd_index != -1:
                lc = lc[lc.mjd == mjd_date_list[mjd_index]]
            xlabel = "MJD"
            ylabel = "magnitude"
            if len(lc) > 0:
                if mday_flag:
                    x = lc.amjd - int(min(lc.amjd))
                else:
                    x = lc.mid_time
                y1 = lc.mag
                y2 = lc.magx
                for i in range(nband):
                    plt.scatter(
                        x[lc.band == band_list[i]],
                        y1[lc.band == band_list[i]],
                        s=2,
                        c=color_list[i],
                        alpha=0.4,
                    )
                    plt.scatter(
                        x[lc.band == band_list[i]],
                        y2[lc.band == band_list[i]],
                        s=4,
                        c=color_list[i],
                        label=band_list[i],
                    )
                    plt.scatter(
                        x[(lc.band == band_list[i]) & (lc.is_bad)],
                        y2[(lc.band == band_list[i]) & (lc.is_bad)],
                        s=20,
                        marker="x",
                        c=color_list[i],
                        label="__nolegend__",
                    )

                plt.legend()
                plt.hlines(
                    ommag[star_index], min(x), max(x), linestyles="dashed", colors="red"
                )
                plt.ylim(
                    np.percentile(y2, 3)
                    - handle_margin(mag_surronding, mag_zoom_index),
                    np.percentile(y2, 97)
                    + handle_margin(mag_surronding, mag_zoom_index),
                )
                plt.gca().invert_yaxis()
        elif plot_flag == "a":
            if mjd_index != -1:
                lc = lc[lc.mjd == mjd_date_list[mjd_index]]
            xlabel = "airmass"
            ylabel = "magnitude"
            if len(lc) > 0:
                x = lc.airmass

                y1 = lc.mag
                y2 = lc.magx
                if not mday_flag:
                    for i in range(nband):
                        plt.scatter(
                            x[lc.band == band_list[i]],
                            y1[lc.band == band_list[i]],
                            s=4,
                            c=color_list[i],
                            label=band_list[i],
                        )
                        plt.scatter(
                            x[(lc.band == band_list[i]) & (lc.is_bad)],
                            y1[(lc.band == band_list[i]) & (lc.is_bad)],
                            s=20,
                            c=color_list[i],
                            marker="x",
                            label="__nolegend__",
                        )
                    plt.legend()
                else:
                    cmap = cm.get_cmap("rainbow")
                    if mjd_index == -1:
                        scatter = plt.scatter(x, y1, s=4, c=lc.mjd, cmap="rainbow")
                        plt.scatter(
                            x[lc.is_bad],
                            y1[lc.is_bad],
                            s=20,
                            marker="x",
                            c=lc[lc.is_bad].mjd,
                            cmap="rainbow",
                        )
                        plt.legend(*scatter.legend_elements())
                    else:
                        for i in range(nband):
                            plt.scatter(
                                x[lc.band == band_list[i]],
                                y1[lc.band == band_list[i]],
                                s=4,
                                c=color_list[i],
                                label=band_list[i],
                            )
                            plt.scatter(
                                x[(lc.band == band_list[i]) & (lc.is_bad)],
                                y1[(lc.band == band_list[i]) & (lc.is_bad)],
                                s=20,
                                c=color_list[i],
                                marker="x",
                                label="__nolegend__",
                            )
                        plt.legend()
                plt.ylim(
                    np.percentile(y2, 3)
                    - handle_margin(mag_surronding, mag_zoom_index),
                    np.percentile(y2, 97)
                    + handle_margin(mag_surronding, mag_zoom_index),
                )
                plt.gca().invert_yaxis()
        elif plot_flag == "f":
            if mjd_index != -1:
                lc = lc[lc.mjd == mjd_date_list[mjd_index]]
            xlabel = "MJD"
            ylabel = "FWHM"
            if len(lc) > 0:
                if mday_flag:
                    x = lc.amjd - int(min(lc.amjd))
                else:
                    x = lc.mid_time
                y1 = lc.fwhm
                for i in range(nband):
                    plt.scatter(
                        x[lc.band == band_list[i]],
                        y1[lc.band == band_list[i]],
                        s=2,
                        c=color_list[i],
                        label=band_list[i],
                    )
                    plt.scatter(
                        x[(lc.band == band_list[i]) & (lc.is_bad)],
                        y1[(lc.band == band_list[i]) & (lc.is_bad)],
                        s=20,
                        marker="x",
                        c=color_list[i],
                        label="__nolegend__",
                    )

                plt.legend()
                plt.ylim(
                    min(y1),
                    max(y1)
                    + handle_margin(fwhm_surronding, fwhm_zoom_index)
                    - fwhm_surronding,
                )
        title = "#: {0:4d}  M = {1:.2f}  RA: {3}  DEC: {4}  MJD+{2:d}".format(
            star_index,
            ommag[star_index],
            int(min(lc.amjd, default=0)),
            *coord_to_str(coord[star_index][0], coord[star_index][1]),
        )
        plt.title(title)
        plt.xlabel(xlabel)
        plt.ylabel(ylabel)

    fig, ax = plt.subplots(1, 1)
    inner_plot(
        frame_info,
        magmatch,
        magx,
        star_index,
        mday_flag,
        mjd_index,
        plot_flag,
        plot_flag_old,
        mag_zoom_index,
        fwhm_zoom_index,
    )
    print(
        "#: {0:4d}  M = {1:5.2f}  RA: {2}  DEC: {3}".format(
            star_index,
            ommag[star_index],
            *coord_to_str(coord[star_index][0], coord[star_index][1]),
        )
    )
    fig.canvas.mpl_connect("key_press_event", on_key)
    plt.show()


if __name__ == "__main__":
    cli()
