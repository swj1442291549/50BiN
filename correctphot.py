import pickle
from tqdm import tqdm
import statsmodels.formula.api as smf
import statsmodels.api as sm
from pathlib import Path
import pandas as pd
from glob import glob
from operator import itemgetter
import click
import numpy as np
import warnings


@click.command()
@click.option("-f", "--file_name", type=str, help="Input pkl file name")
@click.option(
    "--magtype", type=str, default="a", help="magnitude type. a: aperture; p: psf"
)
@click.option("--noc", type=int, default=5, help="Number of selected standard stars")
@click.option(
    "--method", type=str, help="Method of correcting photometry. (Empty, 'xy')"
)
def cli(file_name, magtype, noc, method):
    """Calibrate the instrumental photometry"""

    if file_name is None:
        candidate_file_list = glob("*gcat.pkl")
        if len(candidate_file_list) == 1:
            file_name = candidate_file_list[0]
            print("File: {0}".format(file_name))
        elif len(candidate_file_list) == 0:
            print(
                "WARNING: No *gcat.pkl file is found! Please run command `mergecat` in advance!"
            )
            return
        else:
            print(
                "WARNING: More than one *gcat.pkl is found! Please specify which file to use by `correctphot -f FILE_NAME`"
            )
            return
    if not Path(file_name).is_file():
        print("WARNING: File not found!")
        return

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
        nframe_date_list,
        mjd_date_list,
        ncs,
        posmatch,
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
        "nframe_date_list",
        "mjd_date_list",
        "ncs",
        "posmatch",
    )(
        mergecat_dict
    )


    print("# Star: {0:d}".format(nstar))

    # TODO may change how we select standard stars
    if len(ncs) < noc:
        print("too few std stars selected")
        return
    ncs = ncs[:noc]
    print("# Std Stars: {0:d}".format(len(ncs)))

    if magtype == "a":
        magmatch = apmagmatch
    else:
        magmatch = psfmagmatch
    print("Magtype: {0}".format("Aperture" if magtype == "a" else "PSF"))

    if method is not None:
        if noc < 10:
            print(
                "A mininum number of 10 standard stars is required for least-squares fitting"
            )
            return
        smag_ncs = magmatch[ncs, medframe_index, 0]
        magx, ommag, ommag_err = least_square_correct_phot(
            magmatch, nstar, frame_info, ncs, nframe, posmatch, smag_ncs, method
        )
    else:
        magx, ommag, ommag_err = differential_correct_phot(
            magmatch, nstar, frame_info, ncs, medframe_index, nframe
        )

    # bestframe_index_date_list = get_bestframe_index(ndate, magmatch, ncs, nframe_date_list)
    # magx, ommag, ommag_err, frame_info = airmass_correct_phot(magmatch, nstar, frame_info, ncs, nframe, bestframe_index_date_list, nframe_date_list, ndate, mjd_date_list)

    for i in range(noc):
        print(
            "{0:3d} Std star ID: {1:3d}   mmag: {2:8.5f}   mmag_std: {3:8.5f}".format(
                i, ncs[i], np.nanmean(magx[ncs[i], :, 0]), np.nanstd(magx[ncs[i], :, 0])
            )
        )

    # Save average magnitude for each star
    mmag_catfile_name = "{0}.{1}{2}gcat_mmag".format(
        file_name.split(".")[0], file_name.split(".")[1][0], magtype
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
        file_name.split(".")[0], file_name.split(".")[1][0], magtype
    )

    mergecat_dict = {
        "nframe": nframe,
        "medframe_index": medframe_index,
        "nstar": nstar,
        "ndate": ndate,
        "frame_info": frame_info,
        "nomatch": nomatch,
        "coord": coord,
        "psfmagmatch": psfmagmatch,
        "apmagmatch": apmagmatch,
        "magtype": magtype,
        "magx": magx,
        "ommag": ommag,
        "ommag_err": ommag_err,
        "nframe_date_list": nframe_date_list,
        "mjd_date_list": mjd_date_list,
        "ncs": ncs,
    }
    pickle.dump(mergecat_dict, open(final_catfile_name, "wb"))
    print("Save corrected python pickle data in {0}".format(final_catfile_name))


def airmass_correct_phot(
    magmatch,
    nstar,
    frame_info,
    ncs,
    nframe,
    bestframe_index_date_list,
    nframe_date_list,
    ndate,
    mjd_date_list,
):
    ncs_magmatch_delta = np.copy(magmatch[ncs, :, 0])
    magx_delta = np.zeros((nstar, nframe))
    magx = np.copy(magmatch)
    magmatch_best_noairmass = np.zeros((nstar, ndate))
    bad_frame_index_list = list()
    for i in range(ndate):
        ncs_magmatch_delta[
            :, int(sum(nframe_date_list[:i])) : int(sum(nframe_date_list[: i + 1]))
        ] = np.subtract(
            magmatch[
                ncs,
                int(sum(nframe_date_list[:i])) : int(sum(nframe_date_list[: i + 1])),
                0,
            ].T,
            magmatch[ncs, bestframe_index_date_list[i], 0],
        ).T
        ncs_magmatch_delta_mean_date = np.nanmean(
            ncs_magmatch_delta[
                :, int(sum(nframe_date_list[:i])) : int(sum(nframe_date_list[: i + 1]))
            ],
            axis=0,
        )  # ? maybe not mean
        airmass = frame_info[frame_info.mjd == mjd_date_list[i]].airmass
        mag_delta_date = pd.DataFrame(
            {"airmass": airmass, "delta": ncs_magmatch_delta_mean_date}
        )
        popt, perr, bad_frame_index = fit_airmass_delta_zeropoint(mag_delta_date)
        bad_frame_index_list.append(bad_frame_index)

        magmatch_delta = np.subtract(
            magmatch[
                :,
                int(sum(nframe_date_list[:i])) : int(sum(nframe_date_list[: i + 1])),
                0,
            ].T,
            magmatch[:, bestframe_index_date_list[i], 0],
        ).T
        magx_delta[
            :, int(sum(nframe_date_list[:i])) : int(sum(nframe_date_list[: i + 1]))
        ] = np.subtract(magmatch_delta, np.polyval(popt, airmass))
        magmatch_best = magmatch[:, bestframe_index_date_list[i], 0]
        magmatch_best_noairmass[:, i] = (
            magmatch_best
            - popt[0] * frame_info.loc[bestframe_index_date_list[i]].airmass
        )
    ommag = np.nanmean(magmatch_best_noairmass, axis=1)  # TODO improve
    ommag_err = np.nanstd(magmatch_best_noairmass, axis=1)
    magx[:, :, 0] = np.add(magx_delta.T, ommag).T
    frame_info = frame_info.assign(bad=False)
    frame_info.loc[np.concatenate(bad_frame_index_list), "bad"] = True
    magx[:, np.concatenate(bad_frame_index_list), :] = np.nan
    return magx, ommag, ommag_err, frame_info


def differential_correct_phot(magmatch, nstar, frame_info, ncs, medframe_index, nframe):
    """Correct photometry via differential method

    Args:
        magmatch (array): raw photometry array
        nstar (int): number of star
        frame_info (DataFrame): info
        ncs (list): list of standard stars' index
        medframe_index (int): index of medframe
        nframe (int): number of frame

    Returns:
        magx (array): corrected photometry array
        ommag (float): average magnitude
        ommag_err (float): error of average magnitude
    """
    print("Calibrating stars ...")
    magx = np.copy(magmatch)
    ncs_magmatch_delta = np.subtract(
        magmatch[ncs, :, 0].T, magmatch[ncs, medframe_index, 0]
    ).T
    with warnings.catch_warnings():
        warnings.simplefilter("ignore", category=RuntimeWarning)
        ncs_magmatch_delta_mean = np.nanmean(ncs_magmatch_delta, axis=0)
        magx[:, :, 0] = np.subtract(magmatch[:, :, 0], ncs_magmatch_delta_mean)
    ommag, ommag_err = estimate_ommag(magx, nstar)
    return magx, ommag, ommag_err


def estimate_ommag(magx, nstar):
    with warnings.catch_warnings():
        warnings.simplefilter("ignore", category=RuntimeWarning)
        magx_mean = np.nanmean(magx[:, :, 0], axis=1)
        magx_std = np.nanstd(magx[:, :, 0], axis=1, ddof=3)
        magx_delta = np.subtract(magx[:, :, 0].T, magx_mean).T
        magx_abs_delta_ratio = np.divide(np.abs(magx_delta).T, magx_std).T
        ommag = np.zeros(nstar) * np.nan
        ommag_err = np.zeros(nstar) * np.nan
        for i in range(nstar):
            sel = magx_abs_delta_ratio[i] < 3
            ommag[i] = np.nanmean(magx[i, :, 0][sel])
            ommag_err[i] = np.nanmean(magx[i, :, 1][sel])
    return ommag, ommag_err


def fit_airmass_delta(mag_delta_date):
    mag_delta_date = mag_delta_date[
        0.583 * mag_delta_date.airmass - 0.566 > mag_delta_date.delta
    ]
    popt, pcov = np.polyfit(mag_delta_date.airmass, mag_delta_date.delta, 1, cov=True)
    perr = np.sqrt(np.diag(pcov))
    return popt, perr


def fit_airmass_delta_zeropoint(mag_delta_date):
    # TODO improve
    mag_delta_date_1 = mag_delta_date[
        (0.583 * mag_delta_date.airmass - 0.566 > mag_delta_date.delta)
        & (mag_delta_date.airmass > 1.5)
    ]
    k = np.nanmedian(
        mag_delta_date_1.delta
        / (mag_delta_date_1.airmass - min(mag_delta_date.airmass))
    )
    popt = (k, -min(mag_delta_date.airmass) * k)
    perr = (0, 0)
    # mag_delta_date_1 = mag_delta_date[(0.583 * mag_delta_date.airmass -0.566 > mag_delta_date.delta)]
    # mag_delta_date_1 = mag_delta_date_1.append({"airmass": min(mag_delta_date.airmass), "delta": 0}, ignore_index=True)
    # popt, pcov = np.polyfit(mag_delta_date_1.airmass, mag_delta_date_1.delta, 1, cov=True)
    # perr = np.sqrt(np.diag(pcov))
    bad_frame_index = mag_delta_date[
        (mag_delta_date.delta - np.polyval(popt, mag_delta_date.airmass)) > 0.1
    ].index
    return popt, perr, bad_frame_index


def least_square_correct_phot(
    magmatch, nstar, frame_info, ncs, nframe, posmatch, smag_ncs, method
):
    magx = np.copy(magmatch)
    for i in tqdm(range(nframe)):
        magx[:, i, 0] = np.nan
        try:
            dat = pd.DataFrame(
                {
                    "dmag": smag_ncs - magmatch[ncs, i, 0],
                    "x": posmatch[ncs, i, 0],
                    "y": posmatch[ncs, i, 1],
                }
            )
            if method == "xy":
                est = smf.ols("dmag ~ x + y", data=dat).fit()
            else:
                est = smf.ols("dmag ~ 1", data=dat).fit()
        except:
            continue
        else:
            dat = pd.DataFrame(
                {
                    "mag": magmatch[:, i, 0],
                    "x": posmatch[:, i, 0],
                    "y": posmatch[:, i, 1],
                }
            )
            dmag_pred = est.predict(dat)
            magx[:, i, 0] = dmag_pred + magmatch[:, i, 0]
    ommag, ommag_err = estimate_ommag(magx, nstar)
    return magx, ommag, ommag_err


def get_bestframe_index(ndate, magmatch, ncs, nframe_date_list):
    bestframe_index_date_list = list()
    for i in range(ndate):
        ncs_magmatch_mag_date = magmatch[
            ncs, int(sum(nframe_date_list[:i])) : int(sum(nframe_date_list[: i + 1])), 0
        ]
        ncs_magmatch_magshift_date = np.subtract(
            ncs_magmatch_mag_date.T, np.nanmean(ncs_magmatch_mag_date, axis=1)
        ).T
        ncs_magmatch_magmean_date = np.nanmean(ncs_magmatch_magshift_date, axis=0)
        smoothen = 5
        x_ncs_magmatch_magmean_date = np.pad(
            ncs_magmatch_magmean_date,
            (smoothen // 2, smoothen - smoothen // 2),
            mode="edge",
        )
        x_ncs_magmatch_magmean_date = (
            np.cumsum(
                x_ncs_magmatch_magmean_date[smoothen:]
                - x_ncs_magmatch_magmean_date[:-smoothen]
            )
            / smoothen
        )
        bestframe_index_date_list.append(
            x_ncs_magmatch_magmean_date.argmin() + int(sum(nframe_date_list[:i]))
        )
    return bestframe_index_date_list


if __name__ == "__main__":
    cli()
