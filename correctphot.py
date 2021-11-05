import pickle
from pathlib import Path
import pandas as pd
from operator import itemgetter
from matplotlib import pyplot as plt
import click
import numpy as np
from tqdm import tqdm


@click.command()
@click.argument("input_file_name", type=str)
@click.option(
    "--magtype", type=str, default="a", help="magnitude type. a: aperture; p: psf"
)
@click.option("--noc", type=int, default=10, help="Number of selected standard stars")
def cli(input_file_name, magtype, noc):
    """Input catalog file from catmerge (.pkl)"""
    if not Path(input_file_name).is_file():
        print("File not found!")
        return

    mergecat_dict = pickle.load(open(input_file_name, "rb"))

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
        ncs
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
        "ncs"
    )(
        mergecat_dict
    )

    frame_info = frame_info.assign(
        amjd=frame_info.mjd + frame_info.mid_time / 24 - 0.5
    )  # convert observing time to modified julian day AMJD(1-nf)

    print("total number of stars: {0:d}".format(nstar))

    if len(ncs) < noc:
        print("too few std stars selected")
        return
    else:
        ncs = ncs[:noc]
        print("Std tot: {0:d}".format(noc))
        for i in range(noc):
            print("{0:3d} Std star ID: {1:4d}".format(i, ncs[i]))


    if magtype == "a":
        magmatch = apmagmatch
    else:
        magmatch = psfmagmatch

    magx, ommag, ommag_err = differential_correct_phot(
        magmatch, nstar, frame_info, ncs, medframe_index, nframe
    )

    # bestframe_index_date_list = get_bestframe_index(ndate, magmatch, ncs, nframe_date_list)
    # magx, ommag, ommag_err, frame_info = airmass_correct_phot(magmatch, nstar, frame_info, ncs, medframe_index, nframe, bestframe_index_date_list, nframe_date_list, ndate, mjd_date_list)

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
        "ncs": ncs
    }
    pickle.dump(mergecat_dict, open(final_catfile_name, "wb"))
    print("Save corrected python pickle data in {0}".format(final_catfile_name))

def airmass_correct_phot(magmatch, nstar, frame_info, ncs, medframe_index, nframe, bestframe_index_date_list, nframe_date_list, ndate, mjd_date_list):
    ncs_magmatch_delta = np.copy(magmatch[ncs, :, 0])
    magx_delta = np.zeros((nstar, nframe))
    magx = np.copy(magmatch)
    magmatch_best_noairmass = np.zeros((nstar, ndate))
    bad_frame_index_list = list()
    for i in range(ndate):
        ncs_magmatch_delta[:, int(sum(nframe_date_list[:i])): int(sum(nframe_date_list[:i+1]))] = np.subtract(magmatch[ncs, int(sum(nframe_date_list[:i])): int(sum(nframe_date_list[:i+1])), 0].T, magmatch[ncs, bestframe_index_date_list[i], 0]).T
        ncs_magmatch_delta_mean_date = np.nanmean(ncs_magmatch_delta[:, int(sum(nframe_date_list[:i])): int(sum(nframe_date_list[:i+1]))], axis=0) #? maybe not mean
        airmass = frame_info[frame_info.mjd == mjd_date_list[i]].airmass
        mag_delta_date = pd.DataFrame({'airmass': airmass, "delta": ncs_magmatch_delta_mean_date})
        popt, perr, bad_frame_index = fit_airmass_delta_zeropoint(mag_delta_date)
        bad_frame_index_list.append(bad_frame_index)

        magmatch_delta = np.subtract(magmatch[:, int(sum(nframe_date_list[:i])): int(sum(nframe_date_list[:i+1])), 0].T, magmatch[:, bestframe_index_date_list[i], 0]).T
        magx_delta[:, int(sum(nframe_date_list[:i])): int(sum(nframe_date_list[:i+1]))] = np.subtract(magmatch_delta, np.polyval(popt, airmass))
        magmatch_best = magmatch[:, bestframe_index_date_list[i], 0]
        magmatch_best_noairmass[:, i] = magmatch_best - popt[0] * frame_info.loc[bestframe_index_date_list[i]].airmass
    ommag = np.nanmean(magmatch_best_noairmass, axis=1) # TODO improve
    ommag_err = np.nanstd(magmatch_best_noairmass, axis=1)
    magx[:, :, 0] = np.add(magx_delta.T, ommag).T
    frame_info = frame_info.assign(bad=False)
    frame_info.loc[np.concatenate(bad_frame_index_list), "bad"] = True
    magx[:, np.concatenate(bad_frame_index_list), :] = np.nan
        # x = np.linspace(1.2, 2.2)
        # ax.scatter(mag_delta_date.airmass, mag_delta_date.delta, s=4)
        # ax.plot(x, np.polyval(popt, x), label='date {0}'.format(i + 1))
    # ax.set_xlabel('x')
    # ax.set_ylabel('y')
    # ax.legend()
    # plt.show()
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
    magx = np.zeros_like(magmatch)
    ommag = np.zeros(nstar)
    ommag_err = np.zeros(nstar)
    ut = frame_info.mid_time.values
    magmatch_medframe = magmatch[ncs, medframe_index, 0]
    for ipg in range(nstar):
        magmatch[ipg, magmatch[ipg, :, 0] > 30, 1] = 0
        magmatch[ipg, magmatch[ipg, :, 0] > 30, 0] = 0

    print("Calibrating stars ...")
    for ipg in tqdm(range(nstar)):
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
    return magx, ommag, ommag_err


def fit_airmass_delta(mag_delta_date):
    mag_delta_date = mag_delta_date[0.583 * mag_delta_date.airmass -0.566 > mag_delta_date.delta]
    popt, pcov = np.polyfit(mag_delta_date.airmass, mag_delta_date.delta, 1, cov=True)
    perr = np.sqrt(np.diag(pcov))
    return popt, perr

def fit_airmass_delta_zeropoint(mag_delta_date):
    # TODO improve
    mag_delta_date_1 = mag_delta_date[(0.583 * mag_delta_date.airmass -0.566 > mag_delta_date.delta) & (mag_delta_date.airmass > 1.5)]
    k = np.nanmedian(mag_delta_date_1.delta / (mag_delta_date_1.airmass - min(mag_delta_date.airmass)))
    popt = (k, -min(mag_delta_date.airmass) * k)
    perr = (0, 0)
    # mag_delta_date_1 = mag_delta_date[(0.583 * mag_delta_date.airmass -0.566 > mag_delta_date.delta)]
    # mag_delta_date_1 = mag_delta_date_1.append({"airmass": min(mag_delta_date.airmass), "delta": 0}, ignore_index=True)
    # popt, pcov = np.polyfit(mag_delta_date_1.airmass, mag_delta_date_1.delta, 1, cov=True)
    # perr = np.sqrt(np.diag(pcov))
    bad_frame_index = mag_delta_date[(mag_delta_date.delta - np.polyval(popt, mag_delta_date.airmass)) > 0.1].index
    return popt, perr, bad_frame_index

def get_bestframe_index(ndate, magmatch, ncs, nframe_date_list):
    bestframe_index_date_list = list()
    for i in range(ndate):
        ncs_magmatch_mag_date = magmatch[ncs, int(sum(nframe_date_list[:i])): int(sum(nframe_date_list[:i+1])), 0]
        ncs_magmatch_magshift_date = np.subtract(ncs_magmatch_mag_date.T, np.nanmean(ncs_magmatch_mag_date, axis=1)).T
        ncs_magmatch_magmean_date = np.nanmean(ncs_magmatch_magshift_date, axis=0)
        smoothen = 5
        x_ncs_magmatch_magmean_date = np.pad(ncs_magmatch_magmean_date, (smoothen//2, smoothen-smoothen//2), mode='edge')
        x_ncs_magmatch_magmean_date = np.cumsum(x_ncs_magmatch_magmean_date[smoothen:] - x_ncs_magmatch_magmean_date[:-smoothen]) / smoothen
        bestframe_index_date_list.append(x_ncs_magmatch_magmean_date.argmin() + int(sum(nframe_date_list[:i])))
    return bestframe_index_date_list


if __name__ == "__main__":
    cli()
