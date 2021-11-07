import numpy as np
import pickle
from astropy.time import Time
from datetime import datetime
from glob import glob
import pandas as pd
import click
from astropy.coordinates import SkyCoord, EarthLocation, AltAz
import astropy.units as u
from tqdm import tqdm


@click.command()
@click.option(
    "--phot_flag",
    type=int,
    default=0,
    help="Magnitude type. 0: original aperture photometry; 1: aperture photometry with star center refitted by PSF",
)
@click.option(
    "--dmatch", type=float, default=1.0, help="Position matching distance in arcsec"
)
@click.option(
    "--sdev",
    type=float,
    default=0.006,
    help="Standard deviation for none variable selection",
)
@click.option(
    "--noc",
    type=int,
    help="Minimum number of standard candidates. Not the same as Std Stars output",
)
@click.option(
    "--medframe_factor",
    type=float,
    default=1.2,
    help="A factor on average star number in a frame for reference frame selection",
)
@click.option(
    "--obs_flag",
    type=str,
    default="d",
    help="Observatory flag. 'd': Delingha; 'l': Lenghu",
)
@click.option(
    "-b", "--band", type=str, help="Passband",
)
def cli(phot_flag, dmatch, sdev, medframe_factor, obs_flag, band, noc):
    """Merge the catalogs"""
    if band is None:
        catfile_list = glob("*.allmag{0}".format(phot_flag))
    else:
        catfile_list = glob(
            "*{1}[0-9][0-9][0-9][0-9].allmag{0}".format(phot_flag, band)
        )
    catfile_list.sort()

    # check there is only one filter data
    filter_list = [catfile[-13] for catfile in catfile_list]
    nfilter = len(set(filter_list))
    if nfilter != 1:
        print(
            "The folder contains data of {0} bands: {1}! Please specify one band!".format(
                nfilter, ", ".join(set(filter_list))
            )
        )
        return
    band = filter_list[0]

    # Reading out all individual catalogs into cat_list, info_dict_list, coord_list
    cat_list = list()
    info_dict_list = list()
    coord_list = list()
    print("Reading data ... ")
    for k in tqdm(range(len(catfile_list))):
        cat, info_dict = read_cat_and_info(catfile_list[k])
        if cat is not None:
            cat_list.append(cat)
            info_dict_list.append(info_dict)
            coord_list.append(cat[:, 18:20].astype(float))
    frame_info = pd.DataFrame(info_dict_list)
    nframe = len(frame_info)
    mjd_date_list = np.sort(list(set(frame_info.mjd)))
    nframe_date_list = [
        len(frame_info[frame_info.mjd == mjd_date]) for mjd_date in mjd_date_list
    ]
    ndate = len(mjd_date_list)
    print("Read {0:d} {2} frames of {1:d} nights".format(nframe, ndate, band))

    # Calculate airmass
    mountain = read_obs_location(obs_flag)
    time = Time(frame_info["start_time"])  # should use mid time
    target = SkyCoord(
        np.mean(coord_list[0][:, 0]), np.mean(coord_list[0][:, 1]), unit="deg",
    )
    target_altaz = target.transform_to(AltAz(obstime=time, location=mountain))
    target_airmass = target_altaz.secz
    frame_info = frame_info.assign(airmass=target_airmass)

    medframe_index = find_medframe_index(frame_info, medframe_factor)
    nstar = frame_info.loc[medframe_index]["nstar"]

    # Merge the catalogs
    # Use medframe as a reference, looking for each stars in all other frames by matching coordinates
    apmagmatch = np.zeros((nstar, nframe, 2)) * np.nan
    psfmagmatch = np.zeros((nstar, nframe, 2)) * np.nan
    posmatch = np.zeros((nstar, nframe, 2)) * np.nan
    nomatch = np.zeros(nstar).astype(int)
    print("Matching stars ... ")
    for j in tqdm(range(nstar)):
        ra0, dec0 = coord_list[medframe_index][j]
        for k in range(nframe):
            match_flag = False
            if k != medframe_index:
                sep = np.sqrt(
                    (ra0 - coord_list[k][:, 0]) ** 2 + (dec0 - coord_list[k][:, 1]) ** 2
                )
                if np.min(sep) < dmatch / 3600:
                    match_flag = True
                    cat = cat_list[k]
                    i = np.argmax(sep < dmatch / 3600)
                    apmagmatch[j, k, :] = cat[i, 11:13]
                    psfmagmatch[j, k, :] = cat[i, 7:9]
                    posmatch[j, k, 0] = cat[i, 3]
                    posmatch[j, k, 1] = cat[i, 5]

            else:
                cat = cat_list[k]
                apmagmatch[j, k, 0] = cat[j, 11]
                apmagmatch[j, k, 1] = cat[j, 12]
                psfmagmatch[j, k, 0] = cat[j, 7]
                psfmagmatch[j, k, 1] = cat[j, 8]
                posmatch[j, k, 0] = cat[j, 3]
                posmatch[j, k, 1] = cat[j, 5]
            if not match_flag:
                nomatch[j] += 1

    # Define standard candidate stars for differential photometry
    std = np.arange(nstar)
    nmlim = max(int(nframe * 0.15), 20)  # at most mising in nmlim number of frames
    calib_flag = True
    istd = list()
    while calib_flag:
        istd = std[(nomatch < nmlim)]
        if len(istd) > int(nstar * 0.2):
            calib_flag = False
        else:
            nmlim *= 1.1
        if nmlim >= nstar:
            break
    ic = len(istd)

    # Find non-variable candidate stars for differential photometry
    sigm = np.zeros((ic, ic)) * np.nan
    for k1 in range(ic):
        j1 = istd[k1]
        for k2 in range(k1 + 1, ic - 1):
            j2 = istd[k2]
            m1 = psfmagmatch[j1, :, 0]
            m2 = psfmagmatch[j2, :, 0]
            dm = m1 - m2  # Magnitude difference between j1 and j2
            idm = len(dm[~np.isnan(dm)])  # Number of frame with non-nan records
            sdm = np.nanmean(dm)  # Average magnitude difference
            sig = np.nansum((dm - sdm) ** 2 * np.abs(np.sign(dm)))
            sigm[k1, k2] = np.sqrt(sig / idm) * np.sign(sig)

    if noc is not None:
        sigm_flat = sigm.reshape(-1)
        if len(sigm_flat[sigm_flat < sdev]) < noc:
            print(
                "Less than {1} standard candidates are selected with sdev: {0:.3f}!".format(
                    sdev, noc
                )
            )
            sdev = np.nanpercentile(
                sigm_flat, noc * 100 / len(sigm_flat[~np.isnan(sigm_flat)])
            )
            print("Change sdev to {0:.3f}".format(sdev))
    kstd1 = list()
    kstd2 = list()
    for k1 in range(ic):
        for k2 in range(k1 + 1, ic - 1):
            if sigm[k1, k2] < sdev:
                kstd1.append(k1)
                kstd2.append(k2)

    icc = len(kstd1)

    ncs = list()
    for i in range(ic):
        for j in range(icc):
            if i == kstd2[j]:
                ncs.append(istd[i])
                break
    print("# Std Stars: {0:d}".format(len(ncs)))
    if len(ncs) < 10:
        print(
            "The number of standard stars are too small! To ensure the accuracy of `correctphot`, please consider use `--noc` option"
        )

    with open("stdstar.dat", "w") as f:
        f.write(
            "             ra             dec      apmag  apmag_err     psfmag psfmag_err\n"
        )
        for j in range(len(ncs)):
            f.write(
                "{0:15.8f} {1:15.8f} {2:10.5f} {3:10.5f} {4:10.5f} {5:10.5f}\n".format(
                    coord_list[medframe_index][j, 0],
                    coord_list[medframe_index][j, 1],
                    apmagmatch[j, medframe_index, 0],
                    apmagmatch[j, medframe_index, 1],
                    psfmagmatch[j, medframe_index, 0],
                    psfmagmatch[j, medframe_index, 1],
                )
            )
    print("Save standard stars info in {0}".format("stdstar.dat"))

    # Write merged uncalibrated data into a file
    if ndate > 1:
        mergecat_file_name = "{0}ALL_{1}.{2}gcat.pkl".format(
            info_dict_list[0]["file_name"][1:6],
            info_dict_list[0]["file_name"][12:13],
            phot_flag,
        )
    else:
        mergecat_file_name = "{0}.{1}gcat.pkl".format(
            info_dict_list[0]["file_name"][1:13], phot_flag
        )

    mergecat_dict = {
        "nframe": nframe,  # number of frames
        "medframe_index": medframe_index,  # index of the reference frame (used to match stars)
        "nstar": info_dict_list[medframe_index][
            "nstar"
        ],  # number of stars in the reference frame
        "ndate": ndate,  # number of MJD dates
        "frame_info": frame_info,  # frame into data
        "nomatch": nomatch,  # number of non-match frames
        "coord": coord_list[medframe_index],  # coordinates
        "psfmagmatch": psfmagmatch,  # PSF magnitude
        "apmagmatch": apmagmatch,  # Aperature magnitude
        "nframe_date_list": nframe_date_list,  # number of frames in each date
        "mjd_date_list": mjd_date_list,  # MJD of each date
        "ncs": ncs,  # index of standard stars
        "posmatch": posmatch,  # pos array, same format as magmatch (X, Y)
    }
    pickle.dump(mergecat_dict, open(mergecat_file_name, "wb"))
    print("Save python pickle data in {0}".format(mergecat_file_name))


def read_obs_location(obs_flag):
    """Read Earthlocation for different observatory

    Args:
        obs_flag (str): observatory flag. "d" for Delingha; "l" for Lenghu

    Returns:
        mountain (EarthLocation): cite location on the Earth
    """
    if obs_flag == "d":
        mountain = EarthLocation(
            lat=37.373 * u.deg, lon=97.56 * u.deg, height=3200 * u.m
        )
        return mountain
    elif obs_flag == "l":
        mountain = EarthLocation(
            lat=38.6068 * u.deg, lon=93.8961 * u.deg, height=4200 * u.m
        )
        return mountain


def read_cat_and_info(file_name):
    """Read and sort data and info from cat file

    Sort by apmag2, from the brightest to the faintest

    For those with apmag2 > 30 or apmag2_err > 3 or psfmag > 30 or psfmag_err > 3, change the magnitude to NaN

    Args:
        file_name (str): file name

    Returns:
        cat (array): complete photometry array
        info_dict (dict):
            file_name (str): file_name
            mjd (int): modified Julian date (MJD = JD - 2440000) 
            start_time (datetime): start time
            mid_time (float): UT of observation in hour
            exp (float): exposure time in second
            fwhm (float): fwhm
            aperture (float): aperture in pix
            nstar (int): number of stars in the frame
    """
    cat = pd.read_table(
        file_name,
        delim_whitespace=True,
        skiprows=3,
        names=[
            "sn",
            "RA",
            "DEC",
            "x",
            "dx",
            "y",
            "dy",
            "psfmag",
            "psfmag_err",
            "apmag1",
            "apmag1_err",
            "apmag2",
            "apmag2_err",
            "apmag3",
            "apmag3_err",
            "apmag4",
            "apmag4_err",
            "ID",
        ],
    )
    if len(cat) == 0:
        return None, None
    ra_s = cat.RA.str.split(":", expand=True)
    rah = pd.to_numeric(ra_s[0])
    ram = pd.to_numeric(ra_s[1])
    ras = pd.to_numeric(ra_s[2])
    ra = 15 * (rah + ram / 60 + ras / 3600)
    cat = cat.assign(ra=ra)
    dec_s = cat.DEC.str.split(":", expand=True)
    decd = pd.to_numeric(dec_s[0])
    decm = pd.to_numeric(dec_s[1])
    decs = pd.to_numeric(dec_s[2])
    dec = np.abs(decd) + decm / 60 + decs / 3600
    dec = dec * ((decd == np.abs(decd)) - 0.5) * 2
    cat = cat.assign(dec=dec)
    cat.sort_values("apmag2", inplace=True)
    cat = cat.to_numpy()
    ap_nan_index = (cat[:, 11] > 30) | (cat[:, 12] > 3.0)
    psf_nan_index = (cat[:, 7] > 30) | (cat[:, 8] > 3.0)
    cat[ap_nan_index, 11:13] = np.nan
    cat[psf_nan_index, 7:9] = np.nan

    with open(file_name, "r") as f:
        f.readline()
        info_line = f.readline()
    info_list = list(filter(None, info_line[:-1].split(" ")))
    start_time = datetime.strptime(
        "{0} {1}".format(file_name[6:12], info_list[0]), "%y%m%d %H:%M:%S.%f"
    )
    mid_time = convert_str_to_float(info_list[2])
    exp = convert_str_to_float(info_list[1])
    fwhm = convert_str_to_float(info_list[3])
    aperture = float(info_list[4])
    dt = datetime.strptime(file_name[6:12], "%y%m%d")
    t = Time(dt)
    mjd = int(t.mjd + 1 - 40000)
    nstar = len(cat)

    info_dict = {
        "file_name": file_name,
        "mjd": mjd,
        "start_time": start_time,
        "mid_time": mid_time,
        "exp": exp,
        "fwhm": fwhm,
        "aperture": aperture,
        "nstar": nstar,
    }
    return cat, info_dict


def convert_str_to_float(string):
    """Convert str to float

    To handle the edge case

    Args:
        string (str): string

    Returns:
        f (float): float value
    """
    try:
        f = float(string)
    except:
        f = np.nan
    return f


def find_medframe_index(frame_info, medframe_factor):
    """Find the index of reference frame which has medframe_factor times the mean number of stars

    Args:
        frame_info (DataFrame): info
        medframe_factor (float): number ratio

    Returns:
        medframe_index: index of medframe in catfile_list
    """
    ns = frame_info.nstar
    medframe_index = np.abs(ns - np.sum(ns) / len(ns) * medframe_factor).idxmin()
    print(
        "Reference frame: {0}  # Stars: {1:3d}".format(
            frame_info.loc[medframe_index]["file_name"],
            frame_info.loc[medframe_index]["nstar"],
        )
    )
    return medframe_index


def find_medframe_index_airmass(frame_info):
    """Find the index of reference frame which has the least airmass

    Args:
        frame_info (DataFrame): info

    Returns:
        medframe_index: index of medframe in catfile_list
    """
    medframe_index = frame_info["airmass"].idxmin()
    print(
        "Reference frame: {0}  # Stars: {1:3d}  airmass: {2:.2f}".format(
            frame_info.loc[medframe_index]["file_name"],
            frame_info.loc[medframe_index]["nstar"],
            frame_info.loc[medframe_index]["airmass"],
        )
    )
    return medframe_index


if __name__ == "__main__":
    cli()
