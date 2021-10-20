import numpy as np
from astropy.time import Time
from datetime import datetime
import subprocess
import pandas as pd
import click

@click.command()
@click.argument("phot_flag", type=int)
def main(phot_flag):
    pass


def test(file_name):
 # columns = ((0,10),(11,14),(15,18),(19,22),(23,29),(30,35),
 #               (36,44),(44,45),(46,49),(50,55),(56,61),(62,67),
 #               (68,73),(74,79),(79,80),(81,86),(86,87),(88,108))
 #     string=file.readline()
 #     dataline = [ string[c[0]:c[1]] for c in columns ]
    pass

def read_cat_and_info(file_name):
    """Read and sort data and info from cat file

    Sort by apmag2, from the brightest to the faintest

    Args:
        file_name (str): file name

    Returns:
        cat (DataFrame): data frame of photometry
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
    cat = pd.read_table(file_name, delim_whitespace=True, skiprows=3, names=["sn", 'RA', 'DEC', 'x', 'dx', 'y', 'dy', 'psfmag', 'psfmag_err', 'apmag1', "apmag1_err", "apmag2", "apmag2_err", "apmag3", "apmag3_err", "apmag4", "apmag4_err", "ID"])
    ra_s = cat.RA.str.split(":", expand=True)
    rah=pd.to_numeric(ra_s[0])
    ram=pd.to_numeric(ra_s[1])
    ras=pd.to_numeric(ra_s[2])
    ra = 15 * (rah + ram / 60 + ras / 3600)
    cat = cat.assign(ra=ra)
    dec_s = cat.DEC.str.split(":", expand=True)
    decd=pd.to_numeric(dec_s[0])
    decm=pd.to_numeric(dec_s[1]) 
    decs=pd.to_numeric(dec_s[2])
    dec = np.abs(decd) + decm / 60 + decs / 3600
    dec = dec * ((decd == np.abs(decd)) - 0.5) * 2
    cat = cat.assign(dec=dec)
    cat.sort_values("apmag2", inplace=True)

    with open(file_name, "r") as f:
        header_line = f.readline()
        info_line = f.readline()
    info_list = list(filter(None, info_line[:-1].split(' ')))
    start_time = datetime.strptime("{0} {1}".format(file_name[6:12], info_list[0]), "%y%m%d %H:%M:%S.%f")
    mid_time = float(info_list[2])
    exp = float(info_list[1])
    fwhm = float(info_list[3])
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
            "nstar": nstar
            }
    return cat, info_dict

def find_medframe_index(info_dict_list):
    """Find the index of reference frame which has 1.2 times the mean number of stars

    Args:
        catfile_list (list): list of catfile name

    Returns:
        medframe_index: index of medframe in catfile_list
    """
    nc = list()
    for info_dict in info_dict_list:
        nc.append(info_dict["nstar"])
    nfmean = 1.2 * np.sum(nc) / len(nc)
    medframe_index = np.argmin(np.abs(nc - np.sum(nc)/len(nc) * 1.2))
    print("# frames: {0:3d}  Std frame: {1}  # Stars: {2:3d}".format(len(info_dict_list), info_dict_list[medframe_index]['file_name'], info_dict_list[medframe_index]['nstar']))
    return medframe_index




if __name__ == "__main__":
    # main()
    phot_flag = 0
    dmatch = 1.0 # arcsec
    sdev = 0.006

    file_list_byte = subprocess.check_output("ls *.allmag{0}".format(phot_flag), shell=True)
    catfile_list = list(filter(None, file_list_byte.decode("utf8").split("\n")))
    nframe = len(catfile_list)

    # Reading out all individual catalogs into cat_list, info_dict_list, ra_list and dec_list
    cat_list = list()
    info_dict_list = list()
    ra_list = list()
    dec_list = list()
    for k in range(nframe):
        cat, info_dict = read_cat_and_info(catfile_list[k])
        cat = cat.to_numpy()
        cat_list.append(cat)
        info_dict_list.append(info_dict)
        ra_list.append(cat[:, 18].astype(float))
        dec_list.append(cat[:, 19].astype(float))
    medframe_index = find_medframe_index(info_dict_list)

    cat_ref, info_ref_dict = read_cat_and_info(catfile_list[medframe_index])


    apmagmatch = np.zeros((info_ref_dict["nstar"], nframe, 2))
    psfmagmatch = np.zeros((info_ref_dict["nstar"], nframe, 2))

    nomatch = np.zeros(info_ref_dict["nstar"]).astype(int)
    for j in range(info_ref_dict["nstar"]):
        ra0 = ra_list[medframe_index][j]
        dec0 = dec_list[medframe_index][j]
        for k in range(nframe):
            match_flag = False
            cat = cat_list[k]
            info_dict = info_dict_list[k]
            if k != medframe_index:
                sep = np.sqrt((ra0 - ra_list[k]) ** 2 + (dec0 - dec_list[k]) ** 2)
                if np.min(sep) < dmatch / 3600:
                    match_flag = True
                    i = np.argmax(sep < dmatch / 3600)
                    apmagmatch[j, k, 0] = cat[i, 11]
                    apmagmatch[j, k, 1] = cat[i, 12]
                    psfmagmatch[j, k, 0] = cat[i, 7]
                    psfmagmatch[j, k, 1] = cat[i, 8]

            else:
                apmagmatch[j, k, 0] = cat[j, 11]
                apmagmatch[j, k, 1] = cat[j, 12]
                psfmagmatch[j, k, 0] = cat[j, 7]
                psfmagmatch[j, k, 1] = cat[j, 8]
            if match_flag == False:
                nomatch[j] += 1
        
    mjd_list = [info_dict["mjd"] for info_dict in info_dict_list]
    ndate = len(set(mjd_list))
    print("{0} of nights processed!".format(ndate))

    if ndate > 1:
        mergecat = "{0}ALL_{1}.gcat{2}".format(info_dict_list[0]['file_name'][1:6], info_dict_list[0]["file_name"][12:13], phot_flag)
    else:
        mergecat = "{0}.gcat{1}".format(info_dict_list[0]["file_name"][1:13], phot_flag)


    nmlim = max(int(nframe * 0.15), 20) # at most mising in nmlim number of frames
    calib_flag = True
    while calib_flag:
        ic = 1
        istd = list()
        for j in range(info_dict_list[medframe_index]["nstar"]):
            if nomatch[j] > nmlim:
                continue
            if psfmagmatch[j, medframe_index, 0] < 1 and apmagmatch[j, medframe_index, 0] < 1:
                continue
            istd.append(j)
            ic += 1
            if ic > 100:
                ic -= 1
                calib_flag = False
                print("# Std star : {0}".format(ic))
                break
        if ic < 50:
            nmlim *= 2

    with open("std.dat", "w") as f:
        sigm = np.zeros((ic, ic))
        for k1 in range(ic):
            j1 = istd[k1]
            for k2 in range(k1 + 1, ic - 1):
                j2 = istd[k2]
                idm = 0
                sdm = 0
                sig = 0
                dm = np.zeros(nframe)
                for i in range(nframe):
                    if (psfmagmatch[j1, i, 0] * psfmagmatch[j2, i, 0] != 0):
                        dm[i] = psfmagmatch[j1, i, 0] - psfmagmatch[j2, i, 0]
                        sdm = sdm + dm[i]
                        idm += 1
                sdm = sdm / idm
                for i in range(nframe):
                    if dm[i] != 0:
                        sig = sig + (dm[i] - sdm) ** 2
                if sig != 0:
                    sigm[k1, k2] = np.sqrt(sig / idm)

                if sigm[k1, k2] < sdev:
                    f.write("{0:3d} {1:3d} {2:3d} {3:4d} {4:.10f}\n".format(k1, k2, nomatch[j2], idm, sigm[k1, k2]))

    with open("std.dat", "r") as f:
        lines = f.readlines()
        kstd1 = [int(list(filter(None, line.split(" ")))[0]) for line in lines]
        kstd2 = [int(list(filter(None, line.split(" ")))[1]) for line in lines]
        icc = len(kstd1)

    with open("mdev.dat", "w") as f:
        k = 0
        for i in range(ic):
            for j in range(icc):
                if i == kstd2[j]:
                    f.write("{0:10d} {1:10d} {2:10d}\n".format(k, ic, istd[i]))
                    k += 1
                    break
        kcc = k

    with open("stdstar0n.dat", "w") as f:
        for j in range(kcc):
            f.write("{0:15.8f} {1:15.8f} {2:10.5f} {3:10.5f} {4:10.5f} {5:10.5f}\n".format(ra_list[j][medframe_index], dec_list[j][medframe_index], apmagmatch[j, medframe_index, 0], apmagmatch[j, medframe_index, 1], psfmagmatch[j, medframe_index, 0], psfmagmatch[j, medframe_index, 1]))




