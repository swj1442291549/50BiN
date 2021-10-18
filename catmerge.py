from pathlib import Path
import numpy as np
import jdcal
from astropy.time import Time
from datetime import datetime
import pandas as pd
import click
from util import wc

@click.command()
@click.argument("phot_flag", type=int)
def main(phot_flag):
    pass

def read_file_info(file_name):
    """Read info from file

    from file name and headers

    Args:
        file_name (str): file name

    Returns:
        info (dict):
            mjd (int): modified Julian date (MJD = JD - 2440000) 
            start_time (datetime): start time
            mid_time (float): UT of observation in hour
            exp (float): exposure time in second
            fwhm (float): fwhm
            aperture (float): aperture in pix
    """
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

    info = {
            "mjd": mjd,
            "start_time": start_time,
            "mid_time": mid_time,
            "exp": exp,
            "fwhm": fwhm,
            "aperture": aperture,
            }
    return info

def read_cat(file_name):
    """Read and sort data from cat file

    Sort by PSF mag, from the brightest to the faintest

    Args:
        file_name (str): file name

    Returns:
        df (DataFrame): data frame of photometry
    """
    df = pd.read_table(file_name, delim_whitespace=True, skiprows=3, index_col=0, names=['RA', 'DEC', 'x', 'dx', 'y', 'dy', 'psfmag', 'psfmag_err', 'apmag1', "apmag1_err", "agmag2", "agmag2_err", "agmag3", "agmag3_err", "agmag4", "apmag4_err", "ID"])
    ra_s = df.RA.str.split(":", expand=True)
    df = df.assign(rah=pd.to_numeric(ra_s[0]))
    df = df.assign(ram=pd.to_numeric(ra_s[1]))
    df = df.assign(ras=pd.to_numeric(ra_s[2]))
    ra = 15 * (3600 * df.rah + 60 * df.ram + df.ras)
    df = df.assign(ra=ra)
    dec_s = df.DEC.str.split(":", expand=True)
    df = df.assign(decd=pd.to_numeric(dec_s[0]))
    df = df.assign(decm=pd.to_numeric(dec_s[1]))
    df = df.assign(decs=pd.to_numeric(dec_s[2]))
    dec = 3600 * df.decd + 60 * df.decm + df.decs
    df = df.assign(dec=dec)
    df.sort_values("psfmag", inplace=True)
    return df

def find_medframe_index(catfile_list):
    """Find the index of reference frame which has 1.2 times the mean number of stars

    Use wc to count lines in each catfile
    To account for the headers (2) and column name (1)lines, the actual star number is wc(file_name) - 3

    Args:
        catfile_list (list): list of catfile name

    Returns:
        medframe_index: index of medframe in catfile_list
    """
    nc = list()
    for catfile in catfile_list:
        nc.append(wc(catfile))
    nfmean = 1.2 * np.sum(nc) / len(nc)
    medframe_index = np.argmin(np.abs(nc - np.sum(nc)/len(nc) * 1.2))
    print("# frames: {0:d}  Std frame: {1}  # Stars: {2:d}".format(len(catfile_list), catfile_list[medframe_index], wc(catfile_list[medframe_index]) - 3))
    return medframe_index

    

if __name__ == "__main__":
    # main()
    phot_flag = 0

    cat_paths = Path("./").glob("*.allmag{0}".format(phot_flag))
    catfile_list = [p.name for p in cat_paths]
    medframe_index = find_medframe_index(catfile_list)

    df = read_cat(catfile_list[medframe_index])
    info = read_file_info(catfile_list[medframe_index])
