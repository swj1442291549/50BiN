# 50BiN
Data pipeline for time-series light curve catalogue for 50BiN

Two tools are provided for analysis and visualization:
- `catmerge`
- `plotcurve`


## Usage
### `catmerge`
This program combines all the frames into a single catalogue, which includes three steps:
1. Find the reference frame
2. Merge the frames
3. Find the non-variable candidate stars for differential photometry

The raw data have four types of files, with different extensions:
- `*.alscat`: PSF photometry files
- `*.ap0`: original aperture photometry files
- `*.allmag0`: the combination of the previous two
- `*.allmag1`: similar to `*.allmag0`, but the aperture photometry is calculated with star center refitted by PSF

`phot_flag` (default: 0) controls which aperture photometry (`apmag`) to be used. It has no interference with PSF photometry (`psfmag`).

A reference frame, which will be used to match other frames, is selected based on the number of stars. Note, the size of the catalogue is determined by that of the reference frame. The number of stars for the reference frame is set to be `medframe_factor` (default: 1.2) times the mean number. 

Once the reference frame is decided, each star in all other frames is matched by its coordinates. This is controlled by `dmatch` (default: 1.0 arcsec), setting the maximum matching distance. 

The non-variable candidates are selected if their variance of magnitude difference (`std(m1-m2)`) is smaller than `sdev` (default: 0.006 mag).

The output files include:
- `std.dat`: magnitude difference of non-variable candidates pair
- `mdev.dat`: non-variable candidates index
- `stdstar0.dat`: photometry of non-variable candidates in the reference frame (columns: `ra`, `dec`, `apmag`, `apmag_err`, `psfmag`, `psfmag_err`)
- `*gcat.pkl`: python pickle file that contains all data
```
Usage: catmerge [OPTIONS]

Options:
  --phot_flag INTEGER      Magnitude type. 0: original aperture photometry; 1:
                           aperture photometry with star center refitted by
                           PSF
  --dmatch FLOAT           Position matching distance in arcsec
  --sdev FLOAT             Standard deviation for none variable selection
  --medframe_factor FLOAT  A factor on average star number in a frame for
                           reference frame selection
  --help                   Show this message and exit.
```

### `plotcurve`
This program is to present the data in time-series light curve. It also applies differential stellar photometry. The details of this method will be discussed in the following post. 


```
Usage: plotcurve [OPTIONS] INPUT_FILE_NAME

  Input catalog file from catmerge (.pkl)

Options:
  --magtype TEXT       magnitude type. a: aperture; p: psf
  --noc INTEGER        Number of selected standard stars
  --plot_flag BOOLEAN  Enter interactive plot
  --help               Show this message and exit.
```

The default figure is shown in UT horizontal axis with multiple-day data folded into a single frame. It can be switched to day mode that expands the data into the MJD axis.
You can interact with the plot through the keyboard:
- `n`: next star; next day (day mode)
- `N`: previous star; previous day (day mode)
- `d`: switch day mode
- `s`: save the data of this star 

The output files include:
- `*.orig`: original photometry of one star (`index`, `mjd`, `ut`, `mag`, `mag_err`)
- `*.dat`: corrected photometry of one star (`index`, `mjd`, `ut`, `magx`, `magx_err`)
- `*gcat_mmag`: average brightness of each star (`index`, `ra`, `dec`, `mmag`, `mmag_err`)
- `*gcat_cal.pkl`: python pickle file that contains all data


## Install
To use this script, `virtualenv` is highly recommended (but not required). Please check [this website](https://virtualenv.pypa.io/en/latest/installation.html) for an installation guidance.

Enter the folder that contains all the source files
```bash
$ virtualenv venv
$ . venv/bin/activate
$ pip install --editable .
```

Afterwards, the program should be available:
```bash
$ catmerge
```
