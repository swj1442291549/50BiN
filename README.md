# 50BiN
Data pipeline for time-series light curve catalogue for 50BiN

Three tools are provided for analysis and visualization:
1. `mergecat`
2. `correctphot`
3. `plotcurve`


## Install
To use this script, `virtualenv` is highly recommended (but not required). Please check [this website](https://virtualenv.pypa.io/en/latest/installation.html) for an installation guidance.

```bash
$ virtualenv venv
$ . venv/bin/activate
$ pip install -U git+git://github.com/swj1442291549/50BiN
```

Afterwards, the program should be available:
```bash
$ mergecat
```

## Usage
### `mergecat`
This program combines all the frames into a single catalogue, which includes three steps:
1. Find the reference frame
2. Merge the frames
3. Find the non-variable candidate stars for differential photometry

The raw data have four types of files, with different extensions:
- `*.alscat`: PSF photometry files
- `*.ap0`: original aperture photometry files
- `*.allmag0`: the combination of the previous two
- `*.allmag1`: similar to `*.allmag0`, but the aperture photometry is calculated with star center refitted by PSF

`phot_flag` (default: 0) controls which aperture photometry (`apmag`) to be used. 0: original aperture photometry; 1: aperture photometry with star center refitted by PSF. It has no interference with PSF photometry (`psfmag`).

A reference frame, which will be used to match other frames, is selected based on the number of stars. Note, the size of the catalogue is determined by that of the reference frame. The number of stars for the reference frame is set to be `medframe_factor` (default: 1.2) times the mean number. 

Once the reference frame is decided, each star in all other frames is matched by its coordinates. This is controlled by `dmatch` (default: 1.0 arcsec), setting the maximum matching distance. 

The non-variable candidates are selected if their variance of magnitude difference (`std(m1-m2)`) is smaller than `sdev` (default: 0.006 mag). If there are too few standard stars, an option of `noc` is provided to ensure there are adequate number of standard candidates. This value is not exact the same as the final output of standard stars number.

```
Usage: mergecat.py [OPTIONS]

Options:
  --phot_flag INTEGER      Magnitude type. 0: original aperture photometry; 1:
                           aperture photometry with star center refitted by
                           PSF
  --dmatch FLOAT           Position matching distance in arcsec
  --sdev FLOAT             Standard deviation for none variable selection
  --noc INTEGER            Minimum number of standard candidates. Not the same
                           as Std Stars output
  --medframe_factor FLOAT  A factor on average star number in a frame for
                           reference frame selection
  --obs_flag TEXT          Observatory flag. 'd': Delingha; 'l': Lenghu
  -b, --band TEXT          Passband
  --help                   Show this message and exit.
```

The output files include:
- `stdstar.dat`: photometry of non-variable candidates in the reference frame (columns: `ra`, `dec`, `apmag`, `apmag_err`, `psfmag`, `psfmag_err`)
- `*.{phot_flag}gcat.pkl`: python pickle file that contains all data as a python dictionary

### `correctphot`
This program is to calibrate the instrumental photometry via differential photometry method. 

No external standard stars are required. Instead, the brightest `noc` (default: 5) standard stars from `stdstar.dat` are used to conduct the ensemble photometry, that is, using the average chaning behavior of these standard stars over frames to correct the photometry of the remaining stars. 

`magtype` controls which magnitude to used in this process. `a` for aperture photometry and `p` for psf photometry.

```
Usage: correctphot.py [OPTIONS]

  Calibrate the instrumental photometry

Options:
  -f, --file_name TEXT  Input pkl file name
  --magtype TEXT        magnitude type. a: aperture; p: psf
  --noc INTEGER         Number of selected standard stars
  --help                Show this message and exit.
```
The output files include:
- `*.{phot_flag}{mag_type}gcat_mmag`: average brightness of each star (`index`, `ra`, `dec`, `mmag`, `mmag_err`)
- `*.{phot_flag}{mag_type}gcat_cal.pkl`: python pickle file that contains all data


### `plotcurve`
This program is to present the data in time-series light curve. It can also show the magnitude versus airmass.
```
Usage: plotcurve.py [OPTIONS]

  Plot light curves

Options:
  -f, --file_name TEXT           Input pkl file name
  -i, --init_star_index INTEGER  Index of star to plot
  --help                         Show this message and exit.
```

The default figure is shown in UT horizontal axis with multiple-day data folded into a single frame. It can be switched to day mode that expands the data into the MJD axis.
You can interact with the plot through the keyboard:
- `n`: next star; next day (day mode)
- `N`: previous star; previous day (day mode)
- `d`: switch day mode
- `a`: switch airmass / light curve
- `z`: zoom out
- `Z`: zoom in
- `r`: reset scale
- `s`: save the data of this star 
- `q`: quit

The output files include:
- `*.orig`: original photometry of one star (`index`, `mjd`, `ut`, `mag`, `mag_err`)
- `*.dat`: corrected photometry of one star (`index`, `mjd`, `ut`, `magx`, `magx_err`)

