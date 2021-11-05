from setuptools import setup

setup(
    name="50bin",
    version="0.2.0",
    py_modules=["catmerge", "correctphot", "plotcurve"],
    install_requires=["Click", "astropy", "pandas", "numpy", "matplotlib", "tqdm"],
    entry_points={
        "console_scripts": ["catmerge = catmerge:cli", "plotcurve = plotcurve:cli", "correctphot = correctphot:cli"],
    },
)
