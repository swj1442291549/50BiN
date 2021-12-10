from setuptools import setup

setup(
    name="50bin",
    version="0.5.1",
    py_modules=["mergecat", "correctphot", "plotcurve"],
    install_requires=[
        "Click",
        "astropy",
        "pandas",
        "numpy",
        "matplotlib",
        "tqdm",
        "statsmodels",
    ],
    entry_points={
        "console_scripts": [
            "mergecat = mergecat:cli",
            "plotcurve = plotcurve:cli",
            "correctphot = correctphot:cli",
        ],
    },
)
