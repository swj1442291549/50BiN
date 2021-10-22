from setuptools import setup

setup(
    name='catmerge',
    version='0.1.0',
    py_modules=['catmerge'],
    install_requires=[
        'Click',
        "astropy",
        "pandas",
        "numpy"
    ],
    entry_points={
        'console_scripts': [
            'catmerge = catmerge:cli',
        ],
    },
)
