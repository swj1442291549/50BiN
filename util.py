import numpy as np
from subprocess import check_output

def wc(file_name):
    """Count lines in file

    Args:
        file_name (str): file name

    Returns:
        nc (int): line number
    """
    return int(check_output(["wc", "-l", file_name]).split()[0])




