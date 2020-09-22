# @Author: Andrés Gúrpide <agurpide>
# @Date:   27-08-2020
# @Email:  agurpidelash@irap.omp.eu
# @Last modified by:   agurpide
# @Last modified time: 18-09-2020



import numpy as np
import matplotlib.pyplot as plt
import os
from astropy.time import Time


def plot_obs_file(obs_file, swift_zero_point, time_ax=None):
    """Overlay epochs of XMM, Chandra and Nustar observations on lightcurve.

        Parameters
        ----------
        obs_file: Joint observation file.
        swift_zero_point: Time
    """
    if os.path.isfile(obs_file):

        obs_data = np.genfromtxt(obs_file, usecols=(0, 1, 2, 3), dtype=("U32", "U32", "U32", "U32"), names=True, delimiter="\t")
        for index, epoch in enumerate(obs_data["epoch"]):
            if obs_data[index]["xmmobsid"] != "":
                color = "red"
            elif obs_data[index]["chandraid"] != "":
                color = "green"
            if time_ax is not None:
                time_ax.axvline(x=Time(epoch).mjd - swift_zero_point.mjd, color=color, ls="--")
            else:
                plt.axvline(x=Time(epoch).mjd - swift_zero_point.mjd, color=color, ls="--")
    else:
        print("%s observation file not found" % obs_file)


def create_color_array(data_length, cmap='hsv'):
    """Create an array of colors given the length of a dataset. Useful for plots where a unique color is needed for each dataset.

    The returned colors come from the input map (jet by default).

    Parameters
    ----------
    data_length : The length of your data for the color array creation.

    """
    print("Creating color array for %i datasets" % data_length)
    x = np.arange(data_length)
    ys = [i + x + (i * x)**2 for i in range(data_length)]
    setmap = plt.get_cmap(name=cmap)

    colors = setmap(np.linspace(0, 1, len(ys)))
    return colors


def create_marker_array(data_length):
    """Get an array of markers given the length of a dataset. Useful for plots where a unique marker is needed for each dataset.

    There are 17 different markers and after that they are repeated.
    Parameters
    ----------
    data_length : The length of your data for the marker array creation.

    """
    m = ['o', '^', (12, 1, 50), "s", 'v', 'p', 'P', 'd', '*', 'h', (5, 1, 10), (3, 0, 10), 's', 'v', "8",
         '<', (12, 1, 120), '.', '>', "s", 'p', 'd', '1', '2']
    while data_length > len(m):
        m.extend(m)
    return m
