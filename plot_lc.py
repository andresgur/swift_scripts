import numpy as np
import os
import matplotlib.pyplot as plt
import argparse
from astropy.time import Time


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


def time_todate(x):
    time = Time(x, format='mjd', scale="tt")
    time.format = "iso"
    time.out_subfmt = 'date'
    return time


def draw_arrows(x, y, color, ax=None):
    for i in np.arange(1, len(x)):
        if ax is None:
            plt.annotate("", xy=(x[i - 1], y[i - 1]), xytext=(x[i], y[i]), arrowprops=dict(arrowstyle="<-", shrinkA=10, shrinkB=10, color=color))
        else:
            ax.annotate("", xy=(x[i - 1], y[i - 1]), xytext=(x[i], y[i]), arrowprops=dict(arrowstyle="<-", shrinkA=10, shrinkB=10, color=color))


ap = argparse.ArgumentParser(description='Plot swift hardness-rations')
ap.add_argument("-t", "--threshold", nargs='?', help="Threshold to separate datasets in the plot (weeks). Default does not apply any separation (i.e. all data plotted together", type=int, default="-1")
ap.add_argument("-m", "--minobs", nargs='?', help="Minimum number of observations to accept a chunk of data. Default 0 i.e. all chunks considered", type=int, default="0")
ap.add_argument("-s", "--source", nargs='?', help="Source (for the title)", type=str, default="")
args = ap.parse_args()
if args.threshold > 0:
    time_threshold = args.threshold * 7 * 24 * 3600
    print("Taking chunks of data separated by %d weeks" %args.threshold)
else:
    time_threshold = -1
min_datapoints = args.minobs

input_file = "countrate_hr.txt"
plt.style.use('/home/agurpide/.config/matplotlib/stylelib/paper.mplstyle')
plt.xlabel("HR [1.5 - 10 keV] / [0.3 - 1.5 keV]")
plt.ylabel("Count rate (0.3 - 10 keV) (s$^-1$)")

if os.path.isfile(input_file):
    data = np.genfromtxt("%s" % input_file, names=True, delimiter="\t", dtype=("f8", "f8", "f8", "f8", "f8", "f8", "f8", "f8"),
                         usecols=(0, 1, 2, 3, 4, 5, 6, 7))

    print("Found %d swift observations" % len(data))

    # get source zero point in time (for date calculations)
    f = open("./../zero-point.txt")
    lines = f.readlines()
    f.close()
    zero_point = float(lines[0])

    # first of january of 2001 in mjdr + zero_poin to days
    swift_ref_time = Time(51910, format="mjd", scale="tt")
    swift_zero_point = Time(swift_ref_time.mjd + zero_point / 24 / 3600, format="mjd", scale="tt")
    # create groups
    i = 0
    j = 0
    chunk = 0
    colors = create_color_array(10, "tab10")

    if time_threshold != -1:
        while j < len(data["Time"]):
            if i == len(data["Time"]) - 1:
                print("End reached")
                break
            while ((abs(data["Time"][i] - data["Time"][i + 1])) < time_threshold):
                i += 1
                if i == len(data["Time"]) - 1:
                    i -= 1
                    print("End reached")
                    break

            accumulated_datapoints = (i - j + 1)
            if accumulated_datapoints < min_datapoints:
                print("Less than %d datapoints (%d). Skipping" % (min_datapoints, accumulated_datapoints))
                i += 1
                j = i
            else:
                print("Chunk %d with %d datapoints. i index %d" % (chunk + 1, accumulated_datapoints, i))
                plt.errorbar(data["HR"][j:i + 1], data["Rate"][j:i + 1], xerr=data["HRerr"][j:i + 1],
                             yerr=[-data["Rateneg"][j:i + 1], data["Ratepos"][j:i + 1]], ls="None", marker=".", errorevery=1, markersize=8, elinewidth=0.5,
                             label=" %.1f d (%d) %s" % (float(data["Time"][i] - data["Time"][j]) / 3600 / 24, accumulated_datapoints,
                             time_todate(swift_zero_point.mjd + float(data["Time"][j]) / 3600 / 24).value), color=colors[chunk], markerfacecolor="white")
                draw_arrows(data["HR"][j:i + 1], data["Rate"][j:i + 1], color=colors[chunk])
                i += 1
                j = i
                chunk += 1
            plt.legend()
    else:
        plt.errorbar(data["HR"], data["Rate"], xerr=data["HRerr"], yerr=[-data["Rateneg"], data["Ratepos"]], ls="None", fmt="-", color="black")
    plt.title(args.source)
    plt.savefig("ratio_hr.png")
    plt.show()
else:
    print("Merge curve.qdp and PCHR.qdp into a file named countrate_hr.txt \n #Time\tT_+ve\tT_-ve\tRate\tRatepos\tRateneg\tHR\tHRerr")
