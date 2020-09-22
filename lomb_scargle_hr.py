import numpy as np
import os
import matplotlib.pyplot as plt
import argparse
from astropy.time import Time
from astropy.timeseries import LombScargle
import readingutils as ru


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
            plt.annotate("", xy=(x[i - 1] ,y[i-1]), xytext=(x[i], y[i]), arrowprops=dict(arrowstyle="<-", shrinkA=10, shrinkB=10, color=color))
        else:
            ax.annotate("", xy=(x[i-1], y[i-1]), xytext=(x[i], y[i]), arrowprops=dict(arrowstyle="<-", shrinkA=10, shrinkB=10, color=color))


ap = argparse.ArgumentParser(description='Plot swift hardness-rations')
ap.add_argument("-t", "--threshold", nargs='?', help="Threshold to separate datasets in the plot (weeks). Default does not apply any separation (i.e. all data plotted together", type=int, default="-1")
ap.add_argument("-m", "--maxcountrate", nargs='?', help="Maximum count rate to consider. Default considers all dataset", type=float, default=np.Infinity)
ap.add_argument("-s", "--source", nargs='?', help="Source (for the title)", type=str, default="")
ap.add_argument("--tmin", nargs='?', help="Minimum time", type=float, default=0)
ap.add_argument("--tmax", nargs='?', help="Maximum time", type=float, default=np.Infinity)
ap.add_argument("--hr", nargs='?', help="HR to use, soft or hard", type=str, default="Hard")
args = ap.parse_args()

count_rate_file = "PCHR.qdp"
#input_file = "PCHR.qdp"
plt.style.use('/home/agurpide/.config/matplotlib/stylelib/paper.mplstyle')
plt.xlabel("Period (days)")
plt.ylabel("Power")

if os.path.isfile(count_rate_file):
    data = np.genfromtxt("%s" % count_rate_file, names=True, delimiter="\t", skip_header=2, comments="!",
                                    dtype=("f8", "f8", "f8", "f8", "f8", "f8", "f8", "f8", "f8", "f8", "f8", "f8", "f8", "U22"))
    print("Using %s count rate" % args.hr)
    rate = "HardRate"
    samples_per_peak = 5
    # filter data by time
    data = np.array([row for row in data if row["Time"] >= args.tmin and row["Time"] <= args.tmax and row["%s" % rate] < args.maxcountrate and row["%sSig" % args.hr] > 0])

    errors = data["%serr" % args.hr]

    T_obs = (data["Time"][-1] - data["Time"][0])
    print("Found %d swift observations spanning %.1f days" % (len(data), (T_obs / 3600 / 24)))
    sampling = np.array([(b - a) / 3600 / 24 for a, b in zip(data["Time"], data["Time"][1:])])
    print("Average sampling %.1f days" % (np.median(sampling)))

    ls = LombScargle(data["Time"], data["%s" % rate], dy=errors, fit_mean=True, nterms=1)
    probabilities = [1 - 0.997300, 1 - 0.999936, 1 - 0.95]
    prob_labels = ["3$\sigma$, 4$\sigma$"]
    frequency, power = ls.autopower(minimum_frequency=1 / T_obs, maximum_frequency=1 / (0.5 * 24 * 3600), samples_per_peak=samples_per_peak)
    plt.plot((1 / frequency) / 3600 / 24, power, color="black")
    if ls.nterms == 1:
        power_prob = ls.false_alarm_level(probabilities)
        print(power_prob)
        for prob, prob_label in zip(power_prob, prob_labels):
            plt.axhline(y=prob, ls="--", color="black", label="%s" % prob_label)
        power_max_prob = ls.false_alarm_probability(power.max(), method='bootstrap', samples_per_peak=samples_per_peak)
        print("False alarm probability of maximum power %.2f" % ((1 - power_max_prob) * 100))
    best_frequency = frequency[np.argmax(power)]
    best_period = (1 / best_frequency / 3600 / 24)
    print("Best period (days) %.1f" % best_period)

    theta = ls.model_parameters(best_frequency)

    best_fit_figure, best_fit_ax = plt.subplots(1, 1)
    folded_figure, folded_ax = plt.subplots(1, 1)
    offset = ls.offset()
    time_model = np.linspace(data["Time"][0], data["Time"][-1], 1000)
    y_fit = ls.model(time_model, best_frequency)
    design_matrix = ls.design_matrix(best_frequency, time_model)
    print(theta)
    #phase = best_frequency * data["Time"] - data["Time"][0] * best_frequency
    phase = (data["Time"] - data["Time"][0]) * best_frequency
    phase = phase % 1
    # read zero point in MJD
    zero_point = ru.read_zero_point()
    sorted_indexes = np.argsort(phase)
    bins = np.linspace(0, 1, 10)
    digitized = np.digitize(phase, bins)
    bin_means = [data["%s" % rate][digitized == i].mean() for i in range(1, len(bins))]
    bin_std = [data["%s" % rate][digitized == i].std() for i in range(1, len(bins))]

    sin_phase = np.linspace(0, 1, 100)
    #plt.plot(sin_phase, y_fit, color="black")
    #best_fit_ax.errorbar(phase, data["Rate"], yerr=(-data["Rateneg"] + data["Ratepos"]) / 2, color="green", ls="None", fmt="-", marker=".")
    best_fit_ax.errorbar(data["Time"], data["%s" % rate], yerr=errors, color="green", ls="None", fmt="-", marker=".")
    best_fit_ax.plot(time_model, y_fit)
    plt.xlabel("Time (s) since %.1f" % zero_point)
    folded_ax.errorbar(bins[1:], bin_means, yerr=bin_std, color="green", drawstyle="steps-mid")
    plt.xlim(0, 2)
    plt.show()

else:
    print("% file not found" % (count_rate_file))
