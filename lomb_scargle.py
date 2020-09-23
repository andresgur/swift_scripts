# @Author: Andrés Gúrpide <agurpide>
# @Date:   01-09-2020
# @Email:  agurpidelash@irap.omp.eu
# @Last modified by:   agurpide
# @Last modified time: 22-09-2020



import numpy as np
import os
import matplotlib.pyplot as plt
import argparse
from astropy.time import Time
import math as m
from astropy.timeseries import LombScargle
import readingutils as ru
import plot_utils as pu
from scipy.signal import find_peaks


def phase_folding(time, rate, best_frequency, time_0=0, n_bins=10, cycles=2):
    phases = (time - time_0) * best_frequency
    phases = phases % 1
    phased_bins = np.array([int(np.floor(phase * n_bins)) for phase in phases])
    bins = np.arange(0, n_bins)
    bin_means = [rate[phased_bins == i].mean() for i in bins]
    bin_stds = [rate[phased_bins == i].std() for i in bins]
    bin_means = np.hstack([bin_means, bin_means])
    bin_stds = np.hstack([bin_stds, bin_stds])
    bins = bins / n_bins + 0.05
    two_phase_bins = np.hstack([bins, bins + 1])
    return bin_means, bin_stds, two_phase_bins


def unused_code():

    phase = (data["Time"] - data["Time"][0]) * best_frequency
    phase = phase % 1

    rates = np.hstack([data["Rate"], data["Rate"]])
    phases = np.hstack([phase + 1, phase])

    sorted_indexes = np.argsort(phase)
    bins = np.linspace(0, 1, 11)
    digitized = np.digitize(phase, bins, "left")

    bin_means = [data["Rate"][digitized == i].mean() for i in range(0, len(bins) - 1)]
    n_points_per_bin = [len(data["Rate"][digitized == i]) for i in range(0, len(bins) - 1)]
    print("Mean count rate per phase bin \n")
    print("Phase-bins")
    print(bins)
    print("Count rates")
    print(bin_means)
    print("Points per bin")
    print(n_points_per_bin)

    yerr_pos = [m.sqrt(sum([data**2 for data in data["Ratepos"][digitized == i]])) for i in range(0, len(bins) - 1)]
    yerr_neg = [m.sqrt(sum([data**2 for data in data["Rateneg"][digitized == i]])) for i in range(0, len(bins) - 1)]

    sin_phase = np.linspace(0, 1, 100)


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


ap = argparse.ArgumentParser(description='Compute the generalized periodogram. (see Vanderplass et al 2018 10.3847/1538-4365/aab766.)')
ap.add_argument("-m", "--maxcountrate", nargs='?', help="Maximum count rate to consider. Default considers all dataset", type=float, default=np.Infinity)
ap.add_argument("-s", "--source", nargs='?', help="Source (for the title)", type=str, default="")
ap.add_argument("-fmax", "--fmax", nargs='?', help="Maximum frequency to explore in the lomb scargle periodogram (in 1/days). Default 5 1/days", type=float, default=1 / 5)
ap.add_argument("-fmin", "--fmin", nargs='?', help="Minimum frequency to explore in the lomb scargle periodogram. Default 1 / (2 * Tobs)", type=float, default=-1)
ap.add_argument("--tmin", nargs='?', help="Minimum time in MJD swift days (as in the original Swift file)", type=float, default=0)
ap.add_argument("--tmax", nargs='?', help="Maximum time in MJD swift days (as in the original Swift file)", type=float, default=np.Infinity)
ap.add_argument("-i", "--input_observation_file", nargs='?', help="Path to observation file to mark observations", type=str, default="")
ap.add_argument("-n", "--nterms", nargs='?', help="Number of terms for the periodicity analysis", type=int, default=1)
ap.add_argument("-o", "--outdir", nargs='?', help="Output dir", type=str, default="lomb_scargle")
ap.add_argument("--swift_dir", nargs='?', help="Directory with all swift converted rates", type=str, default="")
args = ap.parse_args()

count_rate_file = "PCCURVE.qdp"
#input_file = "PCHR.qdp"
plt.style.use('/home/agurpide/.config/matplotlib/stylelib/email.mplstyle')

f = open("t0.date")
lines = f.readlines()
f.close()
start_date = lines[1].split("at")[0]
# Periodogram in period log scale
periodogram_figure, periodogram_ax = plt.subplots(1, 1)
plt.xlabel("Period (day)")
plt.ylabel("Lomb-Scargle Power")
# observing window
observing_windonw_figure, obserwing_window_ax = plt.subplots(1, 1)
plt.xlabel("Frequency (days$^{-1}$)")
plt.ylabel("Lomb-Scargle Power")
# periodogram in frequency
periodogram_freq_figure, periodogram_freq_ax = plt.subplots(1, 1)
plt.xlabel("Frequency (days$^{-1}$)")
plt.ylabel("Lomb-Scargle Power")

twin_axis = 0
nterms = args.nterms

if os.path.isfile(count_rate_file):
    data = ru.readPCCURVE("%s" % count_rate_file)
    rate = "Rate"
    samples_per_peak = 5
    # filter data by time
    data = np.array([row for row in data if row["Time"] >= args.tmin and row["Time"] <= args.tmax and row["Rate"] < args.maxcountrate])
    errors = (-data["%sneg" % rate] + data["%spos" % rate]) / 2
    T_obs = (data["Time"][-1] - data["Time"][0])
    print("Found %d swift observations spanning %.1f days" % (len(data), (T_obs / 3600 / 24)))
    sampling = np.array([(b - a) / 3600 / 24 for a, b in zip(data["Time"], data["Time"][1:])])
    print("Median sampling %.1f days" % (np.median(sampling)))
    maximum_frequency = (args.fmax) / 3600 / 24
    minimum_frequency = 1 / (T_obs) if args.fmin == -1 else (args.fmin) / 3600 / 24
    period_range = "%.1f-%.1f" % ((1 / (maximum_frequency * 3600 * 24)), (1 / (minimum_frequency * 3600 * 24)))
    print("Period range explored: %s (days)" % period_range)
    outdir = "%s%s" % (args.outdir, period_range)
    if not os.path.isdir(outdir):
        os.mkdir(outdir)
    # save time range used
    time_range_file = open("%s/time_range.txt" % (outdir), "w+")
    time_range_file.write("%.5f-%.5f" % (data["Time"][0], data["Time"][-1]))
    time_range_file.close()
    # signal transform
    ls = LombScargle(data["Time"], data["%s" % rate], dy=errors, fit_mean=True, nterms=nterms, center_data=False)
    # window transform, do not precenter or fit mean substitute values by ones
    window_ls = LombScargle(data["Time"], np.ones(data["Time"].shape), fit_mean=False, nterms=nterms, center_data=False)

    frequency_window, power_window = window_ls.autopower(minimum_frequency=minimum_frequency,
                                                         maximum_frequency=maximum_frequency, samples_per_peak=samples_per_peak)
    frequency_days = frequency_window * 24 * 3600
    obserwing_window_ax.plot(frequency_days, power_window, color="black")
    peak_power_indexes, properties = find_peaks(power_window, height=0.15)

    results = ""
    for peak_index in peak_power_indexes:
        period = 1 / frequency_days[peak_index]
        results += "%.2f\t%.5f\t%.2f\n" % (power_window[peak_index], frequency_days[peak_index], period)
    window_peaks = open("%s/window_peaks.dat" % (outdir), "w+")
    window_peaks.write("#power\tfrequency\tperiod\n")
    window_peaks.write(results)
    window_peaks.close()

    probabilities = [1 - 0.95, 1 - 0.997300]
    # 1 - 0.999936 4 sigma
    prob_labels = ["2$\sigma$", "3$\sigma$"]

    linestyles = [":", "--"]

    frequency, power = ls.autopower(minimum_frequency=minimum_frequency, maximum_frequency=maximum_frequency,
                                    samples_per_peak=samples_per_peak)
    periodogram_ax.plot((1 / frequency) / 3600 / 24, power, color="black")

    periodogram_freq_ax.plot(frequency_days, power, color="black")
    for ax in [periodogram_freq_ax, obserwing_window_ax]:
        ax.set_xlim(left=min(frequency_days))

    if ls.nterms == 1:
        power_prob = ls.false_alarm_level(probabilities, minimum_frequency=minimum_frequency, maximum_frequency=maximum_frequency, method="bootstrap", samples_per_peak=samples_per_peak, method_kwds={"n_bootstraps": 1000})
        for prob_power, prob_label, linestyle in zip(power_prob, prob_labels, linestyles):
            periodogram_ax.axhline(y=prob_power, ls=linestyle, color="black", label="%s" % prob_label)
            periodogram_freq_ax.axhline(y=prob_power, ls=linestyle, color="black", label="%s" % prob_label)
        power_max_prob = ls.false_alarm_probability(power.max(), minimum_frequency=minimum_frequency, maximum_frequency=maximum_frequency, method="bootstrap", samples_per_peak=samples_per_peak, method_kwds={"n_bootstraps": 1000})
        print("False alarm probability of maximum power %.2f" % ((1 - power_max_prob) * 100))
        print("Period candidates:")
        power_candidates, properties = find_peaks(power, height=power_prob[0].value * 0.9)
        period_candidates = 1 / frequency[power_candidates] / 3600 / 24
        [print("%.2f" % period) for period in period_candidates]
        print("\n")
    results = ""
    for peak_index in power_candidates:
        period = 1 / frequency_days[peak_index]
        prob_peak = ls.false_alarm_probability(power[peak_index], minimum_frequency=minimum_frequency, maximum_frequency=maximum_frequency, method="bootstrap", samples_per_peak=samples_per_peak, method_kwds={"n_bootstraps": 1000})
        results += "%.2f\t%.5f\t%.2f\t%.3f\n" % (power[peak_index], frequency_days[peak_index], period, ((1 - prob_peak) * 100))
    window_peaks = open("%s/power_peaks.dat" % (outdir), "w+")
    window_peaks.write("#power\tfrequency\tperiod\tprobability\n")
    window_peaks.write(results)
    window_peaks.close()
    periodogram_ax.legend()
    periodogram_freq_ax.legend()
    for ax in [periodogram_ax, periodogram_freq_ax, obserwing_window_ax]:
        ax.set_ylim(bottom=0)
    periodogram_figure.savefig("%s/periodogram_period.png" % outdir)
    periodogram_freq_ax.set_title("Periodogram")
    periodogram_freq_figure.savefig("%s/periodogram_freq.png" % outdir)
    obserwing_window_ax.set_title("Observing window")
    observing_windonw_figure.savefig("%s/observing_window_freq.png" % outdir)
    periodogram_ax.plot((1 / frequency_window) / 3600 / 24, power_window, color="grey", ls="--", alpha=0.5)
    periodogram_figure.savefig("%s/periodogram_period_with_window.png" % outdir)
    periodogram_freq_ax.plot(frequency_window * 24 * 3600, power_window, color="grey", ls="--", alpha=0.5)
    plt.savefig("%s/periodogram_frequency_with_window.png" % outdir)
    plt.show()
    best_frequency = frequency[np.argmax(power)]
    #best_frequency = 1 / (145 * 24 * 3600)
    best_period = (1 / best_frequency / 3600 / 24)
    print("Best period (days) %.1f" % best_period)

    theta = ls.model_parameters(best_frequency)

    best_fit_figure, best_fit_ax = plt.subplots(1, 1)
    folded_figure, folded_ax = plt.subplots(1, 1)
    offset = ls.offset()
    time_model = np.linspace(data["Time"][0], data["Time"][-1], 1000)
    y_fit = ls.model(time_model, best_frequency)
    design_matrix = ls.design_matrix(best_frequency, time_model)
    phase_model = np.linspace(0, 2, 20)

    #plt.plot(sin_phase, y_fit, color="black")
    #best_fit_ax.errorbar(phase, data["Rate"], yerr=(-data["Rateneg"] + data["Ratepos"]) / 2, color="green", ls="None", fmt="-", marker=".")
    best_fit_ax.errorbar(data["Time"] / 3600 / 24, data["%s" % rate], yerr=errors, color="black", ls="None", fmt="-", marker=".")
    best_fit_ax.plot(time_model / 3600 / 24, y_fit, color="blue", ls="--")
    best_fit_ax.set_xlabel("Days since %s" % start_date)
    best_fit_ax.set_ylabel("Swift-XRT count rate (0.3 - 10 keV) (ct/s)")
    swift_zero_point = ru.read_zero_point()
    if twin_axis:
        time_ax_2 = best_fit_ax.twiny()
        time_ax_2.set_xlabel('Time (s)')
        xticks2 = best_fit_ax.get_xticks()
        xtick2_labels = [float(time) * 3600 * 24 for time in xticks2]
        time_ax_2.set_xticks(xticks2)
        time_ax_2.set_xticklabels(xtick2_labels)
        time_ax_2.set_xlim(best_fit_ax.get_xlim())

    if args.swift_dir != "":
        if os.path.isdir(args.swift_dir):
            swift_info = ru.read_swift_info("%s/swift_rates.config" % args.swift_dir)
            for row in swift_info:
                swift_rates = ru.read_swift_converted("%s/swift_models/swift_rates_%s.dat" % (args.swift_dir, row["model"].split("_model.xcm")[0]))
                total_band = np.mean(swift_rates["**-0.3_10.0-**"])
                total_band_err = np.std(swift_rates["**-0.3_10.0-**"])
                color = "r"
                if row["chandra"] != "":
                    color = "purple"
                best_fit_ax.errorbar(Time(row["epoch"]).mjd - swift_zero_point.mjd, total_band, yerr=total_band_err, color=color)
                best_fit_ax.scatter(Time(row["epoch"]).mjd - swift_zero_point.mjd, total_band, s=100, zorder=10, marker="*", color=color)
    best_fit_figure.savefig("%s/best_fit_period.png" % outdir)
    bin_means, bin_stds, bins = phase_folding(data["Time"], data["Rate"], best_frequency, 0)
    folded_ax.errorbar(bins, bin_means, color="black", drawstyle="steps-mid", fmt="-", marker="+")
    #folded_ax.plot(phase_model, ls.model(phafse_model / best_frequency, best_frequency))
    phase = (data["Time"] - 0) * best_frequency
    phase = phase % 1
    rates = np.hstack([data["Rate"], data["Rate"]])
    errors = np.hstack([data["Ratepos"], data["Ratepos"]])
    phases = np.hstack([phase + 1, phase])
    #folded_ax.errorbar(phases, rates, yerr=errors, ls="None")
    plt.ylabel("Swift-XRT Count rate (counts / s)")
    plt.xlabel("Phase")
    plt.xlim(0, 2)
    plt.savefig("%s/best_folded.png" % outdir)
    print("Results stored in %s" % outdir)
else:
    print("% file not found" % (count_rate_file))
