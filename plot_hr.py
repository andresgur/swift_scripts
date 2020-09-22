# @Author: Andrés Gúrpide <agurpide>
# @Date:   27-08-2020
# @Email:  agurpidelash@irap.omp.eu
# @Last modified by:   agurpide
# @Last modified time: 18-09-2020



import numpy as np
import os
import matplotlib.pyplot as plt
import argparse
from astropy.time import Time
import plot_utils as pu
import readingutils as ru
import error_utils as eu


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


def update_twin_axis(time_ax, time_ax_2, swift_zero_point):
    xticks2 = time_ax.get_xticks()
    xtick2_labels = [time_todate(swift_zero_point.mjd + float(time)).value for time in xticks2]
    time_ax_2.set_xticks(xticks2)
    time_ax_2.set_xticklabels(xtick2_labels, rotation=45)
    time_ax_2.set_xlim(time_ax.get_xlim())
    time_ax_2.figure.canvas.draw()


def clone_time_axis(time_ax, swift_zero_point=Time(0, format="mjd")):
    """Clone Swift time axis into a date axis.

    Parameters:
    -----------
    time_ax: ax, Time axis
    swift_zero_point: Time, reference zero point time in astropy.time.Time
    """

    time_ax_2 = time_ax.twiny()
    time_ax_2.set_xlabel('Date')
    xticks2 = time_ax.get_xticks()
    xtick2_labels = [time_todate(swift_zero_point.mjd + float(time)).value for time in xticks2]
    time_ax_2.set_xticks(xticks2)
    time_ax_2.set_xticklabels(xtick2_labels)
    time_ax_2.set_xlim(time_ax.get_xlim())
    time_ax.callbacks.connect("ylim_changed", lambda time_ax: update_twin_axis(time_ax, time_ax_2, swift_zero_point))


def time_todate(x):
    time = Time(x, format='mjd', scale="tt")
    time.format = "iso"
    time.out_subfmt = 'date'
    return time


def jd_to_daymonthyear(x, pos):
    '''Format the axis to convert from Julian day to real date.'''
    time = Time(x, format='jd')
    time.format = 'iso'
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
ap.add_argument("-i", "--input_observation_file", nargs='?', help="Path to observation file to mark observations", type=str, default="")
ap.add_argument("--tmin", nargs='?', help="Minimum time in MJD swift days (as in the original Swift file)", type=float, default=0)
ap.add_argument("--tmax", nargs='?', help="Maximum time in MJD swift days (as in the original Swift file)", type=float, default=np.Infinity)
ap.add_argument("--swift_dir", nargs='?', help="Directory with all swift converted rates", type=str, default="")
args = ap.parse_args()
if args.threshold > 0:
    time_threshold = args.threshold * 7 * 24 * 3600
    print("Taking chunks of data separated by %d weeks" % args.threshold)
else:
    time_threshold = -1
min_datapoints = args.minobs

min_snr = 2
hr_error_every = 1

plt.style.use('/home/agurpide/.config/matplotlib/stylelib/paper.mplstyle')

f = open("t0.date")
lines = f.readlines()
f.close()
start_date = lines[1].split("at")[0]

count_rate_file = "PCCURVE.qdp"
hr_file = "PCHR.qdp"
# Figure with lightcurve and hr below
hr_figure, hr_ax = plt.subplots(1, 1)
hr_ax.set_xlabel("HR [1.5 - 10 keV] / [0.3 - 1.5 keV]")
hr_ax.set_ylabel("Count rate (0.3 - 10 keV) (s$^-1$)")
# plot with lightcurve and HR
time_figure, (time_ax, time_hr_ax) = plt.subplots(2, sharex=True, gridspec_kw={'hspace': 0, 'wspace': 0})
time_hr_ax.set_xlabel("Days since %s" % start_date)
time_ax.set_ylabel("Swift-XRT count rate (0.3 - 10 keV) (cts/s)")
time_hr_ax.set_ylabel("HR [1.5 - 10 keV] / [0.3 - 1.5 keV]")
#plt.grid(b=True)
# plot with the lightcurve alone
lightcurve_figure, lightcurve_ax = plt.subplots(1, 1)
lightcurve_ax.set_xlabel("Days since %s" % start_date)
lightcurve_ax.set_ylabel("Swift-XRT count rate (0.3 - 10 keV) (cts/s)")
lightcurve_ax_2 = lightcurve_ax.twiny()
lightcurve_ax_2.set_xlabel('Date')

# Figure with lightcurve and hr next to it
double_panel_figure, (time_double_panel, hr_double_panel) = plt.subplots(1, 2, sharey=True, gridspec_kw={'hspace': 0, 'wspace': 0},
                                                                         figsize=(18, 10))
time_double_panel.set_xlabel("Days since %s" % start_date)
time_double_panel.set_ylabel("Swift-XRT count rate (0.3 - 10 keV) (cts/s)")
hr_double_panel.set_xlabel("HR [1.5 - 10 keV] / [0.3 - 1.5 keV]")


if os.path.isfile(count_rate_file) and os.path.isfile(hr_file):

    data_count_rate = ru.readPCCURVE(count_rate_file)
    data_hr = ru.readPCHR(hr_file)
    data_count_rate = np.array([row for row in data_count_rate if row["Time"] >= args.tmin and row["Time"] <= args.tmax])
    data_hr = np.array([row for row in data_hr if row["Time"] >= args.tmin and row["Time"] <= args.tmax])
    print("Found %d swift observations in %s" % (len(data_count_rate), count_rate_file))
    print("Found %d swift observations in %s" % (len(data_hr), hr_file))
    sampling = np.array([(b - a) / 3600 / 24 for a, b in zip(data_count_rate["Time"], data_count_rate["Time"][1:])])
    print("Median sampling %.1f days" % (np.median(sampling)))
    # remove negative HR
    data_hr = ru.filter_hr(data_hr, 150, 0, True)

    print("Found %d swift observations in %s after filtering" % (len(data_hr), hr_file))

    # select only common obsid
    common_obs, common_1, common_2 = np.intersect1d(data_count_rate["Time"], data_hr["Time"], return_indices=True)
    print("Found %d common observations" % len(common_obs))
    data_count_rate = data_count_rate[common_1]
    data_hr = data_hr[common_2]
    # get source zero point in time (for date calculations)
    swift_zero_point = ru.read_zero_point()

    # create groups
    i = 0
    j = 0
    chunk = 0
    colors = create_color_array(10, "tab20")

    if time_threshold != -1:
        while j < len(data_count_rate["Time"]):
            if i == len(data_count_rate["Time"]) - 1:
                print("End reached")

                break
            while ((abs(data_count_rate["Time"][i] - data_count_rate["Time"][i + 1])) < time_threshold):
                i += 1
                if i == len(data_count_rate["Time"]) - 1:

                    print("End reached")
                    break

            accumulated_datapoints = (i - j + 1)
            if accumulated_datapoints < min_datapoints:
                print("Less than %d datapoints (%d). Skipping" % (min_datapoints, accumulated_datapoints))
                i += 1
                j = i
            else:
                print("Chunk %d with %d datapoints. i index %d" % (chunk + 1, accumulated_datapoints, i))
                hr_ax.errorbar(data_hr["HR"][j:i + 1], data_count_rate["Rate"][j:i + 1], xerr=data_hr["HRerr"][j:i + 1],
                               yerr=[-data_count_rate["Rateneg"][j:i + 1], data_count_rate["Ratepos"][j:i + 1]], ls="None", marker=".", errorevery=hr_error_every, markersize=8, elinewidth=0.5,
                               label=" %.1f d (%d) %s" % (float(data_count_rate["Time"][i] - data_count_rate["Time"][j]) / 3600 / 24, accumulated_datapoints,
                               time_todate(swift_zero_point.mjd + float(data_count_rate["Time"][j]) / 3600 / 24).value), markerfacecolor="none",
                               color=colors[chunk])
                draw_arrows(data_hr["HR"][j:i + 1], data_count_rate["Rate"][j:i + 1], color=colors[chunk], ax=hr_ax)
                time_ax.errorbar(data_count_rate["Time"][j:i + 1] / 3600 / 24, data_count_rate["Rate"][j:i + 1],
                                 yerr=[-data_count_rate["Rateneg"][j:i + 1], data_count_rate["Ratepos"][j:i + 1]],
                                 xerr=[-data_count_rate["T_ve_1"][j:i + 1] / 3600 / 24, data_count_rate["T_ve"][j:i + 1] / 3600 / 24],
                                 ls="None", marker="+", color=colors[chunk])
                time_hr_ax.errorbar(data_count_rate["Time"][j:i + 1] / 3600 / 24, data_hr["HR"][j:i + 1],
                                    xerr=[-data_count_rate["T_ve_1"][j:i + 1] / 3600 / 24, data_count_rate["T_ve"][j:i + 1] / 3600 / 24],
                                    yerr=[-data_count_rate["Rateneg"][j:i + 1], data_count_rate["Ratepos"][j:i + 1]],
                                    ls="None", marker="+", color=colors[chunk])
                lightcurve_ax.errorbar(data_count_rate["Time"][j:i + 1] / 3600 / 24, data_count_rate["Rate"][j:i + 1],
                                       yerr=[-data_count_rate["Rateneg"][j:i + 1], data_count_rate["Ratepos"][j:i + 1]],
                                       xerr=[-data_count_rate["T_ve_1"][j:i + 1] / 3600 / 24, data_count_rate["T_ve"][j:i + 1] / 3600 / 24],
                                       ls="None", marker="+", color=colors[chunk])
                # double panel axes
                time_double_panel.errorbar(data_count_rate["Time"][j:i + 1] / 3600 / 24, data_count_rate["Rate"][j:i + 1],
                                           yerr=[-data_count_rate["Rateneg"][j:i + 1], data_count_rate["Ratepos"][j:i + 1]],
                                           xerr=[-data_count_rate["T_ve_1"][j:i + 1] / 3600 / 24, data_count_rate["T_ve"][j:i + 1] / 3600 / 24],
                                           ls="None", marker=".", color=colors[chunk])
                hr_double_panel.errorbar(data_hr["HR"][j:i + 1], data_count_rate["Rate"][j:i + 1], xerr=data_hr["HRerr"][j:i + 1],
                                          ls="None", marker=".", errorevery=hr_error_every, label=" %.1f d (%d) %s" % (float(data_count_rate["Time"][i] - data_count_rate["Time"][j]) / 3600 / 24, accumulated_datapoints,
                                          time_todate(swift_zero_point.mjd + float(data_count_rate["Time"][j]) / 3600 / 24).value),
                                          color=colors[chunk])
                draw_arrows(data_hr["HR"][j:i + 1], data_count_rate["Rate"][j:i + 1], color=colors[chunk], ax=hr_double_panel)
                i += 1
                j = i
                chunk += 1
            hr_ax.legend()
            hr_double_panel.legend()
    else:
        hr_ax.errorbar(data_hr["HR"], data_count_rate["Rate"], xerr=data_hr["HRerr"], yerr=[-data_count_rate["Rateneg"], data_count_rate["Ratepos"]],
                       ls="None", fmt="-", color="black", marker="+", errorevery=hr_error_every)
        time_ax.errorbar(data_hr["Time"] / 3600 / 24, data_count_rate["Rate"], xerr=[-data_count_rate["T_ve_1"] / 3600 / 24, data_count_rate["T_ve"] / 3600 / 24], yerr=[-data_count_rate["Rateneg"],
                         data_count_rate["Ratepos"]], ls="None", fmt="-", color="black", marker="+")
        lightcurve_ax.errorbar(data_hr["Time"] / 3600 / 24, data_count_rate["Rate"], xerr=[-data_count_rate["T_ve_1"] / 3600 / 24, data_count_rate["T_ve"] / 3600 / 24], yerr=[-data_count_rate["Rateneg"],
                               data_count_rate["Ratepos"]], ls="None", fmt="-", color="black", marker="+")
        # double panel axes
        time_double_panel.errorbar(data_hr["Time"] / 3600 / 24, data_count_rate["Rate"], xerr=[-data_count_rate["T_ve_1"] / 3600 / 24, data_count_rate["T_ve"] / 3600 / 24], yerr=[-data_count_rate["Rateneg"],
                                   data_count_rate["Ratepos"]], ls="None", fmt="-", color="black", marker="+")
        hr_double_panel.errorbar(data_hr["HR"], data_count_rate["Rate"], xerr=data_hr["HRerr"],
                                  ls="None", fmt="-", color="black", marker=".", errorevery=hr_error_every)

    hr_ax.set_title(args.source)
    time_figure.suptitle(args.source)
    lightcurve_ax.set_title(args.source)
    double_panel_figure.suptitle(args.source)

    if args.input_observation_file != "":
        pu.plot_obs_file(args.input_observation_file, swift_zero_point, time_ax)
        pu.plot_obs_file(args.input_observation_file, swift_zero_point, time_hr_ax)
        pu.plot_obs_file(args.input_observation_file, swift_zero_point, lightcurve_ax)
        pu.plot_obs_file(args.input_observation_file, swift_zero_point, time_double_panel)


    # plot observations from other instruments
    time_ax_double_panel_lims = time_double_panel.get_xlim()
    if args.swift_dir != "":
        if os.path.isdir(args.swift_dir):
            swift_info = ru.read_swift_info("%s/fit_goodness.log" % args.swift_dir)
            swift_info.sort(order="epoch")
            colors = pu.create_color_array(len(swift_info), "plasma")
            markers = pu.create_marker_array(len(swift_info))
            marker_edge_color = "black"
            for row, color, marker in zip(swift_info, colors, markers):
                swift_rates = ru.read_swift_converted("%s/swift_models/swift_rates_%s.dat" % (args.swift_dir, row["modelfile"].split("_model.xcm")[0]))
                soft_band = np.mean(swift_rates["**-0.3_1.5-**"])
                soft_band_err = np.std(swift_rates["**-0.3_1.5-**"])
                hard_band = np.mean(swift_rates["**-1.5_10.0-**"])
                hard_band_errr = np.std(swift_rates["**-1.5_10.0-**"])
                total_band = np.mean(swift_rates["**-0.3_10.0-**"])
                total_band_err = np.std(swift_rates["**-0.3_10.0-**"])
                hrratio_errors = eu.error_division(hard_band, soft_band, hard_band_errr, soft_band_err)

                #color = "r"
                hr_ax.scatter(hard_band / soft_band, total_band, s=100, zorder=10, marker="*", color=color, facecolors="none")
                hr_ax.errorbar(hard_band / soft_band, total_band, xerr=hrratio_errors, yerr=total_band_err, ls="None", color=color)
                    #hr_ax.errorbar(swift_rates["**-1.5_10.0-**"] / swift_rates["**-0.3_1.5-**"],
            #              swift_rates["**-0.3_10.0-**"], xerr=hrratio_errors, yerr=swift_rates["**-0.3_10.0-**err"],
            #              marker="None", color="black", ls="None", elinewidth=0.1)
                epoch = ru.get_epoch(row["files"], args.swift_dir).mjd - swift_zero_point.mjd
                print("Counts for epoch %s (%.2f)" % (row["epoch"], total_band))
                time_ax.scatter(epoch, total_band, color=color, s=100, zorder=10, marker=marker)
                time_ax.errorbar(epoch, total_band, yerr=total_band_err, ls="None", color=color)
                lightcurve_ax.scatter(epoch, total_band, color=color, s=100, zorder=10, marker=marker)
                lightcurve_ax.errorbar(epoch, total_band, yerr=total_band_err, ls="None", color=color)
                # double panel plots
                time_double_panel.scatter(epoch, total_band, color=color, s=100, zorder=10, marker=marker, edgecolors=marker_edge_color)
                time_double_panel.errorbar(epoch, total_band, yerr=total_band_err, ls="None", color=color)
                hr_double_panel.scatter(hard_band / soft_band, total_band, s=100, zorder=10, marker=marker, color=color, edgecolors=marker_edge_color)
                hr_double_panel.errorbar(hard_band / soft_band, total_band, xerr=hrratio_errors, yerr=total_band_err,
                                          ls="None", color=color)
    # copy time axis for both double plot and simple lightcurve plot
    clone_time_axis(time_ax, swift_zero_point)
    clone_time_axis(lightcurve_ax, swift_zero_point)
    time_double_panel.set_xlim(time_ax_double_panel_lims)
    double_panel_figure.canvas.draw()


    # double panel plot
    hr_figure.savefig("ratio_hr.png")
    time_figure.savefig("time_countrate.png")
    lightcurve_figure.savefig("lightcurve.png")
    double_panel_figure.savefig("double_panel_light_curve_hr.png")
    #plt.show()

else:
    print("%s or/and %s files not found" % (count_rate_file, hr_file))
