import numpy as np
import os
import matplotlib.pyplot as plt
import readingutils as ru
import argparse
import sys
import error_utils as eu

ap = argparse.ArgumentParser(description='Cluster data points based on HR and count rate')
ap.add_argument("-hr", "--hr", nargs='+', help="List of HR ranges i.e. 0.5:2.0")
ap.add_argument("-c", "--countrate", nargs='+', help="List of count rate ranges i.e. 2.0:2.0 3.0:5.0", type=str)
ap.add_argument("-i", "--input_file", nargs='?', help="Input file with HR and count rates to split data", type=str, default="segregation.txt")
ap.add_argument("--swift_dir", nargs='?', help="Directory with all swift converted rates", type=str, default="")
args = ap.parse_args()


plt.style.use('/home/agurpide/.config/matplotlib/stylelib/paper.mplstyle')

if args.hr is None:
    if os.path.isfile(args.input_file):
        segregation_info = np.genfromtxt(args.input_file, names=True, dtype=("U22, U22"), delimiter="\t", deletechars="")
        print(segregation_info)
        hr_ranges = segregation_info["HR"]
        count_rate_ranges = segregation_info["countrates"]
else:
    hr_ranges = args.hr
    count_rate_ranges = args.countrate
    if len(hr_ranges) != len(count_rate_ranges):
        print("Error number of hr ranges has to be equal to countrate ranges")
        sys.exit()

count_rate_file = "PCCURVE.qdp"
hr_file = "PCHR.qdp"
plt.style.use('/home/agurpide/.config/matplotlib/stylelib/paper.mplstyle')
if os.path.isfile(count_rate_file) and os.path.isfile(hr_file):

    data_count_rate = ru.readPCCURVE(count_rate_file)
    data_hr = ru.readPCHR(hr_file)
    print("Found %d swift observations in %s" % (len(data_count_rate), count_rate_file))
    print("Found %d swift observations in %s" % (len(data_hr), hr_file))
    # remove negative HR
    data_hr = ru.filter_hr(data_hr, 150, 0, True)
    print("Found %d swift observations in %s after filtering" % (len(data_hr), hr_file))
    fig, ax = plt.subplots(1, 1)
    # select only common obsid
    common_obs, common_1, common_2 = np.intersect1d(data_count_rate["Time"], data_hr["Time"], return_indices=True)
    print("Found %d common observations" % len(common_obs))
    data_count_rate = data_count_rate[common_1]
    data_hr = data_hr[common_2]
    chunk_counter = 1
    colors = ["orange", "green", "blue", "yellow", "olive", "purple"]
    #ax.errorbar(data_hr["HR"], data_count_rate["Rate"],
    #            xerr=data_hr["HRerr"], yerr=[-data_count_rate["Rateneg"],
    #            data_count_rate["Ratepos"]], ls="None", fmt="-", color="black", marker="+")

    string_out = "#HR\tcountrate\tobsids\n"
    for hr_range, count_range, color in zip(hr_ranges, count_rate_ranges, colors):
        print("Filtering by %s HR values and %s count rate values" % (hr_range, count_range))
        string_out += "%s\t%s\t" % (hr_range, count_range)
        hr_range = hr_range.split(":")
        count_rate_range = count_range.split(":")
        # filter by count rate
        count_rate_chunk = data_count_rate[(data_count_rate['Rate'] >= float(count_rate_range[0])) & (data_count_rate['Rate'] <= float(count_rate_range[1]))]
        # filter by hr
        hr_chunk = data_hr[(data_hr['HR'] >= float(hr_range[0])) & (data_hr['HR'] <= float(hr_range[1]))]
        common_chunk_obs, indexes_1, indexes_2 = np.intersect1d(hr_chunk["Time"], count_rate_chunk["Time"], return_indices=True)
        final_hr = hr_chunk[indexes_1]
        final_countrate = count_rate_chunk[indexes_2]
        print("Obsids for chunk 1 (total %d/%d)" % (len(common_chunk_obs), len(data_count_rate)))
        obsstring = ""
        for row in final_hr:
            obs = str(row["Obsid"])
            print(obs)
            obsstring += "%s," % obs.split("::ObsID=")[1]
            string_out += "%s," % obs.split("::ObsID=")[1]
        print(obsstring)
        ax.errorbar(final_hr["HR"], final_countrate["Rate"],
                    xerr=final_hr["HRerr"], yerr=[-final_countrate["Rateneg"],
                    final_countrate["Ratepos"]], ls="None", fmt="-", color=color, marker="+")
        chunk_counter += 1
        string_out += "\n"
    file = open("segregation_out.txt", "w")
    file.write(string_out)
    file.close()
    plt.xlabel("HR [1.5 - 10 keV] / [0.3 - 1.5 keV]")
    plt.ylabel("Swift-XRT count rate (0.3 - 10 keV) (cts/s)")
    if args.swift_dir != "":
        if os.path.isdir(args.swift_dir):
            swift_info = ru.read_swift_info("%s/fit_goodness.log" % args.swift_dir)
            for row in swift_info:
                swift_rates = ru.read_swift_converted("%s/swift_models/swift_rates_%s.dat" % (args.swift_dir, row["modelfile"].split("_model.xcm")[0]))
                soft_band = np.mean(swift_rates["**-0.3_1.5-**"])
                soft_band_err = np.std(swift_rates["**-0.3_1.5-**"])
                hard_band = np.mean(swift_rates["**-1.5_10.0-**"])
                hard_band_errr = np.std(swift_rates["**-1.5_10.0-**"])
                total_band = np.mean(swift_rates["**-0.3_10.0-**"])
                total_band_err = np.std(swift_rates["**-0.3_10.0-**"])
                hrratio_errors = eu.error_division(hard_band, soft_band, hard_band_errr, soft_band_err)
                color = "r"
                ax.scatter(hard_band / soft_band, total_band, s=100, zorder=10, marker="*", color=color,
                           facecolors="none")
                ax.errorbar(hard_band / soft_band, total_band, xerr=hrratio_errors, yerr=total_band_err,
                            ls="None", color=color)

    plt.show()
    fig.savefig("segregation.png")
