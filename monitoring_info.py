# @Author: Andrés Gúrpide <agurpide>
# @Date:   16-09-2020
# @Email:  agurpidelash@irap.omp.eu
# @Last modified by:   agurpide
# @Last modified time: 23-09-2020

import matplotlib.pyplot as plt
import argparse
import numpy as np
import readingutils as ru

ap = argparse.ArgumentParser(description='Plot swift monitoring information. To be launched from the Swift data folder')
ap.add_argument("--minSNR", nargs='?', help="Minimum SNR to filter the data", type=int, default="0")

args = ap.parse_args()

plt.style.use('/home/agurpide/.config/matplotlib/stylelib/paper.mplstyle')

count_rate_file = "PCCURVE.qdp"

data_count_rate = ru.readPCCURVE(count_rate_file)
data_count_rate = ru.filter_data(data_count_rate, 0, 0, args.minSNR)
# exposure time
figure, ax = plt.subplots(1, 1)
outfile = "exposure_time.png"
keyword = "Exposure"
ax.hist(data_count_rate["%s" % keyword] / 1000, alpha=1.0, linewidth=3, edgecolor="black", color="#002FA7")
median = np.median(data_count_rate["%s" % keyword])
ax.axvline(x=median / 1000, ls="--", color="black")
ax.set_xlabel("%s (ks)" % keyword)
ax.set_ylabel("N observations")
print("Median %s: %.2f" % (keyword, median))
print("Results saved to %s" % outfile)
figure.savefig("%s" % outfile)

# SN
figure, ax = plt.subplots(1, 1)
outfile = "snr.png"
keyword = "SNR"
ax.hist(data_count_rate["%s" % keyword], alpha=1.0, linewidth=3, edgecolor="black", color="#002FA7")
median = np.median(data_count_rate["%s" % keyword])
ax.axvline(x=median, ls="--", color="black")
print("Median %s: %.2f" % (keyword, median))
ax.set_xlabel("%s" % keyword)
ax.set_ylabel("N observations")
print("Results saved to %s" % outfile)
figure.savefig("%s" % outfile)

# Sigma
figure, ax = plt.subplots(1, 1)
outfile = "sigma.png"
keyword = "Sigma"
ax.hist(data_count_rate["%s" % keyword], alpha=1.0, linewidth=3, edgecolor="black", color="#002FA7")
median = np.median(data_count_rate["%s" % keyword])
ax.axvline(x=median, ls="--", color="black")
print("Median %s: %.2f" % (keyword, median))
ax.set_xlabel("%s" % keyword)
ax.set_ylabel("N observations")
print("Results saved to %s" % outfile)
figure.savefig("%s" % outfile)
