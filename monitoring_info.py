# @Author: Andrés Gúrpide <agurpide>
# @Date:   16-09-2020
# @Email:  agurpidelash@irap.omp.eu
# @Last modified by:   agurpide
# @Last modified time: 16-09-2020

import matplotlib.pyplot as plt
import argparse
import numpy as np
import readingutils as ru

ap = argparse.ArgumentParser(description='Plot swift monitoring information. To be launched from the Swift data folder')
ap.add_argument("-t", "--threshold", nargs='?', help="Threshold to separate datasets in the plot (weeks). Default does not apply any separation (i.e. all data plotted together", type=int, default="-1")
args = ap.parse_args()

plt.style.use('/home/agurpide/.config/matplotlib/stylelib/paper.mplstyle')

count_rate_file = "PCCURVE.qdp"

data_count_rate = ru.readPCCURVE(count_rate_file)
# exposure time
figure, ax = plt.subplots(1, 1)
outfile = "exposure_time.png"
keyword = "Exposure"
ax.hist(data_count_rate["%s" % keyword] / 1000, alpha=1.0, linewidth=3, edgecolor="black", color="#002FA7")
ax.set_xlabel("%s (ks)" % keyword)
ax.set_ylabel("N observations")

print("Results saved to %s" % outfile)
figure.savefig("%s" % outfile)

# SN
figure, ax = plt.subplots(1, 1)
outfile = "snr.png"
keyword = "SNR"
ax.hist(data_count_rate["%s" % keyword], alpha=1.0, linewidth=3, edgecolor="black", color="#002FA7")
ax.set_xlabel("%s" % keyword)
ax.set_ylabel("N observations")
print("Results saved to %s" % outfile)
figure.savefig("%s" % outfile)

# Sigma
figure, ax = plt.subplots(1, 1)
outfile = "sigma.png"
keyword = "Sigma"
ax.hist(data_count_rate["%s" % keyword], alpha=1.0, linewidth=3, edgecolor="black", color="#002FA7")
ax.set_xlabel("%s" % keyword)
ax.set_ylabel("N observations")
print("Results saved to %s" % outfile)
figure.savefig("%s" % outfile)
