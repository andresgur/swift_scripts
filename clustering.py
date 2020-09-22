import numpy as np
import os
import matplotlib.pyplot as plt
from sklearn.neighbors import NearestNeighbors
from sklearn.cluster import DBSCAN
import readingutils as ru
import argparse
from matplotlib.lines import Line2D


def plot_clusters(data, clusters):
    """Plot input data colored by cluster index"""
    colors = ['royalblue', 'maroon', 'forestgreen', 'mediumorchid', 'tan', 'deeppink', 'olive', 'goldenrod', 'lightcyan', 'navy']
    vectorizer = np.vectorize(lambda x: colors[x % len(colors)])
    plt.scatter(data[:, 0], data[:, 1], c=vectorizer(clusters))
    legend_line = [Line2D([0], [0], marker='o', color='w', label='Rejected', markerfacecolor='navy', markersize=15)]
    plt.legend(handles=legend_line)


ap = argparse.ArgumentParser(description='Cluster data points based on HR and count rate')
ap.add_argument("-e", "--eps", nargs='?', help="EPS used in the DBSCAN algorithm. Radius to look for neighbours", type=float, default=0.5)
ap.add_argument("-m", "--min_samples", nargs='?', help="Minimum samples to consider a cluster.", type=int, default=10)
args = ap.parse_args()

eps = args.eps
min_samples = args.min_samples

count_rate_file = "PCCURVE.qdp"
hr_file = "PCHR.qdp"

if os.path.isfile(count_rate_file) and os.path.isfile(hr_file):

    data_count_rate = ru.readPCCURVE(count_rate_file)
    data_hr = ru.readPCHR(hr_file)
    print("Found %d swift observations in %s" % (len(data_count_rate), count_rate_file))
    print("Found %d swift observations in %s" % (len(data_hr), hr_file))
    # remove negative HR
    data_hr = ru.filter_hr(data_hr, 0, 0, True)

    print("Found %d swift observations in %s after filtering" % (len(data_hr), hr_file))

    # select only common obsid
    data_count_rate = np.array([row_data for row_data in data_count_rate if row_data["Obsid"] in data_hr["Obsid"]])

    # convert both arrays into a single 2D array
    input_data = np.vstack((data_hr["HR"], data_count_rate["Rate"])).T
    # get an estimate of the distance between neighbors
    neigh = NearestNeighbors(n_neighbors=2)
    nbrs = neigh.fit(input_data)
    distances, indices = nbrs.kneighbors(input_data)
    distances = np.sort(distances, axis=0)
    plt.figure()
    plt.xlabel("#")
    plt.ylabel("Distance")
    plt.plot(distances[:, 1])
    m = DBSCAN(eps=eps, min_samples=min_samples, n_jobs=5)
    m.fit(input_data)
    clusters = m.labels_
    print(clusters)
    plt.figure()
    plot_clusters(input_data, clusters)
    plt.xlabel("HR")
    plt.ylabel("Count rate")
    plt.figure()
    plt.errorbar(data_hr["HR"], data_count_rate["Rate"],
                 xerr=data_hr["HRerr"], yerr=[-data_count_rate["Rateneg"],
                 data_count_rate["Ratepos"]], ls="None", fmt="-", color="black", marker="+")
    plt.xlabel("HR")
    plt.ylabel("Count rate")
    plt.show()
