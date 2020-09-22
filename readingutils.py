# @Author: Andrés Gúrpide <agurpide>
# @Date:   02-09-2020
# @Email:  agurpidelash@irap.omp.eu
# @Last modified by:   agurpide
# @Last modified time: 16-09-2020



import numpy as np
from astropy.time import Time
import os
from astropy.io import fits


def get_epoch(data_file, parent_path="."):
    """Gets epoch from a data file.xcm from the first file present in the file.

    Returns the epoch in Date format.

    Parameters
    ----------
    data_file: string, The file .xcm containing the datasets loaded.
    parent_path: string, The path to the directory where the file is
    """
    # save current path so we can go back to it later
    current_dir = os.getcwd()
    os.chdir(parent_path)
    parent_dir = parent_path.split("/")[-1]
    os.chdir("../")

    f = open("%s/%s" % (parent_dir, data_file))
    lines = f.readlines()
    f.close()

    for index, line in enumerate(lines):
        if "cd" in line:
            data_dir = line.split("cd")[1].strip()
            os.chdir(data_dir)
            instrument_file = lines[index + 1].split("data 1:1")[1].strip()
            print("Instrument file found: %s for epoch: %s" % (instrument_file, data_file))
            hdu_list = fits.open(instrument_file)
            date_with_time = Time(hdu_list[0].header["DATE-OBS"], format="isot", scale="utc")
            os.chdir(current_dir)
            print("Epoch %s" % date_with_time.value)
            return date_with_time

    print("Error keyword cd not found!!!!!!!!")


def read_zero_point(file="t0.date"):
    """Read zero point value from the input file in MJD

        If the argument `file` isn't passed in, the default is t0.date file from Swift pipeline
        sound is used.

        Parameters
        ----------
        file : str, optional
            The file where to look the zero point (default is t0.date)
        """
    f = open(file)
    lines = f.readlines()
    f.close()
    zero_point = float(lines[3])
    # first of january of 2001 in mjd https://swift.gsfc.nasa.gov/analysis/suppl_uguide/time_guide.html'+ zero_poin to days
    # from this file we have seconds since swift reference time, starting date, modified julian days and julian days (http://www.csgnetwork.com/julianmodifdateconv.html)
    # swift_ref_time = Time(51910, format="mjd", scale="tt")
    swift_zero_point = Time(zero_point, format="mjd", scale="tt")
    return swift_zero_point


def filter_hr(data_hr, minSoftSig=0, minHardSig=0, reject_errors=True):
    """Filters hardness ratio data from Swift data.
        Filtering is done based on minimum soft and hard detection.
        A flag to remove errors higher than the data point value can be given.
        Data points with negative HR are also filtered. Nan values are also removed.

        Parameters
        ----------
        data_hr: ndarray
            The data. Must contain HR, SoftSig, HardSig and HRerrr columns.
        minSoftSig : float, optional
            Minimum soft signal to filer. Default is 0.
        minHardSig : float, optional
            Minimum hard signal to filer. Default is 0.
        reject_errors : boolean, optional
            Whether to reject data points with errors higher than the data point value. Default is True.
        """
    if reject_errors:
        return np.array([row_data for row_data in data_hr if row_data["HR"] > 0 and not
                        np.isnan(row_data["HR"]) and row_data["SoftSig"] > minSoftSig and row_data["HardSig"] > minHardSig and row_data["HRerr"] < row_data["HR"]])
    else:
        return np.array([row_data for row_data in data_hr if row_data["HR"] > 0 and not
                        np.isnan(row_data["HR"]) and row_data["SoftSig"] > minSoftSig and row_data["HardSig"] > minHardSig])


def readPCCURVE(file="PCCURVE.qdp"):
    """Read PCCURVE from Swift data pipeline.

        Parameters
        ----------
        file: str, optional
            The file to be read. Default is PCCURVE.qdp
        """
    print("Reading %s data" % file)
    return np.genfromtxt("%s" % file, names=True, delimiter="\t", skip_header=2, comments="!", dtype=("f8, f8, f8, f8, f8, f8, f8, f8, f8, f8, i8, f8, f8, f8, f8, U30"))


def readPCHR(file="PCHR.qdp"):
    """Read PCHR from Swift data pipeline.

        Parameters
        ----------
        file: str, optional
            The file to be read. Default is PCHR.qdp
        """
    print("Reading %s data" % file)
    return np.genfromtxt("%s" % file, names=True, delimiter="\t", skip_header=2, comments="!", dtype=("f8, f8, f8, f8, f8, f8, f8, f8, f8, f8, f8, f8, f8, U30"))


def read_swift_converted(file="swift_rates.dat"):
    """Read swift file with converted rates from XMM Newton data.

        Parameters
        ----------
        file: str, optional
            Path to the file to be read. Default is swift_rates.dat
        """
    print("Reading %s data" % file)
    return np.genfromtxt("%s" % file, names=True, delimiter="\t", comments="#",
                         dtype=("U32", "f8", "f8", "f8"),
                         filling_values=0, deletechars="%")


def read_swift_info(file="swift_countrates.dat"):
    """Read information file used to convert XMM fluxes to Swift count rates

        Parameters
        ----------
        file: str, optional
            Path to the file to be read. Default is swift_countrates.dat
        """
    print("Reading %s data" % file)
    return np.genfromtxt("%s" % file, names=True, delimiter="\t", comments="#",
                         dtype=("U32", "U32", "U32", "U32", "U32"))
