import numpy as np


def print_errors(value, lower, upper):
    """Compute errors given value lower and upper bounds"""
    err_high = upper - value
    err_lower = value - lower
    print("%.4f\errors{%.4f}{%.4f}" % (value, err_high, err_lower))
    return


def error_division(a, b, a_err, b_err):
    """Compute the error on the division of two values or arrays (a/b)
    ----------
    Parameters
    a: value or np.array  of the numerator
    b: value or np.array of the denominator
    a_err: error on the numerator (or array)
    b_err: error on the denominator (or array)
    """

    return np.sqrt((a_err / b) ** 2 + (a * b_err / b ** 2) ** 2)


def error_multiplication(a, b, a_err, b_err):
    """Compute the error on the multiplication of two values or arrays (a * b)
    ----------
    Parameters
    a: value or np.array  of the first value
    b: value or np.array of the second value
    a_err: error on a (or array)
    b_err: error on b
    """
    return np.sqrt((b * a_err) ** 2 + (b_err * a) ** 2)


def error_power(value, power_index, errorvalue):
    """Compute the error on power of a value (assuming no error on the power index a**b)
    ----------
    Parameters
    value: value or np.array of a
    power_index: value of the power index (b)
    errorvalue: error on the numerator (or array)
    """
    return power_index * errorvalue * value ** (power_index - 1)


def bounds_to_errors(values, lowerbounds, upperbounds):
    ''' Compute errors given the lower and upper bounds of a an array of values.
    Parameters:
    -----------
    value : the central values given by the fit
    lowerbound : its lower bounds
    upperbound : its upper bounds'''

    lower_errors = values - lowerbounds
    upper_errors = upperbounds - values

    for value, lowerbound, upperbound in zip(values, lowerbounds, upperbounds):
        if upperbound < value and upperbound != 0:
            print("Warning upperbound is lower than value!!! %.5f < %.5f" % (upperbound, value))
        if lowerbound > value and lowerbound != 0:
            print("Warning lowerbound is higher than value!!! %.5f > %.5f" % (lowerbound, value))

    uplims = np.zeros(values.shape)

    # lower bound (upper limit)
    uplims[np.where(lowerbounds == 0)] = 1
    lower_errors[np.where(lowerbounds == 0)] = (upperbounds[np.where(lowerbounds == 0)] - values[np.where(lowerbounds == 0)]) * 0.25
    values[np.where(lowerbounds == 0)] = upperbounds[np.where(lowerbounds == 0)]
    upper_errors[np.where(lowerbounds == 0)] = 0

    # upper bound found (lower limit)
    lolims = np.zeros(values.shape)
    lolims[np.where(upperbounds == 0)] = 1
    upper_errors[np.where(upperbounds == 0)] = (values[np.where(upperbounds == 0)] - lowerbounds[np.where(upperbounds == 0)]) * 0.25
    values[np.where(upperbounds == 0)] = lowerbounds[np.where(upperbounds == 0)]
    lower_errors[np.where(upperbounds == 0)] = 0
    return lower_errors, upper_errors, lolims, uplims
