# @Author: Andrés Gúrpide <agurpide>
# @Date:   01-09-2020
# @Email:  agurpidelash@irap.omp.eu
# @Last modified by:   agurpide
# @Last modified time: 22-09-2020



import numpy as np
import os
import matplotlib.pyplot as plt
import argparse
import celerite
import readingutils as ru
from celerite import terms
from scipy.optimize import minimize
import corner
import emcee


def neg_log_like(params, y, gp):
    gp.set_parameter_vector(params)
    return -gp.log_likelihood(y)


def log_probability(params, y, gp):
    gp.set_parameter_vector(params)

    lp = gp.log_prior()
    if not np.isfinite(lp):
        return -np.inf

    ll = gp.log_likelihood(y)
    return ll + lp if np.isfinite(ll) else -np.inf


ap = argparse.ArgumentParser(description='Look for periodicities running a celerite model on the data (see Foreman-Mackey et al 2017. 10.3847/1538-3881/aa9332) ')
ap.add_argument("-m", "--maxcountrate", nargs='?', help="Maximum count rate to consider. Default considers all dataset", type=float, default=np.Infinity)
ap.add_argument("-s", "--source", nargs='?', help="Source (for the title)", type=str, default="")
ap.add_argument("--tmin", nargs='?', help="Minimum time in MJD swift days (as in the original Swift file)", type=float, default=0)
ap.add_argument("--tmax", nargs='?', help="Maximum time in MJD swift days (as in the original Swift file)", type=float, default=np.Infinity)
ap.add_argument("-i", "--input_observation_file", nargs='?', help="Path to observation file to mark observations", type=str, default="")
ap.add_argument("-p", "--period", nargs='?', help="Period (initial guess for the fit)", type=float, default=100)
ap.add_argument("-o", "--outdir", nargs='?', help="Output dir name", type=str, default="celerite")
ap.add_argument("--swift_dir", nargs='?', help="Directory with all swift converted rates", type=str, default="")

args = ap.parse_args()
fit = False

plt.style.use('/home/agurpide/.config/matplotlib/stylelib/paper.mplstyle')
count_rate_file = "PCCURVE.qdp"
data = ru.readPCCURVE("%s" % count_rate_file)
data = np.array([row for row in data if row["Time"] >= args.tmin and row["Time"] <= args.tmax and row["Rate"] < args.maxcountrate])
f = open("t0.date")
lines = f.readlines()
f.close()
start_date = lines[1].split("at")[0]
freeze_omega = False
rate = "Rate"
y = data["%s" % rate]

time = data["Time"]
yerr = (-data["%sneg" % rate] + data["%spos" % rate]) / 2
# A periodic component (initial parameters)
Q = 10000.05
w0 = 2 * np.pi / (args.period * 24 * 3600)
S0 = np.var(y) / (w0 * Q)
w_bounds = (np.log(2 * np.pi / (250 * 24 * 3600)), np.log(2 * np.pi / (20 * 24 * 3600)))
Q_bounds = (np.log(3), 15)
# Q_bounds = (np.log(10000), np.log(10000.1))
bounds = dict(log_S0=(np.log(0.0000001), np.log(100000)), log_Q=Q_bounds, log_omega0=w_bounds)
kernel = terms.SHOTerm(log_S0=np.log(S0), log_Q=np.log(Q), log_omega0=np.log(w0),
                       bounds=bounds)

#kernel.freeze_parameter("log_Q")

gp = celerite.GP(kernel, mean=np.mean(data["Rate"]))
gp.compute(time, yerr)  # You always need to call compute once.
print("Initial log likelihood: {0}".format(gp.log_likelihood(y)))

initial_params = gp.get_parameter_vector()
bounds = gp.get_parameter_bounds()
if fit:
    # solution contains the information about the fit. .x is the best fit parameters
    solution = minimize(neg_log_like, initial_params, method="L-BFGS-B", bounds=bounds, args=(y, gp))
    gp.set_parameter_vector(solution.x)
    print(solution)
    best_params = gp.get_parameter_dict()

    period = (2 * np.pi / (np.exp(best_params["kernel:log_omega0"]) * 3600 * 24))
    q_factor = np.exp(best_params["kernel:log_Q"])
    s_factor = np.exp(best_params["kernel:log_S0"])
    print("Best fit parameters")
    print("Period: %.2f" % period)
    print("Q: %.2f" % q_factor)
    print("S: %.2f" % s_factor)
    celerite_figure, (ax1, ax2) = plt.subplots(2, 1, sharex=True)
    samples = np.linspace(np.min(time), np.max(time), 10000)
    pred_mean, pred_var = gp.predict(y, samples, return_var=True)
    pred_std = np.sqrt(pred_var)
    color = "#ff7f0e"
    ax1.errorbar(time / 3600 / 24, y, yerr=yerr, fmt=".k", capsize=0)
    ax1.set_ylabel("Count-rate ct/s [0.3 - 10 keV]")
    ax1.plot(samples / 3600 / 24, pred_mean, color=color)
    ax1.fill_between(samples / 3600 / 24, pred_mean + pred_std, pred_mean - pred_std, color=color, alpha=0.3,
                     edgecolor="none")
    pred_mean, pred_var = gp.predict(y, time, return_var=True)
    ax2.errorbar(time / 3600 / 24, (y - pred_mean) / yerr, yerr=1, fmt=".k", capsize=0)
    ax2.axhline(y=0, ls="--", color="#afdbf5")
    ax2.set_ylabel("Residuals")
    ax2.set_xlabel("Days since %s" % start_date)
    if not os.path.isdir("%s" % args.outdir):
        os.mkdir("%s" % args.outdir)
    celerite_figure.savefig("%s/fit.png" % args.outdir)
# monte carlo part
# dimension = number of best fit parameters
ndim = len(gp.get_parameter_vector())
# set from the example
nwalkers = 32
# start from the best fit values
initial_samples = gp.get_parameter_vector() + 1e-4 * np.random.randn(nwalkers, ndim)
sampler = emcee.EnsembleSampler(nwalkers, ndim, log_probability, args=(y, gp))
# burn in phase
sampler.run_mcmc(initial_samples, 5000)
final_samples = sampler.get_chain(discard=500, thin=50, flat=True)
# transform names to human
names = list(gp.get_parameter_names())
cols = ["log_omega0", "log_S0", "log_Q"]
inds = [names.index("kernel:%s" % name) for name in cols]
samples = final_samples[:, inds]
samples[:, :-1] = np.exp(samples[:, :-1])
samples[:, 0] = 2 * np.pi / (samples[:, 0] * 3600 * 24)
best_params = gp.get_parameter_dict()
truths = np.array([best_params["kernel:%s" % k] for k in cols])
truths[:-1] = np.exp(truths[:-1])
inds = [names.index("%s" % name) for name in names]
truths[0] = period
corner_fig = corner.corner(samples, truths=truths, smooth=0.5,
                           labels=[r"$P$ (days)", r"$S_0$", r"Q"], show_titles=True)
plt.show()

corner_fig.savefig("%s/corner_fig.png" % args.outdir)
