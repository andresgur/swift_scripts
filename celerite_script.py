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
from shutil import copyfile
import corner
import emcee


def read_config_file(config_file):
    """Read config file with parameter initial values and bounds.

    Parameters
    ----------
    config_file:str, the config file"""

    f = open(config_file)
    lines = f.readlines()
    f.close()
    param_info = lines[1].split("\t")
    return param_info


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


ap = argparse.ArgumentParser(description='Look for periodicities running a celerite model on the data (see Foreman-Mackey et al 2017. 10.3847/1538-3881/aa9332)')
ap.add_argument("-m", "--maxcountrate", nargs='?', help="Maximum count rate to consider. Default considers all dataset", type=float, default=np.Infinity)
ap.add_argument("-c", "--config-file", nargs='?', help="Config file with initial parameter constraints", type=str, default="/home/agurpide/scripts/swift_scripts/celerite_config/parameters.config")
ap.add_argument("-s", "--source", nargs='?', help="Source (for the title)", type=str, default="")
ap.add_argument("--tmin", nargs='?', help="Minimum time in MJD swift days (as in the original Swift file)", type=float, default=0)
ap.add_argument("--tmax", nargs='?', help="Maximum time in MJD swift days (as in the original Swift file)", type=float, default=np.Infinity)
ap.add_argument("-i", "--input_observation_file", nargs='?', help="Path to observation file to mark observations", type=str, default="")
ap.add_argument("-o", "--outdir", nargs='?', help="Output dir name", type=str, default="celerite")
ap.add_argument("-f", "--fit", nargs='?', help="Whether to fit the data and start MCMC with best fit values (0 or 1 for False or True). Default 0", type=int, default=0)
ap.add_argument("--swift_dir", nargs='?', help="Directory with all swift converted rates", type=str, default="")
args = ap.parse_args()

fit = False if args.fit == 0 else True

days_to_seconds = 24 * 3600

plt.style.use('/home/agurpide/.config/matplotlib/stylelib/paper.mplstyle')

count_rate_file = "PCCURVE.qdp"
data = ru.readPCCURVE("%s" % count_rate_file)
# filter data according to input parameters
data = np.array([row for row in data if row["Time"] >= args.tmin and row["Time"] <= args.tmax and row["Rate"] < args.maxcountrate])
f = open("t0.date")
lines = f.readlines()
f.close()
start_date = lines[1].split("at")[0]
rate = "Rate"
y = data["%s" % rate]
time = data["Time"]
yerr = (-data["%sneg" % rate] + data["%spos" % rate]) / 2
# read and set initial parameters and bounds
params = read_config_file(args.config_file)
w = np.log(2 * np.pi / (np.array(params[0].split(":")).astype(float) * days_to_seconds))
S_0 = np.log(np.array(params[1].split(":")).astype(float))
Q = np.log(np.array(params[2].split(":")).astype(float))
bounds = dict(log_S0=(S_0[0], S_0[2]), log_Q=(Q[0], Q[2]), log_omega0=(w[2], w[0]))
kernel = terms.SHOTerm(log_S0=S_0[1], log_Q=Q[1], log_omega0=w[1], bounds=bounds)
gp = celerite.GP(kernel, mean=np.mean(data["Rate"]))
gp.compute(time, yerr)  # You always need to call compute once.
print("Initial log likelihood: {0}".format(gp.log_likelihood(y)))
initial_params = gp.get_parameter_vector()

if fit:
    # solution contains the information about the fit. .x is the best fit parameters
    solution = minimize(neg_log_like, initial_params, method="L-BFGS-B", bounds=gp.get_parameter_bounds(), args=(y, gp))
    gp.set_parameter_vector(solution.x)
    print(solution)
    best_params = gp.get_parameter_dict()
    period = (2 * np.pi / (np.exp(best_params["kernel:log_omega0"]) * days_to_seconds))
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

init_params = gp.get_parameter_vector()
par_names = list(gp.get_parameter_names())

ndim = len(init_params)
# set from the example
nwalkers = 20
# start from the best fit values
initial_samples = gp.get_parameter_vector() + 1e-5 * np.random.randn(nwalkers, ndim)
sampler = emcee.EnsembleSampler(nwalkers, ndim, log_probability, args=(y, gp))
# burn in phase
nsamples = 5000
print("Simulating %s samples" % nsamples)
sampler.run_mcmc(initial_samples, nsamples, progress=True)
acceptance_ratio = sampler.acceptance_fraction
print("Acceptance ratio: (%)")
print(acceptance_ratio)
tau = sampler.get_autocorr_time()
print("Correlation parameters:")
print(tau)
# plot the entire chain
chain = sampler.get_chain(flat=True)
chain_fig, axes = plt.subplots(3, sharex=True, gridspec_kw={'hspace': 0, 'wspace': 0})
means = np.mean(chain, axis=0)
for param, parname, ax, mean in zip(chain.T, par_names, axes, means):
    ax.plot(param, linestyle="None", marker="+", color="black")
    ax.set_ylabel(parname)
    ax.axhline(y=mean)
# discard 1000 samples (found empirically) and plot
discard = 1000
for chain, parname, ax, mean in zip(chain.T, par_names, axes, means):
    ax.axvline(discard * nwalkers, ls="--", color="red")
final_samples = sampler.get_chain(discard=discard, thin=20, flat=True)
chain_fig.savefig("%s/chain_samples.png" % args.outdir)

cols = ["log_omega0", "log_S0", "log_Q"]
inds = [par_names.index("kernel:%s" % c) for c in cols]
samples = final_samples[:, inds]
samples[:, :-1] = np.exp(samples[:, :-1])
2 * np.pi / (np.exp(gp.get_parameter_dict()["kernel:log_omega0"]) * days_to_seconds)
ind_omega = par_names.index("kernel:%s" % "log_omega0")
samples[:, 0] = 2 * np.pi / (samples[:, 0] * days_to_seconds)
medians = np.median(samples, axis=0)
ranges = [(median - 0.9 * median, median + 0.9 * median) for median in medians]
corner_fig = corner.corner(samples, smooth=0.5, range=ranges, labels=[r"$P$ (days)", r"$S_0$", r"Q"], truths=medians)
corner_fig.savefig("%s/corner_fig.png" % args.outdir)

# finally plot final PSD and Model
# MODEL
model_figure, model_ax = plt.subplots()
model_ax.set_xlabel("Days since %s" % start_date)
model_ax.set_ylabel("Count-rate ct/s [0.3 - 10 keV]")
model_ax.errorbar(time / (days_to_seconds), y, yerr=yerr, fmt=".k", capsize=0)
# PSD
psd_figure, psd_ax = plt.subplots()
psd_ax.set_xlabel("Frequency [Days$^{-1}$]")
psd_ax.set_ylabel("Power")
t_samples = 10000
t = np.linspace(np.min(time), np.max(time), t_samples)
min_f = np.exp(w[0])
max_f = np.exp(w[2])
frequencies = np.linspace(min_f, max_f, 5000)
color = "orange"
# draw 1000 samples from the final distributions
n_samples = 1000
psds = np.empty((n_samples, len(frequencies)))
models = np.empty((n_samples, t_samples))

for index, sample in enumerate(final_samples[np.random.randint(len(samples), size=n_samples)]):
    gp.set_parameter_vector(sample)
    psd = gp.kernel.get_psd(frequencies)
    model = gp.predict(y, t, return_cov=False)
    model_ax.plot(t / (days_to_seconds), model, color=color, alpha=0.3)
    psd_ax.plot(frequencies * days_to_seconds / 2 / np.pi, psd, color=color, alpha=0.3)
    models[index] = model
    psds[index] = psd
model_figure.savefig("%s/model_fit_celerite_samples.png" % args.outdir)
psd_figure.savefig("%s/psd_samples.png" % args.outdir)
# Median and standard deviation figures
model_figure, model_ax = plt.subplots()
model_ax.set_xlabel("Days since %s" % start_date)
model_ax.set_ylabel("Count-rate ct/s [0.3 - 10 keV]")
model_ax.errorbar(time / days_to_seconds, y, yerr=yerr, fmt=".k", capsize=0)
# PSD
psd_figure, psd_ax = plt.subplots()
psd_ax.set_xlabel("Frequency [Days$^{-1}$]")
psd_ax.set_ylabel("Power")

p = np.percentile(psds, [32, 50, 68], axis=0)
m = np.percentile(models, [32, 50, 68], axis=0)
model_ax.plot(t / days_to_seconds, m[1], color=color)
model_ax.fill_between(t / days_to_seconds, m[0], m[2], alpha=0.3, color=color)
model_ax.errorbar(time / days_to_seconds, y, yerr=yerr, fmt=".k", capsize=0)
psd_ax.plot(frequencies * days_to_seconds / 2 / np.pi, p[1], color=color)
psd_ax.fill_between(frequencies * days_to_seconds / 2 / np.pi, p[0], p[2], color=color, alpha=0.3)
model_figure.savefig("%s/model_fit_median.png" % args.outdir)
psd_figure.savefig("%s/psd_median.png" % args.outdir)
config_file_name = os.path.basename(args.config_file)
copyfile(args.config_file, "%s/%s/" % (args.outdir, config_file_name))
print("Results stored to %s" % args.outdir)
