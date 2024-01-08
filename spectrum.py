# -*- coding: utf-8 -*-

import logging
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.ticker import MaxNLocator, NullLocator
from matplotlib.colors import LinearSegmentedColormap, colorConverter
from matplotlib.ticker import ScalarFormatter
from scipy.stats import norm
import matplotlib
from input import *
import itertools

try:
    from scipy.ndimage import gaussian_filter
except ImportError:
    gaussian_filter = None

__all__ = ["corner", "hist2d", "quantile"]

matplotlib.rc('text', usetex=True)
matplotlib.rc('font', family='serif')


def spec(xs, xfull, xbin, yfull, ybin, ydata, x_err, y_err, range2, line_color, bins=20, range=None, weights=None, color="k",
         smooth=None, smooth1d=None,
         labels=None, label_kwargs=None,
         show_titles=False, title_fmt=".2f", title_kwargs=None,
         truths=None, truth_color="k",
         #ini_guess = [None]*(len(parameters)-2) + [(rstar,rstar_uncertainty), (g,g_uncertainty)],
         density=True,
         scale_hist=False, quantiles=None, verbose=False, fig=None,
         max_n_ticks=5, top_ticks=False, use_math_text=False, reverse=False,
         hist_kwargs=None, **hist2d_kwargs):
           

    plt.rc('text', usetex=True)
    plt.rc('font', family='serif')

    # Create a new figure if one wasn't provided.
    fig2, ax2 = plt.subplots(figsize=(20, 20))   # fig, axes = plt.subplots(K, K, figsize=(dim, dim))

    # Format the figure.
#    lb = lbdim / dim
#    tr = (lbdim + plotdim) / dim
#    fig.subplots_adjust(left=lb, bottom=lb, right=tr, top=tr, wspace=whspace, hspace=whspace)

    ax2.set_frame_on(True)
    linethick = 1
    line1, = ax2.plot(xfull, yfull, linewidth=linethick, color=line_color, linestyle='-')
    symsize2 = 1
    mew1 = 12
    msize = 12
    elw = 5
    ax2.errorbar(xbin, ybin, xerr=x_err, yerr=y_err, fmt='ks', elinewidth=elw)
    ax2.plot(xbin, ybin, 'ks', mew=mew1, markersize=msize)
    symsize = 4
    ax2.errorbar(xbin, ydata, xerr=x_err, yerr=y_err, fmt='ro', capthick=2, elinewidth=elw, zorder=1000)
    ax2.plot(xbin, ydata, 'ro', mew=mew1, markersize=msize, zorder=1000)
    text_size = 42
    wavelength_min = np.amin(wavelength_bins)
    wavelength_max = np.amax(wavelength_bins)
    transit_min = np.amin(transit_depth)
    transit_max = np.amax(transit_depth)
    if wavelength_min < 1.0:
        ax2.text(0.82*wavelength_min + 0.1*wavelength_max, 2.7*transit_max - 1.7*transit_min, r'\textbf{' + planet_name + r'} \textbf{data}', color='r', fontsize=text_size)
        ax2.text(0.82*wavelength_min + 0.1*wavelength_max, 2.47*transit_max - 1.47*transit_min, r'\textbf{Model (binned)}', color='k', fontsize=text_size)    
    else:
        ax2.text(0.855*wavelength_min + 0.1*wavelength_max, 2.7*transit_max - 1.7*transit_min, r'\textbf{' + planet_name + r'} \textbf{data}', color='r', fontsize=text_size)
        ax2.text(0.855*wavelength_min + 0.1*wavelength_max, 2.47*transit_max - 1.47*transit_min, r'\textbf{Model (binned)}', color='k', fontsize=text_size)
    ax2.set_xlim([wavelength_min - 0.03, wavelength_max + 0.03])
    ax2.set_ylim([1.5*transit_min - 0.5*transit_max, 3*transit_max - 2*transit_min])
    ax2.xaxis.set_major_locator(MaxNLocator(5, prune="lower"))
    ax2.yaxis.set_major_locator(MaxNLocator(5, prune="lower"))
    ax2.xaxis.set_major_formatter(ScalarFormatter(useMathText=use_math_text))
    ax2.yaxis.set_major_formatter(ScalarFormatter(useMathText=use_math_text))
    tick_label_size = 42
    tick_size = 20
    ax2.tick_params(axis='both', which='major', direction='in', length=tick_size, labelsize=tick_label_size, pad=15)
    label_size = 48
    ax2.set_xlabel(r'Wavelength ($\mu$m)', fontsize=label_size, labelpad=15)
    ax2.set_ylabel(r'($R$/$R_{\rm star}$)$^2$ (\%)', fontsize=label_size, labelpad=20)
    #ax2.set_xlabel(r'\textbf{Wavelength (}\boldmath{$\mu$}\textbf{m)}', fontsize=label_size, labelpad=100)
    #ax2.set_ylabel(r'\textbf{(}\boldmath{$R$/$R_{\rm star}$}\textbf{)}\boldmath{$^2$} \textbf{ (\%)}', fontsize=label_size, labelpad=100)
    #ax2.xaxis.set_label_coords(0.5, -0.08)
    #y_label_x = -0.25 + 0.06*K/3
    #ax2.yaxis.set_label_coords(y_label_x, 0.5)

    return fig2
