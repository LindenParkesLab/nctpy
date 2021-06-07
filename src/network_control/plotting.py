from network_control.utils import get_p_val_string
import os, platform
import numpy as np
import scipy as sp

import seaborn as sns
import pkg_resources
import matplotlib as mpl
import matplotlib.pyplot as plt

def set_plotting_params(format='png'):
    if platform.system() == 'Darwin':
        os.system('rm -rf ~/.cache/matplotlib')

    plt.rcParams['pdf.fonttype'] = 42
    plt.rcParams['ps.fonttype'] = 42
    plt.rcParams['savefig.format'] = format

    path = pkg_resources.resource_stream('network_control', 'PublicSans-Thin.ttf')
    prop = mpl.font_manager.FontProperties(fname=path.name)
    plt.rcParams['font.sans-serif'] = prop.get_name()
    plt.rcParams['font.serif'] = prop.get_name()
    plt.rcParams['font.family'] = prop.get_family()
    plt.rcParams['mathtext.fontset'] = 'custom'
    plt.rcParams['mathtext.it'] = 'Public Sans:italic'
    plt.rcParams['mathtext.bf'] = 'Public Sans:bold'
    plt.rcParams['mathtext.cal'] = 'Public Sans'

    plt.rcParams['svg.fonttype'] = 'none'
    sns.set(style='whitegrid', context='paper', font_scale=1, font='Public Sans')


def reg_plot(x, y, xlabel, ylabel, ax, c='gray', add_spearman=False, kdeplot=True, regplot=True):
    if x.shape == y.shape:
        mask_x = ~np.isnan(x)
        mask_y = ~np.isnan(y)
        mask = mask_x * mask_y
        indices = np.where(mask)
    else:
        print('error: input array dimension mismatch.')

    try:
        x = x[indices]
        y = y[indices]
    except:
        pass

    try:
        c = c[indices]
    except:
        pass

    color_blue = sns.color_palette("Set1")[1]
    if kdeplot:
        try:
            sns.kdeplot(x=x, y=y, ax=ax, color='gray', thresh=0.05, alpha=0.25)
        except:
            pass

    if regplot:
        sns.regplot(x=x, y=y, ax=ax, scatter=False, color=color_blue)

    if type(c) == str:
        ax.scatter(x=x, y=y, c=c, s=5, alpha=0.5)
    else:
        ax.scatter(x=x, y=y, c=c, cmap='viridis', s=5, alpha=0.5)
    ax.set_xlabel(xlabel, labelpad=-0.5)
    ax.set_ylabel(ylabel, labelpad=-0.5)
    ax.tick_params(pad=-2.5)
    ax.grid(False)
    sns.despine(right=True, top=True, ax=ax)
    r, r_p = sp.stats.pearsonr(x, y)
    if add_spearman:
        rho, rho_p = sp.stats.spearmanr(x, y)
        textstr = '$\mathit{:}$ = {:.2f}, {:}\n$\\rho$ = {:.2f}, {:}' \
            .format('{r}', r, get_p_val_string(r_p), rho, get_p_val_string(rho_p))
        ax.text(0.05, 0.975, textstr, transform=ax.transAxes,
                verticalalignment='top')
    else:
        textstr = '$\mathit{:}$ = {:.2f}, {:}' \
            .format('{r}', r, get_p_val_string(r_p))
        ax.text(0.05, 0.975, textstr, transform=ax.transAxes,
                verticalalignment='top')
