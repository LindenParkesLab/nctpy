import numpy as np
import scipy as sp
import seaborn as sns

def regplot(x, y, xlabel, ylabel, ax, c='gray'):
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
    try:
        sns.kdeplot(x=x, y=y, ax=ax, color='gray', thresh=0.05, alpha=0.25)
    except:
        pass
    sns.regplot(x=x, y=y, ax=ax, scatter=False, color=color_blue)
    ax.scatter(x=x, y=y, c=c, s=5, alpha=0.5)
    ax.set_xlabel(xlabel)
    ax.set_ylabel(ylabel, labelpad=-1)
    ax.tick_params(pad=-2.5)
    ax.grid(False)
    sns.despine(right=True, top=True, ax=ax)
    pearson_stats = sp.stats.pearsonr(x, y)
    spearman_stats = sp.stats.spearmanr(x, y)
    textstr = 'r = {:.2f}, p = {:.2f} \nrho = {:.2f}, p = {:.2f}'.format(pearson_stats[0], pearson_stats[1],
                                                                         spearman_stats[0], spearman_stats[1])
    ax.text(0.05, 0.975, textstr, transform=ax.transAxes,
            verticalalignment='top')
