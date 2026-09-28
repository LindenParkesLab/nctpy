import os, sys, platform
import inspect
import numpy as np
import scipy as sp

try:
    import nibabel as nib

    import seaborn as sns
    import matplotlib.pyplot as plt
    from nilearn import datasets
    from nilearn import plotting
except ImportError as exc:
    raise ImportError(
        "nctpy.plotting needs the optional plotting dependencies (matplotlib, seaborn, nibabel, nilearn), "
        "which are not installed with nctpy by default. Install them with:\n"
        "    pip install 'nctpy[plot]'\n"
        "or, to run the code printed in the Nature Protocols paper:\n"
        "    pip install 'nctpy[paper]'"
    ) from exc

from nctpy.utils import get_p_val_string


def set_plotting_params(format='png'):
    if platform.system == 'darwin':
        os.system('rm -rf ~/.cache/matplotlib')
    plt.rcParams['pdf.fonttype'] = 42
    plt.rcParams['ps.fonttype'] = 42
    plt.rcParams['savefig.format'] = format
    plt.rcParams['font.size'] = 10

    plt.rcParams['svg.fonttype'] = 'none'
    sns.set_style(style='white')


def reg_plot(x, y, xlabel, ylabel, ax, c='gray', annotate='pearson', regr_line=True, kde=True, fontsize=8):
    if len(x.shape) > 1 and len(y.shape) > 1:
        if x.shape[0] == x.shape[1] and y.shape[0] == y.shape[1]:
            mask_x = ~np.eye(x.shape[0], dtype=bool) * ~np.isnan(x)
            mask_y = ~np.eye(y.shape[0], dtype=bool) * ~np.isnan(y)
            mask = mask_x * mask_y
            indices = np.where(mask)
        else:
            mask_x = ~np.isnan(x)
            mask_y = ~np.isnan(y)
            mask = mask_x * mask_y
            indices = np.where(mask)
    elif len(x.shape) == 1 and len(y.shape) == 1:
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

    # kde plot
    if kde == True:
        try:
            sns.kdeplot(x=x, y=y, ax=ax, color='gray', thresh=0.05, alpha=0.25)
        except:
            pass

    # regression line
    if regr_line == True:
        color_blue = sns.color_palette("Set1")[1]
        sns.regplot(x=x, y=y, ax=ax, scatter=False, color=color_blue)

    # scatter plot
    if type(c) == str:
        ax.scatter(x=x, y=y, c=c, s=5, alpha=0.5)
    else:
        ax.scatter(x=x, y=y, c=c, cmap='viridis', s=5, alpha=0.5)

    # axis options
    ax.set_xlabel(xlabel, labelpad=0)
    ax.set_ylabel(ylabel, labelpad=0)
    # ax.tick_params(pad=-2.5)
    # ax.grid(False)
    # sns.despine(right=True, top=True, ax=ax)
    sns.despine(offset=0, trim=False, left=False, right=True, top=True, bottom=False, ax=ax)
    ax.tick_params(left=True, bottom=True)

    # annotation
    r, r_p = sp.stats.pearsonr(x, y)
    rho, rho_p = sp.stats.spearmanr(x, y)
    if type(annotate) == str:
        if annotate == 'pearson':
            textstr = '$\mathit{:}$ = {:.2f}, {:}'.format('{r}', r, get_p_val_string(r_p))
            ax.text(0.05, 0.975, textstr, transform=ax.transAxes, fontsize=fontsize,
                    verticalalignment='top')
        elif annotate == 'spearman':
            textstr = '$\\rho$ = {:.2f}, {:}'.format(rho, get_p_val_string(rho_p))
            ax.text(0.05, 0.975, textstr, transform=ax.transAxes, fontsize=fontsize,
                    verticalalignment='top')
        elif annotate == 'both':
            textstr = '$\mathit{:}$ = {:.2f}, {:}\n$\\rho$ = {:.2f}, {:}'.format('{r}', r, get_p_val_string(r_p),
                                                                                 rho, get_p_val_string(rho_p))
            ax.text(0.05, 0.975, textstr, transform=ax.transAxes, fontsize=fontsize,
                    verticalalignment='top')
    elif type(annotate) == tuple:
        coef = annotate[0]
        p = annotate[1]
        textstr = 'coef = {:.2f}, {:}'.format(coef, get_p_val_string(p))
        ax.text(0.05, 0.975, textstr, transform=ax.transAxes, fontsize=fontsize, verticalalignment='top')
    else:
        pass


def null_plot(observed, null, xlabel, ax, p_val=None):
    color_blue = sns.color_palette("Set1")[1]
    color_red = sns.color_palette("Set1")[0]
    sns.histplot(x=null, ax=ax, color='gray')
    ax.axvline(x=observed, ymax=1, clip_on=False, linewidth=1, color=color_blue)
    ax.grid(False)
    sns.despine(right=True, top=True, ax=ax)
    ax.set_xlabel(xlabel)
    ax.set_ylabel('counts')

    textstr = 'obs. = {:.0f}'.format(observed)
    ax.text(observed, ax.get_ylim()[1], textstr,
            horizontalalignment='left', verticalalignment='top',
            rotation=270, c=color_blue)

    if p_val:
        textstr = '{:}'.format(get_p_val_string(p_val))
        ax.text(observed - (np.abs(observed)*0.0025), ax.get_ylim()[1], textstr,
                horizontalalignment='right', verticalalignment='top',
                rotation=270, c=color_red)


def roi_to_vtx(roi_data, annot_file):
    labels, ctab, surf_names = nib.freesurfer.read_annot(annot_file)
    vtx_data = np.zeros(labels.shape)

    unique_labels = np.unique(labels)
    if unique_labels[0] == 0:
        unique_labels = unique_labels[1:]

    for i in unique_labels:
        vtx_data[labels == i] = roi_data[i - 1]

    # get min/max for plottin
    x = np.sort(np.unique(vtx_data))

    if x.shape[0] > 1:
        vtx_data_min = x[0]
        vtx_data_max = x[-1]
    else:
        vtx_data_min = 0
        vtx_data_max = 0

    return vtx_data, vtx_data_min, vtx_data_max


def _plot_surf_panel(surf_mesh, surf_map, bg_map, hemi, view, vmin, vmax, cmap, axes):
    # The data are continuous and may be negative, so this uses plot_surf rather than plot_surf_roi, which is
    # for integer label maps and rejects anything else from nilearn 0.13. avg_method='median' is what
    # plot_surf_roi used, so figures are unchanged where the old call worked. darkness was removed in nilearn 0.14.
    kwargs = dict(hemi=hemi, view=view, vmin=vmin, vmax=vmax, bg_map=bg_map, bg_on_data=True, axes=axes,
                  cmap=cmap, colorbar=False, avg_method='median')
    if 'darkness' in inspect.signature(plotting.plot_surf).parameters:
        kwargs['darkness'] = .5
    plotting.plot_surf(surf_mesh, surf_map=surf_map, **kwargs)


def surface_plot(data, lh_annot_file, rh_annot_file,
                 fsaverage=datasets.fetch_surf_fsaverage(mesh='fsaverage5'),
                 order='lr', cmap='viridis', cblim=None):

    # project data to surface
    n_nodes = len(data)
    if order == 'lr':
        vtx_data_lh, _, _ = roi_to_vtx(data[:int(n_nodes/2)], lh_annot_file)
        vtx_data_rh, _, _ = roi_to_vtx(data[int(n_nodes/2):], rh_annot_file)
    elif order == 'rl':
        vtx_data_lh, _, _ = roi_to_vtx(data[int(n_nodes/2):], rh_annot_file)
        vtx_data_rh, _, _ = roi_to_vtx(data[:int(n_nodes/2)], lh_annot_file)

    # get colorbar axes
    if cblim is None:
        if cmap == 'coolwarm':
            vmax = np.round(np.nanmax(np.abs(data)), 1)
            vmin = -vmax
        else:
            vmax = np.nanmax(data)
            vmin = np.nanmin(data)
    else:
        vmax = cblim[0]
        vmin = cblim[1]

    # dummy plot for colorbar
    im = plt.imshow(np.random.random((2, 2)), cmap=cmap, vmin=vmin, vmax=vmax)
    plt.close()

    # main plot
    f, ax = plt.subplots(2, 2, figsize=(2.5, 2.5), subplot_kw={'projection': '3d'})
    _plot_surf_panel(fsaverage['infl_left'], vtx_data_lh, fsaverage['sulc_left'], 'left', 'lateral',
                     vmin, vmax, cmap, ax[0, 0])
    _plot_surf_panel(fsaverage['infl_right'], vtx_data_rh, fsaverage['sulc_right'], 'right', 'lateral',
                     vmin, vmax, cmap, ax[0, 1])
    _plot_surf_panel(fsaverage['infl_left'], vtx_data_lh, fsaverage['sulc_left'], 'left', 'medial',
                     vmin, vmax, cmap, ax[1, 0])
    _plot_surf_panel(fsaverage['infl_right'], vtx_data_rh, fsaverage['sulc_right'], 'right', 'medial',
                     vmin, vmax, cmap, ax[1, 1])

    plt.subplots_adjust(wspace=-0.075, hspace=-0.3)
    cb_ax = f.add_axes([0.9, 0.25, 0.05, 0.5])  # add colorbar
    f.colorbar(im, cax=cb_ax)
    plotting.show()

    return f


def add_module_lines(modules, ax):

    # get unqiue modules
    unique_modules = modules.unique()
    print(unique_modules)

    previous = -1
    for i in np.arange(len(unique_modules)):

        # get box boundaries using first and last occurence of module name
        bool_array = np.asarray(modules == unique_modules[i])
        n = len(bool_array)
        first = -1
        last = -1
        for i in range(0, n):
            if (bool_array[i] != True):
                continue
            if (first == -1):
                first = i
            last = i

        # draw box
        ax.hlines(last + 1, previous + 1, last + 1, colors='w')
        ax.vlines(last + 1, previous + 1, last + 1, colors='w')
        ax.hlines(first, previous + 1, last + 1, colors='w')
        ax.vlines(first, previous + 1, last + 1, colors='w')

        # update previous
        previous = last
