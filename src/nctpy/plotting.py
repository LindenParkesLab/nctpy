"""Figures used in the protocol paper: scatter plots, null distributions, cortical surfaces and module boundaries.

Needs the optional plotting dependencies: ``pip install 'nctpy[plot]'``.
"""

import inspect

import numpy as np
from scipy import stats

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


def set_plotting_params(format="png"):
    """Set the matplotlib and seaborn style used for the protocol paper's figures.

    Embeds TrueType fonts in PDF and PostScript output, keeps SVG text editable, sets the font size to 10 and
    applies seaborn's ``"white"`` style. Changes matplotlib's global ``rcParams``.

    Parameters
    ----------
    format : str, default "png"
        Default file format for ``savefig`` (``rcParams["savefig.format"]``).
    """
    plt.rcParams["pdf.fonttype"] = 42
    plt.rcParams["ps.fonttype"] = 42
    plt.rcParams["savefig.format"] = format
    plt.rcParams["font.size"] = 10

    plt.rcParams["svg.fonttype"] = "none"
    sns.set_style(style="white")


def reg_plot(x, y, xlabel, ylabel, ax, c="gray", annotate="pearson", regr_line=True, kde=True, fontsize=8):
    """Scatter plot of ``y`` against ``x`` with a regression line, a density contour and a correlation.

    Pairs in which either value is NaN are left out. If ``x`` and ``y`` are both square matrices (e.g. two
    node-by-node matrices), their diagonals are left out too and the remaining entries are plotted against each
    other.

    Parameters
    ----------
    x, y : numpy.ndarray
        Values to plot, both 1-D or both 2-D, of the same shape.
    xlabel, ylabel : str
        Axis labels.
    ax : matplotlib.axes.Axes
        Axes to draw on.
    c : str or numpy.ndarray, default "gray"
        A matplotlib colour for every point, or one value per point (same shape as ``x``), coloured with
        ``viridis``.
    annotate : {"pearson", "spearman", "both"} or tuple, default "pearson"
        Text in the top-left corner: Pearson's r, Spearman's rho, or both, each with its p-value. A tuple
        ``(coef, p)`` prints that coefficient and p-value instead. Anything else prints nothing.
    regr_line : bool, default True
        Draw a linear regression line with its confidence band.
    kde : bool, default True
        Draw a kernel density contour behind the points.
    fontsize : float, default 8
        Font size of the annotation.
    """
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
        print("error: input array dimension mismatch.")

    try:
        x = x[indices]
        y = y[indices]
    except Exception:
        pass

    try:
        c = c[indices]
    except Exception:
        pass

    # kde plot (the flags are compared with True, as since 1.0, rather than tested for truth)
    if kde == True:  # noqa: E712
        try:
            sns.kdeplot(x=x, y=y, ax=ax, color="gray", thresh=0.05, alpha=0.25)
        except Exception:
            pass

    # regression line
    if regr_line == True:  # noqa: E712
        color_blue = sns.color_palette("Set1")[1]
        sns.regplot(x=x, y=y, ax=ax, scatter=False, color=color_blue)

    # scatter plot
    if type(c) is str:
        ax.scatter(x=x, y=y, c=c, s=5, alpha=0.5)
    else:
        ax.scatter(x=x, y=y, c=c, cmap="viridis", s=5, alpha=0.5)

    # axis options
    ax.set_xlabel(xlabel, labelpad=0)
    ax.set_ylabel(ylabel, labelpad=0)
    # ax.tick_params(pad=-2.5)
    # ax.grid(False)
    # sns.despine(right=True, top=True, ax=ax)
    sns.despine(offset=0, trim=False, left=False, right=True, top=True, bottom=False, ax=ax)
    ax.tick_params(left=True, bottom=True)

    # annotation
    r, r_p = stats.pearsonr(x, y)
    rho, rho_p = stats.spearmanr(x, y)
    if type(annotate) is str:
        if annotate == "pearson":
            textstr = r"$\mathit{:}$ = {:.2f}, {:}".format("{r}", r, get_p_val_string(r_p))
            ax.text(0.05, 0.975, textstr, transform=ax.transAxes, fontsize=fontsize, verticalalignment="top")
        elif annotate == "spearman":
            textstr = "$\\rho$ = {:.2f}, {:}".format(rho, get_p_val_string(rho_p))
            ax.text(0.05, 0.975, textstr, transform=ax.transAxes, fontsize=fontsize, verticalalignment="top")
        elif annotate == "both":
            textstr = (r"$\mathit{:}$ = {:.2f}, {:}" + "\n" + r"$\rho$ = {:.2f}, {:}").format(
                "{r}", r, get_p_val_string(r_p), rho, get_p_val_string(rho_p)
            )
            ax.text(0.05, 0.975, textstr, transform=ax.transAxes, fontsize=fontsize, verticalalignment="top")
    elif type(annotate) is tuple:
        coef = annotate[0]
        p = annotate[1]
        textstr = "coef = {:.2f}, {:}".format(coef, get_p_val_string(p))
        ax.text(0.05, 0.975, textstr, transform=ax.transAxes, fontsize=fontsize, verticalalignment="top")
    else:
        pass


def null_plot(observed, null, xlabel, ax, p_val=None):
    """Histogram of a null distribution, with the observed value marked.

    Parameters
    ----------
    observed : float
        Observed value of the statistic, drawn as a vertical line and labelled (rounded to an integer).
    null : numpy.ndarray
        Values of the statistic under the null, e.g. one per surrogate network.
    xlabel : str
        Label of the x-axis.
    ax : matplotlib.axes.Axes
        Axes to draw on.
    p_val : float, optional
        p-value to print next to the observed value, e.g. from :func:`nctpy.utils.get_null_p`. Not printed if
        it is None or 0.
    """
    color_blue = sns.color_palette("Set1")[1]
    color_red = sns.color_palette("Set1")[0]
    sns.histplot(x=null, ax=ax, color="gray")
    ax.axvline(x=observed, ymax=1, clip_on=False, linewidth=1, color=color_blue)
    ax.grid(False)
    sns.despine(right=True, top=True, ax=ax)
    ax.set_xlabel(xlabel)
    ax.set_ylabel("counts")

    textstr = "obs. = {:.0f}".format(observed)
    ax.text(
        observed,
        ax.get_ylim()[1],
        textstr,
        horizontalalignment="left",
        verticalalignment="top",
        rotation=270,
        c=color_blue,
    )

    if p_val:
        textstr = "{:}".format(get_p_val_string(p_val))
        ax.text(
            observed - (np.abs(observed) * 0.0025),
            ax.get_ylim()[1],
            textstr,
            horizontalalignment="right",
            verticalalignment="top",
            rotation=270,
            c=color_red,
        )


def roi_to_vtx(roi_data, annot_file):
    """Project one value per parcel onto the vertices of a FreeSurfer annotation.

    Parcel ``k`` (annotation label k >= 1) takes ``roi_data[k - 1]``. Vertices labelled 0 (e.g. the medial
    wall) or -1 (unlabelled) are background and stay 0.

    Parameters
    ----------
    roi_data : numpy.ndarray
        One value per parcel of the hemisphere, in the order of the annotation's labels.
    annot_file : str or os.PathLike
        FreeSurfer annotation file (``.annot``) for the hemisphere.

    Returns
    -------
    vtx_data : numpy.ndarray
        One value per vertex.
    vtx_data_min, vtx_data_max : float
        Minimum and maximum of ``vtx_data`` (both 0 if it is constant).
    """
    labels = nib.freesurfer.read_annot(annot_file)[0]
    vtx_data = np.zeros(labels.shape)
    for i in np.unique(labels[labels > 0]):
        vtx_data[labels == i] = roi_data[i - 1]

    # get min/max for plotting
    x = np.unique(vtx_data)
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
    kwargs = dict(
        hemi=hemi,
        view=view,
        vmin=vmin,
        vmax=vmax,
        bg_map=bg_map,
        bg_on_data=True,
        axes=axes,
        cmap=cmap,
        colorbar=False,
        avg_method="median",
    )
    if "darkness" in inspect.signature(plotting.plot_surf).parameters:
        kwargs["darkness"] = 0.5
    plotting.plot_surf(surf_mesh, surf_map=surf_map, **kwargs)


def surface_plot(data, lh_annot_file, rh_annot_file, fsaverage=None, order="lr", cmap="viridis", cblim=None):
    """Plot one value per parcel on the inflated cortical surface, in lateral and medial views of each hemisphere.

    The first half of ``data`` is drawn on one hemisphere and the second half on the other, through
    :func:`roi_to_vtx`. The figure is shown and returned.

    Parameters
    ----------
    data : numpy.ndarray
        One value per parcel, for both hemispheres.
    lh_annot_file, rh_annot_file : str or os.PathLike
        FreeSurfer annotation files of the left and right hemispheres, on the mesh of ``fsaverage``.
    fsaverage : dict-like, optional
        Surface meshes and sulcal depth maps (``infl_left``, ``infl_right``, ``sulc_left``, ``sulc_right``), as
        returned by :func:`nilearn.datasets.fetch_surf_fsaverage`. Defaults to fsaverage5, loaded when the plot
        is drawn.
    order : {"lr", "rl"}, default "lr"
        Whether ``data`` lists the left hemisphere's parcels first ("lr") or the right's ("rl").
    cmap : str, default "viridis"
        Matplotlib colormap. With ``"coolwarm"`` the colour limits are symmetric about zero.
    cblim : tuple of float, optional
        Colour limits as ``(vmax, vmin)``. Default: the range of ``data`` (see ``cmap``).

    Returns
    -------
    matplotlib.figure.Figure
        A 2 x 2 grid of views with a colorbar.
    """
    # fsaverage5 is loaded when the plot is drawn; until 1.1.0 the default was evaluated on importing nctpy.plotting
    if fsaverage is None:
        fsaverage = datasets.fetch_surf_fsaverage(mesh="fsaverage5")

    # project data to surface
    n_nodes = len(data)
    if order == "lr":
        vtx_data_lh, _, _ = roi_to_vtx(data[: int(n_nodes / 2)], lh_annot_file)
        vtx_data_rh, _, _ = roi_to_vtx(data[int(n_nodes / 2) :], rh_annot_file)
    elif order == "rl":
        vtx_data_lh, _, _ = roi_to_vtx(data[int(n_nodes / 2) :], rh_annot_file)
        vtx_data_rh, _, _ = roi_to_vtx(data[: int(n_nodes / 2)], lh_annot_file)

    # get colorbar axes
    if cblim is None:
        if cmap == "coolwarm":
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
    f, ax = plt.subplots(2, 2, figsize=(2.5, 2.5), subplot_kw={"projection": "3d"})
    _plot_surf_panel(
        fsaverage["infl_left"], vtx_data_lh, fsaverage["sulc_left"], "left", "lateral", vmin, vmax, cmap, ax[0, 0]
    )
    _plot_surf_panel(
        fsaverage["infl_right"], vtx_data_rh, fsaverage["sulc_right"], "right", "lateral", vmin, vmax, cmap, ax[0, 1]
    )
    _plot_surf_panel(
        fsaverage["infl_left"], vtx_data_lh, fsaverage["sulc_left"], "left", "medial", vmin, vmax, cmap, ax[1, 0]
    )
    _plot_surf_panel(
        fsaverage["infl_right"], vtx_data_rh, fsaverage["sulc_right"], "right", "medial", vmin, vmax, cmap, ax[1, 1]
    )

    plt.subplots_adjust(wspace=-0.075, hspace=-0.3)
    cb_ax = f.add_axes([0.9, 0.25, 0.05, 0.5])  # add colorbar
    f.colorbar(im, cax=cb_ax)
    plotting.show()

    return f


def add_module_lines(modules, ax):
    """Draw white boxes around the blocks of a matrix plot that belong to the same module.

    For a node-by-node matrix whose rows and columns are sorted by module, each module's block on the diagonal
    is outlined. Prints the modules found.

    Parameters
    ----------
    modules : pandas.Series
        Module (e.g. functional system) of each node, in the order of the matrix's rows; nodes of the same
        module must be contiguous.
    ax : matplotlib.axes.Axes
        Axes holding the matrix plot, e.g. from :func:`seaborn.heatmap`.
    """
    # get unqiue modules
    unique_modules = modules.unique()
    print(unique_modules)

    previous = -1
    for module in unique_modules:
        # get box boundaries using first and last occurence of module name
        where = np.flatnonzero(np.asarray(modules == module))
        first, last = int(where[0]), int(where[-1])

        # draw box
        ax.hlines(last + 1, previous + 1, last + 1, colors="w")
        ax.vlines(last + 1, previous + 1, last + 1, colors="w")
        ax.hlines(first, previous + 1, last + 1, colors="w")
        ax.vlines(first, previous + 1, last + 1, colors="w")

        # update previous
        previous = last
