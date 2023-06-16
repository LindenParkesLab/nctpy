# %% import packages
import os
import numpy as np

# %% plotting
import matplotlib.pyplot as plt
plt.rcParams.update({'font.size': 10})
plt.rcParams['svg.fonttype'] = 'none'

# %% directories
projdir = '/Users/lindenmp/Google-Drive-Rutgers/work/research_projects/nctpy/'
datadir = os.path.join(projdir, 'data')
resultsdir = os.path.join(projdir, 'results')

# %%
matlab_outputs = os.path.join(projdir, 'matlab_scripts')

t = np.loadtxt(os.path.join(matlab_outputs, 't.csv'), delimiter=',')

X = np.loadtxt(os.path.join(matlab_outputs, 'X.csv'), delimiter=',')
Y = np.loadtxt(os.path.join(matlab_outputs, 'Y.csv'), delimiter=',')
xdot = np.loadtxt(os.path.join(matlab_outputs, 'xdot.csv'), delimiter=',')
ydot = np.loadtxt(os.path.join(matlab_outputs, 'ydot.csv'), delimiter=',')

sol = np.loadtxt(os.path.join(matlab_outputs, 'sol.csv'), delimiter=',')
solI = np.loadtxt(os.path.join(matlab_outputs, 'solI.csv'), delimiter=',')
solC = np.loadtxt(os.path.join(matlab_outputs, 'solC.csv'), delimiter=',')

u = np.loadtxt(os.path.join(matlab_outputs, 'u.csv'), delimiter=',')

# %% plot: (un)controlled trajectory
for which_plot in ['uncontrolled', 'controlled']:
    orange = [245/255, 140/255, 100/255]
    blue = [171/255, 201/255, 234/255]
    gray = [100/255, 100/255, 100/255]

    # vector field
    f, ax = plt.subplots(figsize=(1.25, 1.25))
    ax.quiver(X, Y, xdot, ydot, color=gray)
    # uncontrolled trajectory
    if which_plot == 'uncontrolled':
        ax.plot(sol[0, :], sol[1, :], 'k', linewidth=1)
        ax.plot(sol[0, 0], sol[1, 0], 'kx', markersize=5)
        ax.plot(sol[0, -1], sol[1, -1], 'ko', markersize=5)
    if which_plot == 'controlled':
        ax.plot(sol[0, :], sol[1, :], 'k', linewidth=1, alpha=0.25)
        ax.plot(sol[0, 0], sol[1, 0], 'kx', markersize=5, alpha=0.25)
        ax.plot(sol[0, -1], sol[1, -1], 'ko', markersize=5, alpha=0.25)

        ax.plot(solC[0, :], solC[1, :], 'k', linewidth=1)
        ax.plot(solC[0, 0], solC[1, 0], 'kx', markersize=5)
        ax.plot(solC[0, -1], solC[1, -1], 'ko', markersize=5)
    ax.set_xlabel('x1')
    ax.set_ylabel('x2')
    ax.set_xlim([-1, 1])
    ax.set_ylim([-1, 1])
    f.savefig(os.path.join(resultsdir, which_plot+'_vector_field.svg'), dpi=300, bbox_inches='tight', pad_inches=0)
    plt.close()

    # time series, x1
    f, ax = plt.subplots(figsize=(1, 0.75))
    if which_plot == 'uncontrolled':
        ax.plot(t, sol[0, :], color=orange, linewidth=1)
        ax.plot(t[0], sol[0, 0], 'x', markersize=2.5, color=orange)
        ax.plot(t[-1], sol[0, -1], 'o', markersize=2.5, color=orange)
    if which_plot == 'controlled':
        ax.plot(t, sol[0, :], color=orange, linewidth=1, alpha=0.25)
        ax.plot(t[0], sol[0, 0], 'x', markersize=2.5, color=orange, alpha=0.25)
        ax.plot(t[-1], sol[0, -1], 'o', markersize=2.5, color=orange, alpha=0.25)

        ax.plot(t, solC[0, :], color=orange, linewidth=1)
        ax.plot(t[0], solC[0, 0], 'x', markersize=2.5, color=orange)
        ax.plot(t[-1], solC[0, -1], 'o', markersize=2.5, color=orange)
        # ax.plot(t, u, color=blue, linewidth=1, linestyle="--")
    ax.set_xlabel('t')
    # ax.set_ylabel('activity')
    ax.set_ylim([-1, 1])
    f.savefig(os.path.join(resultsdir, which_plot+'_x1.svg'), dpi=300, bbox_inches='tight', pad_inches=0)
    plt.close()

    # time series, x2
    f, ax = plt.subplots(figsize=(1, 0.75))
    if which_plot == 'uncontrolled':
        ax.plot(t, sol[1, :], color=orange, linewidth=1)
        ax.plot(t[0], sol[1, 0], 'x', markersize=2.5, color=orange)
        ax.plot(t[-1], sol[1, -1], 'o', markersize=2.5, color=orange)
    if which_plot == 'controlled':
        ax.plot(t, sol[1, :], color=orange, linewidth=1, alpha=0.25)
        ax.plot(t[0], sol[1, 0], 'x', markersize=2.5, color=orange, alpha=0.25)
        ax.plot(t[-1], sol[1, -1], 'o', markersize=2.5, color=orange, alpha=0.25)

        ax.plot(t, solC[1, :], color=orange, linewidth=1)
        ax.plot(t[0], solC[1, 0], 'x', markersize=2.5, color=orange)
        ax.plot(t[-1], solC[1, -1], 'o', markersize=2.5, color=orange)
        ax.plot(t, u, color=blue, linewidth=1, linestyle=":")
    ax.set_xlabel('t')
    # ax.set_ylabel('activity')
    ax.set_ylim([-1, 1])
    f.savefig(os.path.join(resultsdir, which_plot+'_x2.svg'), dpi=300, bbox_inches='tight', pad_inches=0)
    plt.close()
