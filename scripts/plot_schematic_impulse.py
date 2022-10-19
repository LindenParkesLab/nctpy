# %% import packages
import os
import numpy as np

# %% plotting
import matplotlib.pyplot as plt
plt.rcParams.update({'font.size': 10})
plt.rcParams['svg.fonttype'] = 'none'

# %% directories
projdir = '/Users/lindenmp/Google-Drive-Penn/work/research_projects/control_package/'
datadir = os.path.join(projdir, 'data')
resultsdir = os.path.join(projdir, 'results')

# %%
matlab_outputs = os.path.join(projdir, 'matlab_outputs')

impulse = 'x1'
# impulse = 'x2'
t = np.loadtxt(os.path.join(matlab_outputs, 't_{0}.csv'.format(impulse)), delimiter=',')

X = np.loadtxt(os.path.join(matlab_outputs, 'X_{0}.csv'.format(impulse)), delimiter=',')
Y = np.loadtxt(os.path.join(matlab_outputs, 'Y_{0}.csv'.format(impulse)), delimiter=',')
xdot = np.loadtxt(os.path.join(matlab_outputs, 'xdot_{0}.csv'.format(impulse)), delimiter=',')
ydot = np.loadtxt(os.path.join(matlab_outputs, 'ydot_{0}.csv'.format(impulse)), delimiter=',')

sol = np.loadtxt(os.path.join(matlab_outputs, 'sol_{0}.csv'.format(impulse)), delimiter=',')
solI = np.loadtxt(os.path.join(matlab_outputs, 'solI_{0}.csv'.format(impulse)), delimiter=',')
solC = np.loadtxt(os.path.join(matlab_outputs, 'solC_{0}.csv'.format(impulse)), delimiter=',')

u = np.loadtxt(os.path.join(matlab_outputs, 'u_{0}.csv'.format(impulse)), delimiter=',')

# %% plot: impulse response
# which_plot = 'uncontrolled'
which_plot = 'controlled'

orange = [245/255, 140/255, 100/255]
blue = [171/255, 201/255, 234/255]
gray = [100/255, 100/255, 100/255]

# vector field
f, ax = plt.subplots(figsize=(1.25, 1.25))
ax.quiver(X, Y, xdot, ydot, color=gray)
if which_plot == 'uncontrolled':
    ax.plot(sol[0, :], sol[1, :], 'k', linewidth=1)
    ax.plot(sol[0, 0], sol[1, 0], 'kx', markersize=5)
    ax.plot(sol[0, -1], sol[1, -1], 'ko', markersize=5)
if which_plot == 'controlled':
    ax.plot(sol[0, :], sol[1, :], 'k', linewidth=1, alpha=0.25)
    ax.plot(sol[0, 0], sol[1, 0], 'kx', markersize=5, alpha=0.25)
    ax.plot(sol[0, -1], sol[1, -1], 'ko', markersize=5, alpha=0.25)

    ax.plot(solI[0, :], solI[1, :], 'k', linewidth=1)
    ax.plot(solI[0, 0], solI[1, 0], 'x', markersize=5, color=blue)
    ax.plot(solI[0, -1], solI[1, -1], 'ko', markersize=5)

    ax.plot([sol[0, 0], solI[0, 0]], [sol[1, 0], solI[1, 0]], color=blue, linewidth=1.5, linestyle=":")
ax.set_xlabel('x1')
ax.set_ylabel('x2')
f.savefig(os.path.join(resultsdir, 'impulse_{0}_{1}_vector_field.svg'.format(impulse, which_plot)),
          dpi=300, bbox_inches='tight', pad_inches=0)
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
    ax.plot(t, solI[0, :], color=orange, linewidth=1)
    if impulse == 'x1':
        ax.plot(t[0], solI[0, 0], 'x', markersize=2.5, color=blue)
        ax.plot([t[0], t[0]], [sol[0, 0], solI[0, 0]], color=blue, linewidth=1.5, linestyle=":")
    else:
        ax.plot(t[0], solI[0, 0], 'x', markersize=2.5, color=orange)
    ax.plot(t[-1], solI[0, -1], 'o', markersize=2.5, color=orange)
ax.set_xlabel('t')
# ax.set_ylabel('activity')
ax.set_ylim([-0.55, 0.55])
f.savefig(os.path.join(resultsdir, 'impulse_{0}_{1}_x1.svg'.format(impulse, which_plot)),
          dpi=300, bbox_inches='tight', pad_inches=0)
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

    ax.plot(t, solI[1, :], color=orange, linewidth=1)
    if impulse == 'x2':
        ax.plot(t[0], solI[1, 0], 'x', markersize=2.5, color=blue)
        ax.plot([t[0], t[0]], [sol[1, 0], solI[1, 0]], color=blue, linewidth=1.5, linestyle=":")
    else:
        ax.plot(t[0], solI[1, 0], 'x', markersize=2.5, color=orange)
    ax.plot(t[-1], solI[1, -1], 'o', markersize=2.5, color=orange)
ax.set_xlabel('t')
# ax.set_ylabel('activity')
ax.set_ylim([-0.55, 0.55])
f.savefig(os.path.join(resultsdir, 'impulse_{0}_{1}_x2.svg'.format(impulse, which_plot)),
          dpi=300, bbox_inches='tight', pad_inches=0)
plt.close()
