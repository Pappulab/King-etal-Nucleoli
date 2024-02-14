import os
import numpy as np
from typing import Union, Dict, Tuple, Any, List, Iterable
from mpl_toolkits.axes_grid1 import make_axes_locatable
import matplotlib.pyplot as plt
from matplotlib.cm import get_cmap
from matplotlib.colors import LinearSegmentedColormap
import matplotlib
import scipy.stats as st
import pandas as pd


matplotlib.rcParams['pdf.fonttype'] = 42
matplotlib.rcParams['ps.fonttype'] = 42



def plot_calculated_pH(all_calculated_pHs: Iterable[float],
                       savename: str,
                       curves_used: str,
                       vmin: float = None,
                       vmax: float = None) -> None:
    filename = savename.replace('_pH.pdf', '')
    basename = os.path.basename(filename)

    fig = plt.figure(figsize=(10, 10))

    fontsize=16
    fig.suptitle(basename)
    ax = fig.add_subplot(121)
    ax.set_title(f'Median calculated pH ({curves_used})', fontsize=fontsize)
    ax.set_xlabel('X', fontsize=fontsize)
    ax.set_ylabel('Y', fontsize=fontsize)

    median = np.median(all_calculated_pHs, axis=0)
    variance = np.var(all_calculated_pHs, axis=0)

    cmap = 'viridis'
    im = ax.imshow(median, cmap=cmap, vmin=vmin, vmax=vmax)
    divider = make_axes_locatable(ax)
    cax = divider.append_axes('right', size='5%', pad=0.1)
    colorbar = fig.colorbar(im, cax=cax)
    colorbar.ax.tick_params(axis='both', which='major', labelsize=fontsize * 0.8)
    colorbar.ax.tick_params(axis='both', which='minor', labelsize=fontsize * 0.8)
    colorbar.ax.set_ylabel('Calculated pH', fontsize=fontsize)
    [x.set_linewidth(3) for x in ax.spines.values()]
    ax.tick_params(axis='both', which='major', labelsize=fontsize)
    ax.tick_params(axis='both', which='minor', labelsize=fontsize)

    ax = fig.add_subplot(122)
    ax.set_xlabel('X', fontsize=fontsize)
    ax.set_ylabel('Y', fontsize=fontsize)
    ax.set_title('Variance of Calculated pH', fontsize=fontsize)
    im = ax.imshow(variance, cmap=cmap)
    divider = make_axes_locatable(ax)
    cax = divider.append_axes('right', size='5%', pad=0.1)
    colorbar = fig.colorbar(im, cax=cax)
    colorbar.ax.tick_params(axis='both', which='major', labelsize=fontsize * 0.8)
    colorbar.ax.tick_params(axis='both', which='minor', labelsize=fontsize * 0.8)
    colorbar.ax.set_ylabel('$\sigma^2$', fontsize=fontsize)
    [x.set_linewidth(3) for x in ax.spines.values()]
    ax.tick_params(axis='both', which='major', labelsize=fontsize)
    ax.tick_params(axis='both', which='minor', labelsize=fontsize)
    plt.subplots_adjust(wspace=0.5)
    fig.savefig(savename)

    plt.show()


def plot_filtered_calibration_curves(actual_pHs,
                                     raw_data,
                                     fits_x,
                                     fits_y,
                                     savename,
                                     colors):
    fig = plt.figure(figsize=(12, 10))
    ax = fig.add_subplot(111)
    ax.set_xlabel('pH')
    ax.set_ylabel('$Em_i / Em_j$')
    for idx, ratio in enumerate(fits_y):
        fitted_curve = fits_y[ratio]
        raw_data_curve = raw_data[ratio]
        ax.plot(actual_pHs, raw_data_curve, 'o', markersize=5, color=colors[idx])
        ax.plot(fits_x, fitted_curve, linestyle='-', linewidth=2, color=colors[idx], label=f'{ratio}')

    legend = fig.legend(loc='center right', bbox_to_anchor=(1, 0.5))
    fig.subplots_adjust(right=0.85)

    [x.set_linewidth(3) for x in ax.spines.values()]
    fig.savefig(savename)
    plt.show()


def plot_identified_fits(fits, all_data_wavelengths, savename, cmap_name = 'rainbow'):
    # https://stackoverflow.com/questions/56654952/how-to-draw-cell-borders-in-imshow
    def highlight_cell(x,y, ax=None, **kwargs):
        rect = plt.Rectangle((x-.5, y-.5), 1, 1, fill=False, **kwargs)
        ax = ax or plt.gca()
        ax.add_patch(rect)
        return ax

    cmap = get_cmap(cmap_name)  # tab20b
    N = (np.ceil(len(all_data_wavelengths) / 10) * 10).astype(int)
    colors = cmap(np.linspace(0, 1, N))
    new_cmap = LinearSegmentedColormap.from_list('custom_colormap', colors)

    num_wavelengths = len(all_data_wavelengths)
    blank = np.zeros((num_wavelengths, num_wavelengths), dtype=float)
    blank += np.nan
    fig = plt.figure(figsize=(8, 6))
    ax = fig.add_subplot(111)
    indices_i = list()
    indices_j = list()
    for index, ratio in enumerate(fits):
        a, b = ratio
        i = all_data_wavelengths.index(a)
        j = all_data_wavelengths.index(b)
        indices_i.append(i)
        indices_j.append(j)
        blank[i, j] = index

    labels = ['%.1f' % w for w in all_data_wavelengths]
    ticks = np.arange(len(all_data_wavelengths))

    im = ax.imshow(blank, cmap=new_cmap, zorder=-1)
    ax.set_xlabel(r'Intensity measurement at wavelength $\lambda$')
    ax.set_ylabel(r'Intensity measurement at wavelength $\lambda$')
    ax.set_xticks(ticks, labels, rotation=45)
    ax.set_yticks(ticks, labels, rotation=45)

    for i, j in zip(indices_i, indices_j):
        ax = highlight_cell(j, i, ax=ax, color='black', linewidth=2)
    [x.set_linewidth(3) for x in ax.spines.values()]
    fig.savefig(savename)
    plt.show()

    colors = im.cmap(im.norm(np.arange(len(fits))))
    return colors



def violin_plot_of_dilute_patches(df,
                                  all_patches,
                                  violin_labels,
                                  savename):
    import seaborn as sns
    sns.set_style('whitegrid')

    # Wrap around index function
    def wrap_around_index(lst, index):
        length = len(lst)
        if length == 0:
            return None
        else:
            return lst[index % length]

    filename = savename.replace('_patches.pdf', '')
    basename = os.path.basename(filename)

    fig = plt.figure(figsize=(12, 10))
    ax = fig.add_subplot(111)
    ax.set_title(f'Calculated pH of dilute patches ({basename})')
    ax.set_ylabel('Calculated pH')
    ax.set_xlabel('Calibration curve')

    colors = 'darkblue crimson'.split()
    palette = {label:wrap_around_index(colors, i // len(all_patches)) for i, label in enumerate(violin_labels,start=0)}
    g = sns.stripplot(df, color='grey', palette=palette)
    for points in g.collections:
        points.set_alpha(0.75)

    g = sns.violinplot(df, ax=g, palette=palette)
    for violin in g.collections:
        violin.set_alpha(0.6)
    [x.set_linewidth(3) for x in ax.spines.values()]
    fig.savefig(savename)
    plt.show()


def vilion_plot_of_condensates(df,
                              savename,
                              confidence_interval = 0.8):
    import seaborn as sns
    sns.set_style('whitegrid')
    colors = 'darkblue crimson'.split()
    filename = savename.replace('_condensates.pdf', '')
    basename = os.path.basename(filename)

    fig = plt.figure(figsize=(10, 8))
    ax = fig.add_subplot(111)
    ax.set_title(f'Calculated pH of condensates ({basename})')
    ax.set_ylabel('Median calculated pH')
    ax.set_xlabel('Calibration curve')
    g = sns.stripplot(df)
    g = sns.violinplot(df, ax=g)
    for violin in g.collections:
        violin.set_alpha(0.6)
    [x.set_linewidth(3) for x in ax.spines.values()]
    fig.savefig(savename)
    plt.show()
