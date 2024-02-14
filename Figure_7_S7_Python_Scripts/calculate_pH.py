#!/usr/bin/env python
import os
import pandas as pd
import numpy as np
from pathlib import Path
from typing import Union, Dict, Tuple, Any, List, Iterable
import matplotlib.pyplot as plt
from scipy.stats import linregress
from scipy.optimize import curve_fit, fsolve
from mpl_toolkits.axes_grid1 import make_axes_locatable
from argparse import ArgumentParser, FileType, Namespace
from tifffile import TiffFile
from matplotlib.cm import get_cmap
from czifile import imread
from pprint import pprint
from skimage import measure, filters
from skimage.measure import label, regionprops, regionprops_table
import matplotlib
from scipy.ndimage import label
from skimage.util import view_as_windows
from scipy.ndimage import gaussian_filter
from mpl_toolkits.mplot3d import Axes3D
from tabulate import tabulate
from plotting import (
    vilion_plot_of_condensates,
    violin_plot_of_dilute_patches,
    plot_identified_fits,
    plot_filtered_calibration_curves,
    plot_calculated_pH
)


# Supremely useful answer to convert a DataFrame to a TSV file: https://stackoverflow.com/a/35974742/866930
def to_fwf(df, fname):
    content = tabulate(df.values.tolist(), list(df.columns), tablefmt="plain", floatfmt='.7f')
    open(fname, 'w').write(content)
pd.DataFrame.to_fwf = to_fwf  # attach the `to_fwf` as an instance method to Pandas Dataframes.


# Necessary for editable PDF and PostScript files in Illustrator
matplotlib.rcParams['pdf.fonttype']     = 42
matplotlib.rcParams['ps.fonttype']      = 42
matplotlib.rcParams['font.sans-serif']  = 'Arial'
matplotlib.rcParams['font.family']      = 'sans-serif'


# Type annotations for easier
FloatList = List[float]
FloatTuple = Tuple[float]
ArrayOrList = Union[FloatList, np.array]


def validate_arguments(parsed_args: Namespace) -> None:
    args = parsed_args
    if args.calibration_measurements_file is None:
        raise RuntimeError('A filename containing measurements is required. Please provide one. Exiting.')

    if args.tifffile is None and args.czifile is None:
        raise RuntimeError('No TIFF or CZI filepath / filename provided for analysis. Please provide one. Exiting.')

    if args.all_wavelengths is None:
        raise RuntimeError('No file containing all wavelengths found. Please provide one. Exiting.')


def get_input_image_name(args) -> Path:
    filepath = None
    if args.tifffile is not None and args.czifile is None:
        filepath = Path(args.tifffile.name).expanduser()

    elif args.tifffile is None and args.czifile is not None:
        filepath = Path(args.czifile.name).expanduser()

    return filepath


# This is an ad-hoc algorithm that tries to add a new NxM mask to an existing 2D array containing X masks of size NxM.
# Depending on the chosen size, convergence isn't possible and manual field selections are recommended. It was originally
# intended to help with the random selections of the dilute phase within an image of the condensate. Given its brittle
# nature which necessitates a heuristic approach that is a special case of the rectangular packing problem, use of this
# function / algorithm should be avoided where possible and manual selections should be prioritized.
def find_non_intersecting_subsections(channels, mask, N, M):
    # Label the mask
    labeled_mask, num_labels = label(mask)

    # Determine the dimensions of the mask
    height, width = mask.shape

    num_tries = 0
    all_coords = list()
    selections = list()
    temp = list()
    x = None
    while True:
        if num_tries == 10000:
            print(f'Error: could not find non intersecting subsections after {num_tries} attempts. Please use manual selections.')
            break

        j = np.random.randint(0, height - M + 1, 1)[0]
        i = np.random.randint(0, width - N + 1, 1)[0]

        # check for collisions in the labelled
        selection = labeled_mask[j:j+M, i:i+N].copy()
        if np.any(selection):
            num_tries += 1
            continue

        coords = (j, j+M, i, i+N)
        all_channel_selections = list()
        for channel in channels:
            s = channels[channel][j:j+M, i:i+N]
            all_channel_selections.append(s)
        median = np.max(all_channel_selections, axis=0)

        norm = (median - np.min(median)) / (np.max(median) - np.min(median))
        norm_flat = norm.flatten()
        indices = np.where(norm_flat > 0.5)
        if len(indices[0]) / len(norm_flat) > 0.0075:
            num_tries += 1
            continue

        all_coords.append(coords)
        labeled_mask[j:j+M, i:i+N] = 10

        patch = dict()
        for channel in channels:
            patch[channel] = channels[channel][j:j+M, i:i+N]
        selections.append(patch)
    return selections, all_coords


def read_calibration_measurements(filename_or_obj: Union[FileType, str, Path]) -> pd.DataFrame:
    data = pd.read_csv(filename_or_obj, sep='\s+')

    # determine that we've read the data properly
    if len(data.columns) <= 1:
        raise RuntimeError('Only 1 column read from calibration data measurements. Cannot continue. Exiting.')
    return data


def extract_wavelengths_and_pHs(data: pd.DataFrame) -> (pd.Series, pd.Series):
    all_wavelengths = data['wavelength'].to_list()
    all_pHs = data.columns.to_list()[1:]
    actual_pHs = [float(col.lower().replace('ph', '')) for col in data.columns if 'ph' in col.lower()]
    return all_wavelengths, all_pHs, actual_pHs


def calculate_measurement_ratios(data: pd.DataFrame) -> (Dict[FloatTuple, FloatList], Dict[float, List[FloatList]]):
    all_wavelengths, all_pHs, actual_pHs = extract_wavelengths_and_pHs(data)
    all_ratios = dict()
    all_data_for_plotting = dict()
    for pH, actual_pH in zip(all_pHs, actual_pHs):
        to_plot = list()
        for i, a in enumerate(all_wavelengths):
            row = list()
            for j, b in enumerate(all_wavelengths):
                ratio = data[pH][i]/data[pH][j]
                row.append(ratio)
                key = (a, b)
                if key not in all_ratios:
                    all_ratios[key] = [ratio]
                else:
                    all_ratios[key].append(ratio)
            to_plot.append(row)
        all_data_for_plotting[actual_pH] = to_plot
    return all_ratios, all_data_for_plotting


def calculate_rsquared(real_data: ArrayOrList,
                       fit_data: ArrayOrList) -> float:
    residuals = real_data - fit_data
    ss_res = np.sum(residuals**2.0)
    ss_tot = np.sum((real_data - np.mean(real_data))**2.0)
    r_squared = 1.0 - (ss_res / ss_tot)
    return r_squared


def calculate_fit(actual_pHs,
                  raw_data_curve,
                  sign=1.0) -> Tuple[np.array, np.array, float]:
    fit         = curve_fit(lambda t, a, b: a * np.exp(sign * b * t), actual_pHs, raw_data_curve)
    popt, pcov  = fit
    a, b        = popt
    y_fit       = a * np.exp(sign * b * np.array(actual_pHs))
    rsquared    = calculate_rsquared(raw_data_curve, y_fit)
    return y_fit, fit, rsquared


def filter_ratios(input_filename: Path,
                  data: pd.DataFrame,
                  analyze_area: str = 'lower',
                  cutoff_min_slope: float = 0.5,
                  minimum_rsquared: float = 0.97) -> Tuple[Dict[FloatTuple, ArrayOrList],
                                                           Dict[FloatTuple, FloatTuple],
                                                           Dict[FloatTuple, float]]:
    def highlight_cell(x,y, ax=None, **kwargs):
        rect = plt.Rectangle((x-.5, y-.5), 1, 1, fill=False, **kwargs)
        ax = ax or plt.gca()
        ax.add_patch(rect)
        return ax

    if analyze_area.lower() not in 'lower upper'.split():
        raise RuntimeError(f'Unsupported direction "{direction}".')

    all_wavelengths, all_pHs, actual_pHs = extract_wavelengths_and_pHs(data)
    all_ratios, all_data_for_plotting = calculate_measurement_ratios(data)

    all_heatmaps = None
    lower_triangle_indices = np.tril_indices(len(all_wavelengths))
    upper_triangle_indices = np.triu_indices(len(all_wavelengths))
    indices = lower_triangle_indices
    sign = 1.0
    if analyze_area.lower() == 'upper':
        indices = upper_triangle_indices
        sign = -1.0

    to_print = dict()
    to_print['lambda1'] = list()
    to_print['lambda2'] = list()
    to_print['a'] = list()
    to_print['b'] = list()
    to_print['raw_min'] =  list()
    to_print['raw_max'] = list()
    to_print['fit_min'] = list()
    to_print['fit_max'] = list()

    fits = dict()
    fits_params = dict()
    fits_quality = dict()
    count = 1
    for x, y in zip(*indices):
        a = all_wavelengths[x]
        b = all_wavelengths[y]
        ratio = (a, b)
        raw_data_curve = all_ratios[ratio]

        # check the slope and first of the data
        regression = linregress(actual_pHs, raw_data_curve)
        if abs(regression.slope) >= cutoff_min_slope:
            y_fit, fit, rsquared = calculate_fit(actual_pHs, raw_data_curve, sign)
            if rsquared >= minimum_rsquared:
                fits[ratio]         = y_fit
                fits_params[ratio]  = fit
                fits_quality[ratio] = rsquared


                scale = 1.75
                fontsize=16 * scale
                extension = 'tiff'
                ###fig = plt.figure(figsize=(12, 8))
                ###ax = fig.add_subplot(111)
                print(len(actual_pHs), len(raw_data_curve))

                a, b = fits_params[ratio][0]
                x = np.linspace(actual_pHs[0], actual_pHs[-1], 50)
                y = a * np.exp(sign * b * x)

                colormap = np.array('red orange cyan lightgreen magenta'.split())
                categories = np.arange(len(actual_pHs))
                ###ax.scatter(actual_pHs, raw_data_curve, s=150, edgecolor='k', linewidth=3, alpha=1.0, c=colormap[categories])
                rs = f'$(R^2={rsquared:.3f})$'
                label = 'Fit to curve $a*e^{-bx}$ %s' % rs
                ###ax.plot(x, y, linestyle=':', linewidth=3, label=label, zorder=-1)

                ###ax.set_title(f'Filtered Calibration Curve for $\\lambda_i={ratio[1]}$, $\\lambda_j={ratio[0]}$', fontsize=fontsize)
                ###ax.set_xlabel('pH', fontsize=fontsize)
                ###ax.set_ylabel(r'$\frac{Em_j}{Em_i}$', fontsize=fontsize)

                ###ax.tick_params(axis='both', which='major', labelsize=fontsize)
                ###ax.tick_params(axis='both', which='minor', labelsize=fontsize)
                ###[x.set_linewidth(3) for x in ax.spines.values()]

                ###ax.legend(fontsize=fontsize)
                ###fig.savefig(f'example_curve_{scale}.{extension}')
                ###plt.show()

                # exit()
                print(ratio[0], ratio[1], fit[0][0], fit[0][1], min(raw_data_curve), max(raw_data_curve), y_fit.min(), y_fit.max())
                to_print['lambda1'].append(ratio[0])
                to_print['lambda2'].append(ratio[1])
                to_print['a'].append(fit[0][0])
                to_print['b'].append(fit[0][1])
                to_print['raw_min'].append(min(raw_data_curve))
                to_print['raw_max'].append(max(raw_data_curve))
                to_print['fit_min'].append(y_fit.min())
                to_print['fit_max'].append(y_fit.max())


                count += 1
    to_print = pd.DataFrame(to_print)
    print(to_print.to_fwf(f'{input_filename.stem}_fits.txt'))

    return fits, fits_params, fits_quality


def calculate_pH(fit_coeff: FloatTuple,
                 y: float) -> float:
    # `a` and `b` correspond to the variables in `y_fit = a * np.exp(b * np.array(x))`.
    a, b = fit_coeff
    return np.log(y/a)/b


def calculate_pH_scaled_y(fit_coeff: ArrayOrList,
                          y_fit: ArrayOrList,
                          y_scaled: float,
                          sign=1.0) -> float:
    a, b = fit_coeff
    return sign * (1.0/b) * np.log(((y_scaled * (y_fit.max() - y_fit.min())) + y_fit.min())/a)


def read_TiffFile(filename: Union[str, Path], all_wavelengths) -> Dict[float, np.ndarray]:
    channels = dict()
    with TiffFile(filename) as tiff:
        tiff_channels       = tiff.asarray()
        num_tiff_channels   = tiff_channels.shape[0]
        num_wavelengths     = len(all_wavelengths)
        if num_tiff_channels != num_wavelengths:
            raise RuntimeError(f'The number of tiff channels ({num_tiff_channels}) '
                               f'must match the number of wavelengths ({num_wavelengths}).')
        for channel, tdata in enumerate(tiff_channels):
            channels[all_wavelengths[channel]] = tdata
    return channels


def cleanup_mask(mask, min_area = 50):
    cleaned = np.zeros_like(mask)
    labelled = measure.label(mask)
    condensates_indices = dict()

    areas = dict()
    num_objects = int(np.max(labelled))
    count = 1
    for obj in range(1, num_objects+1):
        indices = np.where(labelled == obj)
        selection = labelled[indices]
        if len(selection.flatten()) >= min_area:
            cleaned[indices] = count
            areas[count] = len(indices[0])
            condensates_indices[count] = indices
            count += 1
    return cleaned, areas, condensates_indices


def select_data_for_bootstrapping(median, cleaned_mask, condensate=True):
    name = 'Dilute Phase'
    indices = np.where(cleaned_mask == 0)
    if condensate:
        name = 'Condensate'
        indices = np.where(cleaned_mask > 0)

    coordinates = [(x, y) for x, y in zip(*indices)]
    blank = np.zeros_like(median)
    data_to_sample = list()
    for coord in coordinates:
        x, y = coord
        blank[x, y] = median[x, y]
        data_to_sample.append(median[x, y])
    return data_to_sample


def bootstrap_condensates(median, areas, condensates_mask, sample_fraction = 0.9):
    analysis = dict()
    for index in areas:
        indices = np.where(condensates_mask == index)
        condensate = median[indices].flatten()
        c_all_mins, c_all_maxes, c_all_means, c_all_medians = bootstrap(condensate, int(sample_fraction * len(condensate)))
        analysis[index] = (c_all_mins, c_all_maxes, c_all_means, c_all_medians)
    return analysis


def calculate_pHs_from_measurements_and_fits(channels, fits, fits_params, sign):
    first_channel = channels[list(channels.keys())[0]]
    first_channel_num_cells = first_channel.shape[0] * first_channel.shape[1]

    all_pH_fields = list()
    all_pH_ratios = list()
    all_c_ratios = np.zeros(len(fits) * first_channel_num_cells)
    all_pHs = np.zeros(len(fits) * first_channel_num_cells)

    for index, ratio in enumerate(fits):
        a, b = ratio
        channel_a = channels[a]
        channel_b = channels[b]
        c_ratio = channel_a / channel_b

        scaled_measurement = (c_ratio - c_ratio.min())/(c_ratio.max() - c_ratio.min())
        calc_pH = calculate_pH_scaled_y(fits_params[ratio][0], fits[ratio], scaled_measurement, sign)
        i = index * first_channel_num_cells
        j = (index * first_channel_num_cells) + first_channel_num_cells
        all_c_ratios[i:j] = scaled_measurement.flatten()
        all_pHs[i:j] = calc_pH.flatten()

        all_pH_fields.append(calc_pH)
        all_pH_ratios.append(ratio)
    combined = np.array(all_pH_fields)
    return combined, all_c_ratios, all_pHs, all_pH_ratios


def bootstrap(data_to_sample, sample_size, num_iterations=100):
    all_mins = list()
    all_maxes = list()
    all_means = list()
    all_medians = list()
    for iteration in range(num_iterations):
        sampled = np.random.choice(data_to_sample, sample_size)
        all_mins.append(sampled.min())
        all_maxes.append(sampled.max())
        all_means.append(sampled.mean())
        all_medians.append(np.median(sampled))
    return all_mins, all_maxes, all_means, all_medians


def output_data(data):
    print(f'Min: {np.min(data)}, Max: {np.max(data)}, Mean: {np.mean(data)}, Median: {np.median(data)}')


def read_czifile(filename, new_shape=None):
    raw_image = imread(filename)
    shape = tuple(s for s in raw_image.shape if s > 1)
    image = imread(filename).reshape(shape)
    if new_shape is not None:
        image = image.reshape(new_shape)
    fig = plt.figure(111)
    ax = fig.add_subplot(111)
    ax.imshow(image[0], cmap='jet')
    plt.show()


def identify_and_calibrate_channels_measurements(args):
    data                                        = read_calibration_measurements(args.calibration_measurements_file)
    all_possible_wavelengths                    = np.loadtxt(args.all_wavelengths)
    all_ratios, all_data_for_plotting           = calculate_measurement_ratios(data)
    input_filename                              = get_input_image_name(args)
    fits, fits_params, fits_quality             = filter_ratios(input_filename, data, analyze_area=args.calibration_curves_to_use)

    if args.tifffile is not None and args.czifile is None:
        channels = read_TiffFile(args.tifffile, all_possible_wavelengths)
        channel_shape = channels[all_possible_wavelengths[0]].shape
        reshaped = np.zeros((len(channels), channel_shape[0], channel_shape[1]))
        for idx, channel in enumerate(channels):
            reshaped[idx,:,:] = channels[channel]

    elif args.tifffile is None and args.czifile is not None:
        raw_channels = imread(args.czifile)

        shape = tuple(s for s in raw_channels.shape if s > 1)
        reshaped = raw_channels.reshape(shape)
        channels = {k:reshaped[i,:,:] for i, k in enumerate(all_possible_wavelengths)}
    return data, channels, reshaped


def calculate_dilute_pH(patch):
    flat_patch = patch.flatten()
    return bootstrap(flat_patch, int(0.9 * len(flat_patch)))


def select_random_dilute_samples(channels, num_patches, patch_size):
    all_patches = list()
    for patch_num in range(1, num_patches+1):
        all_patch_channels = dict()
        for channel_index, wavelength in enumerate(channels):
            channel_intensity = channels[wavelength]
            csize = channel_intensity.shape[0]

            i = np.random.randint(0, csize - patch_size, 1)[0]
            j = np.random.randint(0, csize - patch_size, 1)[0]

            patch = channel_intensity[j:j+patch_size, i:i+patch_size]
            all_patch_channels[channel_index] = patch
        all_patches.append(all_patch_channels)
    return all_patches


def select_dilute_samples(channels, all_patch_coords):
    all_patches = list()
    for patch_coords in all_patch_coords:
        yl, yu, xl, xu = patch_coords  # coords are always y-x in numpy (i.e. row-major)
        all_patch_channels = dict()
        for channel_index, wavelength in enumerate(channels):
            channel_intensity = channels[wavelength]
            patch = channel_intensity[yl:yu, xl:xu]
            all_patch_channels[channel_index] = patch
        all_patches.append(all_patch_channels)
    return all_patches


# Note that the
def determine_patches(args, channels, no_patches=False):
    all_patches = list()
    if no_patches:
        return all_patches

    if args.patch_coordinates is not None:
        if len(args.patch_coordinates) > 0:
            for group in args.patch_coordinates:
                if len(group) % 4 != 0:
                    raise RuntimeError(f'Please pass coordinates in groups of 4 (currently: {len(group)}).')
        all_patches = select_dilute_samples(channels, args.patch_coordinates)
        return all_patches

    if args.num_random_dilute_samples is not None and args.random_dilute_sample_size is not None:
       all_patches = select_random_dilute_samples(channels, args.num_random_dilute_samples, args.random_dilute_sample_size)
    return all_patches


def plot_mask(channels, mask, median, variance, area_to_analyze):
    fontsize = 16

    cmap = 'viridis'
    fig, axes = plt.subplots(nrows=2, ncols=2, figsize=(10, 10))

    ax = axes[0, 0]
    ax.set_xlabel('X', fontsize=fontsize)
    ax.set_ylabel('Y', fontsize=fontsize)
    im = ax.imshow(median, cmap=cmap)
    divider = make_axes_locatable(ax)
    cax = divider.append_axes('right', size='5%', pad=0.1)
    colorbar = fig.colorbar(im, cax=cax)
    colorbar.ax.set_ylabel('Calculated pH', fontsize=fontsize)
    colorbar.ax.tick_params(axis='both', which='major', labelsize=fontsize * 0.8)
    colorbar.ax.tick_params(axis='both', which='minor', labelsize=fontsize * 0.8)
    [x.set_linewidth(3) for x in ax.spines.values()]
    ax.tick_params(axis='both', which='major', labelsize=fontsize)
    ax.tick_params(axis='both', which='minor', labelsize=fontsize)

    ax = axes[0, 1]
    ax.set_xlabel('X', fontsize=fontsize)
    ax.set_ylabel('Y', fontsize=fontsize)
    im = ax.imshow(variance, cmap=cmap)
    divider = make_axes_locatable(ax)
    cax = divider.append_axes('right', size='5%', pad=0.1)
    colorbar = fig.colorbar(im, cax=cax)
    colorbar.ax.set_ylabel('$\sigma^2$', fontsize=fontsize)
    colorbar.ax.tick_params(axis='both', which='major', labelsize=fontsize * 0.8)
    colorbar.ax.tick_params(axis='both', which='minor', labelsize=fontsize * 0.8)
    [x.set_linewidth(3) for x in ax.spines.values()]
    ax.tick_params(axis='both', which='major', labelsize=fontsize)
    ax.tick_params(axis='both', which='minor', labelsize=fontsize)

    ax = axes[1, 0]
    ax.set_xlabel('X', fontsize=fontsize)
    ax.set_ylabel('Y', fontsize=fontsize)

    im = ax.imshow(mask, cmap='gray')
    divider = make_axes_locatable(ax)
    cax = divider.append_axes('right', size='5%', pad=0.1)
    colorbar = fig.colorbar(im, cax=cax)
    colorbar.ax.tick_params(axis='both', which='major', labelsize=fontsize * 0.8)
    colorbar.ax.tick_params(axis='both', which='minor', labelsize=fontsize * 0.8)
    [x.set_linewidth(3) for x in ax.spines.values()]
    ax.tick_params(axis='both', which='major', labelsize=fontsize)
    ax.tick_params(axis='both', which='minor', labelsize=fontsize)

    ax = axes[1, 1]
    ax.set_xlabel('X', fontsize=fontsize)
    ax.set_ylabel('Y', fontsize=fontsize)

    dilute_indices = np.where(mask == 0)
    condensate_indices = np.where(mask > 0)
    first_channel = list(channels.keys())[0]
    canvas = np.zeros_like(channels[first_channel])
    canvas[condensate_indices] = channels[first_channel][condensate_indices].flatten().mean()
    canvas[dilute_indices] = channels[first_channel][dilute_indices].flatten().mean()

    im = ax.imshow(canvas, cmap=cmap)
    divider = make_axes_locatable(ax)
    cax = divider.append_axes('right', size='5%', pad=0.1)
    colorbar = fig.colorbar(im, cax=cax)
    colorbar.ax.set_ylabel('Average pH', fontsize=fontsize)
    colorbar.ax.tick_params(axis='both', which='major', labelsize=fontsize * 0.8)
    colorbar.ax.tick_params(axis='both', which='minor', labelsize=fontsize * 0.8)
    [x.set_linewidth(3) for x in ax.spines.values()]
    ax.tick_params(axis='both', which='major', labelsize=fontsize)
    ax.tick_params(axis='both', which='minor', labelsize=fontsize)

    fig.subplots_adjust(wspace=0.5, hspace=0.2)
    fig.savefig(f'all_{area_to_analyze}.pdf')
    plt.show()


def plot_pH(channels, median, variance, savename):
    fontsize = 16
    cmap = 'viridis'
    fig, axes = plt.subplots(nrows=1, ncols=2, figsize=(16, 6))

    ax = axes[0]
    ax.set_xlabel('X', fontsize=fontsize)
    ax.set_ylabel('Y', fontsize=fontsize)
    im = ax.imshow(median, cmap=cmap)
    divider = make_axes_locatable(ax)
    cax = divider.append_axes('right', size='5%', pad=0.1)
    colorbar = fig.colorbar(im, cax=cax)
    colorbar.ax.set_ylabel('Calculated pH', fontsize=fontsize)
    colorbar.ax.tick_params(axis='both', which='major', labelsize=fontsize * 0.8)
    colorbar.ax.tick_params(axis='both', which='minor', labelsize=fontsize * 0.8)
    [x.set_linewidth(3) for x in ax.spines.values()]
    ax.tick_params(axis='both', which='major', labelsize=fontsize)
    ax.tick_params(axis='both', which='minor', labelsize=fontsize)

    ax = axes[1]
    ax.set_xlabel('X', fontsize=fontsize)
    ax.set_ylabel('Y', fontsize=fontsize)
    im = ax.imshow(variance, cmap=cmap)
    divider = make_axes_locatable(ax)
    cax = divider.append_axes('right', size='5%', pad=0.1)
    colorbar = fig.colorbar(im, cax=cax)
    colorbar.ax.set_ylabel('$\sigma^2$', fontsize=fontsize)
    colorbar.ax.tick_params(axis='both', which='major', labelsize=fontsize * 0.8)
    colorbar.ax.tick_params(axis='both', which='minor', labelsize=fontsize * 0.8)
    [x.set_linewidth(3) for x in ax.spines.values()]
    ax.tick_params(axis='both', which='major', labelsize=fontsize)
    ax.tick_params(axis='both', which='minor', labelsize=fontsize)
    fig.savefig(savename)



def main():
    parser = ArgumentParser()

    parser.add_argument('--bootstrap-ratio',
                        '-b',
                        help='The ratio by which to select the number of points for bootstrapping.',
                        type=float,
                        default=0.9)

    parser.add_argument('--calibration-measurements-file',
                        '-c',
                        type=FileType('r'),
                        help='The file containing the calibration measurements.',
                        default=None)

    parser.add_argument('--czifile',
                        '-z',
                        help='The CZI image whose pHs are unknown.',
                        type=FileType('rb'),
                        default=None)

    parser.add_argument('--tifffile',
                        '-i',
                        help='The image whose pHs are unknown.',
                        type=FileType('rb'),
                        default=None)

    parser.add_argument('--all-wavelengths',
                        '-w',
                        help='A text file containing a list of all the wavelengths.',
                        type=FileType('r'),
                        default=None)

    parser.add_argument('--condensate-mask',
                        '-m',
                        type=FileType('r'),
                        help='The condensate mask',
                        default=None)

    parser.add_argument('--calibration-curves-to-use',
                        '-u',
                        help='',
                        type=str,
                        choices='upper lower'.split(),
                        default='upper')

    parser.add_argument('--num-random-dilute-samples',
                        '-nd',
                        type=int,
                        help='The number of dilute samples to perform.',
                        default=10)

    parser.add_argument('--random-dilute-sample-size',
                        '-ds',
                        help='The symmetric size of the dilute sample size (i.e. NxN).',
                        default=100)

    parser.add_argument('--patch-coordinates',
                        '-pc',
                        help='The coordinates corresponding to the location of patches in row-major form (yl,yu,xl,xu).',
                        type=int,
                        nargs=4,
                        action='append',
                        default=list())

    parser.add_argument('--crop-coordinates',
                        '-cc',
                        help='The coordinates corresponding to the location of crops in row-major form (yl,yu,xl,xu).',
                        type=int,
                        nargs=4,
                        action='append',
                        default=None)

    args = parser.parse_args()
    validate_arguments(args)

    data, channels, raw_channels = identify_and_calibrate_channels_measurements(args)
    all_data_wavelengths, all_pHs, actual_pHs   = extract_wavelengths_and_pHs(data)
    all_possible_wavelengths = np.loadtxt(args.all_wavelengths.name).tolist()
    skip_channels_indices = [idx for idx, c in enumerate(all_possible_wavelengths) if c not in all_data_wavelengths]

    no_patches = False
    if args.patch_coordinates is None:
        no_patches = True

    all_patches = determine_patches(args, channels, len(args.patch_coordinates) == 0)
    cleaned_mask = None
    areas = dict()
    putative_patches = list()
    if args.condensate_mask is not None:
        with TiffFile(args.condensate_mask.name) as mtiff:
            mask = mtiff.asarray()
        cleaned_mask, areas, condensates_indices = cleanup_mask(mask)

    patches = select_dilute_samples(channels, args.patch_coordinates)
    for patch in patches:
        for channel in patch:
            plt.imshow(patch[channel], cmap='jet')
            plt.show()

    if args.condensate_mask is not None and no_patches:
        # find patches with low means
        putative_patches, putative_patches_indices = find_non_intersecting_subsections(channels, mask, 100, 100)
        if len(putative_patches_indices) < 5:
            raise RuntimeError(f'Error with patches selection. Only {len(putative_patches_indices)} found.')

        patches = select_dilute_samples(channels, putative_patches_indices)

        medians = dict()
        for idx, patch in enumerate(patches):
            averages = list()
            for channel in patch:
                average = patch[channel].flatten().mean()
                averages.append(average)
            median_average = np.median(averages)
            medians[idx] = median_average

        sorted_medians = {k:v for k, v in sorted(medians.items(), key=lambda item: item[1])}
        all_patches = list()
        count = 0
        for patch in sorted_medians:
            all_patches.append(patches[patch])
            count += 1
            if count == 5:
                break
        print()

    violin_plot_series = list()
    violin_labels = list()
    condensate_medians = dict()
    area_to_analyze = args.calibration_curves_to_use
    sign = 1.0
    if area_to_analyze.lower() == 'upper':
        sign = -1.0

    input_filename = get_input_image_name(args)
    all_ratios, all_data_for_plotting = calculate_measurement_ratios(data)
    fits, fits_params, fits_quality = filter_ratios(input_filename, data, area_to_analyze, sign)
    combined, all_c_ratios, all_pH_calculations, all_pH_ratios = calculate_pHs_from_measurements_and_fits(channels, fits, fits_params, sign)
    savename = os.path.basename(input_filename).replace('.czi', '.pdf')
    plot_mask(channels, cleaned_mask, np.median(combined, axis=0), np.var(combined, axis=0), area_to_analyze)

    all_patches_averages = list()
    for patch in all_patches:
        channels_patch_averages = dict()
        for channel in patch:
            if channel in skip_channels_indices:
                continue

            patch_average = patch[channel][:,:].flatten().mean()
            channels_patch_averages[channel] = patch_average
        all_patches_averages.append(channels_patch_averages)

    for idx, patches_mean_intensity in enumerate(all_patches_averages, start=1):
        all_calc_pH = list()
        print(f'PATCH {channel} {idx}')
        for filtered_curve in range(len(all_pH_ratios)):
            ratio = all_pH_ratios[filtered_curve]
            a, b = ratio
            a_index = all_possible_wavelengths.index(a)
            b_index = all_possible_wavelengths.index(b)

            a_intensity = patches_mean_intensity[a_index]
            b_intensity = patches_mean_intensity[b_index]
            ratio_ab = a_intensity / b_intensity

            scaled_measurement = (ratio_ab - fits[ratio].min())/(fits[ratio].max() - fits[ratio].min())
            calc_pH = calculate_pH_scaled_y(fits_params[ratio][0], fits[ratio], scaled_measurement, sign)
            print(area_to_analyze, ratio, calc_pH)

            all_calc_pH.append(calc_pH)
        print(area_to_analyze, 'PATCH Median: ', np.median(all_calc_pH))
        print()

        label = f'Patch-{idx} ({area_to_analyze})'
        violin_labels.append(label)
        violin_plot_series.append(pd.Series(all_calc_pH, name=label))

        # # PLOT CALCULATED PH
        calculated_pH_savename = f'{input_filename.stem}_{area_to_analyze}_pH.pdf'
        plot_calculated_pH(combined, calculated_pH_savename, area_to_analyze)

        identified_fits_savename = f'identified_fits_{area_to_analyze}.pdf'
        fits_colors = plot_identified_fits(fits, all_data_wavelengths, identified_fits_savename)
        fits_savename = f'fits_{area_to_analyze}.pdf'
        fits_moredata = dict()
        for ratio in fits_params:
            a, b = fits_params[ratio][0]
            x = np.linspace(actual_pHs[0], actual_pHs[-1], 50)
            fits_moredata[ratio] = a * np.exp(sign * b * x)
        plot_filtered_calibration_curves(actual_pHs, all_ratios, x, fits_moredata, fits_savename, fits_colors)


        # condensates
        all_condensate_averages = list()
        for condensate_index in condensates_indices:
            indices = condensates_indices[condensate_index]

            condensate_averages = dict()
            for channel_index in range(raw_channels.shape[0]):
                if channel_index in skip_channels_indices:
                    continue

                channel_intensity = raw_channels[channel_index,:,:]
                channel_condensate_mask = channel_intensity[indices[0], indices[1]]
                condensate_sum = channel_condensate_mask.sum()
                average = condensate_sum / len(indices[0])
                condensate_averages[channel_index] = average
                print(channel_intensity)
            all_condensate_averages.append(condensate_averages)

        all_condensates_calc_pH_medians = list()
        for idx, condensate_averages in enumerate(all_condensate_averages, start=1):
            all_calc_pH = list()
            print(f'CONDENSATE {idx}')
            for filtered_curve in range(len(all_pH_ratios)):
                ratio = all_pH_ratios[filtered_curve]
                a, b = ratio
                a_index = all_possible_wavelengths.index(a)
                b_index = all_possible_wavelengths.index(b)

                a_intensity = condensate_averages[a_index]
                b_intensity = condensate_averages[b_index]
                ratio_ab = a_intensity / b_intensity

                scaled_measurement = (ratio_ab - fits[ratio].min())/(fits[ratio].max() - fits[ratio].min())
                calc_pH = calculate_pH_scaled_y(fits_params[ratio][0], fits[ratio], scaled_measurement, sign)
                print(area_to_analyze, ratio, calc_pH)

                all_calc_pH.append(calc_pH)
            all_condensates_calc_pH_medians.append(np.median(all_calc_pH))
            print(area_to_analyze, 'CONDENSATE Median: ', np.median(all_calc_pH))
            print()
        condensate_medians[area_to_analyze] = all_condensates_calc_pH_medians[:]


    print(len(violin_plot_series), 1111, len(all_patches))
    if len(violin_plot_series) > 0:
        patches_df = pd.concat(violin_plot_series, axis=1)
        with open(f'{input_filename.stem}_patches.txt', 'w') as ifile:
            ifile.write(patches_df.to_string(index=False))
        patches_savename = f'{input_filename.stem}_patches.pdf'
        violin_plot_of_dilute_patches(patches_df, all_patches, violin_labels, patches_savename)

    condensates_df = pd.DataFrame(condensate_medians)
    with open(f'{input_filename.stem}_condensates.txt', 'w') as ifile:
        ifile.write(condensates_df.to_string(index=False))
    condensates_savename = f'{input_filename.stem}_condensates.pdf'
    vilion_plot_of_condensates(condensates_df, condensates_savename)

    for filtered_curve in range(len(all_pH_ratios)):
        ratio = all_pH_ratios[filtered_curve]
        print(ratio)


if __name__ == '__main__':
    main()
