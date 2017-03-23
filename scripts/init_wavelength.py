# coding: utf-8

""" Generate a rough, initial wavelength solution for a 1D spectrum. """

from __future__ import division, print_function

# Standard library
import os
from os.path import abspath, expanduser, exists, join
import glob

# Third-party
from astropy.io import fits
import astropy.units as u
import ccdproc
from matplotlib.widgets import SpanSelector
import matplotlib.pyplot as plt
import numpy as np
from PyQt5 import QtCore, QtGui, QtWidgets
from PyQt5.QtCore import Qt
import yaml
from scipy.interpolate import InterpolatedUnivariateSpline

# Package
from longslit.wavelength import fit_emission_line
from longslit.log import logger

def get_line_props(pixels, flux, xmin, xmax, sigma0=4.):
    i1 = int(np.floor(xmin))
    i2 = int(np.ceil(xmax))+1
    pix = pixels[i1:i2]
    _flux = flux[i1:i2]
    norm_flux = (_flux - _flux.min()) / (_flux.max() - _flux.min())

    try:
        line_props = fit_emission_line(pix, norm_flux,
                                       centroid0=pix[norm_flux.argmax()],
                                       sigma0=sigma0,
                                       amp0=np.sqrt(2*np.pi*sigma0**2),
                                       offset0=0.)
    except ValueError as e:
        msg = "Failed to fit line! See terminal for more information."
        logger.error(msg)
        logger.error(str(e))
        # textbox.setText(msg) # TODO: turn all of this into a class so I can get that
        return None

    line_props['amp'] *= (_flux.max() - _flux.min())
    line_props['const'] += _flux.min()
    return line_props

def draw_line_marker(line_props, wavelength, ax):
    peak = line_props['amp']/(np.sqrt(2*np.pi)*line_props['stddev'])
    print(peak)
    centroid = line_props['centroid']

    space = 0.05*peak
    ax.plot([centroid,centroid], [peak+space,peak+3*space],
            lw=1., linestyle='-', marker='', alpha=0.5, c='#2166AC')
    ax.text(centroid, peak+4*space, "{:.3f} $\AA$".format(wavelength),
            ha='center', va='bottom', rotation='vertical')

def gui_solution(pixels, flux, fig, ax, line_list=None):
    # map_dict = dict()
    # map_dict['wavelength'] = []
    # map_dict['pixel'] = []
    # line_widths = [] # for auto-solving lines

    # FOR TESTING:
    map_dict = {'wavelength': [5460.7399999999998, 7245.1665999999996, 6717.0429999999997, 6678.2762000000002, 6598.9529000000002, 6532.8822], 'pixel': [1226.9646056472734, 349.38080535972756, 610.93127457855871, 630.09556101866531, 668.9871368080278, 701.444940640303]}
    line_widths = [1.5] # for auto-solving lines

    # Add a side menu for specifying the wavelength of the selected line
    panel = QtWidgets.QWidget()

    label = QtWidgets.QLabel("Enter wavelength [Å]:", panel)
    help_label = QtWidgets.QLabel("Enter a wavelength value, zoom in, then select the region "
                                  "containing the emission line at the specified wavelength.",
                                  panel)
    help_label.setStyleSheet("font-style: italic;")
    textbox = QtWidgets.QLineEdit(parent=panel)

    panel.layout = QtWidgets.QGridLayout(panel)
    panel.layout.addWidget(label, 0, 0, 1, 1)
    panel.layout.addWidget(textbox, 0, 1, 1, 1)
    panel.layout.addWidget(help_label, 1, 0, 1, 2, Qt.AlignCenter)

    main_window = fig.canvas.manager.window
    dock = QtWidgets.QDockWidget("Enter wavelength:", main_window)
    main_window.addDockWidget(Qt.BottomDockWidgetArea, dock)
    dock.setWidget(panel)
    # -------------------------------------------------------

    def on_select(xmin, xmax):
        wvln = textbox.text().strip()

        if wvln == '':
            textbox.setText('Error: please enter a wavelength value before selecting')
            return

        wave_val = float(wvln)

        # if line_list specified, find closest line from list:
        if line_list is not None:
            absdiff = np.abs(np.array(line_list) - wave_val)
            idx = absdiff.argmin()
            if absdiff[idx] > 1.:
                logger.error("Couldn't find precise line corresponding to "
                             "input {:.3f}".format(wave_val))
                return

            logger.info("Snapping input wavelength {:.3f} to line list "
                        "value {:.3f}".format(wave_val, line_list[idx]))
            wave_val = line_list[idx]

        line_props = get_line_props(pixels, flux, xmin, xmax)
        draw_line_marker(line_props, wave_val, ax)

        fig.suptitle('')
        plt.draw()
        fig.canvas.draw()

        map_dict['wavelength'].append(wave_val)
        map_dict['pixel'].append(line_props['centroid'])
        line_widths.append(line_props['stddev'])

    # A 1D span selector to highlight a given line
    span = SpanSelector(ax, on_select, 'horizontal', useblit=True,
                        rectprops=dict(alpha=0.5, facecolor='red'))
    span.set_active(False)

    def enable_span():
        tb = fig.canvas.manager.toolbar

        if span.active:
            span.set_active(False)
        else:
            span.set_active(True)

        if tb._active == 'PAN':
            tb.pan()
            tb._actions['pan'].setChecked(False)

        if tb._active == 'ZOOM':
            tb.zoom()
            tb._actions['zoom'].setChecked(False)

    span_control = QtWidgets.QPushButton('λ')
    span_control.setCheckable(True)
    span_control.clicked.connect(enable_span)
    span_control.setStyleSheet("color: #de2d26; font-weight: bold;")
    fig.canvas.manager.toolbar.addWidget(span_control)

    if line_list is not None:
        def auto_identify():
            _idx = np.argsort(map_dict['wavelength'])
            wvln = np.array(map_dict['wavelength'])[_idx]
            pixl = np.array(map_dict['pixel'])[_idx]

            # build an approximate wavelength solution to predict where lines are
            spl = InterpolatedUnivariateSpline(wvln, pixl, k=3)

            predicted_pixels = spl(line_list)

            new_waves = []
            new_pixels = []

            lw = np.median(line_widths)
            for pix_ctr,xmin,xmax,wave in zip(predicted_pixels,
                                              predicted_pixels-3*lw,
                                              predicted_pixels+3*lw,
                                              line_list):
                lp = get_line_props(pixels, flux, xmin, xmax,
                                    sigma0=lw)

                if lp is None or lp['amp'] < 1000.: # HACK
                    # logger.error("Failed to fit predicted line at {:.3f}A, {:.2f}pix"
                    #              .format(wave, pix_ctr))
                    continue

                draw_line_marker(lp, wave, ax)
                new_waves.append(wave)
                new_pixels.append(pix_ctr)

            fig.canvas.draw()

            fig2,axes2 = plt.subplots(2,1,figsize=(6,10))
            axes2[0].plot(new_pixels, new_waves, linestyle='none', marker='o')

            coef = np.polynomial.chebyshev.chebfit(new_pixels, new_waves, deg=5)
            pred = np.polynomial.chebyshev.chebval(new_pixels, coef)
            axes2[1].plot(new_pixels, new_waves-pred,
                          linestyle='none', marker='o')

            plt.show()

        autoid_control = QtWidgets.QPushButton('auto-identify')
        autoid_control.clicked.connect(auto_identify)
        fig.canvas.manager.toolbar.addWidget(autoid_control)

    plt.show()

    return map_dict

def main(config_file, linelist_file, wavelength_data_file, aperture1d=None):
    """ """

    for fname in [config_file, linelist_file, wavelength_data_file]:
        if fname is not None and not exists(fname):
            raise IOError("File '{}' does not exist!".format(fname))

    # read in global configuration
    config_file = os.path.abspath(config_file)
    with open(config_file, 'r') as f:
        config = yaml.load(f.read())

    # read linelist if specified
    if linelist_file is not None:
        linelist = np.genfromtxt(linelist_file, usecols=[0], dtype=float)

    else:
        linelist = None

    # now get root data directory:
    parts = config_file.split(os.sep)
    root_path = os.sep.join(parts[:-3])

    if wavelength_data_file is None: # find a COMP lamp:
        for root,dirs,files in os.walk(root_path):
            if len(glob.glob(join(root, '*.fit*'))) == 0:
                continue

            ic = ccdproc.ImageFileCollection(root)
            if len(ic.files) == 0:
                continue

            hdu = None
            # TODO: configurable imagetyp name
            for hdu,wavelength_data_file in ic.hdus(return_fname=True, imagetyp='COMP'):
                break

            if hdu is not None:
                break
        else:
            raise RuntimeError('No comp lamp spectra were found starting from '
                               'top-level path: {}'.format(root_path))

    else:
        hdu = fits.open(wavelength_data_file)[0]

    # read 2D CCD data
    ccd_data = hdu.data
    if ccd_data.ndim != 2:
        raise ValueError("CCD frame data in file '{}' has {} dimensions. Expected 2."
                         .format(wavelength_data_file, ccd_data.ndim))

    if config['dispersion_axis'] == 1:
        ccd_data = ccd_data.T # make the dispersion axis = 0

    # create 1D pixel and flux grids
    pixels = np.arange(ccd_data.shape[0])

    # remove overscan region
    trimmed = ccd_data[:,:config['overscan']] # TODO: customize overscan region

    if aperture1d is not None:
        # TODO: add a way for the user to specify the columns to extract
        raise NotImplementedError('sorry')

    else:
        # 32 pixels around center
        ctr = int(np.median(trimmed.shape[1]))
        c1 = ctr - 16
        c2 = ctr + 16 + 1

    flux = np.median(trimmed[:,c1:c2], axis=-1)

    # make the plot
    fig,ax = plt.subplots(1, 1)
    ax.plot(pixels, flux, marker='', drawstyle='steps-mid')
    ax.set_xlim(pixels.min(), pixels.max())

    # since we require QT anyways...
    figManager = plt.get_current_fig_manager()
    figManager.window.showMaximized()

    wave_pix_map = gui_solution(pixels, flux, fig=fig, ax=ax,
                                line_list=linelist)
    print(wave_pix_map)

    return

    # sort by pixel and write to file
    _ix = np.argsort(wave_to_pix['pixel'])
    pix_wvl = zip(np.array(wave_to_pix['pixel'])[_ix],
                  np.array(wave_to_pix['wavelength'])[_ix])
    with open(outputpath, 'w') as f:
        txt = ["# pixel wavelength"]
        for row in pix_wvl:
            txt.append("{:.5f} {:.5f}".format(*row))
        f.write("\n".join(txt))

if __name__ == "__main__":
    from argparse import ArgumentParser
    import logging

    # Define parser object
    parser = ArgumentParser(description="")

    vq_group = parser.add_mutually_exclusive_group()
    vq_group.add_argument('-v', '--verbose', action='count', default=0, dest='verbosity')
    vq_group.add_argument('-q', '--quiet', action='count', default=0, dest='quietness')

    parser.add_argument("-c", "--config", dest="config_file", required=True, type=str,
                        help="Path to configuration file.")
    parser.add_argument("--linelist", dest="linelist_file", type=str, default=None,
                        help="Path to a text file where the 0th column is a list of "
                             "emission lines for the comparison lamp. Default is to "
                             "require the user to enter exact wavelengths.")
    parser.add_argument("--comp", dest="wavelength_data_file", type=str, default=None,
                        help="Path to a specific comp lamp file. Default is to find "
                             "the first one within the data directory structure.")

    args = parser.parse_args()

    # Set logger level based on verbose flags
    if args.verbosity != 0:
        if args.verbosity == 1:
            logger.setLevel(logging.DEBUG)
        else: # anything >= 2
            logger.setLevel(1)

    elif args.quietness != 0:
        if args.quietness == 1:
            logger.setLevel(logging.WARNING)
        else: # anything >= 2
            logger.setLevel(logging.ERROR)

    else: # default
        logger.setLevel(logging.INFO)

    kwargs = vars(args)
    kwargs.pop('verbosity')
    kwargs.pop('quietness')
    main(**kwargs)
