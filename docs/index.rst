.. |br| raw:: html

   <br />

######################
Long Slit Spectroscopy
######################

This package contains tools for reducing 2D images to 1D, wavelength
calibrated spectra. The package contains class and function-level access
to most of the procedures needed to do this, and also provides a 'pipeline'
interface to the reduction scripts that is customizable through a set of
configuration files.

==================
Pipeline interface
==================

Terminology:

* **collection**:
    A set of data that can be reduced together. This doesn't necessarily
    correspond to all data from a single observing run. For example, if on night
    3 of a 4 night observing run you change the tilt of the diffraction grating,
    this could significantly change the wavelength solutions of spectra on
    nights 3 and 4 relative to 1 and 2. These would be considered two separate
    *collections* and the pipeline would have to be run on these subsets of the
    full data separately.

Using the package as a pipeline assumes that there is some *collection* of
data that can be reduced in parallel following the same procedure, all with
the same rough wavelength solution. This is because the same rough wavelength
solution is used as an initialization for the per-source wavelength solutions.

For the examples below, we'll assume the simpler case in which all of the data
can be reduced identically, and that they all live under some root directory
``/path/to/observing_run``.

With that in mind, the first thing to do in using the pipeline interface is to
initialize the pipeline using the ``init_pipeline.py`` script. We need to
specify the path to the root of the data, and a name for the pipeline run:

.. code-block:: bash

    cd /path/to/longslit
    python scripts/init_pipeline.py --rootpath=/path/to/observing_run --name=mdm-rrlyrae

#. **Assumed file structure**
    The pipeline assumes that all of the data are stored under the same
    top-level directory. Within that directory, the data may organized in
    further subdirectories. For example, if the top-level path is called
    ``observing_run``, both ``observing_run/night1/<data>``,
    ``observing_run/night2/<data>``, etc., and ``observing_run/<all data>``
    are supported. After running ``init_pipeline.py`` (as above), the pipeline
    will create a directory in ``observing_run/longslit/`` to store metadata
    about the reduction pipeline status as the data are processed.

#. **Creating the global configuration file**
    Initializing the pipeline will also create a template global configuration
    file. The globally configurable parameters are:

    * ``name``: The name of the pipeline run. The metadata for a given pipeline
        run will be stored in ``observing_run/longslit/<name>`` (e.g.,
        information about what files have already been processed). This will be
        filled with the name you specify when running ``init_pipeline.py``
    * ``dispersion_axis``: The index of the axis along the CCD that corresponds
        to the dispersion (e.g., wavelength) axis. This follows `numpy`
        convention, so ``axis=0`` means the 'row' direction is the dispersion
        axis.
    * ``overscan``: The index of the row/column where the overscan region
        begins. This is always assumed to be in the direction orthogonal to the
        dispersion axis. So, if ``dispersion_axis=0``, setting ``overscan=256``
        means that the overscan region begins at the pixel column at index 256.
    * ``path_exclude``: A pattern or list of patterns that are checked against
        all subdirectories of the top-level data path. When a subdirectory name
        matches one of the patterns, it is excluded from the processing queue.
        If this is specified, all other subdirectories are included by default.
    * ``path_include``: A pattern or list of patterns that are checked against
        all subdirectories of the top-level data path. When a subdirectory name
        matches one of the patterns, it is included in the processing queue.
        If this is specified, all other subdirectories are excluded by default.

#. **Initialize the wavelength solution**
    The next thing that needs to be done is to specify an initial, rough
    identification of spectral lines in a wavelength calibration (comparison
    lamp) spectrum. This is either done interactively through a `matplotlib`
    GUI, or by specifying (pixel, wavelength) pairs through the command line
    interface (you can specify one or the other when running the
    ``init_wavelength.py`` -- use the flag ``--cli`` for the command-line
    interface or ``--gui`` for the GUI).

    a. **Using the GUI to initialize the wavelength solution**
        For example:

        # TODO: what will the default filename be?

        .. code-block:: bash

            cd /path/to/longslit
            python scripts/init_wavelength.py \
            --config=/path/to/observing_run/longslit/mdm-rrlyrae/mdm-rrlyrae.yml \
            --gui

        Running ``init_wavelength.py`` with the GUI interface will open two
        interactive windows.

        - Main window should be the plot of the 1D spectrum, median'd over
            some region.
        - Second window should be a table view of (pixel, wavelength guess)
            for each of the 3 1D traces

#. **The next thing**
    More stuff
