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
        Continuing with our example, execute:

        # TODO: what will the default filename be?

        .. code-block:: bash

            cd /path/to/longslit
            python scripts/init_wavelength.py \
            --config=/path/to/observing_run/longslit/mdm-rrlyrae/mdm-rrlyrae.yml \
            --gui

        This will start the GUI interface.

        The first window will ask you to first select a range of columns or rows
        to extract a 1D spectrum from (by taking the median of this region to
        remove cosmic rays).

        The first window will close, and the next window will then display the
        1D comparison lamp spectrum vs. pixel. Here you'll need to identify as
        many lines as you can using a reference spectrum from your observatory.
        When you see a line you know, use the 'scrubbing' tool to select a range
        of pixels around the line that contains just the line you are
        identifying. When you select a line, you'll see a text box appear --
        enter the exact wavelength of the contained line in Angstroms. If the
        line is blended or part of a multi-plet, pick a different one! Make sure
        the lines you identify span a suitable range of pixels and cover the
        spectrum fairly well.

        As you identify lines, the script will output (to a cache file and to
        the command line) the best-fit wavelengths of the lines within the
        selected regions by fitting a Gaussian to the pixels and flux values.

    b. **Using the CLI to initialize the wavelength solution**
        TODO: Fill this in.


    .. warning::

        This might fail if there is significant tilt or rotation of the CCD with
        respect to the disperser. In such a case, the emission lines will curve
        or tilt across the spatial axis of the CCD and may make it hard to
        fit for the 2D wavelength solution. It would be nice to support these
        cases too!

#. ****
    More stuff
