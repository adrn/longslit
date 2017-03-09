# Standard library
import os
import pathlib

# Third-party
import yaml

# Package
from .log import logger

global_config_filename = 'spec_global_config.yml'
config_filename = 'spec_config.yml'

# names of FITS OBJECT keywords to skip

class GlobalConfig(object):

    def __init__(self, name, overscan,
                 dispersion_axis=0, # TODO: can only be 0 or 1 - validate!
                 path_exclude=None, path_include=None,
                 **kwargs):
        """
        Global configuration settings.

        Parameters
        ----------
        name : str
            Name of the pipeline run.
        overscan : int
            Index at which the overscan region begins.
            TODO: assumes it is always greater than index - modify to allow other
            sepcifidcations of the overscan region.
        dispersion_axis : int (optional)
            The dispersion axis index along the CCD (default = 0).
        path_exclude : list (optional)
            Path names to exclude from the reduction process, include all others.
        path_include : list (optional)
            Path names to include in the reduction process, ignore all others.
        """

        for key in kwargs:
            logger.warn("GlobalConfig received unrecognized keyword argument '{}'. "
                        "Ignoring...".format(key))

        self.dispersion_axis = int(dispersion_axis)

        # TODO HACK: overscan region index must be an integer specifying the row/column
        #   *after which* is the overscan region. This should support specifying a
        #   rectangular region with indices
        self.overscan = int(overscan) # TODO: HACK, only support index

        # if skip_objects is None:
        #     skip_objects = []
        # elif isinstance(skip_objects, str):
        #     skip_objects = [skip_objects]
        # self.skip_objects = list(skip_objects)

        # TODO NOTE: this might not be necessary if everything is handled by ccdproc
        self.path_include = path_include
        self.path_exclude = path_exclude

        if init_wavelength_filename is None:
            raise ValueError("You must specify the wavelength filename.")
        self.init_wavelength_filename = str(init_wavelength_filename)

    def write(self, f, overwrite=False):
        """
        Write current configuration settings to a file-like object.
        TODO:
        """

        if f.hasattr('write'):
            f.write(yaml.dumps(vars(self)))

        else:
            if os.path.exists(f) and not overwrite:
                raise IOError("File exists! Use overwrite=True to overwrite.")

            with open(f, 'w') as f:
                f.write(yaml.dumps(vars(self)))

    @classmethod
    def from_path(cls, path):
        """
        Check for a config file in the specified path, or any parent path.
        """

        path = pathlib.Path(path).expanduser().absolute().resolve()
        init_path = path / '' # make a copy

        if path.is_file():
            raise ValueError("Path '{}' is a file - specify a directory."
                             .format(str(path)))

        # Find all config files starting from the input path, crawling upwards
        #   through the directory structure
        root = pathlib.Path('/')
        config_files = []
        while path != root:
            if (path / config_filename).is_file():
                config_files.append(path)
            path = (path / '..').resolve()

        if len(config_files) == 0:
            raise IOError("Couldn't find any valid config files '{}' in '{}' or parent directories."
                          .format(config_filename, str(init_path)))

        # Now open each config file and load config settings. Deeper config files get
        #   precedence, so we loop over the config filenames in reverse order, replacing
        #   as we go.
        config_dict = dict()
        for filename in config_filename[::-1]:
            with open(str(filename)) as f:
                this_config_dict = yaml.load(f)
            config_dict.update(**this_config_dict)

        return cls(**config_dict)
