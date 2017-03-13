# coding: utf-8

""" Initialize the 1D spectral reduction pipeline. """

# Standard library
import os
from os.path import abspath, expanduser, exists, join
import sys

# Third-party
import yaml

# Package
from longslit.log import logger

def main(name, rootpath):

    rootpath = abspath(expanduser(rootpath))
    if not exists(rootpath):
        raise IOError("Path '{}' doesn't exist!".format(rootpath))

    pipeline_path = join(rootpath, 'longslit', name)
    os.makedirs(pipeline_path, exist_ok=True)
    logger.debug('Pipeline output path: {}'.format(pipeline_path))

    global_config_filename = join(pipeline_path, '{}-config.yml'.format(name))

    if exists(global_config_filename):
        logger.error("Config file already exists at '{}'\n ignoring..."
                     .format(global_config_filename))
        sys.exit(1)

    defaults = dict()
    defaults['name'] = name
    defaults['dispersion_axis'] = 0
    defaults['overscan'] = 0
    defaults['path_exclude'] = ['']
    defaults['path_include'] = ['']
    defaults['gain'] = '2.7 electron/adu'
    defaults['read_noise'] = '7.9 electron'

    with open(global_config_filename, 'w') as f:
        yaml.dump(defaults, f)

    logger.info('Created template pipeline global config file at: {}'
                .format(global_config_filename))

    for k,v in defaults.items():
        logger.debug('{} = {}'.format(k, v))

if __name__ == "__main__":
    from argparse import ArgumentParser
    import logging

    # Define parser object
    parser = ArgumentParser(description="")
    parser.add_argument("-v", "--verbose", action="store_true", dest="verbose",
                        default=False, help="Be chatty! (default = False)")
    parser.add_argument("-q", "--quiet", action="store_true", dest="quiet",
                        default=False, help="Be quiet! (default = False)")

    parser.add_argument("-n", "--name", dest="name", required=True, type=str,
                        help="The name of this reduction pipeline run.")
    parser.add_argument("--rootpath", dest="rootpath", required=True,
                        type=str, help="Path to root directory containing data files.")

    args = parser.parse_args()

    # Set logger level based on verbose flags
    if args.verbose:
        logger.setLevel(logging.DEBUG)
    elif args.quiet:
        logger.setLevel(logging.ERROR)
    else:
        logger.setLevel(logging.INFO)

    main(name=args.name, rootpath=args.rootpath)
