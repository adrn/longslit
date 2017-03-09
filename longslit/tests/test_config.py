# Standard library
import os
import sys

# Third-party
from astropy import log as logger
import astropy.units as u
import numpy as np
import pytest

# Project
from ..config import Config

def test_config_one(tmpdir):
    """ only global config """

    # first set up the directory structure
    root_path = tmpdir / 'run'
    n1 = root_path / 'n1'
    n2 = root_path / 'n2'

    # make a config file

def test_config_two():
    # global config, overwrite some settings with subdir config
    pass
