TODO:
* Need a way to keep track of which files / objects are done
* Maybe initially, need to search over all files and make a database?
* Store stuff like that in <top-level>/longslit/?
* For a single observing run with two instrument calibrations, could have
  two separate "global" configuration files - just name them differently in
  <top-level>/longslit/... and specify different path include/exclude
  patterns
* With no filename, ``init_wavelength.py`` should just find a COMP FITS file
  automatically and use that. Can be overridden by specifying a path to a
  specific filename.
  * GUI mode should be like Glue - select in pixel direction (snap to pixels)
    over a line, then type in the peak wavelength. Under the hood, fit a
    Gaussian to find centroid, then display Gaussian over line??
  * If there are multiple maxima, reject the selection and raise an error?
  * For now, user specifies range of columns and just roll with that - deal
    with the full 2D solution later?
  * If no line list specified, add a text box to let the user add lines
  * Add a command-line option to specify rows/columns to median for the 1D trace
    (maybe --aperture1d??)
* To be robust to cosmic rays, need to operate on masked arrays? Sigma clipping?
