Changelog
=========

<!--next-version-placeholder-->

0.4.4 \[04.June.2021\]
----------------------

### Fixed

-   HDF5 to ACQ converter now works.

0.4.0
-----

First version on PyPI!

### Fixed

-   Frequency array now has the correct channel width (it was ever so
    slightly off).

### Changed

-   `decode_file` now only decodes the file, rather than writing. Use
    `convert_file` to convert a file to a different format.

0.3.2
-----

### Fixed

-   Errors in writing out files (now have tests)

0.3.1
-----

### Fixed

-   Error reading fastspec header due to first VERSION line.

0.3.0 \[27.Dec.2019\]
---------------------

### Added

-   Ability to read in new fastspec header
-   New format: npz
-   Ability to save all metadata

### Changed

-   Structure of the code is now class-based and extendible
-   **No longer outputs uncalibrated temperature, but rather the
    dimensionless ratio, Q**

### Fixed

-   Line-incrementing problem when reading

0.2.0
-----

### Added

-   HDF5 output format
-   Ability to specify output format via CLI
-   Metadata to outputs

### Fixed

-   Removed warnings when dividing by zero
-   Installation without numpy
-   Progress bar iterations

0.1.0
-----

-   Initial version that reads ACQ into memory.
