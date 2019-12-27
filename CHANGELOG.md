# Changelog

# 0.3.0 [27.Dec.2019]

## Added
- Ability to read in new fastspec header
- New format: npz
- Ability to save all metadata

## Changed
- Structure of the code is now class-based and extendible
- **No longer outputs uncalibrated temperature, but rather the dimensionless ratio, Q**

## Fixed
- Line-incrementing problem when reading


# 0.2.0

## Added
- HDF5 output format
- Ability to specify output format via CLI
- Metadata to outputs

## Fixed
- Removed warnings when dividing by zero
- Installation without numpy
- Progress bar iterations

# 0.1.0

- Initial version that reads ACQ into memory.
