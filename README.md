# read_acq

Read ACQ files for EDGES.


## Installation

At this time, you will first need to ensure that `numpy` is installed (`pip install numpy`)
or (`conda install numpy`).

Then, you can simply run `pip install git+https://github.com/edges-collab/read_acq`.

If you wish to develop `read_acq`, do the following::

    conda install numpy
    git clone https://github.com/edges-collab/read_acq
    cd read_acq
    pip install -e .
    
## Usage

`read_acq` can be used in a Python interpreter or directly via the command line. 

To use the CLI, a single command is provided:

    acq convert <data.acq> [<data2.acq> ...]
    
This will convert the file provided to a `.mat` format, and place the resulting `.mat`
file in the same location as the original datafile (but with different extension). 
The command can be run from anywhere on the system, and the file given can be a 
relative or absolute path. 

Multiple data files will be given, and each will be converted. Wildcards may also be
used in any of the filenames, eg.:

    acq convert data/*.acq
    
The main point of entry for the Python interface is `decode_file`:

    >>> from read_acq import decode_file
    >>> data = decode_file("my_data.acq", write_formats=[])
    
By default, this function will also write the file in `.mat` format, but you can turn
that off by providing `write_formats=[]`. The output is a `numpy` array of the data.
Several more options are provided, use `help(decode_file)` in an interpreter to see
all the options. 