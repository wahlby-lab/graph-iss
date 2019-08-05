# Python preprocessing library for 2D data

Python script for tile registration, normalization, signal candidate detection, and merging.

Usage:

` $ python <LIBRARY-PATH>/main.py <INPUT-CSV1> <OUT-FOLDER> <OUT-FILE1> <OUT-FILE2> <NUM-THREADS> <H-THRESH> <INPUT-CSV2> <OUT-FILE3> <OUT-FILE4> <REG-TYPE>`

  - `<LIBRARY-PATH>`: Path to pre-processing library
  - `<INPUT-CSV1>`: Input csv file contaning an Anduril file array of input image paths.
  - `<OUT-FOLDER>`: Output folder where to store aligned images
  - `<OUT-FILE1>`: Output dataframe of signal candidate detection after merging
  - `<OUT-FILE2>`: Output dataframe of signal candidate deteCtion before merging
  - `<NUM-THREADS>`: Number of thread used for parallelizing sequencing rounds
  - `<H-THRESH>`: h-maxima threshold
  - `<INPUT-CSV2>`: csv file containing image normalization intervals
  - `<OUT-FILE3>`: Output hdf5 object storing normalized images
  - `<OUT-FILE4>`: Output hdf5 object storing signal candidates binary masks
  - `<REG-TYPE>`: Input string setting the type of registration procedure. (valid arguments: "DO1" : if only the general stain of the first sequencing round is available, "DO" : if a general stain image is available for each sequencing round


Python script for signal probability prediction.

Usage:

` $ python <NETWORK-LIB-PATH>/main.py <INPUT-CSV1> <OUT-FOLDER> <OUT-FILE1> <OUT-FILE2> <NUM-THREADS> <H-THRESH> <INPUT-CSV2> <OUT-FILE3> <OUT-FILE4>

  - `<NETWORK-LIB-PATH>`: Path to network library
  - `<INPUT-FILE1>`: Input dataframe of signal candidate detection after merging
  - `<OUT-FILE1>`: Output dataframe of signal candidate predictions
  - `<P-TH>`: (Optional) probability threshold for signal candidate predictions
