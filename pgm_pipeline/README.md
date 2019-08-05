# Python preprocessing library for graph-based decoding.

Python script for signal candidate decoding through a graphical model and quality assesment.

Usage:

` $ python <LIBRARY-PATH>/main.py <INPUT-CSV1> <OUT-FOLDER> <OUT-FILE1> <OUT-FILE2> <NUM-THREADS> <H-THRESH> <INPUT-CSV2> <OUT-FILE3> <OUT-FILE4> <REG-TYPE>`

  - `<LIBRARY-PATH>`: Path to decoding library
  - `<INPUT-FILE1>`: Dataframe of signal probability predictions
  - `<OUT-FILE1>`: Output dataframe of decoded barcodes
  - `<OUT-CSV1>`: Output CSV of decoded barcodes
  - `<NUM-THREADS>`: Number of thread used for parallelization
  - `<D-TH>`: d<sub>th</sub> parameter for building graph connected components
  - `<D-MAX>`: d<sub>max</sub> parameter for building graph connected components
  - `<OUT-FILE2>`: Output hdf5 storing a list of sets of signals detections contributing to the decoded barcodes
  - `<SEARCH-TYPE>`: Enable barcodes decoding guided by the list of targeted barcodes (valid arguments, "prior": feature enabled, other string: feature disabled)
  - `<INPUT-CSV>`: Input CSV file of targeted barcode sequences (codebook)
  - `<INPUT-FILE2>`: Input dataframe of signal candidate detection before merging
  - `<INPUT-FILE3>`: Input hdf5 object storing normalized images
  - `<NETWORK-LIB-PATH>`: Path to signal probability prediction library
