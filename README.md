# graph-ISS
Graph-based Decoding for In Situ Sequencing (ISS).

This repository contains the primary source code implementation of the graph-based image analysis pipeline for processing in situ sequencing data and ipython notebooks for reproducing publication analysis results and figures [ref].
The image decoding pipeline consists in three Python 3 library packages for 2D and 3D data proccesing and a Anduril2 [2] pipeline implementing the decoding workflow.

[2] *Cervera, Alejandra, et al. "Anduril 2: upgraded large-scale data integration framework." Bioinformatics (2019).*

### Decoding Pipeline Install Requirements
#### Anduril2
Anduril 2 is a workflow platform for high-throughput analysis of biomedical data. Workflows are constructed using Scala 2.11 and executed in parallel locally or on Linux clusters using a workflow engine optimized for iterative development. Documentation and installation instructions are available at: http://anduril.org.

#### Bio-Format Command Line Tools
Bio-format command line tools are necessary for dividing whole slide images in smaller tiles for faster computation. Bio-format command line tools can be downloaded from https://www.openmicroscopy.org/bio-formats.

#### Python Library Requirements
Create a conda evironment named "pgm_pipeline":

``` $ conda create --name pgm_pipeline```  

Activate environment:

``` $ conda activate pgm_pipeline```

Install the following python packages:  
  - `joblib>=0.13.2`
  - `keras>=2.2.4`
  - `networkx>=2.3`
  - `numpy>=1.13.1`
  - `pandas>=0.23.4`
  - `scikit-image>=0.13.0`
  - `scikit-learn>=0.21.3`
  - `scipy>=0.19.1`
  - `pytables>=3.4.2`
  - `tqdm>=4.32.2`
  
Install `SimpleElastix` inside the virtual environment following installation instructions and documentation available at https://simpleelastix.readthedocs.io.

To deactivate the conda enviroment:

``` $ conda deactivate```

### Analysis Notebook Install Requirements
The following python packages are required for running the notebooks:
  - `joblib>=0.13.2`
  - `matplotlib>=2.2.2`
  - `networkx>=2.3`
  - `nimfa>=1.3.4`
  - `numpy>=1.13.1`
  - `opencv-python>=3.4.1.15`
  - `pandas>=0.23.4`
  - `scikit-image>=0.13.0`
  - `scikit-learn>=0.21.3`
  - `scipy>=0.19.1`
  - `seaborn>=0.9.0`
  - `SpatialDE>=1.1.3`
  - `tqdm>=4.32.2`
  - `umap-learn>=0.3.9`

### Data Download
An example ISS data [3] for testing Anduril decoding pipeline and decoding results for reproducing publication analyses can be downloaded from: https://doi.org/10.5281/zenodo.3357950.

[3] *Ke, Rongqin, et al. "In situ sequencing for RNA analysis in preserved tissue and cells." Nature methods 10.9 (2013): 857.*

### Anduril Pipeline Example Usage
An example of 2D Anduril decoding pipeline is availabel for testing (`ISS_Anduril_Pipeline_Example.scala`). For running the test example, `<GRAPH-ISS-FOLDER>` and `<DATA-FOLDER>` strings in the scala file should be replaced respectively with graph-iss and downloaded data folder local paths.

To lunch the execution, run the command:

` $ ./ISS_Anduril_Pipeline_Example.scala`
