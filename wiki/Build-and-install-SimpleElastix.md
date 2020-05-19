```
# Install SimpleElastix requirements
sudo apt-get install cmake r-base r-base-dev ruby ruby-dev tcl tcl-dev tk tk-dev

# Clone repository
git clone https://github.com/SuperElastix/SimpleElastix.git

# Create build folder
cd SimpleElastix
mkdir build
cd build

# Configure build
# Replace <path-to-miniconda-pgm_pipeline-env> with your pgm_pipeline environment path, and <x> with your python 3 version
cmake -DPYTHON_EXECUTABLE:FILEPATH=<path-to-miniconda-pgm_pipeline-env>/bin/python3 \
-DPYTHON_LIBRARY:FILEPATH=<path-to-miniconda-pgm_pipeline-env>/lib/libpython3.<x>m.so \
-DPYTHON_INCLUDE_DIR:PATH=<path-to-miniconda-pgm_pipeline-env>/include/python3.<x>m \
-DBUILD_EXAMPLES=OFF \
-DBUILD_TESTING=OFF \
../SuperBuild

# Build
make -j4

# Install python package
cd SimpleITK-build/Wrapping/Python/
python Packaging/setup.py install

```
