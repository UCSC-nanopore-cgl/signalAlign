#!/bin/sh
# setup environment
cwd=`pwd`
export PATH="$PATH:$cwd/src/cmake/bin"

# ubuntu install
sudo apt-get update
sudo apt-get install -y --no-install-recommends autoconf wget git curl build-essential libbz2-dev zlib1g-dev liblzma-dev libeigen3-dev libreadline-gplv2-dev libncursesw5-dev libssl-dev libsqlite3-dev tk-dev libgdbm-dev libc6-dev libcurl4-openssl-dev ca-certificates python3.7-dev python3-pip python3.7-venv samtools

# create src dir
mkdir -p src
cd src

# install cmake
mkdir -p cmake && cd cmake
wget https://cmake.org/files/v3.17/cmake-3.17.0-Linux-x86_64.sh --no-check-certificate
sh cmake-3.17.0-Linux-x86_64.sh --skip-license

# install htslib
cd ..
wget https://github.com/samtools/htslib/releases/download/1.9/htslib-1.9.tar.bz2 --no-check-certificate
tar -vxjf htslib-1.9.tar.bz2
rm htslib-1.9.tar.bz2
cd htslib-1.9
./configure --prefix /usr/local --enable-plugins CPPFLAGS="-fPIC" CFLAGS="-fPIC"
make
sudo make install

# boost install
cd ..
wget -O boost_1_69_0.tar.gz https://sourceforge.net/projects/boost/files/boost/1.69.0/boost_1_69_0.tar.gz/download --no-check-certificate
tar -xzf boost_1_69_0.tar.gz >/dev/null
rm boost_1_69_0.tar.gz
cd boost_1_69_0
./bootstrap.sh --with-libraries=system,date_time,filesystem,iostreams,coroutine,context,regex,thread,atomic >/dev/null
./b2 -d0 cxxflags="-fPIC" cflags="-fPIC" link=static -a
sudo ./b2 -d0 install

# hdf5 install
cd ..
wget -q https://support.hdfgroup.org/ftp/HDF5/releases/hdf5-1.10/hdf5-1.10.4/src/hdf5-1.10.4.tar.gz
tar -xzf hdf5-1.10.4.tar.gz
cd hdf5-1.10.4
./configure --enable-threadsafe --disable-hl --prefix=/usr/local/
make
sudo make install

# setup python
cd ..
curl https://bootstrap.pypa.io/get-pip.py -o get-pip.py
python3.7 get-pip.py
rm get-pip.py
python3.7 -m pip install --upgrade pip
python3.7 -m pip -q install setuptools cython virtualenv pytest

# install embed
git clone --recursive https://github.com/adbailey4/embed_fast5.git --branch bailey-dev
cd embed_fast5
python3.7 -m pip install .
python3.7 -m pytest

# install bwa
cd ..
git clone https://github.com/lh3/bwa.git
cd bwa
make

# install signalAlign
cd ..
git clone --recursive https://github.com/UCSC-nanopore-cgl/signalAlign.git --branch bailey-dev
cd signalAlign && mkdir -p build && cd build
cmake .. -DCMAKE_INSTALL_PREFIX=. -DCMAKE_VERBOSE_MAKEFILE=ON -DBUILD_SHARED_LIBS=ON -DCMAKE_BUILD_TYPE=RELEASE
make -j 4
cd ..
python3.7 -m pip install .
python3.7 -m pytest

# install vbz_compression
cd .. && git clone https://github.com/nanoporetech/vbz_compression.git
cd vbz_compression
git submodule update --init
mkdir build
cd build
cmake -D CMAKE_BUILD_TYPE=Release -D ENABLE_CONAN=OFF -D ENABLE_PERF_TESTING=OFF -D ENABLE_PYTHON=OFF ..
make -j 4
sudo make install
cd ../../..

