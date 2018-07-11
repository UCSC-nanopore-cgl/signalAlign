rootPath = ./
include ./include.mk

libSources = impl/*.c
libHeaders = inc/*.h
libTests = tests/*.c

signalAlignDependencies =  ${basicLibsDependencies}
signalAlignLib = ${basicLibs}

test_directory = ${rootPath}/src/signalalign/tests/
scrappie_build = ${rootPath}/scrappie/build

htsLib = -L././htslib -lhts
LIBS= -lsz -lz -lm

HDF5?=install
# Default to automatically installing hdf5
ifeq ($(HDF5), install)
    H5_LIB=./lib/libhdf5_hl.a ./lib/libhdf5.a
    H5_INCLUDE=-I./include
    LIBS += -ldl
else
    # Use system-wide hdf5
    H5_LIB=
    H5_INCLUDE=
    LIBS += -lhdf5
endif


all : sL bD hs python-utils ${libPath}/signalAlignLib.a ${signalAlignBin}/signalAlignLibTests \
	  ${signalAlignBin}/compareDistributions \
	  ${signalAlignBin}/signalMachine ${signalAlignBin}/runSignalAlign \
	  ${signalAlignBin}/variantCallingLib.py ${signalAlignBin}/alignmentAnalysisLib.py \
	  ${signalAlignBin}/buildHdpUtil ${signalAlignBin}/trainModels all_tests \
	  externals nanoporeParams python_setup ${scrappie_build}/scrappie \

python-utils :
	cd python_utils && python3 setup.py install


${rootPath}/lib/libhdf5.a:
	if [ ! -e hdf5-1.8.14.tar.gz ]; then wget https://support.hdfgroup.org/ftp/HDF5/releases/hdf5-1.8/hdf5-1.8.14/src/hdf5-1.8.14.tar.gz; fi
	tar -xzf hdf5-1.8.14.tar.gz || exit 255
	cd hdf5-1.8.14 && ./configure --enable-threadsafe --prefix=`pwd`/.. && make && make install


${scrappie_build}/scrappie :
	cd scrappie && \
    if ! [ -d ${scrappie_build} ]; then \
        mkdir build; \
    fi; \
	cd build && \
	cmake .. && \
	make

core : sL bD ${libPath}/signalAlignLib.a ${signalAlignBin}/signalAlignLibTests ${signalAlignBin}/signalMachine

install: all pip_install

clean_light:
	if [ -d ${signalAlignBin} ]; then rm -r ${signalAlignBin}; fi
	rm -f ${libPath}/signalAlignLib.a

clean :
	if [ -d ${signalAlignBin} ]; then rm -r ${signalAlignBin}; fi
	rm -f ${libPath}/signalAlignLib.a
	cd externalTools && make clean
	cd scrappie && make clean

python_setup :
	python3 setup.py install

pip_install : .FORCE
	pip3 install -e .

signalAlignLib : ${libPath}/signalAlignLib.a

sL :
	cd sonLib && CFLAGS="${CFLAGS} -fPIC" make

bD :
	mkdir -v -p ${rootPath}bin

externals :
	cd externalTools && make all

test_files := $(shell find $(test_directory) -name '*.py')

test :
	for i in ${test_files}; do \
		python $$i; \
		if [ $$? -ne 0 ]; then\
		exit -1;\
		fi;\
	done
#	cd ${binPath} && ./sonLibTests
	cd python_utils && pytest
	cd scrappie && make test

# //		exit "$$?"; \

${signalAlignBin}/compareDistributions : compareDistributions.c ${libPath}/signalAlignLib.a ${signalAlignDependencies}
	${cxx} ${cflags}  -I inc -I${libPath} -o ${signalAlignBin}/compareDistributions compareDistributions.c ${libPath}/signalAlignLib.a ${signalAlignLib}

${signalAlignBin}/signalAlignLibTests : ${libTests} tests/*.h ${libPath}/signalAlignLib.a ${signalAlignDependencies} ${rootPath}/lib/libhdf5.a .FORCE
	${cxx} ${cflags}  -I inc -I${libPath} ${H5_INCLUDE} -I${htsLibRootPath} -I${htsLibPath} -Wno-error -o ${signalAlignBin}/signalAlignLibTests ${libTests} ${libPath}/signalAlignLib.a ${signalAlignLib} ${H5_LIB} ${LIBS} ${htsLib}

${signalAlignBin}/signalMachine : signalMachine.c ${libPath}/signalAlignLib.a ${signalAlignDependencies}
	${cxx} ${cflags}  -I inc -I${libPath} -I${htsLibRootPath} -I${htsLibPath} -o ${signalAlignBin}/signalMachine signalMachine.c ${libPath}/signalAlignLib.a ${signalAlignLib}  ${htsLib}

nanoporeParams : estimateNanoporeParams.c ${libPath}/signalAlignLib.a ${signalAlignDependencies}
	${cxx} ${cflags}  -I inc -I${libPath} -o ${signalAlignBin}/estimateNanoporeParams estimateNanoporeParams.c ${libPath}/signalAlignLib.a ${signalAlignLib}
	cp ${rootPath}src/signalalign/scripts/nanoporeParamRunner.py ${signalAlignBin}/nanoporeParamRunner
	chmod +x ${signalAlignBin}/nanoporeParamRunner

${signalAlignBin}/buildHdpUtil : buildHdpUtil.c ${libPath}/signalAlignLib.a ${signalAlignDependencies}
	${cxx} ${cflags}   -I inc -I${libPath} -o ${signalAlignBin}/buildHdpUtil buildHdpUtil.c ${libPath}/signalAlignLib.a ${signalAlignLib}

${signalAlignBin}/runSignalAlign : ${rootPath}src/signalalign/scripts/runSignalAlign.py
	cp ${rootPath}src/signalalign/scripts/runSignalAlign.py ${signalAlignBin}/runSignalAlign
	chmod +x ${signalAlignBin}/runSignalAlign

${signalAlignBin}/trainModels : ${rootPath}src/signalalign/train/trainModels.py
	cp ${rootPath}src/signalalign/train/trainModels.py ${signalAlignBin}/trainModels
	chmod +x ${signalAlignBin}/trainModels

all_tests : .FORCE

	chmod +x ${test_directory}/*

${signalAlignBin}/zayante : ${rootPath}src/signalalign/scripts/zayante.py
	cp ${rootPath}src/signalalign/scripts/zayante.py ${signalAlignBin}/zayante
	chmod +x ${signalAlignBin}/zayante

${signalAlignBin}/empire : ${rootPath}src/signalalign/scripts/empire.py
	cp ${rootPath}src/signalalign/scripts/empire.py ${signalAlignBin}/empire
	chmod +x ${signalAlignBin}/empire

${signalAlignBin}/variantCallingLib.py : ${rootPath}src/signalalign/scripts/variantCallingLib.py
	cp ${rootPath}src/signalalign/scripts/variantCallingLib.py ${signalAlignBin}/variantCallingLib.py

${signalAlignBin}/alignmentAnalysisLib.py : ${rootPath}src/signalalign/scripts/alignmentAnalysisLib.py
	cp ${rootPath}src/signalalign/scripts/alignmentAnalysisLib.py ${signalAlignBin}/alignmentAnalysisLib.py

${libPath}/signalAlignLib.a : ${libSources} ${libHeaders} ${stBarDependencies} ${rootPath}/lib/libhdf5.a
	${cxx} ${cflags} -fPIC -I inc -I ${libPath}/ ${H5_INCLUDE} -I ${htsLibRootPath} -I ${htsLibPath}  ${htsLib} -c ${libSources} ${H5_LIB} ${LIBS}
	ar rc signalAlignLib.a *.o
	ranlib signalAlignLib.a
	rm *.o
	mv signalAlignLib.a ${libPath}/
	cp ${libHeaders} ${libPath}/

hs :
	cd htslib && make


.FORCE:
