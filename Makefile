rootPath = ./
include ./include.mk

libSources = impl/*.c
libHeaders = inc/*.h
libTests = tests/*.c

signalAlignDependencies =  ${basicLibsDependencies}
signalAlignLib = ${basicLibs}

test_directory = ${rootPath}/src/signalalign/tests/
scrappie_build = ${rootPath}/scrappie/build

LIBS=-lz -llzma -lbz2 -lcurl -lpthread -lcrypto
#-ldeflate

HDF5 ?= install
HTS ?= install

HDF5_VERSION ?= 1.10.4
tmp_LIBS =

# Default to automatically installing hdf5
ifeq ($(HDF5), install)
    H5_LIB += ./lib/libhdf5.a
    H5_INCLUDE = -I./include
    LIBS += -ldl
else
    # Use system-wide hdf5
    H5_LIB =
    H5_INCLUDE ?=
    LIBS += -lhdf5
endif

# Default to build and link the libhts submodule
ifeq ($(HTS), install)
    HTS_LIB += ./htslib/libhts.a
    HTS_INCLUDE = -I./htslib
else
    # Use system-wide htslib
    HTS_LIB =
    HTS_INCLUDE =
    LIBS += -lhts
endif



#
#
# If this library is a dependency the user wants HDF5 to be downloaded and built.
#

.PHONY: all
all : bD lib/libhdf5.a htslib/libhts.a sL ${libPath}/signalAlignLib.a ${signalAlignBin}/signalAlignLibTests \
	  ${signalAlignBin}/compareDistributions ${signalAlignBin}/kmerEventAlign \
	  ${signalAlignBin}/signalMachine ${signalAlignBin}/runSignalAlign \
	  ${signalAlignBin}/variantCallingLib.py ${signalAlignBin}/alignmentAnalysisLib.py \
	  ${signalAlignBin}/buildHdpUtil ${signalAlignBin}/trainModels all_tests \
	  externals ${signalAlignBin}/filterReads ${signalAlignBin}/extract \
	  ${signalAlignBin}/sequencing_summary ${signalAlignBin}/plot_kmer_distributions \
	  ${signalAlignBin}/plot_variant_accuracy ${signalAlignBin}/compare_trained_models \
	  ${signalAlignBin}/remove_sa_analyses ${signalAlignBin}/plot_labelled_read \
	  ${signalAlignBin}/plot_em_model_distributions

lib/libhdf5.a:
	if [ ! -e hdf5-$(HDF5_VERSION).tar.gz ]; then \
		version_major_minor=`echo "$(HDF5_VERSION)" | sed -E 's/\.[0-9]+$$//'`; \
		wget https://support.hdfgroup.org/ftp/HDF5/releases/hdf5-$${version_major_minor}/hdf5-$(HDF5_VERSION)/src/hdf5-$(HDF5_VERSION).tar.gz; \
	fi

	tar -xzf hdf5-$(HDF5_VERSION).tar.gz || exit 255
	cd hdf5-$(HDF5_VERSION) && \
		./configure --enable-threadsafe --disable-hl --libdir=`pwd`/../lib --includedir=`pwd`/../include --prefix=`pwd`/.. && \
		make && make install

# Build libhts
#
htslib/libhts.a:
	cd htslib && make || exit 255


cflags += $(H5_INCLUDE) $(HTS_INCLUDE)


#
#${scrappie_build}/scrappie :
#	cd scrappie && \
#    if ! [ -d ${scrappie_build} ]; then \
#        mkdir build; \
#    fi; \
#	cd build && \
#	cmake .. && \
#	make

core : sL bD ${libPath}/signalAlignLib.a ${signalAlignBin}/signalAlignLibTests ${signalAlignBin}/signalMachine

install: all pip_install

clean_light:
	if [ -d ${signalAlignBin} ]; then rm -r ${signalAlignBin}; fi
	rm -f ${libPath}/signalAlignLib.a

clean :
	if [ -d ${signalAlignBin} ]; then rm -r ${signalAlignBin}; fi
	if [ -d build/ ]; then rm -r build/; fi
	if [ -d lib/ ]; then rm -r lib/; fi
	if [ -d dist/ ]; then rm -r dist/; fi
	if [ -d eventdetection/lib/ ]; then rm -r eventdetection/lib/; fi
	if [ -d eventdetection/build/ ]; then rm -r eventdetection/build/; fi
	if [ -d eventdetection/dist/ ]; then rm -r eventdetection/dist/; fi

	rm -f ${libPath}/signalAlignLib.a
	cd externalTools && make clean
	cd sonLib && make clean

python_develop :
	which python3
	python3 setup.py develop

python_install :
	which python3
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
	export FAIL=0 ; \
	for i in ${test_files}; do \
		python3 $$i; \
		if [ $$? -ne 0 ]; then\
			echo "\nTEST FAIL $$i\n";\
			export FAIL=1;\
		fi;\
	done; \
	if [ $$FAIL -ne 0 ]; then exit -1; fi;


${signalAlignBin}/compareDistributions : compareDistributions.c ${libPath}/signalAlignLib.a ${signalAlignDependencies}
	${cxx} ${cflags}  -I inc -I${libPath} -o ${signalAlignBin}/compareDistributions compareDistributions.c ${libPath}/signalAlignLib.a ${signalAlignLib} ${HTS_LIB} ${H5_LIB} ${LIBS}

${signalAlignBin}/signalAlignLibTests : ${libTests} tests/*.h ${libPath}/signalAlignLib.a ${signalAlignDependencies} .FORCE
	${cxx} ${cflags}  -I inc -I${libPath}   -Wno-error -o ${signalAlignBin}/signalAlignLibTests ${libTests} ${libPath}/signalAlignLib.a ${signalAlignLib} ${HTS_LIB} ${H5_LIB} ${LIBS}

${signalAlignBin}/signalMachine : signalMachine.c ${libPath}/signalAlignLib.a ${signalAlignDependencies}
	${cxx} ${cflags}  -I inc -I${libPath}   -o ${signalAlignBin}/signalMachine signalMachine.c ${libPath}/signalAlignLib.a ${signalAlignLib}  ${HTS_LIB} ${H5_LIB} ${LIBS}

${signalAlignBin}/kmerEventAlign : kmerEventAlign.c ${libPath}/signalAlignLib.a ${signalAlignDependencies}
	${cxx} ${cflags}  -I inc -I${libPath}   -o ${signalAlignBin}/kmerEventAlign kmerEventAlign.c ${libPath}/signalAlignLib.a ${signalAlignLib}  ${HTS_LIB} ${H5_LIB}  ${LIBS}

${signalAlignBin}/extract : extract.c ${libPath}/signalAlignLib.a ${signalAlignDependencies}
	${cxx} ${cflags}  -I inc -I${libPath}   -o ${signalAlignBin}/extract extract.c ${libPath}/signalAlignLib.a ${signalAlignLib}   ${HTS_LIB} ${H5_LIB} ${LIBS}

#nanoporeParams : estimateNanoporeParams.c ${libPath}/signalAlignLib.a ${signalAlignDependencies}
#	${cxx} ${cflags}  -I inc -I${libPath} -o ${signalAlignBin}/estimateNanoporeParams estimateNanoporeParams.c ${libPath}/signalAlignLib.a ${signalAlignLib} ${LIBS}
#	cp ${rootPath}src/signalalign/scripts/nanoporeParamRunner.py ${signalAlignBin}/nanoporeParamRunner
#	chmod +x ${signalAlignBin}/nanoporeParamRunner

${signalAlignBin}/buildHdpUtil : buildHdpUtil.c ${libPath}/signalAlignLib.a ${signalAlignDependencies}
	${cxx} ${cflags}   -I inc -I${libPath} -o ${signalAlignBin}/buildHdpUtil buildHdpUtil.c ${libPath}/signalAlignLib.a ${signalAlignLib} ${HTS_LIB} ${H5_LIB} ${LIBS}

${signalAlignBin}/runSignalAlign : ${rootPath}src/signalalign/scripts/runSignalAlign.py
	cp ${rootPath}src/signalalign/scripts/runSignalAlign.py ${signalAlignBin}/runSignalAlign
	chmod +x ${signalAlignBin}/runSignalAlign

${signalAlignBin}/trainModels : ${rootPath}src/signalalign/train/trainModels.py
	cp ${rootPath}src/signalalign/train/trainModels.py ${signalAlignBin}/trainModels
	chmod +x ${signalAlignBin}/trainModels

${signalAlignBin}/sequencing_summary : ${rootPath}src/signalalign/visualization/sequencing_summary.py
	cp ${rootPath}src/signalalign/visualization/sequencing_summary.py ${signalAlignBin}/sequencing_summary
	chmod +x ${signalAlignBin}/sequencing_summary

${signalAlignBin}/plot_kmer_distributions : ${rootPath}src/signalalign/visualization/plot_kmer_distributions.py
	cp ${rootPath}src/signalalign/visualization/plot_kmer_distributions.py ${signalAlignBin}/plot_kmer_distributions
	chmod +x ${signalAlignBin}/plot_kmer_distributions

${signalAlignBin}/plot_em_model_distributions : ${rootPath}src/signalalign/visualization/plot_em_model_distributions.py
	cp ${rootPath}src/signalalign/visualization/plot_em_model_distributions.py ${signalAlignBin}/plot_em_model_distributions
	chmod +x ${signalAlignBin}/plot_em_model_distributions

${signalAlignBin}/plot_variant_accuracy : ${rootPath}src/signalalign/visualization/plot_variant_accuracy.py
	cp ${rootPath}src/signalalign/visualization/plot_variant_accuracy.py ${signalAlignBin}/plot_variant_accuracy
	chmod +x ${signalAlignBin}/plot_variant_accuracy

${signalAlignBin}/plot_labelled_read : ${rootPath}src/signalalign/visualization/plot_labelled_read.py
	cp ${rootPath}src/signalalign/visualization/plot_labelled_read.py ${signalAlignBin}/plot_labelled_read
	chmod +x ${signalAlignBin}/plot_labelled_read

${signalAlignBin}/compare_trained_models : ${rootPath}src/signalalign/visualization/compare_trained_models.py
	cp ${rootPath}src/signalalign/visualization/compare_trained_models.py ${signalAlignBin}/compare_trained_models
	chmod +x ${signalAlignBin}/compare_trained_models

${signalAlignBin}/remove_sa_analyses : ${rootPath}src/signalalign/remove_sa_analyses.py
	cp ${rootPath}src/signalalign/remove_sa_analyses.py ${signalAlignBin}/remove_sa_analyses
	chmod +x ${signalAlignBin}/remove_sa_analyses

all_tests : .FORCE

	chmod +x ${test_directory}/*

${signalAlignBin}/zayante : ${rootPath}src/signalalign/scripts/zayante.py
	cp ${rootPath}src/signalalign/scripts/zayante.py ${signalAlignBin}/zayante
	chmod +x ${signalAlignBin}/zayante

${signalAlignBin}/empire : ${rootPath}src/signalalign/scripts/empire.py
	cp ${rootPath}src/signalalign/scripts/empire.py ${signalAlignBin}/empire
	chmod +x ${signalAlignBin}/empire

${signalAlignBin}/variantCallingLib.py : ${rootPath}src/signalalign/scripts/variantCallingLib.py
	cp ${rootPath}src/signalalign/scripts/variantCallingLib.py ${signalAlignBin}/variantCallingLib.pyq

${signalAlignBin}/filterReads : ${rootPath}src/signalalign/filter_reads.py
	cp ${rootPath}src/signalalign/filter_reads.py ${signalAlignBin}/filterReads
	chmod +x ${signalAlignBin}/filterReads


${signalAlignBin}/alignmentAnalysisLib.py : ${rootPath}src/signalalign/scripts/alignmentAnalysisLib.py
	cp ${rootPath}src/signalalign/scripts/alignmentAnalysisLib.py ${signalAlignBin}/alignmentAnalysisLib.py

${libPath}/signalAlignLib.a : ${libSources} ${libHeaders} ${stBarDependencies}
	${cxx} ${cflags} -fPIC -Iinc/ -I${libPath}/  -c ${libSources} ${HTS_LIB} ${H5_LIB} ${LIBS}
	ar rc signalAlignLib.a *.o
	ranlib signalAlignLib.a
	rm *.o
	mv signalAlignLib.a ${libPath}/
	cp ${libHeaders} ${libPath}/


.FORCE: 
