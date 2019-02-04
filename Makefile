rootPath = ./
include ./include.mk

libSources = impl/*.c
libHeaders = inc/*.h
libTests = tests/*.c

signalAlignDependencies =  ${basicLibsDependencies}
signalAlignLib = ${basicLibs}

test_directory = ${rootPath}/src/signalalign/tests/
scrappie_build = ${rootPath}/scrappie/build

LIBS= -lhts -lhdf5

all : sL bD ${libPath}/signalAlignLib.a ${signalAlignBin}/signalAlignLibTests \
	  ${signalAlignBin}/compareDistributions ${signalAlignBin}/kmerEventAlign \
	  ${signalAlignBin}/signalMachine ${signalAlignBin}/runSignalAlign \
	  ${signalAlignBin}/variantCallingLib.py ${signalAlignBin}/alignmentAnalysisLib.py \
	  ${signalAlignBin}/buildHdpUtil ${signalAlignBin}/trainModels all_tests \
	  externals python_setup ${signalAlignBin}/filterReads ${signalAlignBin}/extract \
	  ${signalAlignBin}/sequencing_summary ${signalAlignBin}/plot_kmer_distributions \
	  ${signalAlignBin}/plot_variant_accuracy ${signalAlignBin}/compare_trained_models


#${rootPath}/lib/libhdf5.a:
#	if [ ! -e hdf5-1.10.2.tar.gz ]; then wget https://support.hdfgroup.org/ftp/HDF5/releases/hdf5-1.10/hdf5-1.10.2/src/hdf5-1.10.2.tar.gz; fi
#	tar -xzf hdf5-1.10.2.tar.gz || exit 255
#	cd hdf5-1.10.2 && ./configure --enable-threadsafe --disable-hl --prefix=`pwd`/.. && make && make install

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

python_setup :
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
		python $$i; \
		if [ $$? -ne 0 ]; then\
			echo "\nTEST FAIL $$i\n";\
			export FAIL=1;\
		fi;\
	done; \
	if [ $$FAIL -ne 0 ]; then exit -1; fi;


${signalAlignBin}/compareDistributions : compareDistributions.c ${libPath}/signalAlignLib.a ${signalAlignDependencies}
	${cxx} ${cflags}  -I inc -I${libPath} -o ${signalAlignBin}/compareDistributions compareDistributions.c ${libPath}/signalAlignLib.a ${signalAlignLib} ${LIBS}

${signalAlignBin}/signalAlignLibTests : ${libTests} tests/*.h ${libPath}/signalAlignLib.a ${signalAlignDependencies} .FORCE
	${cxx} ${cflags}  -I inc -I${libPath}   -Wno-error -o ${signalAlignBin}/signalAlignLibTests ${libTests} ${libPath}/signalAlignLib.a ${signalAlignLib}  ${LIBS}

${signalAlignBin}/signalMachine : signalMachine.c ${libPath}/signalAlignLib.a ${signalAlignDependencies}
	${cxx} ${cflags}  -I inc -I${libPath}   -o ${signalAlignBin}/signalMachine signalMachine.c ${libPath}/signalAlignLib.a ${signalAlignLib} ${LIBS}

${signalAlignBin}/kmerEventAlign : kmerEventAlign.c ${libPath}/signalAlignLib.a ${signalAlignDependencies}
	${cxx} ${cflags}  -I inc -I${libPath}   -o ${signalAlignBin}/kmerEventAlign kmerEventAlign.c ${libPath}/signalAlignLib.a ${signalAlignLib}    ${LIBS}

${signalAlignBin}/extract : extract.c ${libPath}/signalAlignLib.a ${signalAlignDependencies}
	${cxx} ${cflags}  -I inc -I${libPath}   -o ${signalAlignBin}/extract extract.c ${libPath}/signalAlignLib.a ${signalAlignLib}    ${LIBS}

#nanoporeParams : estimateNanoporeParams.c ${libPath}/signalAlignLib.a ${signalAlignDependencies}
#	${cxx} ${cflags}  -I inc -I${libPath} -o ${signalAlignBin}/estimateNanoporeParams estimateNanoporeParams.c ${libPath}/signalAlignLib.a ${signalAlignLib} ${LIBS}
#	cp ${rootPath}src/signalalign/scripts/nanoporeParamRunner.py ${signalAlignBin}/nanoporeParamRunner
#	chmod +x ${signalAlignBin}/nanoporeParamRunner

${signalAlignBin}/buildHdpUtil : buildHdpUtil.c ${libPath}/signalAlignLib.a ${signalAlignDependencies}
	${cxx} ${cflags}   -I inc -I${libPath} -o ${signalAlignBin}/buildHdpUtil buildHdpUtil.c ${libPath}/signalAlignLib.a ${signalAlignLib} ${LIBS}

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

${signalAlignBin}/plot_variant_accuracy : ${rootPath}src/signalalign/visualization/plot_variant_accuracy.py
	cp ${rootPath}src/signalalign/visualization/plot_variant_accuracy.py ${signalAlignBin}/plot_variant_accuracy
	chmod +x ${signalAlignBin}/plot_variant_accuracy

${signalAlignBin}/compare_trained_models : ${rootPath}src/signalalign/visualization/compare_trained_models.py
	cp ${rootPath}src/signalalign/visualization/compare_trained_models.py ${signalAlignBin}/compare_trained_models
	chmod +x ${signalAlignBin}/compare_trained_models

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
	${cxx} ${cflags} -fPIC -Iinc/ -I${libPath}/  -c ${libSources}  ${LIBS}
	ar rc signalAlignLib.a *.o
	ranlib signalAlignLib.a
	rm *.o
	mv signalAlignLib.a ${libPath}/
	cp ${libHeaders} ${libPath}/


.FORCE: python_setup
