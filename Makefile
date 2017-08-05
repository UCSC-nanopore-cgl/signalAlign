rootPath = ./
include ./include.mk

libSources = impl/*.c
libHeaders = inc/*.h
libTests = tests/*.c

signalAlignDependencies =  ${basicLibsDependencies}
signalAlignLib = ${basicLibs}

all : sL bD mA ${libPath}/signalAlignLib.a ${signalAlignBin}/signalAlignLibTests ${signalAlignBin}/compareDistributions \
	  ${signalAlignBin}/signalMachine ${signalAlignBin}/runSignalAlign \
	  ${signalAlignBin}/signalAlignLib.py ${signalAlignBin}/variantCallingLib.py ${signalAlignBin}/alignmentAnalysisLib.py \
	  ${signalAlignBin}/buildHdpUtil ${signalAlignBin}/trainModels ${signalAlignBin}/hdp_pipeline ${signalAlignBin}/testSignalAlign \
	  externals nanoporeParams
	  #${signalAlignBin}/zayante ${signalAlignBin}/bonnyDoon \
	  #${signalAlignBin}/empire ${signalAlignBin}/jamison \

core : sL bD ${libPath}/signalAlignLib.a ${signalAlignBin}/signalAlignLibTests ${signalAlignBin}/signalMachine

install: all pip_install

clean_light:
	if [ -d ${signalAlignBin} ]; then rm -r ${signalAlignBin}; fi
	rm -f ${libPath}/signalAlignLib.a

clean :
	if [ -d ${signalAlignBin} ]; then rm -r ${signalAlignBin}; fi
	#rm -r ${signalAlignBin}
	rm -f ${libPath}/signalAlignLib.a
	cd externalTools && make clean

pip_install:
	pip install -e .

signalAlignLib : ${libPath}/signalAlignLib.a

sL :
	cd sonLib && make

mA :
	cd marginAlign && make


bD :
	mkdir -v -p ${rootPath}bin

externals :
	cd externalTools && make all

test :
	#cd ${signalAlignBin} && ./signalAlignLibTests
	cd ${signalAlignBin} && ./testSignalAlign
	#cd ${binPath} && ./sonLibTests

${signalAlignBin}/compareDistributions : compareDistributions.c ${libPath}/signalAlignLib.a ${signalAlignDependencies}
	${cxx} ${cflags} -I inc -I${libPath} -o ${signalAlignBin}/compareDistributions compareDistributions.c ${libPath}/signalAlignLib.a ${signalAlignLib}

${signalAlignBin}/signalAlignLibTests : ${libTests} tests/*.h ${libPath}/signalAlignLib.a ${signalAlignDependencies}
	${cxx} ${cflags} -I inc -I${libPath} -Wno-error -o ${signalAlignBin}/signalAlignLibTests ${libTests} ${libPath}/signalAlignLib.a ${signalAlignLib}

${signalAlignBin}/signalMachine : signalMachine.c ${libPath}/signalAlignLib.a ${signalAlignDependencies}
	${cxx} ${cflags} -I inc -I${libPath} -o ${signalAlignBin}/signalMachine signalMachine.c ${libPath}/signalAlignLib.a ${signalAlignLib}

nanoporeParams : estimateNanoporeParams.c ${libPath}/signalAlignLib.a ${signalAlignDependencies}
	${cxx} ${cflags} -I inc -I${libPath} -o ${signalAlignBin}/estimateNanoporeParams estimateNanoporeParams.c ${libPath}/signalAlignLib.a ${signalAlignLib}
	cp ${rootPath}src/signalalign/scripts/nanoporeParamRunner.py ${signalAlignBin}/nanoporeParamRunner
	chmod +x ${signalAlignBin}/nanoporeParamRunner

${signalAlignBin}/buildHdpUtil : buildHdpUtil.c ${libPath}/signalAlignLib.a ${signalAlignDependencies}
	${cxx} ${cflags} -I inc -I${libPath} -o ${signalAlignBin}/buildHdpUtil buildHdpUtil.c ${libPath}/signalAlignLib.a ${signalAlignLib}

${signalAlignBin}/runSignalAlign : ${rootPath}src/signalalign/scripts/runSignalAlign.py
	cp ${rootPath}src/signalalign/scripts/runSignalAlign.py ${signalAlignBin}/runSignalAlign
	chmod +x ${signalAlignBin}/runSignalAlign

${signalAlignBin}/trainModels : ${rootPath}src/signalalign/scripts/trainModels.py
	cp ${rootPath}src/signalalign/scripts/trainModels.py ${signalAlignBin}/trainModels
	chmod +x ${signalAlignBin}/trainModels

${signalAlignBin}/hdp_pipeline : ${rootPath}src/signalalign/scripts/hdp_pipeline.py
	cp ${rootPath}src/signalalign/scripts/hdp_pipeline.py ${signalAlignBin}/hdp_pipeline
	chmod +x ${signalAlignBin}/hdp_pipeline

${signalAlignBin}/testSignalAlign : ${rootPath}src/signalalign/scripts/testSignalAlign.py
	cp ${rootPath}src/signalalign/scripts/testSignalAlign.py ${signalAlignBin}/testSignalAlign
	chmod +x ${signalAlignBin}/testSignalAlign

${signalAlignBin}/zayante : ${rootPath}src/signalalign/scripts/zayante.py
	cp ${rootPath}src/signalalign/scripts/zayante.py ${signalAlignBin}/zayante
	chmod +x ${signalAlignBin}/zayante

${signalAlignBin}/bonnyDoon : ${rootPath}src/signalalign/scripts/bonnyDoon.py
	cp ${rootPath}src/signalalign/scripts/bonnyDoon.py ${signalAlignBin}/bonnyDoon
	chmod +x ${signalAlignBin}/bonnyDoon

${signalAlignBin}/empire : ${rootPath}src/signalalign/scripts/empire.py
	cp ${rootPath}src/signalalign/scripts/empire.py ${signalAlignBin}/empire
	chmod +x ${signalAlignBin}/empire

${signalAlignBin}/jamison : ${rootPath}src/signalalign/scripts/jamison.py
	cp ${rootPath}src/signalalign/scripts/jamison.py ${signalAlignBin}/jamison
	chmod +x ${signalAlignBin}/jamison

${signalAlignBin}/signalAlignLib.py : ${rootPath}src/signalalign/scripts/signalAlignLib.py
	cp ${rootPath}src/signalalign/scripts/signalAlignLib.py ${signalAlignBin}/signalAlignLib.py

${signalAlignBin}/variantCallingLib.py : ${rootPath}src/signalalign/scripts/variantCallingLib.py
	cp ${rootPath}src/signalalign/scripts/variantCallingLib.py ${signalAlignBin}/variantCallingLib.py

${signalAlignBin}/alignmentAnalysisLib.py : ${rootPath}src/signalalign/scripts/alignmentAnalysisLib.py
	cp ${rootPath}src/signalalign/scripts/alignmentAnalysisLib.py ${signalAlignBin}/alignmentAnalysisLib.py

${libPath}/signalAlignLib.a : ${libSources} ${libHeaders} ${stBarDependencies}
	${cxx} ${cflags} -I inc -I ${libPath}/ -c ${libSources}
	ar rc signalAlignLib.a *.o
	ranlib signalAlignLib.a
	rm *.o
	mv signalAlignLib.a ${libPath}/
	cp ${libHeaders} ${libPath}/
