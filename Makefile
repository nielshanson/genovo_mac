LIBS =-lgsl -lgslcblas

LIBTOOL=/opt/local/bin/glibtool
LIBDIRS=/opt/local/lib/
CPP = $(LIBTOOL) --mode=compile --tag=cc g++ -c -I${INCLUDE_PATH} ${INCDIRS} ${LIBDIRS} -O3
LINK = $(LIBTOOL) --mode=link --tag=cc g++ -I${INCLUDE_PATH} ${MPIFLAGS}  ${INCDIRS} -static ${LIBDIRS} -O3

BINARIES = finalize compute_score_denovo assemble

.SILENT:

all: $(BINARIES)

%:	%.cc libjv.la
	$(LINK) -o $@ $@.cc libjv.la -I $(LIBDIRS) $(INCDIRS) $(LIBS)

deploy: all
	mkdir -p ../workDir
	cp $(BINARIES) ../workDir
tidy:
	rm -f *.lo *.o *.la *.a .libs/* 

test_%: test_%.cc libjv.la
	$(LINK) -o $@ $@.cc libjv.la -I $(LIBDIRS) $(INCDIRS) $(LIBS)

clean: 
	rm -f *.lo *.o *.la *.a .libs/* $(TESTS)

jvDataGenerator.lo: jvDataGenerator.h jvDataGenerator.cc jvCommon.h
	$(CPP)  jvDataGenerator.cc 

jvSampler.lo: jvSampler.h jvSampler.cc jvCommon.h
	$(CPP)  jvSampler.cc

jvState.lo: jvState.h jvState.cc jvCommon.h
	$(CPP)  jvState.cc

jvFastaParser.lo: jvFastaParser.h jvFastaParser.cc jvCommon.h
	$(CPP)  jvFastaParser.cc

jvStateSerializer.lo: jvStateSerializer.cc jvStateSerializer.h jvCommon.h
	$(CPP)  jvStateSerializer.cc

jvAlign.lo: jvAlign.cc jvAlign.h jvCommon.h
	$(CPP)  jvAlign.cc

jvAlign2.lo: jvAlign2.cc jvAlign2.h jvCommon.h
	$(CPP)  jvAlign2.cc

libjv.la: jvState.lo jvSampler.lo jvDataGenerator.lo jvFastaParser.lo jvStateSerializer.lo jvAlign.lo jvAlign2.lo
	$(LINK)   -o libjv.la jvState.lo jvSampler.lo jvDataGenerator.lo jvFastaParser.lo jvStateSerializer.lo jvAlign.lo jvAlign2.lo
