CC = g++
CCOPTS = -Wall -ansi -O2
CCLIBS = -L. -lprofilehmm -lboost_program_options -lgsl -lgslcblas
OBJS = phmmtools-common.o findcluster.o multialign.o
PHMMLIB = libprofilehmm.a
PHMMSRC = profilehmm/backward.cpp profilehmm/baum-welch.cpp profilehmm/extendedexp.hpp profilehmm/extendedexp.cpp profilehmm/forward.cpp profilehmm/maintenance.cpp profilehmm/profilehmm.hpp profilehmm/scoring.cpp profilehmm/settings.hpp profilehmm/settings.cpp profilehmm/state.cpp profilehmm/viterbi.cpp
TAR = phmmtools.tar.bz2

all: findcluster multialign

findcluster: findcluster.o phmmtools-common.o $(PHMMLIB)
	$(CC) $(CCOPTS) -o $@ $< phmmtools-common.o $(CCLIBS)

multialign: multialign.o phmmtools-common.o $(PHMMLIB)
	$(CC) $(CCOPTS) -o $@ $< phmmtools-common.o $(CCLIBS)

findcluster.o: findcluster.cpp
	$(CC) $(CCOPTS) -c -o $@ $<

phmmtools-common.o: phmmtools-common.cpp phmmtools-common.hpp
	$(CC) $(CCOPTS) -c -o $@ $<

multialign.o: multialign.cpp
	$(CC) $(CCOPTS) -c -o $@ $<

$(PHMMLIB): $(PHMMSRC)
	$(MAKE) -C profilehmm clean
	$(MAKE) -C profilehmm
	mv profilehmm/$(PHMMLIB) ./

tar: findcluster.cpp gpl makefile multialign.cpp phmmtools-common.cpp phmmtools-common.hpp
	tar cfj $(TAR) *.cpp *.hpp makefile gpl

clean:
	rm -f findcluster multialign $(PHMMLIB) core $(OBJS) $(TAR)
	$(MAKE) -C profilehmm clean

tidy:
	rm -f $(OBJS) core
	$(MAKE) -C profilehmm tidy

