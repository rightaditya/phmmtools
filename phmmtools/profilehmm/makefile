TARGET = libprofilehmm.a
CC = g++
CCOPTS = -Wall -ansi -c -O2
OBJS = backward.o baum-welch.o extendedexp.o forward.o maintenance.o scoring.o settings.o state.o viterbi.o
TAR = libprofilehmm0.9.tar.bz2

$(TARGET): $(OBJS)
	ar -rcs $(TARGET) $(OBJS)

backward.o: backward.cpp profilehmm.hpp settings.hpp
	$(CC) $(CCOPTS) -o $@ $<

baum-welch.o: baum-welch.cpp profilehmm.hpp settings.hpp
	$(CC) $(CCOPTS) -o $@ $<

extendedexp.o: extendedexp.cpp extendedexp.hpp
	$(CC) $(CCOPTS) -o $@ $<

forward.o: forward.cpp profilehmm.hpp settings.hpp
	$(CC) $(CCOPTS) -o $@ $<

maintenance.o: maintenance.cpp profilehmm.hpp settings.hpp
	$(CC) $(CCOPTS) -o $@ $<

scoring.o: scoring.cpp profilehmm.hpp settings.hpp
	$(CC) $(CCOPTS) -o $@ $<

settings.o: settings.cpp settings.hpp
	$(CC) $(CCOPTS) -o $@ $<

state.o: state.cpp profilehmm.hpp settings.hpp
	$(CC) $(CCOPTS) -o $@ $<

viterbi.o: viterbi.cpp profilehmm.hpp settings.hpp
	$(CC) $(CCOPTS) -o $@ $<

tar: backward.cpp baum-welch.cpp extendedexp.cpp extendedexp.hpp forward.cpp gpl maintenance.cpp makefile profilehmm.hpp scoring.cpp settings.cpp settings.hpp state.cpp viterbi.cpp
	tar cfj $(TAR) *.hpp *.cpp gpl makefile

clean:
	rm -f $(OBJS) $(TARGET) $(TAR) core

tidy:
	rm -f $(OBJS) core

