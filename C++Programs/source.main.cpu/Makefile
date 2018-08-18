CC=g++
TIMER=/usr/bin/time
OPT=
NOPT=-Ofast  $(OPT)

all: main.o wavesim.o arrayman.o  abcon.o rnw.o hider.o
	$(CC) -std=c++11 *.o -o ./wve_nacdgpu $(NOPT)

arrayman.o : arrayman.cc
	$(CC) -std=c++11 -c $(NOPT) arrayman.cc

rnw.o : rnw.cc
	$(CC) -std=c++11 -c $(NOPT) rnw.cc

abcon.o : abcon.cc
	$(CC) -std=c++11 -c $(NOPT) abcon.cc

hider.o	: hider.cc
	$(CC) -std=c++11 -c $(NOPT) hider.cc

wavesim.o : wavesim.cc
	$(CC) -std=c++11 -c $(NOPT) wavesim.cc

main.o : main.cc
	$(CC) -std=c++11 -c $(NOPT) main.cc

clean :
	rm *.o
