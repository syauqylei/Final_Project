CC=g++
NOPT=-Ofast

all: main.o wavesim.o arrayman.o  abcon.o rnw.o 
	$(CC) -std=c++11 -pg *.o -o ./wve_nacdgpu $(NOPT)

arrayman.o : arrayman.cc
	$(CC) -std=c++11 -c -pg  $(NOPT) arrayman.cc

rnw.o : rnw.cc
	$(CC) -std=c++11 -c -pg  $(NOPT) rnw.cc

abcon.o : abcon.cc
	$(CC) -std=c++11 -c -pg  $(NOPT) abcon.cc

wavesim.o : wavesim.cc
	$(CC) -std=c++11 -c -pg  $(NOPT)  wavesim.cc
	
main.o : main.cc
	$(CC) -std=c++11 -c -pg  $(NOPT) main.cc

clean :
	rm *.o
