CC=pgc++
TIMER=/usr/bin/time
OPT= -ta=tesla:cc60,managed -Minfo=accel
NOPT=-fast -Minfo=opt  $(OPT)

all: main.o wavesim.o arrayman.o  abcon.o rnw.o hider.o
	$(CC) -acc -std=c++11 *.o -o ./wve_fdparout $(NOPT)

arrayman.o : arrayman.cc
	$(CC) -acc -std=c++11 -c $(NOPT) arrayman.cc

rnw.o : rnw.cc
	$(CC) -acc -std=c++11 -c $(NOPT) rnw.cc

abcon.o : abcon.cc
	$(CC) -acc -std=c++11 -c $(NOPT) abcon.cc

hider.o	: hider.cc
	$(CC) -acc -std=c++11 -c $(NOPT) hider.cc

wavesim.o : wavesim.cc
	$(CC) -acc -std=c++11 -c $(NOPT) wavesim.cc

main.o : main.cc
	$(CC) -acc  -std=c++11 -c $(NOPT) main.cc

clean :
	rm *.o
