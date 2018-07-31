CC=pgc++
TIMER=/usr/bin/time
OPT= -Minfo=accel
NOPT=-fast -Minfo=opt  $(OPT)

all: main.o rnw.o arrayman.o
	$(CC) -std=c++11 -pg *.o -o ./mute $(NOPT)

rnw.o : rnw.cc
	$(CC) -std=c++11 -c -pg -acc $(NOPT) rnw.cc

arrayman.o : arrayman.cc
	$(CC) -std=c++11 -c -pg -acc $(NOPT) arrayman.cc

main.o : main.cc
	$(CC) -std=c++11 -c -pg -acc $(NOPT) main.cc

clean :
	rm *.o
