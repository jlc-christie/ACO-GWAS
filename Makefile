CFLAGS = -O -Wall
FinalYearProject: main.o input.o stats.o
	g++ $(CFLAGS) -o FinalYearProject main.o input.o stats.o
main.o: main.cpp
	g++ $(CFLAGS) -c main.cpp
input.o: input.cpp
	g++ $(CFLAGS) -c input.cpp
stats.o: stats.cpp
	g++ $(CFLAGS) -c stats.cpp
optimized:
	g++ -O3 -Wall -o FinalYearProject main.o input.o stats.o
run:
	./FinalYearProject
clean:
	rm -f core *.o
