test.exe: test.o BiggerInts.o
	g++ -O3 -o test.exe -Wall -Wextra -Wshadow -Wpedantic *.o -std=c++17

clean:
	rm -f *.o
	rm -f test.exe

BiggerInts.o: BiggerInts.cpp BiggerInts.h
	g++ -O3 -c -Wall -Wextra -Wshadow -Wpedantic -std=c++17 BiggerInts.cpp
test.o: test.cpp BiggerInts.h
	g++ -O3 -c -Wall -Wextra -Wshadow -Wpedantic -std=c++17 test.cpp