test.exe: *.cpp *.h
	g++ -O3 -o test.exe -Wall -Wextra -Wshadow -Wpedantic *.cpp -std=c++17

