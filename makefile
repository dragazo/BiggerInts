test.exe: *.cpp *.h
	g++ -O4 -o test.exe -Wall -Wextra -Wshadow -Wpedantic *.cpp -std=c++17

