if [ $1 ]; then
	g++ -std=c++11 -g -Wall main.cpp
else
	g++ -std=c++11 -O3 -Ofast main.cpp
fi
