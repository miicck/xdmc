if [ $1 ]; then
	g++ -Wall main.cpp
else
	g++ -O3 -Ofast main.cpp
fi
