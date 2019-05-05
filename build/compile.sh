g++ -Wall -I/usr/local/include -c ../src/*.cpp
g++ *.o -lgsl -lgslcblas -lm -o ../DVR
