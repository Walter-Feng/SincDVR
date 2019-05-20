g++ -Wall -I /usr/local/include -g -c ../src/*.cpp
g++ *.o -lgsl -lgslcblas -lm -g -o ../DVR
