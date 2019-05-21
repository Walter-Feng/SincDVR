g++ -Wall -I /usr/local/include -g -c ../src/main.cpp ../src/FunctionLib.cpp ../src/fileprint.cpp
g++ main.o FunctionLib.o fileprint.o -static -lgsl -lgslcblas -lm -g -o ../DVR
