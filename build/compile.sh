#!/bin/sh
source ~/.bashrc
g++ ../src/*.cpp -Wall -I ${GSL_HOME}/include -c -g
g++ *.o -L ${GSL_HOME}/lib -lgsl -lgslcblas -lm -g -o ../DVR
