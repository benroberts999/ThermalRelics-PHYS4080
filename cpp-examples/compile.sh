#!/bin/bash

FLAGS='-O3 -lgsl -lblas -Wall -Wpedantic -Wextra -Wconversion -Wshadow'

#I need these on my mac (to include GSL). Not required on linux
MAC_FLAGS='-L/usr/local/opt/gnu-scientific-library/lib/ -I/usr/local/opt/gnu-scientific-library/include/'

g++ ode.cpp $FLAGS $MAC_FLAGS -o ode