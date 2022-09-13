#!/bin/bash

# Always enable all the warnings (these are the bare minimum)
WARNS='-Wall -Wpedantic -Wextra -Wconversion -Wshadow'

# We need to link against GSL libraries:
LIBS='-lgsl -lblas'

#I need these on my mac (to include GSL). Not required on linux
# (Mac seems to install these in non-standard location, which is why we need to
# specifically tell the compiler where to look for these libraries)
MAC_FLAGS='-L/usr/local/opt/gnu-scientific-library/lib/ -I/usr/local/opt/gnu-scientific-library/include/'

g++ -O3 ode.cpp $LIBS $WARNS $MAC_FLAGS -o ode &&
g++ -O3 integrate.cpp $LIBS $WARNS $MAC_FLAGS -o integrate