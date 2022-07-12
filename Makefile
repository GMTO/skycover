SRCS = point.cpp probe.cpp star.cpp stargroup.cpp polygon.cpp collisions.cpp prod.cpp shadow.cpp

# O3 gave incorrect results on Linux
INCLUDE = -I/Users/bxin/simulations/boost_1_79_0/

CFLAGS = $(INCLUDE) -O3
#CFLAGS = $(INCLUDE) -g -gdwarf-2

CC   = g++

# CfA compiler location
#CC   = /opt/stow/gcc-6.2.0/bin/g++ -Wl,-rpath=/opt/stow/gcc-6.2.0/lib64 

all: probegeom skycov agwsvalid

# we now use agwsprobe.ipynb to create the txt files. do not overwrite that here.
#probegeom:
#	python agwsprobe.py 15 76

skycov: skycov.cpp $(SRCS)
	$(CC) $(CFLAGS) -std=c++14 $(SRCS) skycov.cpp -o skycov 

agwsvalid: agwsvalid.cpp $(SRCS)
	$(CC) $(CFLAGS) -std=c++14 $(SRCS) agwsvalid.cpp -o agwsvalid

clean:
	rm skycov agwsvalid

dgwfshadow:
	$(CC) $(CFLAGS) -I/Users/bxin/simulations/boost_1_79_0/ dgwfshadow.cpp

shadow.o:	shadow.cpp
	$(CC) $(CFLAGS) -c -I/Users/bxin/simulations/boost_1_79_0/ shadow.cpp
