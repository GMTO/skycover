SRCS = point.cpp probe.cpp star.cpp stargroup.cpp polygon.cpp collisions.cpp prod.cpp shadow.cpp

# O3 gave incorrect results on Linux
INCLUDE = -I/home/bmcleod/include

CFLAGS = $(INCLUDE) -O3
#CFLAGS = $(INCLUDE) -g -gdwarf-2

CC   = g++

# CfA compiler location
#CC   = /opt/stow/gcc-6.2.0/bin/g++ -Wl,-rpath=/opt/stow/gcc-6.2.0/lib64 

all: probegeom skycov agwsvalid

probegeom:
	python agwsprobe.py 15 76

skycov: skycov.cpp $(SRCS)
	$(CC) $(CFLAGS) -std=c++11 $(SRCS) skycov.cpp -o skycov 

agwsvalid: agwsvalid.cpp $(SRCS)
	$(CC) $(CFLAGS) -std=c++11 $(SRCS) agwsvalid.cpp -o agwsvalid

clean:
	rm skycov agwsvalid

dgwfshadow:
	$(CC) $(CFLAGS) -I/home/bmcleod/include dgwfshadow.cpp

shadow.o:	shadow.cpp
	$(CC) $(CFLAGS) -c -I/home/bmcleod/include shadow.cpp
