EXE = wfg
OBJ = avl.o 

OPT = -O3 
DEBUG = -g -Werror -pedantic

march_error = $(error please define an architecture, e.g., 'make march=pentium')
ifdef ARCH
	OPT=$(ARCH) $(COPT)
else
	ifdef march
		OPT += -march=$(march)
	endif
endif


CC = gcc
CFLAGS = -std=c99 -Wall

$(EXE): $(OBJ) wfg.c wfg.h wfgmain.c architecture
	$(CC) $(CFLAGS) $(OPT) -Dopt=0 -o wfg0 wfg.c wfgmain.c $(OBJ)
	$(CC) $(CFLAGS) $(OPT) -Dopt=1 -o wfg1 wfg.c wfgmain.c $(OBJ)
	$(CC) $(CFLAGS) $(OPT) -Dopt=2 -o wfg2 wfg.c wfgmain.c $(OBJ)
	$(CC) $(CFLAGS) $(OPT) -Dopt=3 -o wfg3 wfg.c wfgmain.c $(OBJ)

architecture:
	$(if $(march), , $(march_error))

wfg2.so: wfg.c wfg.h
	$(CC) $(CFLAGS) -Dopt=2 -fPIC -shared wfg.c -o wfg2.so

debug: wfg0debug wfg1debug wfg2debug wfg3debug

wfg0debug: wfg.c wfg.h wfgmain.c avl.c
	$(CC) $(CFLAGS) $(DEBUG)  -Dopt=0 -o wfg0debug wfg.c wfgmain.c avl.c
wfg1debug: wfg.c wfg.h wfgmain.c avl.c
	$(CC) $(CFLAGS) $(DEBUG)  -Dopt=1 -o wfg1debug wfg.c wfgmain.c avl.c
wfg2debug: wfg.c wfg.h wfgmain.c avl.c
	$(CC) $(CFLAGS) $(DEBUG)  -Dopt=2 -o wfg2debug wfg.c wfgmain.c avl.c
wfg3debug: wfg.c wfg.h wfgmain.c avl.c
	$(CC) $(CFLAGS) $(DEBUG)  -Dopt=3 -o wfg3debug wfg.c wfgmain.c avl.c

valgrind: debug
	valgrind -v --leak-check=full --show-reachable=yes ./wfg2debug ../random/3d/ran.50pts.3d.20 

%.o: %.c
	$(CC) $(CFLAGS) -c $<

clean: 
	rm -f wfg[0-4] exc[0-4] *.o 
	rm -rf *.dSYM

VERSION=$(shell git tag)

release: clean
	rm -rf '/tmp/WFG_'$(VERSION)
	mkdir -p '/tmp/WFG_'$(VERSION)
	cp -rf * '/tmp/WFG_'$(VERSION)
	cd /tmp/; tar czvf $(PWD)/../'WFG_'$(VERSION).tar.gz 'WFG_'$(VERSION)
