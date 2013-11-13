#if your machine uses cc instead of gcc, edit the following:
export CC=gcc
export FC=g77
export MAKE=make

#set the base directory
export BASE_DIR = $(shell pwd)

#uncomment and set PGPLOT files if installed
#export PGHEAD = \"/usr/local/include/cpgplot.h\"
#export PGPLOT = -L/usr/local/lib -lcpgplot -lpgplot -lpng -L/usr/X11R6/lib -lX11 -L/usr/lib -lgcc

#compilation flags
export CFLAGS= -O3 -Wall

#uncomment if you will be using the WFPC2 module
#export USEWFPC2=1

#uncomment if you will be using the ACS module
#export USEACS=1

#uncomment if you will be using the WFC3 module
#export USEWFC3=1

#set the directory where Tiny Tim binaries are located (ACS/WFPC2/WFC3 PSF generation):
#export USETT=1
#export CFLAGS+= -DTTDIR='"/data0/tinytim-7.0"'

#set maximum number of images, and maximum number of stars
export CFLAGS+= -DMAXNIMG=100
export CFLAGS+= -DMAXNSTARS=2000000

#optionally print and/or crash if a NaN occurs (intended for development only)
export CFLAGS+= -DNAN_PRINT
#export CFLAGS+= -DNAN_CRASH

##### NO FURTHER ADJUSTMENTS BELOW THIS LINE ####

#add base directory settings to flags
export CFLAGS += -DBASEDIR='"$(BASE_DIR)"' -I$(BASE_DIR)/include

#definitions
ALLEXE = dolphot_lib.o fits_lib.o bin/addstars bin/apphot bin/calcsky bin/dolphot bin/synthimg bin/fakelist
ALL = $(ALLEXE) ccdred
DOLPHOT_DEP = dolphot_lib.o fits_lib.o dolphot_common.h dolphot_defs.h
DOLPHOT_OBJ = dolphot_lib.o fits_lib.o
DOLPHOT_FOBJ =

#linked libraries
export LIBS = -lm

#add WFPC2 modules, if selected
ifdef USEWFPC2
export CFLAGS+= -DUSEWFPC2
ALL+= wfpc2
DOLPHOT_DEP += wfpc2
DOLPHOT_OBJ += wfpc2/wfpc2phot.o wfpc2/wfpc2filters.o wfpc2/wfpc2distort.o
DOLPHOT_FOBJ += wfpc2/wfpc2filters.o
endif

#add ACS modules, if selected
ifdef USEACS
export CFLAGS+= -DUSEACS
ALL+= acs
DOLPHOT_DEP += acs
DOLPHOT_OBJ += acs/acsphot.o acs/acsfilters.o acs/acsdistort.o
DOLPHOT_FOBJ += acs/acsfilters.o
endif

#add WFC3 modules, if selected
ifdef USEWFC3
export CFLAGS+= -DUSEWFC3
ALL+= wfc3
DOLPHOT_DEP += wfc3
DOLPHOT_OBJ += wfc3/wfc3phot.o wfc3/wfc3filters.o wfc3/wfc3distort.o
DOLPHOT_FOBJ += wfc3/wfc3filters.o
endif

#add display module, if selected
ifdef PGPLOT
ALL+= display.o
ALLEXE+= display.o
endif

all: $(ALL)

ifdef USEWFPC2
wfpc2::
	cd wfpc2; $(MAKE)
endif

ifdef USEACS
acs::
	cd acs; $(MAKE)
endif

ifdef USEWFC3
wfc3::
	cd wfc3; $(MAKE)
endif

ccdred::
	cd ccdred; $(MAKE)

dolphot_lib.o: dolphot_lib.c include/dolphot.h
	$(CC) $(CFLAGS) -c -o $@ $<

fits_lib.o: fits_lib.c include/dolphot.h include/fits.h
	$(CC) $(CFLAGS) -c -o $@ $<

bin/addstars: addstars.c $(DOLPHOT_DEP)
	$(CC) $(CFLAGS) -o $@ addstars.c $(DOLPHOT_OBJ) $(LIBS)

ifdef PGPLOT
display.o: display.c dolphot_lib.o fits_lib.o
	$(CC) $(CFLAGS) -c -o $@ $< -DPGHEAD=$(PGHEAD)

bin/apphot: apphot.c dolphot_lib.o fits_lib.o display.o
	$(CC) $(CFLAGS) -c apphot.c -DPGPLOT -DPGHEAD=$(PGHEAD)
	$(FC) -o $@ apphot.o dolphot_lib.o fits_lib.o display.o $(LIBS) $(PGPLOT)
	rm apphot.o

bin/dolphot: dolphot.c $(DOLPHOT_DEP)
	$(CC) $(CFLAGS) -c dolphot.c -DPGPLOT -DPGHEAD=$(PGHEAD)
	$(FC) -o $@ dolphot.o $(DOLPHOT_OBJ) $(LIBS) $(PGPLOT)
	rm dolphot.o
#	$(CC) $(CFLAGS) -o $@ dolphot.c -DPGPLOT -DPGHEAD=$(PGHEAD) $(DOLPHOT_OBJ) $(LIBS) $(PGPLOT) -lcpgplot -lpgplot -lX11 -L/usr/lib -lgcc
else
bin/apphot: apphot.c dolphot_lib.o fits_lib.o
	$(CC) $(CFLAGS) -o $@ apphot.c dolphot_lib.o fits_lib.o $(LIBS)

bin/dolphot: dolphot.c $(DOLPHOT_DEP)
	$(CC) $(CFLAGS) -o $@ dolphot.c $(DOLPHOT_OBJ) $(LIBS)
endif

bin/calcsky: calcsky.c dolphot_lib.o fits_lib.o
	$(CC) $(CFLAGS) -o $@ calcsky.c dolphot_lib.o fits_lib.o $(LIBS)

bin/synthimg: synthimg.c dolphot_lib.o fits_lib.o
	$(CC) $(CFLAGS) -o $@ synthimg.c dolphot_lib.o fits_lib.o $(LIBS)

bin/fakelist: fakelist.c fakeproc.o fakeproc.h $(DOLPHOT_FOBJ)
	$(CC) $(CFLAGS) -o $@ fakelist.c fakeproc.o $(DOLPHOT_FOBJ) $(LIBS)

clean:
	rm $(ALLEXE)
	cd ccdred ; $(MAKE) clean
ifdef USEWFPC2
	cd wfpc2 ; $(MAKE) clean
endif
ifdef USEACS
	cd acs ; $(MAKE) clean
endif
ifdef USEWFC3
	cd wfc3 ; $(MAKE) clean
endif
