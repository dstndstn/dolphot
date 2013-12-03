ALL = ../bin/wfc3IRandersonpsf ../bin/wfc3UVISandersonpsf ../bin/wfc3distort ../bin/wfc3fitdistort ../bin/wfc3mask wfc3phot.o ../bin/wfc3showpsf ../bin/wfc3fakelist wfc3filters.o wfc3distort.o ../bin/wfc3undistort ../bin/wfc3unmask

ifdef USETT
ALL += ../bin/wfc3ttpar ../bin/wfc3makepsf
endif

all: $(ALL)

wfc3phot.o: wfc3phot.c wfc3psfdata.h ../include/fits.h wfc3distort.h wfc3filters.h ../dolphot_defs.h
	$(CC) $(CFLAGS) -c -o $@ $<

../bin/wfc3IRandersonpsf: wfc3IRandersonpsf.c ../dolphot_lib.o ../fits_lib.o wfc3psfdata.h ../include/fits.h ../include/dolphot.h
	$(CC) $(CFLAGS) -o $@ $< ../dolphot_lib.o ../fits_lib.o

../bin/wfc3UVISandersonpsf: wfc3UVISandersonpsf.c ../dolphot_lib.o ../fits_lib.o wfc3psfdata.h ../include/fits.h ../include/dolphot.h
	$(CC) $(CFLAGS) -o $@ $< ../dolphot_lib.o ../fits_lib.o

../bin/wfc3distort: wfc3distort_main.c wfc3distort.o wfc3distort.h wfc3filters.o wfc3filters.h ../include/dolphot.h
	$(CC) $(CFLAGS) -o $@ $< ../dolphot_lib.o wfc3distort.o wfc3filters.o $(LIBS)

../bin/wfc3fitdistort: wfc3fitdistort.c wfc3distort.o wfc3distort.h wfc3filters.o wfc3filters.h ../include/dolphot.h
	$(CC) $(CFLAGS) -o $@ $< ../dolphot_lib.o wfc3distort.o wfc3filters.o $(LIBS)

../bin/wfc3undistort: wfc3undistort.c wfc3distort.o wfc3distort.h wfc3filters.o wfc3filters.h ../include/dolphot.h
	$(CC) $(CFLAGS) -o $@ $< ../dolphot_lib.o wfc3distort.o wfc3filters.o $(LIBS)

../bin/wfc3mask: wfc3mask.c ../dolphot_lib.o ../fits_lib.o ../include/fits.h ../include/dolphot.h wfc3psfdata.h
	$(CC) $(CFLAGS) -o $@ $< ../dolphot_lib.o ../fits_lib.o $(LIBS)

ifdef USETT
../bin/wfc3ttpar: wfc3ttpar.c
	$(CC) $(CFLAGS) -o $@ $<

../bin/wfc3makepsf: wfc3makepsf.c ../dolphot_lib.o ../fits_lib.o wfc3psfdata.h ../include/fits.h ../include/dolphot.h
	$(CC) $(CFLAGS) -o $@ $< ../dolphot_lib.o ../fits_lib.o
endif

../bin/wfc3showpsf: wfc3showpsf.c ../dolphot_lib.o wfc3psfdata.h ../include/dolphot.h
	$(CC) $(CFLAGS) -o $@ $< ../dolphot_lib.o $(LIBS)

wfc3filters.o: wfc3filters.c
	$(CC) $(CFLAGS) -c -o $@ $<

wfc3distort.o: wfc3distort.c
	$(CC) $(CFLAGS) -c -o $@ $<

../bin/wfc3fakelist: wfc3fakelist.c wfc3filters.h wfc3filters.o
	$(CC) $(CFLAGS) -o $@ $< wfc3filters.o $(LIBS)

../bin/wfc3unmask: wfc3unmask.c ../dolphot_lib.o ../fits_lib.o ../include/fits.h ../include/dolphot.h wfc3psfdata.h
	$(CC) $(CFLAGS) -o $@ $< ../dolphot_lib.o ../fits_lib.o $(LIBS)

clean:
	rm $(ALL)

