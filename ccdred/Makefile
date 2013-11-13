ALL = ccdproc_lib.o ../bin/applymask ../bin/bintab2txt ../bin/calib ../bin/ccdproc ../bin/crrej ../bin/darkcombine ../bin/fillhead ../bin/flatcombine ../bin/fringecombine ../bin/imarith ../bin/imcombine ../bin/imhead ../bin/imstat ../bin/makergb ../bin/setbadpix ../bin/splitgroups ../bin/subfringe ../bin/zerocombine

all: $(ALL)

ccdproc_lib.o: ccdproc_lib.c ccdproc.h ../include/dolphot.h ../include/fits.h
	$(CC) $(CFLAGS) -c -o $@ $<

../bin/applymask: applymask.c ../dolphot_lib.o ../fits_lib.o
	$(CC) $(CFLAGS) -o $@ $< ../dolphot_lib.o ../fits_lib.o $(LIBS)

../bin/bintab2txt: bintab2txt.c ../dolphot_lib.o ../fits_lib.o
	$(CC) $(CFLAGS) -o $@ $< ../dolphot_lib.o ../fits_lib.o $(LIBS)

../bin/calib: calib.c
	$(CC) $(CFLAGS) -o $@ $< $(LIBS)

../bin/ccdproc: ccdproc.c ../dolphot_lib.o ../fits_lib.o ccdproc.h ccdproc_lib.o
	$(CC) $(CFLAGS) -o $@ $< ../dolphot_lib.o ../fits_lib.o ccdproc_lib.o $(LIBS)

../bin/crrej: crrej.c ../dolphot_lib.o ../fits_lib.o ccdproc.h ccdproc_lib.o
	$(CC) $(CFLAGS) -o $@ $< ../dolphot_lib.o ../fits_lib.o ccdproc_lib.o $(LIBS)

../bin/darkcombine: darkcombine.c ../dolphot_lib.o ../fits_lib.o ccdproc.h ccdproc_lib.o
	$(CC) $(CFLAGS) -o $@ $< ../dolphot_lib.o ../fits_lib.o ccdproc_lib.o $(LIBS)

../bin/fillhead: fillhead.c ../dolphot_lib.o ../fits_lib.o
	$(CC) $(CFLAGS) -o $@ $< ../dolphot_lib.o ../fits_lib.o $(LIBS)

../bin/flatcombine: flatcombine.c ../dolphot_lib.o ../fits_lib.o ccdproc.h ccdproc_lib.o
	$(CC) $(CFLAGS) -o $@ $< ../dolphot_lib.o ../fits_lib.o ccdproc_lib.o $(LIBS)

../bin/fringecombine: fringecombine.c ../dolphot_lib.o ../fits_lib.o ccdproc.h ccdproc_lib.o
	$(CC) $(CFLAGS) -o $@ $< ../dolphot_lib.o ../fits_lib.o ccdproc_lib.o $(LIBS)

../bin/imarith: imarith.c ../dolphot_lib.o ../fits_lib.o
	$(CC) $(CFLAGS) -o $@ $< ../dolphot_lib.o ../fits_lib.o $(LIBS)

../bin/imcombine: imcombine.c ../dolphot_lib.o ../fits_lib.o ccdproc.h ccdproc_lib.o
	$(CC) $(CFLAGS) -o $@ $< ../dolphot_lib.o ../fits_lib.o ccdproc_lib.o $(LIBS)

../bin/imhead: imhead.c ../dolphot_lib.o ../fits_lib.o
	$(CC) $(CFLAGS) -o $@ $< ../dolphot_lib.o ../fits_lib.o $(LIBS)

../bin/imstat: imstat.c ../dolphot_lib.o ../fits_lib.o ccdproc.h ccdproc_lib.o
	$(CC) $(CFLAGS) -o $@ $< ../dolphot_lib.o ../fits_lib.o ccdproc_lib.o $(LIBS)

../bin/makergb: makergb.c ../dolphot_lib.o ../fits_lib.o
	$(CC) $(CFLAGS) -o $@ $< ../dolphot_lib.o ../fits_lib.o $(LIBS)

../bin/setbadpix: setbadpix.c ../dolphot_lib.o ../fits_lib.o
	$(CC) $(CFLAGS) -o $@ $< ../dolphot_lib.o ../fits_lib.o $(LIBS)

../bin/splitgroups: splitgroups.c ../dolphot_lib.o ../fits_lib.o
	$(CC) $(CFLAGS) -o $@ $< ../dolphot_lib.o ../fits_lib.o $(LIBS)

../bin/subfringe: subfringe.c ../dolphot_lib.o ../fits_lib.o ccdproc.h ccdproc_lib.o
	$(CC) $(CFLAGS) -o $@ $< ../dolphot_lib.o ../fits_lib.o ccdproc_lib.o $(LIBS)

../bin/zerocombine: zerocombine.c ../dolphot_lib.o ../fits_lib.o ccdproc.h ccdproc_lib.o
	$(CC) $(CFLAGS) -o $@ $< ../dolphot_lib.o ../fits_lib.o ccdproc_lib.o $(LIBS)

clean:
	rm $(ALL)
