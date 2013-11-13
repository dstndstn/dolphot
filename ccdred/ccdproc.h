#ifndef CCDPROC_H
#define CCDPROC_H

#include <fits.h>

#ifdef CCDPROC_C
int SPECIAL=1,MASK=1,OVERSCAN=2,ZERO=1,DARK=0,FLAT=1,TRIM=1,MERGEAMP=1;
char imager[81],oversec[81],zerobase[161],darkbase[161],flatbase[161],trimsec[81],datasec[81],ccdsize[81],ccdsec[81],class[81],detsec[81],typekw[161],filtkw[161],darktime[161];
#else
extern int SPECIAL,MASK,OVERSCAN,ZERO,DARK,FLAT,TRIM,MERGEAMP;
extern char imager[81],oversec[81],zerobase[161],darkbase[161],flatbase[161],trimsec[81],datasec[81],ccdsize[81],ccdsec[81],class[81],detsec[81],typekw[161],filtkw[161],darktime[161];

extern int ccdprocparam(char*,char*);
extern int parsesec(imtype*,char*,int*,int*,int*,int*,int*,int*);
extern int getsec(imtype*,cardtype,int*,int*,int*,int*,int*,int*);
#endif
#endif
