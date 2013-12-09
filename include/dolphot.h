#ifndef DOLPHOT_H
#define DOLPHOT_H
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

typedef float ***imgtype;
typedef float **chiptype;

typedef float ** const constchiptype;

#ifndef DOLPHOT_C
extern int ferr(FILE*,char*);
extern void merr(void);
extern chiptype allocchip(int,int);
extern chiptype allocchipchar(int,int);
extern void freechip(chiptype,int,int);
extern imgtype allocimg(int,int,int);
extern imgtype allocimgchar(int,int,int);
extern void freeimg(imgtype,int,int,int);
extern int parseparam(char*,int(*)(char*,char*));
extern void paramfile(char*,int(*)(char*,char*));
extern void paramfile1(char*,int(*)(char*,char*));

extern void endian2(char*ptr,int N);
extern void endian4(char*ptr,int N);
extern void endian8(char*ptr,int N);

extern size_t ffread(void*,size_t,size_t,FILE*);
extern size_t ffwrite(void*,size_t,size_t,FILE*);
extern size_t sfread(void*,size_t,size_t,FILE*);
extern size_t sfwrite(void*,size_t,size_t,FILE*);

extern char camera[81];
#endif

#endif
