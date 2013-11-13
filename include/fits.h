#ifndef FITS_H
#define FITS_H
#include "dolphot.h"
typedef char cardtype[81];

typedef struct {
   cardtype name;
   int length;
   char type;
   size_t offset;
} bintabledatatype;

typedef struct {
   int Ncards,Nmax,X,Y,Z,bits,pcount,tfields;
   double bzero,bscale;
   char xtension[81];
   cardtype*cards;
   imgtype data;
   bintabledatatype *bintabdata;
} imtype;

typedef struct {
   int Next;
   imtype img;
   imtype *ext;
} ftype;

#ifdef FITS_C
char fitsstr[8][81]={"GAIN    =","EXPTIME =","RNOISE  =","MINVAL  =","MAXVAL  =","EPOCH   =","AIRMASS =","EXPTIME0="};
double fitsval[8];
#else
extern char fitsstr[8][81];
extern double fitsval[8];

extern void getparam(void);
extern void addcard(imtype*,char[81]);
//extern void readcards(FILE*,ftype*,imtype*);
extern FILE* readfitsh(char*,ftype*,int);
extern void readexth(FILE*,imtype*,int);
extern void readchip(FILE*,chiptype,imtype*);
extern void readimage(FILE*,imtype*);
extern void readbody(FILE*,ftype*,int);
extern void readfits(char*,ftype*,int);
extern void skipimage(FILE*,imtype*);
extern void readfitsinfo(char*,ftype*,int);
extern void parsecards(imtype*,double*,double*,double*,float*,float*,double*,double*,double*,int,int);
extern char* getcard(int,cardtype*,char*);
extern char* getcardval(imtype*,char*,int);
extern int gettablevalint(imtype*,char*,int,int,int);
extern float gettablevalfloat(imtype*,char*,int,int,int);
extern double gettablevaldouble(imtype*,char*,int,int,int);
extern char*gettablevalstring(imtype*,char*,int,int);
extern void insertcards(imtype*,double,double,double,float,float,double,double,double);
//extern void writecards(FILE*,ftype*,imtype*);
extern FILE* writefitsh(char*,ftype*,int);
extern void writeexth(FILE*,imtype*,int);
extern void writechip(FILE*,chiptype,imtype*);
extern void writeimage(FILE*,imtype*);
extern void writefits(char*,ftype*,int);
extern void freeim(imtype*);
extern void imcopy(imtype*,imtype*);
extern void fitscopy(ftype*,ftype*);
extern void freefits(ftype*);
extern void freeinfo(ftype*);
extern int isimage(imtype*);
#endif
#endif
