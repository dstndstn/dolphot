#define CCDPROC_C
#include "ccdproc.h"

void perr(char *str) {
   printf("%s\n",str);
   exit(1);
}

int ccdprocparam(char*var,char*val) {
   int i;
   char *ptr;

   if (!strcasecmp(var,"imager")) {strcpy(imager,val); return 1;}
   if (!strcasecmp(var,"typekw")) {strcpy(typekw,val); return 1;}
   if (!strcasecmp(var,"filtkw")) {strcpy(filtkw,val); return 1;}
   if (!strcasecmp(var,"oversec")) {strcpy(oversec,val); return 1;}
   if (!strcasecmp(var,"zerobase")) {strcpy(zerobase,val); return 1;}
   if (!strcasecmp(var,"darkbase")) {strcpy(darkbase,val); return 1;}
   if (!strcasecmp(var,"darktime")) {strcpy(darktime,val); return 1;}
   if (!strcasecmp(var,"flatbase")) {strcpy(flatbase,val); return 1;}   
   if (!strcasecmp(var,"trimsec")) {strcpy(trimsec,val); return 1;}
   if (!strcasecmp(var,"datasec")) {strcpy(datasec,val); return 1;}
   if (!strcasecmp(var,"ccdsize")) {strcpy(ccdsize,val); return 1;}
   if (!strcasecmp(var,"ccdsec")) {strcpy(ccdsec,val); return 1;}
   if (!strcasecmp(var,"detsec")) {strcpy(detsec,val); return 1;}
   i=strtol(val,&ptr,10);
   if (!*ptr) {
      if (!strcasecmp(var,"SPECIAL")) {SPECIAL=i; if (SPECIAL<0 || SPECIAL>1) perr("Mask image must be 0 or 1"); return 1;}
      if (!strcasecmp(var,"MASK")) {MASK=i; if (MASK<0 || MASK>1) perr("Mask image must be 0 or 1"); return 1;}
      if (!strcasecmp(var,"OVERSCAN")) {OVERSCAN=i; if (OVERSCAN<0 || OVERSCAN>2) perr("Overscan correction must be 0, 1, or 2"); return 1;}
      if (!strcasecmp(var,"ZERO")) {ZERO=i; if (ZERO<0 || ZERO>1) perr("Zero correction must be 0 or 1"); return 1;}
      if (!strcasecmp(var,"DARK")) {DARK=i; if (DARK<0 || DARK>1) perr("Dark correction must be 0 or 1"); return 1;}
      if (!strcasecmp(var,"FLAT")) {FLAT=i; if (FLAT<0 || FLAT>1) perr("Flat correction must be 0 or 1"); return 1;}
      if (!strcasecmp(var,"TRIM")) {TRIM=i; if (TRIM<0 || TRIM>1) perr("Trim image must be 0 or 1"); return 1;}
      if (!strcasecmp(var,"MERGEAMP")) {MERGEAMP=i; if (MERGEAMP<0 || MERGEAMP>1) perr("Merge amplifiers must be 0 or 1"); return 1;}
   }
   return 0;
}

/*
  imager - imager name (str)
  typekw - header keyword for image type (str)
  filtkw - header keyword for filter (str)
  SPECIAL - apply crosstalk correction (int, 0=n/1=y)
  MASK - use bad pixel mask (int, 0=n/1=y)
  OVERSCAN - make overscan correction (int, 0=n/1=const/2=each row)
  oversec - overscan section (str)
  ZERO - make zero correction (int, 0=n/1=y)
  zerobase - zero base name (str)
  DARK - make dark correction (int, 0=n/1=y)
  darkbase - dark base name (str)
  darktime - header keyword for dark time (str)
  FLAT - make flat correction (int, 0=n/1=y)
  flatbase - flat base name (str)
  TRIM - trim image (int, 0=n/1=y)
  trimsec - trim section (str)
  datasec - data section (str)
  MERGEAMP - merge amplifiers (int, 0=n/1=y)
  ccdsize - actual CCD size (str)
  ccdsec - section of CCD used by image (str)
  detsec - section of detecter used by image (str)
*/

int parsesec(imtype*img,char*str,int*x1,int*x2,int*y1,int*y2,int*z1,int*z2) {
   char*ptr;

   while (*str==32) str++;
   if (*(str++)!='[') return 1;
   while (*str==32) str++;
   if (*str=='*') {
      *x1=0;
      *x2=img->X;
      str++;
   }
   else {
      *x1=strtol(str,&ptr,10)-1;
      if (ptr==str) return 1;
      str=ptr;
      while (*str==32) str++;
      if (*(str++)!=':') return 1;
      *x2=strtol(str,&ptr,10);
      if (ptr==str) return 1;
      str=ptr;
   }
   while (*str==32) str++;
   if (*(str++)!=',') return 1;
   while (*str==32) str++;
   if (*str=='*') {
      *y1=0;
      *y2=img->Y;
      str++;
   }
   else {
      *y1=strtol(str,&ptr,10)-1;
      if (ptr==str) return 1;
      str=ptr;
      while (*str==32) str++;
      if (*(str++)!=':') return 1;
      *y2=strtol(str,&ptr,10);
      if (ptr==str) return 1;
      str=ptr;
   }
   while (*str==32) str++;
   if (*str==']') {
      *z1=0;
      *z2=img->Z;
      return 0;
   }
   if (*(str++)!=',') return 1;
   while (*str==32) str++;
   if (*str=='*') {
      *z1=0;
      *z2=img->Z;
      str++;
   }
   else {
      *z1=strtol(str,&ptr,10)-1;
      if (ptr==str) return 1;
      str=ptr;
      while (*str==32) str++;
      if (*(str++)!=':') return 1;
      *z2=strtol(str,&ptr,10);
      if (ptr==str) return 1;
      str=ptr;
   }
   while (*str==32) str++;
   if (*(str++)!=']') return 1;
   return 0;
}

int getsec(imtype*img,cardtype str,int*x1,int*x2,int*y1,int*y2,int*z1,int*z2) {
   if (!parsesec(img,str,x1,x2,y1,y2,z1,z2)) return 1;
   if (!parsesec(img,getcardval(img,str,0),x1,x2,y1,y2,z1,z2)) return 1;
   *x1=*y1=*z1=0;
   *x2=img->X;
   *y2=img->Y;
   *z2=img->Z;
   return 0;
}
