#include <fits.h>
#include "wfc3psfdata.h"

ftype fits,pam;
double RN,EXP,EXP0,EPOCH,FIXED_ET=-1.;
float DMIN,DMAX;
int MASKCR=1,NCOMBINE=1,FIXED_NC=-1,USE_WHT=0,USE_FLAT_FLAG=0;
int offsetx,offsety;
int Next0;
int FORCE_SUBARRAY=0,FORCE_OFFSETX=0,FORCE_OFFSETY=0; // 3=IR, 5=UVIS1, 6=UVIS2

char* eitherstring(char*val,int err) {
   // astrodrizzle
   if (Next0==4) {
      return gettablevalstring(fits.ext+3,val,0,err);
   }
   // FLT or multidrizzle
   return getcardval(&(fits.img),val,err);
}

double eitherdouble(char*val,int err) {
   // astrodrizzle
   if (Next0==4) {
      return gettablevaldouble(fits.ext+3,val,0,0,err);
   }
   // FLT or multidrizzle
   return atof(getcardval(&(fits.img),val,err));
}

int WFC3type(ftype *f) {
   int i;
   char detector[81],aperture[81];

   if (strcmp(getcardval(&(fits.img),"FILETYPE",1),"SCI") || strcmp(getcardval(&(fits.img),"TELESCOP",1),"HST") || strcmp(eitherstring("INSTRUME",1),"WFC3")) {printf("**Format error (filetype,telescop,instrume)\n"); exit(-1);}

   strcpy(detector,getcardval(&(fits.img),"DETECTOR",1));
   strcpy(aperture,getcardval(&(fits.img),"APERTURE",1));

   if (FORCE_SUBARRAY!=0) {
      int XMAX=4096,YMAX=2051;
      if (FORCE_SUBARRAY>3) {
	 if (strcmp(detector,"UVIS")) printf("Warning: detector is %s; UVIS was commanded\n",detector);
	 if (Next0!=3) printf("Warning: expect three-extension image for UVIS subarray\n");
      }
      else {
	 XMAX=YMAX=1014;
	 if (strcmp(detector,"IR")) printf("Warning: detector is %s; IR was commanded\n",detector);
	 if (Next0!=5) printf("Warning: expect five-extension image for IR subarray\n");
      }
      if (Next0<3) {
	 printf("Error: at least three extensions required\n");
	 exit(-1);
      }
      else if (f->ext[0].Z!=1) {
	 printf("Error: primary image has third dimension\n");
	 exit(-1);
      }
      else if (FORCE_OFFSETX<0 || FORCE_OFFSETY<0 || f->ext[0].X+FORCE_OFFSETX>XMAX || f->ext[0].Y+FORCE_OFFSETY>YMAX) {
	 printf("Error: subarray would fall outside %dx%d image\n",XMAX,YMAX);
	 exit(-1);
      }
      else if (f->ext[0].X!=f->ext[1].X || f->ext[0].Y!=f->ext[1].Y || f->ext[1].Z!=1 || f->ext[0].X!=f->ext[2].X || f->ext[0].Y!=f->ext[2].Y || f->ext[2].Z!=1) {
	 printf("Error: data quality image does not match primary image size\n");
	 exit(-1);
      }

      offsetx = FORCE_OFFSETX;
      offsety = FORCE_OFFSETY;
      return FORCE_SUBARRAY;
   }

   offsetx = 0;
   offsety = 0;
   // UVIS image types
   if (!strcmp(detector,"UVIS")) {
      //undrizzled UVIS
      if (Next0==6) {
	 for (i=0;i<6 && f->ext[i].X==4096 && f->ext[i].Y==2051 && f->ext[i].Z==1;i++);
	 if (i==6) return 1;
	 return 0;
      }
      if (Next0==3) {
	 if (!strcmp(aperture,"UVIS1-2K2A-SUB")) {
	    for (i=0;i<3 && f->ext[i].X==2047 && f->ext[i].Y==2050 && f->ext[i].Z==1;i++);
	    if (i==3) return 5;
	 }
	 if (!strcmp(aperture,"UVIS1-2K2B-SUB")) {
	    for (i=0;i<3 && f->ext[i].X==2047 && f->ext[i].Y==2050 && f->ext[i].Z==1;i++);
	    if (i==3) {offsetx=2049; return 5;}
	 }
	 if (!strcmp(aperture,"UVIS2-2K2C-SUB")) {
	    for (i=0;i<3 && f->ext[i].X==2047 && f->ext[i].Y==2050 && f->ext[i].Z==1;i++);
	    if (i==3) {offsety=1; return 6;}
	 }
	 if (!strcmp(aperture,"UVIS2-2K2D-SUB")) {
	    for (i=0;i<3 && f->ext[i].X==2047 && f->ext[i].Y==2050 && f->ext[i].Z==1;i++);
	    if (i==3) {offsetx=2049; offsety=1; return 6;}
	 }
	 if (!strcmp(aperture,"UVIS2-M1K1C-SUB")) {
	    for (i=0;i<3 && f->ext[i].X==1024 && f->ext[i].Y==1024 && f->ext[i].Z==1;i++);
	    if (i==3) {offsetx=1024; offsety=1027; return 6;}
	 }
	 if (!strcmp(aperture,"UVIS2-C1K1C-SUB")) {
	    for (i=0;i<3 && f->ext[i].X==1025 && f->ext[i].Y==1024 && f->ext[i].Z==1;i++);
	    if (i==3) {offsety=1; return 6;}
	 }
	 if (!strcmp(aperture,"UVIS1-C512A-SUB")) {
	    for (i=0;i<3 && f->ext[i].X==513 && f->ext[i].Y==512 && f->ext[i].Z==1;i++);
	    if (i==3) {offsety=1538; return 5;}
	 }
	 if (!strcmp(aperture,"UVIS2-M512C-SUB")) {
	    for (i=0;i<3 && f->ext[i].X==512 && f->ext[i].Y==512 && f->ext[i].Z==1;i++);
	    if (i==3) {offsetx=1536; offsety=1539; return 6;}
	 }
	 if (!strcmp(aperture,"UVIS2-C512C-SUB")) {
	    for (i=0;i<3 && f->ext[i].X==513 && f->ext[i].Y==512 && f->ext[i].Z==1;i++);
	    if (i==3) {offsety=1; return 6;}
	 }
      }
      //drizzled UVIS
      if (Next0==3 || Next0==4) {
	 if (f->ext[0].Z==1 && f->ext[1].X==f->ext[0].X && f->ext[1].Y==f->ext[0].Y && f->ext[1].Z==f->ext[0].Z && (USE_WHT!=0 || (f->ext[2].X==f->ext[0].X && f->ext[2].Y==f->ext[0].Y))) {
	    printf("Irregular size; assuming drizzled\n");
	    return 2;
	 }
      }
   }

   // IR image types
   if (!strcmp(detector,"IR")) {
      //undrizzled IR
      if (Next0==5) {
	 for (i=0;i<5 && f->ext[i].X==1014 && f->ext[i].Y==1014 && f->ext[i].Z==1;i++);
	 if (i==5) return 3;
	 if (!strcmp(aperture,"IRSUB64") || !strcmp(aperture,"IRSUB64-FIX")) {
	    for (i=0;i<5 && (i>=3 || (f->ext[i].X==64 && f->ext[i].Y==64 && f->ext[i].Z==1));i++);
	    if (i==5) {offsetx=offsety=480; return 3;}
	 }
	 if (!strcmp(aperture,"IRSUB128") || !strcmp(aperture,"IRSUB128-FIX")) {
	    for (i=0;i<5 && (i>=3 || (f->ext[i].X==128 && f->ext[i].Y==128 && f->ext[i].Z==1));i++);
	    if (i==5) {offsetx=offsety=448; return 3;}
	 }
	 if (!strcmp(aperture,"IRSUB256") || !strcmp(aperture,"IRSUB256-FIX")) {
	    for (i=0;i<5 && (i>=3 || (f->ext[i].X==256 && f->ext[i].Y==256 && f->ext[i].Z==1));i++);
	    if (i==5) {offsetx=offsety=384; return 3;}
	 }
	 if (!strcmp(aperture,"IRSUB512") || !strcmp(aperture,"IRSUB512-FIX")) {
	    for (i=0;i<5 && (i>=3 || (f->ext[i].X==512 && f->ext[i].Y==512 && f->ext[i].Z==1));i++);
	    if (i==5) {offsetx=offsety=256; return 3;}
	 }
      }

      // drizzled
      if (Next0==3 || Next0==4) {
	 if (f->ext[0].Z==1 && f->ext[1].X==f->ext[0].X && f->ext[1].Y==f->ext[0].Y && f->ext[1].Z==f->ext[0].Z && (USE_WHT!=0 || (f->ext[2].X==f->ext[0].X && f->ext[2].Y==f->ext[0].Y))) {
	    printf("Irregular size; assuming drizzled\n");
	    return 4;
	 }
      }
   }
   return 0;
}

void WFC3getcards(int ext,int MODE,int drz) {
   int x,y;
   char str[81];

   if (!strcmp(getcardval(&(fits.img),"CCDAMP",1),"ABCD")) {
      if (MODE==1) {
	 RN=0.5*(eitherdouble("READNSEC",1)+eitherdouble("READNSED",1));
      }
      else if (MODE==2) {
	 RN=0.5*(eitherdouble("READNSEA",1)+eitherdouble("READNSEB",1));
      }
      else RN=0.25*(eitherdouble("READNSEA",1)+eitherdouble("READNSEB",1)+eitherdouble("READNSEC",1)+eitherdouble("READNSED",1));
   }
   else if (!strcmp(getcardval(&(fits.img),"CCDAMP",1),"A")) RN=eitherdouble("READNSEA",1);
   else if (!strcmp(getcardval(&(fits.img),"CCDAMP",1),"B")) RN=eitherdouble("READNSEB",1);
   else if (!strcmp(getcardval(&(fits.img),"CCDAMP",1),"C")) RN=eitherdouble("READNSEC",1);
   else if (!strcmp(getcardval(&(fits.img),"CCDAMP",1),"D")) RN=eitherdouble("READNSED",1);
   else {printf("**Format error (ccdamp)\n"); exit(-1);}

   if (FIXED_ET>0) EXP=FIXED_ET;
   else EXP=atof(getcardval(&(fits.img),"EXPTIME",1));
   if (drz) {
      DMIN=DMAX=fits.ext[ext].data[0][0][0];
      for (y=0;y<fits.ext[ext].Y;y++) for (x=0;x<fits.ext[ext].X;x++) {
	 if (fits.ext[ext].data[0][y][x]<DMIN) DMIN=fits.ext[ext].data[0][y][x];
	 else if (fits.ext[ext].data[0][y][x]>DMAX) DMAX=fits.ext[ext].data[0][y][x];
      }
   }
   else {
      DMIN=(float)atof(getcardval(fits.ext+ext,"GOODMIN",1))-1;
      DMAX=(float)atof(getcardval(fits.ext+ext,"GOODMAX",1))+1;
   }
   EPOCH=0.5*(atof(getcardval(&(fits.img),"EXPSTART",1))+atof(getcardval(&(fits.img),"EXPEND",1)));
   if (FIXED_NC>0) {
      NCOMBINE=FIXED_NC;
      printf("Setting number of images to %d\n",NCOMBINE);
   }
   else {
      strcpy(str,getcardval(&(fits.img),"NCOMBINE",0));
      if (str[0]) NCOMBINE=atoi(getcardval(&(fits.img),"NCOMBINE",0));
      if (!str[0]) printf("NCOMBINE keyword not found; assuming data came from single image\n");
      else if (NCOMBINE<1) {
	 printf("Illegal NCOMBINE setting; assuming data came from single image\n");
	 NCOMBINE=1;
      }
      else printf("NCOMBINE = %d\n",NCOMBINE);
   }
   EXP0=EXP/NCOMBINE;
   return;
}

// ignoring type 64 (warm pixel) on assumption it is fixed OK
// ignoring type 512 in IR by default on assumption those pixels are fine
// ignoring type 16384 (ghost/crosstalk)
void WFC3_DQ_mask(ftype*f,int ext,int drz,int ir) {
   int x,y,dq,i;
   for (y=0;y<f->ext[ext].Y;y++) for (x=0;x<f->ext[ext].X;x++) if (!drz) {
      dq=(int)(f->ext[ext+2].data[0][y][x]+0.5);
      if (ir) {
	 // 1024 unused at present
	 if (dq&256) f->ext[ext].data[0][y][x]=DMAX+1; // 256
	 else if ((dq&512) && USE_FLAT_FLAG) f->ext[ext].data[0][y][x]=DMIN-1;
	 else if (dq&6335) f->ext[ext].data[0][y][x]=DMIN-1; // 4096+2048+128+63
	 else if ((dq&8192) && MASKCR) f->ext[ext].data[0][y][x]=DMIN-1;
      }
      else {
	 if (dq&2304) f->ext[ext].data[0][y][x]=DMAX+1; // 2048+256
	 else if (dq&5823) f->ext[ext].data[0][y][x]=DMIN-1; // 4096+1024+512+128+63
	 else if ((dq&8192) && MASKCR) f->ext[ext].data[0][y][x]=DMIN-1; // 8192
      }
   }
   else {
      dq=0;
      if (USE_WHT) {
	 if (f->ext[ext+1].data[0][y][x]>0.0) dq=1;
      }
      else {
	 for (i=0;i<f->ext[ext].Z && !dq;i++) if ((int)(f->ext[ext+2].data[i][y][x]+0.5)!=0) dq=1;
      }
      if (dq==0) f->ext[ext].data[0][y][x]=DMIN-1;
   }
   return;
}

void WFC3_PAM_mult(int ext,char*pamfn,double mult) {
   int x,y;
   float DMIN1,DMAX1;
   char fn[321];

   if (DMAX>0.) DMAX1=1.5*DMAX;
   else DMAX1=0.;
   if (DMIN<0.) DMIN1=1.5*DMIN;
   else DMIN1=0.;
   sprintf(fn,"%s/wfc3/data/%s",BASEDIR,pamfn);
   readfits(fn,&pam,0);
   if (pam.Next!=1 || pam.ext[0].X<fits.ext[ext].X+offsetx || pam.ext[0].Y<fits.ext[ext].Y+offsety || pam.ext[0].Z!=fits.ext[ext].Z) {
      printf("Error in PAM file %s\n",pamfn);
      exit(-1);
      return;
   }
   for (y=0;y<fits.ext[ext].Y;y++) for (x=0;x<fits.ext[ext].X;x++) if (fits.ext[ext].data[0][y][x]>DMIN && fits.ext[ext].data[0][y][x]<DMAX) {
      fits.ext[ext].data[0][y][x]*=pam.ext[0].data[0][y+offsety][x+offsetx]/mult;
      fits.ext[ext+1].data[0][y][x]*=pam.ext[0].data[0][y+offsety][x+offsetx]/mult;
   }
   else if (fits.ext[ext].data[0][y][x]>=DMAX) fits.ext[ext].data[0][y][x]=DMAX1+1;
   else fits.ext[ext].data[0][y][x]=DMIN1-1;
   freefits(&pam);
   DMIN=DMIN1;
   DMAX=DMAX1;
   return;
}

void WFC3_exptime_mult(int ext,double mult) {
   int x,y;

   DMAX*=EXP/mult;
   DMIN*=EXP/mult;
   for (y=0;y<fits.ext[ext].Y;y++) for (x=0;x<fits.ext[ext].X;x++) {
      fits.ext[ext].data[0][y][x]*=EXP/mult;
      fits.ext[ext+1].data[0][y][x]*=EXP/mult;
   }
   return;
}

/*
  auto-solve for noise parameters:
  y = noise^2 = rdnoise^2/gain^2 + data/gain = a + b*data
  chi^2 = SUM[(a+b*x-y)^2/sig^2] = SUM[(a/y+b*x/y-1)^2]
  dchi^2/da = 0 = SUM[(a/y+b*x/y-1)/y]
  dchi^2/db = 0 = SUM[(a/y+b*x/y-1)*x/y]
  I=SUM(1/yy); X=SUM(x/yy); Y=SUM(1/y); XX=SUM(xx/yy); XY=SUM(x/y);
  0 = a*I+b*X-Y   = a*I*X + b*X*X - X*Y   = a*XX*I + b*X*XX - XX*Y
  0 = a*X+b*XX-XY = a*I*X + b*I*XX - I*XY = a*X*X + b*X*XX - X*XY
  b = (I*XY-X*Y)/(I*XX-X*X) = 1/gain
  a = (XX*Y-X*XY)/(I*XX-X*X) = rdnoise^2/gain^2
  rdnoise = gain*sqrt((XX*Y-X*XY)/(I*XX-X*X))

  assuming gain=1 (should be the case):
  rdnoise^2 is simply the average of noise^2-data
  0 = a*I+X-Y  -->  a=(Y-X)/I
*/
void WFC3setcards(int ext,int ir) {
   int i;
   /*
   if (ir)
   {
      int x,y;
      chiptype data,noise;
      double I=0,X=0,Y=0,XX=0,XY=0,tGN,tRN;

      data=fits.ext[ext].data[0];
      noise=fits.ext[ext+1].data[0];
      for (y=0;y<fits.ext[ext].Y;y++) for (x=0;x<fits.ext[ext].X;x++) if (data[y][x]>DMIN && data[y][x]<DMAX) {
	 int dq=(int)(fits.ext[ext+2].data[0][y][x]+0.5);
	 if ((dq&256)==0 && (dq&8192)==0) {
	    double isig2;
	    isig2=1./(noise[y][x]*noise[y][x]*noise[y][x]*noise[y][x]);
	    printf("%f %f ^^^^\n",noise[y][x],data[y][x]);
	    I+=isig2;
	    Y+=noise[y][x]*noise[y][x]*isig2;
	    if (data[y][x]>0.) {
	       X+=data[y][x]*isig2;
	       XX+=data[y][x]*data[y][x]*isig2;
	       XY+=data[y][x]*noise[y][x]*noise[y][x]*isig2;
	    }
	 }
      }
      tGN=(XX*I-X*X)/(XY*I-X*Y);
      tRN=sqrt((XX*Y-X*XY)/(XX*I-X*X))*tGN;
      printf("%f %f\n",tGN,tRN);
   }
   */
   insertcards(fits.ext+ext,1.,RN*sqrt(NCOMBINE),EXP,DMIN,DMAX,EPOCH,0.0,EXP0);
   freeim(fits.ext+ext+1);
   freeim(fits.ext+ext+2);
   if (ir && Next0==5) {
      freeim(fits.ext+ext+3);
      freeim(fits.ext+ext+4);
      fits.Next-=4;
      for (i=ext+1;i<fits.Next;i++) memcpy(fits.ext+i,fits.ext+i+4,sizeof(imtype));
   }
   else {
      fits.Next-=2;
      for (i=ext+1;i<fits.Next;i++) memcpy(fits.ext+i,fits.ext+i+2,sizeof(imtype));
   }
   return;
}

int main(int argc,char**argv) {
   int i,tp;
   char card[81],*ptr;

   if (argc<2) {
      printf("Usage: %s <<-flags>> <fits files>\n",*argv);
      printf(" -keepcr      leaves fixed cosmic ray pixels in image\n");
      printf(" -exptime=#   overrides exposure time; needed for _drz images\n");
      printf(" -ncombine=#  overrides NCOMBINE keyword\n");
      printf(" -maskflat    to use Bad or uncertain flat value mask for WFC3/IR\n");
      printf(" -usewht      to use weight extension instead of context for _drz\n");
      printf(" -uvis1sub=#,# to force UVIS1 subarray shifted by #,# from full image\n");
      printf(" -uvis2sub=#,# to force UVIS2 subarray shifted by #,# from full image\n");
      printf(" -irsub=#,#   to force IR subarray shifted by #,# from full image\n");
      return 1;
   }
   for (i=1;i<argc;i++) if (!strcasecmp(argv[i],"-KEEPCR")) MASKCR=0;
   else if (!strncasecmp(argv[i],"-EXPTIME=",9)) FIXED_ET=atof(argv[i]+9);
   else if (!strncasecmp(argv[i],"-NCOMBINE=",10)) FIXED_NC=atof(argv[i]+10);
   else if (!strcasecmp(argv[i],"-MASKFLAT")) USE_FLAT_FLAG=1;
   else if (!strcasecmp(argv[i],"-USEWHT")) USE_WHT=1;
   else if (!strncasecmp(argv[i],"-irsub=",7)) {
      FORCE_SUBARRAY = 3;
      FORCE_OFFSETX = strtol(argv[i]+7,&ptr,10);
      if (*ptr!=',') {
	 printf("Cannot parse subarray definition \"%s\"\n",argv[i]);
	 return 1;
      }
      FORCE_OFFSETY = atoi(ptr+1);
   }
   else if (!strncasecmp(argv[i],"-uvis1sub=",10)) {
      FORCE_SUBARRAY = 5;
      FORCE_OFFSETX = strtol(argv[i]+10,&ptr,10);
      if (*ptr!=',') {
	 printf("Cannot parse subarray definition \"%s\"\n",argv[i]);
	 return 1;
      }
      FORCE_OFFSETY = atoi(ptr+1);
   }
   else if (!strncasecmp(argv[i],"-uvis2sub=",10)) {
      FORCE_SUBARRAY = 6;
      FORCE_OFFSETX = strtol(argv[i]+10,&ptr,10);
      if (*ptr!=',') {
	 printf("Cannot parse subarray definition \"%s\"\n",argv[i]);
	 return 1;
      }
      FORCE_OFFSETY = atoi(ptr+1);
   }
   else {
      readfits(argv[i],&fits,1);
      Next0=fits.Next;
      tp=WFC3type(&fits);
      if (Next0>0 && strcmp(getcardval(fits.ext,"DOL_WFC3",0),"")) printf("%s already run through wfc3mask\n",argv[i]);
      else if (tp<1 || tp>6) printf("%s is not a WFC3 fits file\n",argv[i]);
      else {
	 if (tp==1) { // FLT or CRJ UVIS image
	    WFC3getcards(0,2,0);
	    WFC3_DQ_mask(&fits,0,0,0);
	    WFC3_PAM_mult(0,"UVIS2wfc3_map.fits",wfc3_ctmult[2]);
	    WFC3setcards(0,0);
	    WFC3getcards(1,1,0);
	    WFC3_DQ_mask(&fits,1,0,0);
	    WFC3_PAM_mult(1,"UVIS1wfc3_map.fits",wfc3_ctmult[1]);
	    WFC3setcards(1,0);
	    addcard(fits.ext,"DOL_WFC3=                    2 / DOLPHOT WFC3 tag                               ");
	    addcard(fits.ext+1,"DOL_WFC3=                    1 / DOLPHOT WFC3 tag                               ");
	    fits.Next=2;
	 }
	 else if (tp==5) { // FLT or CRJ UVIS1 only
	    WFC3getcards(0,1,0);
	    WFC3_DQ_mask(&fits,0,0,0);
	    WFC3_PAM_mult(0,"UVIS1wfc3_map.fits",wfc3_ctmult[1]);
	    WFC3setcards(0,0);
	    addcard(fits.ext,"DOL_WFC3=                    1 / DOLPHOT WFC3 tag                               ");
	    if (offsetx!=0) {
	       sprintf(card,"DOL_OFFX=                 %4d / Origin of subarray relative to full chip       ",offsetx);
	       addcard(fits.ext,card);
	    }
	    if (offsety!=0) {
	       sprintf(card,"DOL_OFFY=                 %4d / Origin of subarray relative to full chip       ",offsety);
	       addcard(fits.ext,card);
	    }
	    fits.Next=1;
	 }
	 else if (tp==6) { // FLT or CRJ UVIS2 only
	    WFC3getcards(0,2,0);
	    WFC3_DQ_mask(&fits,0,0,0);
	    WFC3_PAM_mult(0,"UVIS2wfc3_map.fits",wfc3_ctmult[2]);
	    WFC3setcards(0,0);
	    addcard(fits.ext,"DOL_WFC3=                    2 / DOLPHOT WFC3 tag                               ");
	    if (offsetx!=0) {
	       sprintf(card,"DOL_OFFX=                 %4d / Origin of subarray relative to full chip       ",offsetx);
	       addcard(fits.ext,card);
	    }
	    if (offsety!=0) {
	       sprintf(card,"DOL_OFFY=                 %4d / Origin of subarray relative to full chip       ",offsety);
	       addcard(fits.ext,card);
	    }
	    fits.Next=1;
	 }
	 else if (tp==2) { // DRZ UVIS image
	    WFC3getcards(0,-1,1);
	    WFC3_DQ_mask(&fits,0,1,0);
	    WFC3_exptime_mult(0,wfc3_ctmult[1]);
	    WFC3setcards(0,0);
	    fits.Next=1;
	    addcard(fits.ext,"DOL_WFC3=                   -2 / DOLPHOT WFC3 tag                               ");
	 }
	 else if (tp==3) { // FLT or CRJ IR image
	    WFC3getcards(0,0,0);
	    WFC3_DQ_mask(&fits,0,0,1);
	    WFC3_PAM_mult(0,"ir_wfc3_map.fits",wfc3_ctmult[0]);
	    WFC3_exptime_mult(0,1.0);
	    WFC3setcards(0,1);
	    fits.Next=1;
	    addcard(fits.ext,"DOL_WFC3=                    0 / DOLPHOT WFC3 tag                               ");
	    if (offsetx!=0) {
	       sprintf(card,"DOL_OFFX=                 %4d / Origin of subarray relative to full chip       ",offsetx);
	       addcard(fits.ext,card);
	    }
	    if (offsety!=0) {
	       sprintf(card,"DOL_OFFY=                 %4d / Origin of subarray relative to full chip       ",offsety);
	       addcard(fits.ext,card);
	    }
	 }
	 else if (tp==4) { // DRZ IR image
	    WFC3getcards(0,-2,1);
	    WFC3_DQ_mask(&fits,0,1,1);
	    WFC3_exptime_mult(0,wfc3_ctmult[0]);
	    WFC3setcards(0,1);
	    fits.Next=1;
	    addcard(fits.ext,"DOL_WFC3=                   -1 / DOLPHOT WFC3 tag                               ");
	 }
	 writefits(argv[i],&fits,1);
      }
      freefits(&fits);
   }
   return 0;
}
