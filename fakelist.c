#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <assert.h>
#include "fakeproc.h"

typedef enum {
   NONE,
   WFPC2,
   ACS,
   WFC3
} instrument;

typedef struct {
   instrument inst;
   int filt;
   int UBVRI;
} imdatatype;

typedef struct {
   instrument inst;
   int cm;
   int filt;
} hstmodetype;

#ifdef USEWFPC2
#include "wfpc2/wfpc2filters.h"
int *wfpc2fused;
#endif
#ifdef USEACS
#include "acs/acsfilters.h"
int *acsfused,acscm=1;
#endif
#ifdef USEWFC3
#include "wfc3/wfc3filters.h"
int *wfc3fused;
#endif

imdatatype imdata[2];
hstmodetype *hstmode;
int MODE=0,Nimg;
FILE *finfo;

void readinfo(char*basefn) {
   char str[321],*ptr;
   int i;

   sprintf(str,"%s.info",basefn);
   if ((finfo=fopen(str,"r"))==NULL) {
      fprintf(stderr,"Cannot read \"%s\"\n",str);
      exit(-1);
   }
   fgets(str,321,finfo);
   Nimg=strtol(str,&ptr,10);
   if (Nimg<1 || strcmp(ptr," sets of output data\n")) {
      fprintf(stderr,"%s.info is not a dolphot info file\n",basefn);
      exit(-1);
   }
   hstmode=(hstmodetype*)calloc(Nimg,sizeof(hstmodetype));
#ifdef USEWFPC2
   wfpc2fused=(int*)calloc(WFPC2_NFILTERS,sizeof(int)); assert(wfpc2fused!=NULL);
#endif
#ifdef USEACS
   acsfused=(int*)calloc(ACS_NFILTERS,sizeof(int)); assert(acsfused!=NULL);
#endif
#ifdef USEWFC3
   wfc3fused=(int*)calloc(WFC3_NFILTERS,sizeof(int)); assert(wfc3fused!=NULL);
#endif
   for (i=0;i<Nimg;i++) {
      fgets(str,321,finfo); // filename
      fgets(str,321,finfo); // epoch
   }
   return;
}

void setsmag(double*smag) {
   int f0=-1,i,sOK=0;
   float cui=0,mlt=1,dcol[7]={2.677834,1.930288,1.00000,0.484694,0.000000,-0.642059,-1.03668};

   for (i=0;i<7;i++) if (smag[i]<99.5) {
      if (f0==-1) {
	 f0=i;
	 cui=smag[i];
	 mlt=dcol[i];
      }
      else {
	 sOK=1;
	 cui-=smag[i];
	 mlt-=dcol[i];
      }
   }
   if (!sOK) {
      for (i=0;i<7;i++) smag[i]=99.999;
      return;
   }
   cui/=mlt;
   for (i=0;i<7;i++) smag[i]=smag[f0]+(dcol[i]-dcol[f0])*cui;
   return;
}

void outstar(double m,double c) {
   int i;
   double smag[7]={99.999,99.999,99.999,99.999,99.999,99.999,99.999};
   static double *VMAG=NULL,*dVMAG,*TMAG;

   if (MODE==1) {
      printf(" %6.3f %6.3f\n",m,m-c);
      return;
   }
   if (imdata[0].inst==NONE && imdata[1].inst==NONE) {
      smag[imdata[0].filt]=m;
      smag[imdata[1].filt]=smag[imdata[0].filt]-c;
   }
   else if (imdata[0].inst==imdata[1].inst) {
      int NFILTERS;
      switch(imdata[0].inst) {
#ifdef USEWFPC2
      case WFPC2:
	 NFILTERS = WFPC2_NFILTERS;
	 break;
#endif
#ifdef USEACS
      case ACS:
	 NFILTERS = ACS_NFILTERS;
	 break;
#endif
#ifdef USEWFC3
      case WFC3:
	 NFILTERS = WFC3_NFILTERS;
	 break;
#endif
      default:
	 fprintf(stderr,"Stupid error; report to developer\n");
	 exit(-1);
      }
      if (VMAG==NULL) {
	 VMAG=(double*)calloc(NFILTERS,sizeof(double)); assert(VMAG!=NULL);
	 dVMAG=(double*)calloc(NFILTERS,sizeof(double)); assert(dVMAG!=NULL);
	 TMAG=(double*)calloc(NFILTERS,sizeof(double)); assert(TMAG!=NULL);
      }
      for (i=0;i<NFILTERS;i++) {
	 VMAG[i]=TMAG[i]=99.999;
	 dVMAG[i]=9.999;
      }
      VMAG[imdata[0].filt]=m;
      VMAG[imdata[1].filt]=m-c;
      dVMAG[imdata[0].filt]=0.01;
      dVMAG[imdata[1].filt]=0.01;
      switch(imdata[0].inst) {
#ifdef USEWFPC2
      case WFPC2:
	 WFPC2transform(VMAG,dVMAG,TMAG);
	 break;
#endif
#ifdef USEACS
      case ACS:
	 ACStransform(acscm,VMAG,dVMAG,TMAG);
	 break;
#endif
#ifdef USEWFC3
      case WFC3:
	 WFC3transform(VMAG,dVMAG,TMAG);
	 break;
#endif
      default:
	 fprintf(stderr,"Stupid error; report to developer\n");
	 exit(-1);
      }
      if (imdata[0].UBVRI>=0) smag[imdata[0].UBVRI]=TMAG[imdata[0].filt];
      if (imdata[1].UBVRI>=0) smag[imdata[1].UBVRI]=TMAG[imdata[1].filt];
   }
   else {
      fprintf(stderr,"Cannot mix transformed and untransformed mags\n");
      exit(-1);
   }
   setsmag(smag);
#ifdef USEWFPC2
   for (i=0;i<WFPC2_NFILTERS;i++) if (wfpc2fused[i]) {
      if (i==imdata[0].filt) printf(" %6.3f",m);
      else if (i==imdata[1].filt) printf(" %6.3f",m-c);
      else printf(" %6.3f",WFPC2untransform(i,smag));
   }
#endif
#ifdef USEACS
   for (i=0;i<ACS_NFILTERS;i++) if (acsfused[i]) {
      if (i==imdata[0].filt) printf(" %6.3f",m);
      else if (i==imdata[1].filt) printf(" %6.3f",m-c);
      else if (acscm==0 || ACSfilters[i].zp[1]<0) printf(" %6.3f",ACSuntransform(0,i,smag));
      else printf(" %6.3f",ACSuntransform(1,i,smag));
   }
#endif
#ifdef USEWFC3
   for (i=0;i<WFC3_NFILTERS;i++) if (wfc3fused[i]) {
      if (i==imdata[0].filt) printf(" %6.3f",m);
      else if (i==imdata[1].filt) printf(" %6.3f",m-c);
      else printf(" %6.3f",WFC3untransform(i,smag));
   }
#endif
   printf("\n");
   return;
}

int procchip(void) {
   char str[321],*ptr;
#if defined(USEWFPC2) || defined(USEACS) || defined(USEWFC3)
   char *ptr2;
#endif
   int i,ext,chip,X1,Y1,Nfake;
   double4 *fake;

   while (fgets(str,321,finfo) && strncmp(str,"EXTENSION",9));
   if (strncmp(str,"EXTENSION",9)) return 0;
   ext=strtol(str+9,&ptr,10);
   if (strncmp(ptr," CHIP",5)) {
      fprintf(stderr,"Format error in .info\n");
      exit(-1);
   }
   chip=atoi(ptr+5);
   fgets(str,321,finfo);
   if (strcmp(str,"Limits\n")) {
      fprintf(stderr,"No limits information in .info file\n");
      exit(-1);
   }
   fscanf(finfo,"%d %d %d %d",&X0,&X1,&Y0,&Y1);
   fgets(str,321,finfo);
   nxy[0]=1+(Y1-Y0-1)/64;
   nxy[1]=1+(X1-X0-1)/64;
   xystep[0]=(double)(Y1-Y0)/nxy[0];
   xystep[1]=(double)(X1-X0)/nxy[1];
   for (i=0;i<Nimg;i++) hstmode[i].inst=NONE;
   fgets(str,321,finfo); // first "camera-specific info" line
#ifdef USEWFPC2
   if (!strncmp(str,"* WFPC2-specific info",21)) {
      fgets(str,321,finfo);
      while (!strncmp(str,"* image ",8)) {
	 i=strtol(str+8,&ptr,10);
	 assert(i>0 && i<=Nimg && ptr[0]==':' && ptr[1]==' ');
	 ptr+=2;
	 for (ptr2=ptr;*ptr2 && *ptr2!=' ';ptr2++);
	 if (!*ptr2 || ptr2[1]<'0' || ptr2[1]>'9') {fprintf(stderr,"Format error in .info\n"); exit(-1);}
	 *ptr2=0;
	 hstmode[i-1].filt=WFPC2findfilt(ptr);
	 hstmode[i-1].cm=strtol(ptr2+1,&ptr,10);
	 wfpc2fused[hstmode[i-1].filt]=1;
	 fgets(str,321,finfo);
      }
   }
#endif
#ifdef USEACS
   if (!strncmp(str,"* ACS-specific info",19)) {
      int ncm[2]={0,0};
      fgets(str,321,finfo);
      while (!strncmp(str,"* image ",8)) {
	 i=strtol(str+8,&ptr,10);
	 assert(i>0 && i<=Nimg && ptr[0]==':' && ptr[1]==' ');
	 ptr+=2;
	 for (ptr2=ptr;*ptr2 && *ptr2!=' ';ptr2++);
	 if (!*ptr2 || ptr2[1]<'0' || ptr2[1]>'9') {fprintf(stderr,"Format error in .info\n"); exit(-1);}
	 *ptr2=0;
	 hstmode[i-1].filt=ACSfindfilt(ptr);
	 hstmode[i-1].cm=strtol(ptr2+1,&ptr,10);
	 ncm[(1+hstmode[i-1].cm)/2]++;
	 acsfused[hstmode[i-1].filt]=1;
	 fgets(str,321,finfo);
      }
      if (ncm[0]>ncm[1]) acscm=0;
      else acscm=1;
   }
#endif
#ifdef USEWFC3
   if (!strncmp(str,"* WFC3-specific info",20)) {
      fgets(str,321,finfo);
      while (!strncmp(str,"* image ",8)) {
	 i=strtol(str+8,&ptr,10);
	 assert(i>0 && i<=Nimg && ptr[0]==':' && ptr[1]==' ');
	 ptr+=2;
	 for (ptr2=ptr;*ptr2 && *ptr2!=' ';ptr2++);
	 if (!*ptr2 || ptr2[1]<'0' || ptr2[1]>'9') {fprintf(stderr,"Format error in .info\n"); exit(-1);}
	 *ptr2=0;
	 hstmode[i-1].filt=WFC3findfilt(ptr);
	 hstmode[i-1].cm=strtol(ptr2+1,&ptr,10);
	 wfc3fused[hstmode[i-1].filt]=1;
	 fgets(str,321,finfo);
      }
   }
#endif
   process(ext,chip,&Nfake,&fake);
   for (i=0;i<Nfake;i++) {
      printf("%d %d %7.2f %7.2f",ext,chip,fake[i][0],fake[i][1]);
      outstar(fake[i][2],fake[i][3]);
   }
   return 1;
}

void usage(char*exe) {
   fprintf(stderr,"Usage: %s <phot> <Filt1> <Filt2> <F1min> <F1max> <Cmin> <Cmax>\n",exe);
   fprintf(stderr,"Flags\n");
   fprintf(stderr,"  -USECMD=<filename> reads CMD distribution of stars\n");
   fprintf(stderr,"  -USEXY=<filename>  reads position distribution of stars\n");
   fprintf(stderr,"  -NSTAR=#           set number of stars per chip (default 50000)\n");
   exit(-1);
}

int findUBVRI(char*str) {
   char UBVRI[7][2]={"U","B","V","R","I","J","H"};
   int i;

   for (i=0;i<7;i++) if (!strcmp(str,UBVRI[i])) return i;
   fprintf(stderr,"Error: cannot identify filter %s\n",str);
   exit(-1);
}

int findUBVRIc(char c) {
   char UBVRI[8]="UBVRIJH";
   int i;

   for (i=0;i<7;i++) if (c==UBVRI[i]) return i;
   return -1;
}

void findfilt(int i,char *str) {
#ifdef USEWFPC2
   if (!strncasecmp(str,"WFPC2_",6)) {
      imdata[i].inst = WFPC2;
      imdata[i].filt = WFPC2findfilt(str+6);
      imdata[i].UBVRI = findUBVRIc(WFPC2filters[imdata[i].filt].color);
   }
   else
#endif
#ifdef USEACS
   if (!strncasecmp(str,"ACS_",4)) {
      imdata[i].inst = ACS;
      imdata[i].filt = ACSfindfilt(str+4);
      imdata[i].UBVRI = findUBVRIc(ACSfilters[imdata[i].filt].color);
   }
   else
#endif
#ifdef USEWFC3
   if (!strncasecmp(str,"WFC3_",5)) {
      imdata[i].inst = WFC3;
      imdata[i].filt = WFC3findfilt(str+5);
      imdata[i].UBVRI = findUBVRIc(WFC3filters[imdata[i].filt].color);
   }
   else
#endif
   {
      imdata[i].inst = NONE;
      imdata[i].filt = findUBVRI(str);
      imdata[i].UBVRI = imdata[i].filt;
   }
}

int main(int argc,char**argv) {
   int i;
   double mmax,cmax;

   if (argc<8) usage(*argv);
#ifdef USEWFPC2
   WFPC2initfilters();
#endif
#ifdef USEACS
   ACSinitfilters();
#endif
#ifdef USEWFC3
   WFC3initfilters();
#endif
   findfilt(0,argv[2]);
   findfilt(1,argv[3]);
   mmin=CMDSTEP*((int)(atof(argv[4])/CMDSTEP+100.5)-100);
   mmax=CMDSTEP*((int)(atof(argv[5])/CMDSTEP+100.5)-100);
   cmin=CMDSTEP*((int)(atof(argv[6])/CMDSTEP+100.5)-100);
   cmax=CMDSTEP*((int)(atof(argv[7])/CMDSTEP+100.5)-100);
   if ((mmax<mmin) || cmax<cmin) usage(*argv);
   ncmd[0]=(int)((mmax-mmin)/CMDSTEP+0.5);
   ncmd[1]=(int)((cmax-cmin)/CMDSTEP+0.5);
   for (i=8;i<argc;i++) {
      if (!strncasecmp(argv[i],"-usecmd=",8)) strcpy(cmdfn,argv[i]+8);
      else if (!strncasecmp(argv[i],"-usexy=",7)) strcpy(xyfn,argv[i]+7);
      else if (!strncasecmp(argv[i],"-nstar=",7)) NSTAR=atof(argv[i]+7);
      else usage(*argv);
   }
   readinfo(argv[1]);
   while (procchip());
   fclose(finfo);
   return 0;
}
