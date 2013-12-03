#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <assert.h>
#include <time.h>
#include "wfc3filters.h"

#define CMDSTEP 0.125

int cmdfilt[2],nxy[2],ncmd[2],ubvri[2],Nimg,*filt,*cm,X0,X1,Y0,Y1,*fused;
double **NXY,xystep[2],**NCMD,lim[2][2],*exptime,NSTAR=50000;
double *MJD;
FILE *finfo;
char xyfn[321]="",cmdfn[321]="";

void readinfo(char*basefn) {
   char str[321],*ptr;
   int i;

   sprintf(str,"%s.info",basefn);
   if ((finfo=fopen(str,"r"))==NULL) {
      printf("Cannot read \"%s\"\n",str);
      exit(-1);
   }
   fgets(str,321,finfo);
   Nimg=strtol(str,&ptr,10);
   if (Nimg<1 || strcmp(ptr," sets of output data\n")) {
      printf("%s.info is not a dolphot info file\n",basefn);
      exit(-1);
   }
   filt=(int*)calloc(Nimg,sizeof(int)); assert(filt!=NULL);
   fused=(int*)calloc(WFC3_NFILTERS,sizeof(int)); assert(fused!=NULL);
   cm=(int*)calloc(Nimg,sizeof(int)); assert(cm!=NULL);
   exptime=(double*)calloc(Nimg,sizeof(double)); assert(exptime!=NULL);
   MJD=(double*)calloc(Nimg,sizeof(double)); assert(MJD!=NULL);
   for (i=0;i<Nimg;i++) {
      fgets(str,321,finfo);
      fgets(str,321,finfo);
      MJD[i]=atof(str);
   }
   return;
}

void readxy(int ext0,int chip0) {
   FILE *f;
   int ext,chip,i,j,TOT=0;
   double x,y;
   char str[321];

   if ((f=fopen(xyfn,"r"))==NULL) {
      printf("Cannot read %s\n",xyfn);
      exit(-1);
   }
   while (fscanf(f,"%d %d %lf %lf",&ext,&chip,&x,&y)==4) {
      if (ext==ext0 && chip==chip0 && x>=X0 && y>=Y0) {
	 i=(int)((y-Y0)/xystep[0]);
	 j=(int)((x-X0)/xystep[1]);
	 if (i<nxy[0] && j<nxy[1]) {
	    NXY[i][j]++;
	    TOT++;
	 }
      }
      fgets(str,321,f);
   }
   fclose(f);
   for (i=0;i<nxy[0];i++) for (j=0;j<nxy[1];j++) NXY[i][j]/=(double)TOT;
   return;
}

void addcmd(int i,int j,double m) {
   if (i>=0 && i<ncmd[0] && j>=0 && j<ncmd[1]) NCMD[i][j]+=m;
}

void readcmd(int ext0,int chip0) {
   FILE *f;
   int ext,chip,i,ii,j;
   double x,y,m,c,TOT=0.;
   char str[321];

   if ((f=fopen(cmdfn,"r"))==NULL) {
      printf("Cannot read %s\n",cmdfn);
      exit(-1);
   }
   while (fscanf(f,"%d %d %lf %lf %lf %lf",&ext,&chip,&x,&y,&m,&c)==6) {
      if (ext==ext0 && chip==chip0) {
	 i=(int)((m-lim[0][0])/CMDSTEP+10)-10;
	 j=(int)((c-lim[1][0])/CMDSTEP+10)-10;
	 addcmd(i-1,j-1,0.25);
	 addcmd(i-1,j,0.5);
	 addcmd(i-1,j+1,0.25);
	 addcmd(i,j-1,0.5);
	 addcmd(i,j,1.);
	 addcmd(i,j+1,0.5);
	 addcmd(i+1,j-1,0.25);
	 addcmd(i+1,j,0.5);
	 addcmd(i+1,j+1,0.25);
      }
      fgets(str,321,f);
   }
   fclose(f);
   for (i=0;i<ncmd[0]-1;i++) for (ii=i+1;ii<ncmd[0];ii++) for (j=0;j<ncmd[1];j++) if (NCMD[i][j]>NCMD[ii][j]) NCMD[ii][j]=NCMD[i][j];
   for (i=0;i<ncmd[0];i++) for (j=0;j<ncmd[1];j++) TOT+=NCMD[i][j];
   for (i=0;i<ncmd[0];i++) for (j=0;j<ncmd[1];j++) NCMD[i][j]/=(double)TOT;
   return;
}

// from NR;
#define IM1 2147483563
#define IM2 2147483399
#define AM (1.0/IM1)
#define IMM1 (IM1-1)
#define IA1 40014
#define IA2 40692
#define IQ1 53668
#define IQ2 52774
#define IR1 12211
#define IR2 3791
#define NTAB 32
#define NDIV (1+IMM1/NTAB)
#define EPS 1.2e-7
#define RNMX (1.0-EPS)
double ran2(void) {
   int j;
   long k;
   static long idum2=123456789;
   static long iy=0;
   static long iv[NTAB];
   double temp;
   static long seed=-1;

   if (seed <= 0) {
      seed=-1-time(NULL);
      if (-(seed) < 1) seed=1;
      else seed = -(seed);
      idum2=(seed);
      for (j=NTAB+7;j>=0;j--) {
	 k=(seed)/IQ1;
	 seed=IA1*(seed-k*IQ1)-k*IR1;
	 if (seed < 0) seed += IM1;
	 if (j < NTAB) iv[j] = seed;
      }
      iy=iv[0];
   }
   k=(seed)/IQ1;
   seed=IA1*(seed-k*IQ1)-k*IR1;
   if (seed < 0) seed += IM1;
   k=idum2/IQ2;
   idum2=IA2*(idum2-k*IQ2)-k*IR2;
   if (idum2 < 0) idum2 += IM2;
   j=iy/NDIV;
   iy=iv[j]-idum2;
   iv[j] = seed;
   if (iy < 1) iy += IMM1;
   if ((temp=AM*iy) > RNMX) return RNMX;
   else return temp;
}
#undef IM1
#undef IM2
#undef AM
#undef IMM1
#undef IA1
#undef IA2
#undef IQ1
#undef IQ2
#undef IR1
#undef IR2
#undef NTAB
#undef NDIV
#undef EPS
#undef RNMX

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

   if (cmdfilt[0]<0 && cmdfilt[1]<0) {
      smag[-cmdfilt[0]-1]=m;
      smag[-cmdfilt[1]-1]=smag[-cmdfilt[0]-1]-c;
   }
   else if (cmdfilt[0]>=0 && cmdfilt[1]>=0) {
      if (VMAG==NULL) {
	 VMAG=(double*)calloc(WFC3_NFILTERS,sizeof(double)); assert(VMAG!=NULL);
	 dVMAG=(double*)calloc(WFC3_NFILTERS,sizeof(double)); assert(dVMAG!=NULL);
	 TMAG=(double*)calloc(WFC3_NFILTERS,sizeof(double)); assert(TMAG!=NULL);
      }
      for (i=0;i<WFC3_NFILTERS;i++) {
	 VMAG[i]=TMAG[i]=99.999;
	 dVMAG[i]=9.999;
      }
      VMAG[cmdfilt[0]]=m;
      VMAG[cmdfilt[1]]=m-c;
      dVMAG[cmdfilt[0]]=0.01;
      dVMAG[cmdfilt[1]]=0.01;
      WFC3transform(VMAG,dVMAG,TMAG);
      if (ubvri[0]>=0) smag[ubvri[0]]=TMAG[cmdfilt[0]];
      if (ubvri[1]>=0) smag[ubvri[1]]=TMAG[cmdfilt[1]];
   }
   else {
      printf("Cannot mix transformed and untransformed mags\n");
      exit(-1);
   }
   setsmag(smag);
   for (i=0;i<WFC3_NFILTERS;i++) if (fused[i]) {
      if (i==cmdfilt[0]) printf(" %6.3f",m);
      else if (i==cmdfilt[1]) printf(" %6.3f",m-c);
      else printf(" %6.3f",WFC3untransform(i,smag));
   }
   printf("\n");
   return;
}

int procchip(void) {
   char str[321],*ptr,*ptr2;
   double *dptr,Nd;
   int i,j,ii,jj,N,n,ext,chip;

   while (fgets(str,321,finfo) && strncmp(str,"EXTENSION",9));
   if (strncmp(str,"EXTENSION",9)) return 0;
   ext=strtol(str+9,&ptr,10);
   if (strncmp(ptr," CHIP",5)) {printf("Format error in .info\n"); exit(-1);}
   chip=atoi(ptr+5);
   fgets(str,321,finfo);
   if (strcmp(str,"Limits\n")) {
      printf("No limits information in .info file\n");
      exit(-1);
   }
   fscanf(finfo,"%d %d %d %d",&X0,&X1,&Y0,&Y1);
   nxy[0]=1+(Y1-Y0-1)/64;
   nxy[1]=1+(X1-X0-1)/64;
   xystep[0]=(double)(Y1-Y0)/nxy[0];
   xystep[1]=(double)(X1-X0)/nxy[1];
   NXY=(double**)calloc(nxy[0],sizeof(double*)); assert(NXY!=NULL);
   dptr=(double*)calloc(nxy[0]*nxy[1],sizeof(double)); assert(dptr!=NULL);
   for (i=0;i<nxy[0];i++) NXY[i]=dptr+i*nxy[1];
   fgets(str,321,finfo);
   fgets(str,321,finfo);
   while (!feof(finfo) && strcmp(str,"* WFC3-specific info\n")) fgets(str,321,finfo);
   if (strcmp(str,"* WFC3-specific info\n")) {
      printf("No WFC3 filter/exptime information in .info file\n");
      exit(-1);
   }
   for (i=0;i<Nimg;i++) {
      fgets(str,321,finfo);
      if (strncmp(str,"* image ",8)) {printf("Format error in .info\n"); exit(-1);}
      ptr=str+8;
      for (;*ptr && *ptr!=' ';ptr++);
      for (;*ptr && *ptr==' ';ptr++);
      for (ptr2=ptr;*ptr2 && *ptr2!=' ';ptr2++);
      if (!*ptr2 || ptr2[1]<'0' || ptr2[1]>'9') {printf("Format error in .info\n"); exit(-1);}
      *ptr2=0;
      filt[i]=WFC3findfilt(ptr);
      cm[i]=strtol(ptr2+1,&ptr,10);
      exptime[i]=atof(ptr);
      fused[filt[i]]=1;
   }
   for (i=0;i<ncmd[0];i++) for (j=0;j<ncmd[1];j++) NCMD[i][j]=0;
   if (xyfn[0]) readxy(ext,chip);
   else for (i=0;i<nxy[0];i++) for (j=0;j<nxy[1];j++) NXY[i][j]=1./(double)(nxy[0]*nxy[1]);
   if (cmdfn[0]) readcmd(ext,chip);
   else for (i=0;i<ncmd[0];i++) for (j=0;j<ncmd[1];j++) NCMD[i][j]=1./(double)(ncmd[0]*ncmd[1]);
   for (i=0;i<nxy[0];i++) for (j=0;j<nxy[1];j++) for (ii=0;ii<ncmd[0];ii++) for (jj=0;jj<ncmd[1];jj++) {
      Nd=NXY[i][j]*NCMD[ii][jj]*NSTAR;
      N=(int)Nd;
      if (ran2()<Nd-N) N++;
      for (n=0;n<N;n++) {
	 printf("%d %d %7.2f %7.2f",ext,chip,X0+(j+ran2())*xystep[1],Y0+(i+ran2())*xystep[0]);
	 outstar(lim[0][0]+(ii+ran2())*CMDSTEP,lim[1][0]+(jj+ran2())*CMDSTEP);
      }
   }
   free(NXY[0]);
   free(NXY);
   return 1;
}

void usage(char*exe) {
   printf("Usage: %s <phot> <Filt1> <Filt2> <F1min> <F1max> <Cmin> <Cmax>\n",exe);
   printf("Flags\n");
   printf("  -USECMD=<filename> reads CMD distribution of stars\n");
   printf("  -USEXY=<filename>  reads position distribution of stars\n");
   printf("  -NSTAR=#           set number of stars per chip (default 50000)\n");
   exit(-1);
}

int findUBVRIJH(char*str) {
   char UBVRIJH[7][2]={"U","B","V","R","I","J","H"};
   int i;

   for (i=0;i<7;i++) if (!strcmp(str,UBVRIJH[i])) return -1-i;
   return WFC3findfilt(str);
}

int findUBVRIJHc(char c) {
   char UBVRIJH[8]="UBVRIJH";
   int i;

   for (i=0;i<7;i++) if (c==UBVRIJH[i]) return i;
   return -1;
}

int main(int argc,char**argv) {
   int i;
   double*ptr;

   fprintf(stderr,"Warning: wfc3fakelist is deprecated and will be removed in a future release.\nConsider using fakelist instead.\n");
   if (argc<8) usage(*argv);
   WFC3initfilters();
   cmdfilt[0]=findUBVRIJH(argv[2]);
   if (cmdfilt[0]>=0) ubvri[0]=findUBVRIJHc(WFC3filters[cmdfilt[0]].color);
   cmdfilt[1]=findUBVRIJH(argv[3]);
   if (cmdfilt[1]>=0) ubvri[1]=findUBVRIJHc(WFC3filters[cmdfilt[1]].color);
   lim[0][0]=CMDSTEP*((int)(atof(argv[4])/CMDSTEP+100.5)-100);
   lim[0][1]=CMDSTEP*((int)(atof(argv[5])/CMDSTEP+100.5)-100);
   lim[1][0]=CMDSTEP*((int)(atof(argv[6])/CMDSTEP+100.5)-100);
   lim[1][1]=CMDSTEP*((int)(atof(argv[7])/CMDSTEP+100.5)-100);
   if ((lim[0][1]<lim[0][0]) || lim[1][1]<lim[1][0]) usage(*argv);
   ncmd[0]=(int)((lim[0][1]-lim[0][0])/CMDSTEP+0.5);
   ncmd[1]=(int)((lim[1][1]-lim[1][0])/CMDSTEP+0.5);
   NCMD=(double**)calloc(ncmd[0],sizeof(double*)); assert(NCMD!=NULL);
   ptr=(double*)calloc(ncmd[0]*ncmd[1],sizeof(double)); assert(ptr!=NULL);
   for (i=0;i<ncmd[0];i++) NCMD[i]=ptr+i*ncmd[1];
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
