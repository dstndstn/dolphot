#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#include <time.h>

// Configurable/shared inputs
typedef double double4[4];
int X0,Y0,nxy[2],ncmd[2];
double NSTAR=50000,xystep[2],CMDSTEP=0.125,mmin,cmin;
char xyfn[321]="",cmdfn[321]="";

// Local variables
static double **NXY,**NCMD;

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
	 i=(int)((m-mmin)/CMDSTEP+10)-10;
	 j=(int)((c-cmin)/CMDSTEP+10)-10;
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

void process(int ext,int chip,int*Nfake,double4**fake)
{
   int i,j,ii,jj,N,n;
   double *dptr;

   NCMD=(double**)calloc(ncmd[0],sizeof(double*)); assert(NCMD!=NULL);
   dptr=(double*)calloc(ncmd[0]*ncmd[1],sizeof(double)); assert(dptr!=NULL);
   for (i=0;i<ncmd[0];i++) NCMD[i]=dptr+i*ncmd[1];
   NXY=(double**)calloc(nxy[0],sizeof(double*)); assert(NXY!=NULL);
   dptr=(double*)calloc(nxy[0]*nxy[1],sizeof(double)); assert(dptr!=NULL);
   for (i=0;i<nxy[0];i++) NXY[i]=dptr+i*nxy[1];
   for (i=0;i<ncmd[0];i++) for (j=0;j<ncmd[1];j++) NCMD[i][j]=0;
   if (xyfn[0]) readxy(ext,chip);
   else for (i=0;i<nxy[0];i++) for (j=0;j<nxy[1];j++) NXY[i][j]=1./(double)(nxy[0]*nxy[1]);
   if (cmdfn[0]) readcmd(ext,chip);
   else for (i=0;i<ncmd[0];i++) for (j=0;j<ncmd[1];j++) NCMD[i][j]=1./(double)(ncmd[0]*ncmd[1]);
   *Nfake=0;
   for (i=0;i<nxy[0];i++) for (j=0;j<nxy[1];j++) for (ii=0;ii<ncmd[0];ii++) for (jj=0;jj<ncmd[1];jj++) (*Nfake)+=(int)(NXY[i][j]*NCMD[ii][jj]*NSTAR+1);
   *fake=(double4*)calloc(*Nfake,sizeof(double4)); assert(*fake!=NULL);
   *Nfake=0;
   for (i=0;i<nxy[0];i++) for (j=0;j<nxy[1];j++) for (ii=0;ii<ncmd[0];ii++) for (jj=0;jj<ncmd[1];jj++) {
      N=(int)(NXY[i][j]*NCMD[ii][jj]*NSTAR+ran2());
      for (n=0;n<N;n++) {
	 (*fake)[*Nfake][0] = X0+(j+ran2())*xystep[1];
	 (*fake)[*Nfake][1] = Y0+(i+ran2())*xystep[0];
	 (*fake)[*Nfake][2] = mmin+(ii+ran2())*CMDSTEP;
	 (*fake)[*Nfake][3] = cmin+(jj+ran2())*CMDSTEP;
	 (*Nfake)++;
      }
   }
   free(NXY[0]);
   free(NXY);
}
