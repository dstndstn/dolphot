#include "dolphot_common.h"
#include <assert.h>

typedef float float3[3];

// uncomment out to use single aperture correction
//#define SINGLE_APCOR

//int use_mmap = 0;

// Standard photometry flags
int IMDIFF=0,NOSKY34=0;
fntype outfn;
int XMIN,XMAX,YMIN,YMAX;
int Nstars;

/*
#include <stdint.h>
#include <mach/mach_time.h>
void tictoc(int i) {
   static uint64_t start;
   if (i==0) {
      start = mach_absolute_time();
   }
   else {
      uint64_t elapsed;
      elapsed = mach_absolute_time() - start;
      printf("%f sec elapsed\n",elapsed/1.0e9);
   }
}
*/

//NR Utilities

//NR Routines for FFT;

#define NR_END 1
#define FREE_ARG char*
#define SIGN(a,b) ((b) >= 0.0 ? fabs(a) : -fabs(a))
static float maxarg1,maxarg2;
#define FMAX(a,b) (maxarg1=(a),maxarg2=(b),(maxarg1) > (maxarg2) ?\
        (maxarg1) : (maxarg2))

void nrerror(char error_text[]) {
   fprintf(stderr,"Numerical Recipes run-time error...\n");
   fprintf(stderr,"%s\n",error_text);
   fprintf(stderr,"...now exiting to system...\n");
   exit(1);
}

float *vector(int nl, int nh) {
   float *v;

   v=(float *)malloc((size_t) ((nh-nl+1+NR_END)*FLOATSIZE));
   if (!v) nrerror("allocation failure in vector()");
   return v-nl+NR_END;
}

void free_vector(float *v, int nl, int nh) {
   free((FREE_ARG) (v+nl-NR_END));
   return;
}


float **matrix(int nrl, int nrh, int ncl, int nch) {
   int i, nrow=nrh-nrl+1,ncol=nch-ncl+1;
   float **m;

   m=(float **) malloc((size_t)((nrow+NR_END)*PTRSIZE));
   if (!m) nrerror("allocation failure 1 in matrix()");
   m += NR_END;
   m -= nrl;

   m[nrl]=(float *) malloc((size_t)((nrow*ncol+NR_END)*FLOATSIZE));
   if (!m[nrl]) nrerror("allocation failure 2 in matrix()");
   m[nrl] += NR_END;
   m[nrl] -= ncl;

   for(i=nrl+1;i<=nrh;i++) m[i]=m[i-1]+ncol;

   return m;
}

void free_matrix(float **m, int nrl, int nrh, int ncl, int nch) {
   free((FREE_ARG) (m[nrl]+ncl-NR_END));
   free((FREE_ARG) (m+nrl-NR_END));
   return;
}

//photometry variables
int fitpsf,NO_SKY=0,NOSUBTRACT=0;
int **indx,*tindx;
float **snmap,**snmmap;
char **ran,**snmapflag;

//file and frame variables;
reopenableFile *fdata,*fsky,*fpsf,*fres;
FILE *of,*fapcor=NULL,*fpsfs=NULL,*ffakeout,*fwarn,*fverb;
ftype refpsf;
chiptype refpsfimg=NULL;

#ifdef PGPLOT
#include PGHEAD
typedef struct {
   float x,y,dx,dy;
} alignDiagDataType;
typedef struct {
   int N,XMIN,XMAX,YMIN,YMAX;
   float sig;
   alignDiagDataType *data;
} alignDiagType;
alignDiagType ***alignDiagData=0;
typedef struct {
   int N;
   float dmmean,mmin,mmax,dmmin,dmmax;
   float *dm,*m;
   int inst,filt;
   char imgstr[161];
   char filtstr[41];
} apcorDiagType;
apcorDiagType ***apcorDiagData=0;
typedef struct {
   int N;
   float mid;
   float*data;
} psfDiagType;
psfDiagType ***psfDiagData=0;

static int DIAG_NEXT,*DIAG_NZ;
static int DIAG_EXT,DIAG_Z;

void allocDiagData1(int Next) {
   int i;
   if (DiagPlotType[0]==0) return;
   DIAG_NEXT = Next;
   alignDiagData = (alignDiagType***)calloc(Nimg,PTRSIZE);
   apcorDiagData = (apcorDiagType***)calloc(Nimg,PTRSIZE);
   psfDiagData = (psfDiagType***)calloc(Nimg,PTRSIZE);
   DIAG_NZ = (int*)calloc(Next+1,INTSIZE);
   if (!alignDiagData || !apcorDiagData || !psfDiagData || !DIAG_NZ) merr();
   for (i=0;i<Nimg;i++) {
      alignDiagData[i] = (alignDiagType**)calloc(Next+1,PTRSIZE);
      apcorDiagData[i] = (apcorDiagType**)calloc(Next+1,PTRSIZE);
      psfDiagData[i] = (psfDiagType**)calloc(Next+1,PTRSIZE);
      if (!alignDiagData[i] || !apcorDiagData[i] || !psfDiagData[i]) merr();
   }
}

void allocDiagData2(int ext) {
   int i;
   if (DiagPlotType[0]==0) return;
   DIAG_NZ[ext] = dataim[0].Z;
   for (i=0;i<Nimg;i++) {
      alignDiagData[i][ext] = (alignDiagType*)calloc(dataim[0].Z,sizeof(alignDiagType));
      apcorDiagData[i][ext] = (apcorDiagType*)calloc(dataim[0].Z,sizeof(apcorDiagType));
      psfDiagData[i][ext] = (psfDiagType*)calloc(dataim[0].Z,sizeof(psfDiagType));
      if (!alignDiagData[i][ext] || !apcorDiagData[i][ext] || !psfDiagData[i][ext]) merr();
   }
}

void setDiagFrame(int ext,int z) {
   DIAG_EXT = ext;
   DIAG_Z = z;
}

void plotDiagData(void) {
   int img,ext,z,i,j,nplot;
   int insts[100],filts[100],Nif=0,Nifs[100];
   float dmmean[100][100],midpix[100][100];
   char str[161],filtstrs[100][41];
   if (DiagPlotType[0]==0) return;
   if (!strcasecmp(DiagPlotType,"gif")) sprintf(str,"%s.diag.gif/gif",outfn);
   else if (!strcasecmp(DiagPlotType,"ps")) sprintf(str,"%s.diag.ps/vcps",outfn);
   else if (!strcasecmp(DiagPlotType,"png")) sprintf(str,"%s.diag.png/png",outfn);
   else {
      printf("Error: unknown DiagPlotType %s\n",DiagPlotType);
      return;
   }
   cpgopen(str);
   cpgpap(0,1);
   cpgsubp(2,2);
   nplot=4;
   for (img=0;img<Nimg;img++) {
      for (ext=0;ext<=DIAG_NEXT;ext++) {
	 for (z=0;z<DIAG_NZ[ext];z++) {
	    int first=1;
	    if (apcorDiagData[img][ext][z].N>0) {
	       float x[2],y[2];
	       while (nplot<4) {cpgpage(); nplot++;}
	       first=0; nplot=1;
	       x[0] = apcorDiagData[img][ext][z].mmin;
	       x[1] = apcorDiagData[img][ext][z].mmax;
	       y[0] = apcorDiagData[img][ext][z].dmmean;
	       y[1] = apcorDiagData[img][ext][z].dmmean;
	       cpgenv(apcorDiagData[img][ext][z].mmin-0.1,apcorDiagData[img][ext][z].mmax+0.1,apcorDiagData[img][ext][z].dmmax+0.05,apcorDiagData[img][ext][z].dmmin-0.05,0,0);
	       cpgpt(apcorDiagData[img][ext][z].N,apcorDiagData[img][ext][z].m,apcorDiagData[img][ext][z].dm,-4);
	       cpgsci(2);
	       cpgline(2,x,y);
	       cpgsci(1);
	       cpglab("Instrumental mag","Aperture Correction",apcorDiagData[img][ext][z].imgstr);
	       if (psfDiagData[img][ext][z].N>0) {
		  for (i=0;i<Nif && (apcorDiagData[img][ext][z].inst!=insts[i] || apcorDiagData[img][ext][z].filt!=filts[i]);i++);
		  if (i==Nif) {
		     strcpy(filtstrs[Nif],apcorDiagData[img][ext][z].filtstr);
		     insts[Nif] = apcorDiagData[img][ext][z].inst;
		     filts[Nif++] = apcorDiagData[img][ext][z].filt;
		  }
		  dmmean[i][Nifs[i]] = apcorDiagData[img][ext][z].dmmean;
		  midpix[i][Nifs[i]] = psfDiagData[img][ext][z].mid;
		  Nifs[i]++;
	       }
	    }
	    if (psfDiagData[img][ext][z].N>0) {
	       float tr[6]={0,0,0,0,0,0};
	       float min=0,max=0;
	       if (first) {
		  while (nplot<4) {cpgpage(); nplot++;}
		  first=0; nplot=1;
	       }
	       else nplot++;
	       for (i=0;i<psfDiagData[img][ext][z].N*psfDiagData[img][ext][z].N;i++) {
		  if (psfDiagData[img][ext][z].data[i]<min) min=psfDiagData[img][ext][z].data[i];
		  else if (psfDiagData[img][ext][z].data[i]>max) max=psfDiagData[img][ext][z].data[i];
	       }
	       tr[0]=tr[3]=-(psfDiagData[img][ext][z].N+1)/2;
	       tr[1]=tr[5]=1;
	       cpgenv(tr[0]+0.5,tr[0]+psfDiagData[img][ext][z].N+0.5,tr[0]+0.5,tr[0]+psfDiagData[img][ext][z].N+0.5,1,0);
	       cpggray(psfDiagData[img][ext][z].data,psfDiagData[img][ext][z].N,psfDiagData[img][ext][z].N,1,psfDiagData[img][ext][z].N,1,psfDiagData[img][ext][z].N,min,max,tr);
	       sprintf(str,"PSF Correction Image (max=%5.3f, min=%5.3f)",max,min);
	       cpglab("X","Y",str);
	    }
	    if (alignDiagData[img][ext][z].N>0) {
	       float *x,*y,sig;
	       int N;
	       if (first) {
		  while (nplot<4) {cpgpage(); nplot++;}
		  first=0; nplot=2;
	       }
	       else nplot+=2;
	       x = (float*)calloc(alignDiagData[img][ext][z].N,FLOATSIZE);
	       y = (float*)calloc(alignDiagData[img][ext][z].N,FLOATSIZE);
	       if (!x || !y) merr();
	       sig = alignDiagData[img][ext][z].sig*0.707;
	       if (sig<0.1) sig=0.1;
	       cpgenv(alignDiagData[img][ext][z].XMIN,alignDiagData[img][ext][z].XMAX,alignDiagData[img][ext][z].YMIN,alignDiagData[img][ext][z].YMAX,1,0);
	       N=0;
	       for (i=0;i<alignDiagData[img][ext][z].N;i++) {
		  if (fabs(alignDiagData[img][ext][z].data[i].dx)<=sig) {
		     x[N]=alignDiagData[img][ext][z].data[i].x;
		     y[N++]=alignDiagData[img][ext][z].data[i].y;
		  }
	       }
	       cpgpt(N,x,y,-4);
	       N=0;
	       for (i=0;i<alignDiagData[img][ext][z].N;i++) {
		  if (alignDiagData[img][ext][z].data[i].dx<-sig) {
		     x[N]=alignDiagData[img][ext][z].data[i].x;
		     y[N++]=alignDiagData[img][ext][z].data[i].y;
		  }
	       }
	       cpgsci(4);
	       cpgpt(N,x,y,-4);
	       N=0;
	       for (i=0;i<alignDiagData[img][ext][z].N;i++) {
		  if (alignDiagData[img][ext][z].data[i].dx>sig) {
		     x[N]=alignDiagData[img][ext][z].data[i].x;
		     y[N++]=alignDiagData[img][ext][z].data[i].y;
		  }
	       }
	       cpgsci(2);
	       cpgpt(N,x,y,-4);
	       cpgsci(1);
	       sprintf(str,"X Alignment Residuals (threshold=%5.3f)",sig);
	       cpglab("Xref","Yref",str);
	       cpgenv(alignDiagData[img][ext][z].XMIN,alignDiagData[img][ext][z].XMAX,alignDiagData[img][ext][z].YMIN,alignDiagData[img][ext][z].YMAX,1,0);
	       N=0;
	       for (i=0;i<alignDiagData[img][ext][z].N;i++) {
		  if (fabs(alignDiagData[img][ext][z].data[i].dy)<=sig) {
		     x[N]=alignDiagData[img][ext][z].data[i].x;
		     y[N++]=alignDiagData[img][ext][z].data[i].y;
		  }
	       }
	       cpgpt(N,x,y,-4);
	       N=0;
	       for (i=0;i<alignDiagData[img][ext][z].N;i++) {
		  if (alignDiagData[img][ext][z].data[i].dy<-sig) {
		     x[N]=alignDiagData[img][ext][z].data[i].x;
		     y[N++]=alignDiagData[img][ext][z].data[i].y;
		  }
	       }
	       cpgsci(4);
	       cpgpt(N,x,y,-4);
	       N=0;
	       for (i=0;i<alignDiagData[img][ext][z].N;i++) {
		  if (alignDiagData[img][ext][z].data[i].dy>sig) {
		     x[N]=alignDiagData[img][ext][z].data[i].x;
		     y[N++]=alignDiagData[img][ext][z].data[i].y;
		  }
	       }
	       cpgsci(2);
	       cpgpt(N,x,y,-4);
	       cpgsci(1);
	       sprintf(str,"Y Alignment Residuals (threshold=%5.3f)",sig);
	       cpglab("Xref","Yref",str);
	       free(x);
	       free(y);
	    }
	 }
      }
   }
   for (i=0;i<Nif;i++) {
      float xmin,xmax,ymin,ymax;
      xmin=xmax=midpix[i][0];
      ymin=ymax=dmmean[i][0];
      for (j=1;j<Nifs[i];j++) {
	 if (midpix[i][j]<xmin) xmin=midpix[i][j];
	 else if (midpix[i][j]>xmax) xmax=midpix[i][j];
	 if (dmmean[i][j]<ymin) ymin=dmmean[i][j];
	 else if (dmmean[i][j]>ymax) ymax=dmmean[i][j];
      }
      if (i==0) {
	 while (nplot<4) {cpgpage(); nplot++;}
      }
      cpgenv(xmin-0.001,xmax+0.001,ymin-0.01,ymax+0.01,0,0);
      cpgpt(Nifs[i],midpix[i],dmmean[i],-4);
      cpglab("Central PSF Adjustment","Aperture Correction",filtstrs[i]);
   }
   cpgend();
}
#endif

Inline float skyval(int img,int x,int y) {
   if (fsky[img].lastoffset<0) return 0.;
   if (x<0) x=0; else if (x>=dataim[img].X) x=dataim[img].X-1;
   if (y<0) y=0; else if (y>=dataim[img].Y) y=dataim[img].Y-1;
   return sky[img][y][x];
}

Inline float noise(int img,int x,int y,float ssky) {
   float n;
   n=data[img][y][x]-ssky;
   if (n<0) n=0.;
   if (ssky>0.) n+=ssky;
   return n/iGAIN[img]+iRN[img];
}

Inline float modelnoise(int img,int x,int y,int dx,int dy,float s,float ssky) {
   float n;

   n=data[img][y][x]-res[img][y][x];
   if (n<0.) n=0.;
   if (psf[dy][dx]>0. && s>0.) n+=s*psf[dy][dx];
   if (n<data[img][y][x]-ssky) n=data[img][y][x]-ssky;
   if (ssky>0.) n+=ssky;
   return n/iGAIN[img]+iRN[img];
}

Inline int peak2(int img,int x,int y) {
   if (x>0 && data[img][y][x-1]>data[img][y][x] && data[img][y][x-1]<iDMAX[img]) return 0;
   if (x<dataim[img].X-1 && data[img][y][x+1]>data[img][y][x] && data[img][y][x+1]<iDMAX[img]) return 0;
   if (y>0 && data[img][y-1][x]>data[img][y][x] && data[img][y-1][x]<iDMAX[img]) return 0;
   if (y<dataim[img].Y-1 && data[img][y+1][x]>data[img][y][x] && data[img][y+1][x]<iDMAX[img]) return 0;
   return 1;
}

void setindx(int i) {
   int j;
   if (i==0) memset(indx[0],0,FLOATSIZE*X*Y*SubResRef*SubResRef);
   else {
      for (j=0;j<X*SubResRef;j++) tindx[j]=i;
      for (j=0;j<Y*SubResRef;j++) memcpy(indx[j],tindx,FLOATSIZE*X*SubResRef);
   }
   return;
}

void markstars(int x,int y) {
   int xx,yy;
   x*=SubResRef;
   y*=SubResRef;
   int dx = rMark*SubResRef;
   for (xx=-dx;xx<=dx;xx++) if (x+xx>=XMIN*SubResRef && x+xx<XMAX*SubResRef) for (yy=-dx;yy<=dx;yy++) if (y+yy>=YMIN*SubResRef && y+yy<YMAX*SubResRef) if (indx[y+yy][x+xx]>=0) stars[indx[y+yy][x+xx]].flag=1;
   return;
}

float calcskyval(int n0,float*list,int img,float*unc) {
   int i,n,CONT=1;
   double R=0,RR=0,smlt=2.25,tsig;

   if (smlt<SkySig) smlt=SkySig;
   n=n0;
   while (CONT && n) {
      R=RR=0.;
      for (i=0;i<n;i++) {
	 R+=list[i];
	 RR+=list[i]*list[i];
      }
      R/=n;
      if (n>1) RR=smlt*sqrt(1.+(RR-R*R*n)/(n-1));
      else RR=smlt;
      if (R<0) tsig=sqrt(1.+iRN[img])*3.5;
      else tsig=sqrt(1.+iRN[img]+R/iGAIN[img])*3.5;
      CONT=0;
      for (i=0;i<n;i++) if (fabs(list[i]-R)>RR) {
	 list[i--]=list[--n];
	 CONT=1;
      }
      if (!CONT && RR>tsig && smlt>SkySig) {
	 smlt*=0.9;
	 if (smlt<SkySig) smlt=SkySig;
	 CONT=1;
      }
   }
   if (!n) return 2.e30;
   if (unc!=NULL) *unc=RR/sqrt(n);
   return R;
}

void clear_sky(void) {
   int img;

   for (img=0;img<Timg;img++) sky_set[img]=0;
   return;
}

float getsky_norm(int img,int x,int y,float*unc) {
   static float *list=NULL;
   double s=0.;
#ifdef USESKYSIG
   double sig=0.;
#endif
   int n=0,x1,y1,dx,dy,rsky,r2;

   if (NO_SKY && fsky[img].lastoffset>=0) return skyval(img,x,y);
   if (sky_set[img]) {
#ifdef USESKYSIG
      if (unc) *unc=sky_unc[img];
#endif
      return sky_val[img];
   }
   rsky=(int)(RSky1[img]+0.999);
   if (list==NULL) {
      list=(float*)calloc((2*rsky+1)*(2*rsky+1),FLOATSIZE);
      if (!list) merr();
   }
   for (x1=x-rsky;x1<=x+rsky;x1+=SkipSky) for (y1=y-rsky;y1<=y+rsky;y1+=SkipSky) if (ppixOK(img,x1,y1) && (dx=abs(x-x1))>=RSky0[img] && (dy=abs(y-y1))>=RSky0[img] && (r2=dx*dx+dy*dy)>=RSky0[img]*RSky0[img] && r2<=RSky1[img]*RSky1[img]) {
      list[n++]=res[img][y1][x1];
#ifdef USESKYSIG
      sig+=noise(img,x1,y1,0.);
#endif
   }
#ifdef USESKYSIG
   if (n) sig=sqrt(sig)/n;
#endif
   s=calcskyval(n,list,img,NULL);
   if (s>1.e30) return 0.;
   sky_set[img]=1;
   sky_val[img]=s;
#ifdef USESKYSIG
   sky_unc[img]=sig;
   if (unc) *unc=sig;
#endif
   return s;
}

float getsky_small(int img,int x,int y,float*unc) {
   int x1,y1,i,N=0,cont;
   float r2,av,sd;
   static float *list;
   static int first=1;

   if (first) {
      list=(float*)calloc(4*rpsfMax*rpsfMax,FLOATSIZE);
      if (!list) merr();
      first=0;
   }
#ifdef USESKYSIG
   if (unc) *unc=0.;
#endif
   for (y1=y-RPSF[img];y1<=y+RPSF[img];y1++) for (x1=x-RPSF[img];x1<=x+RPSF[img];x1++) if (ppixOK(img,x1,y1)) {
      r2=sqrt((x1-x)*(x1-x)+(y1-y)*(y1-y));
      if (r2>RAper[img]+1 && r2<=RPSF[img]) {
	 list[N++]=res[img][y1][x1];
#ifdef USESKYSIG
	 if (unc) (*unc)+=noise(img,x1,y1,0.);
#endif
      }
   }
   if (unc && N) *unc=sqrt(*unc)/N;
   cont=1;
   av=sd=0;
   if (N==1) av=list[0];
   while (cont && N>1) {
      av=sd=0;
      for (i=0;i<N;i++) av+=list[i];
      av/=N;
      for (i=0;i<N;i++) sd+=(list[i]-av)*(list[i]-av);
      sd=sqrt(sd/(N-1))*2.5;
      cont=0;
      for (i=0;i<N;) {
	 if (fabs(list[i]-av)>sd) {
	    list[i]=list[--N];
	    cont=1;
	 }
	 else i++;
      }
   }
   return av;
}

float getsky(int img,int x,int y,float*unc) {
   if (FitSky==0 && fsky[img].lastoffset>=0) {
#ifdef USESKYSIG
      if (unc) *unc=0.;
#endif
      return skyval(img,x,y);
   }
   if (FitSky==2) return getsky_small(img,x,y,unc);
   return getsky_norm(img,x,y,unc);
}

void centroid(int IMG,int x0,int y0,double *X,double *Y) {
   int img,x,y,xx,yy,csat,nsat,nok;
   float min,xc1,yc1,xc=0,yc=0,wt,twt=0;
   for (img=0;img<Timg;img++) if (img==IMG || (IMG<0 && img<Nimg)) {
      double tx,ty;
      wt=xc1=yc1=0;
      shift(img,x0,y0,&tx,&ty,1);
      x=(int)(tx+0.5);
      y=(int)(ty+0.5);
      min=0;
      csat=nsat=nok=0;
      for (yy=-RCentroid;yy<=RCentroid;yy++) for (xx=-RCentroid;xx<=RCentroid;xx++) if (abs(xx)+abs(yy)<RCentroid*2 && ppixOK(img,x+xx,y+yy) && res[img][y+yy][x+xx]<min) min=res[img][y+yy][x+xx];
      for (yy=-RCentroid;yy<=RCentroid;yy++) for (xx=-RCentroid;xx<=RCentroid;xx++) if (abs(xx)+abs(yy)<RCentroid*2) {
	 if (ppixOK(img,x+xx,y+yy)) {
	    if (res[img][y+yy][x+xx]>min) {
	       xc1+=xx*(res[img][y+yy][x+xx]-min);
	       yc1+=yy*(res[img][y+yy][x+xx]-min);
	       wt+=(res[img][y+yy][x+xx]-min);
	    }
	    nok++;
	 }
	 else if (posOK(img,x+xx,y+yy) && data[img][y+yy][x+xx]>=iDMAX[img]) {
	    xc1+=xx*(iDMAX[img]-min);
	    yc1+=yy*(iDMAX[img]-min);
	    wt+=(iDMAX[img]-min);
	    if (xx==0 && yy==0) csat=1;
	    else csat=2;
	    nsat++;
	 }
	 else if (ppixOK(img,x-xx,y-yy)) {
	    xc1+=xx*(res[img][y-yy][x-xx]-min);
	    yc1+=yy*(res[img][y-yy][x-xx]-min);
	    wt+=(res[img][y-yy][x-xx]-min);
	 }
      }
      if (nsat>=nok) {xc1=yc1=wt=0.;}
      else if (csat==1) {xc1*=0.05; yc1*=0.05; wt*=0.05;}
      else if (csat) {xc1*=0.005; yc1*=0.005; wt*=0.005;}
      twt+=wt;
      xc+=xc1+wt*(float)(tx-x);
      yc+=yc1+wt*(float)(ty-y);
   }
   if (twt>0) {
      *X=xc/twt+x0+0.5;
      *Y=yc/twt+y0+0.5;
   }
   else {
      *X=x0+0.5;
      *Y=y0+0.5;
   }
   return;
}

int star_flag(double x0,double y0,int IMG) {
   int img,nsat,nbad,ntot,n1,x1,y1,ix,iy,cx,cy,best=15,flag;
   double x,y;

   for (img=0;img<Nimg;img++) if (IMG<0 || img==IMG) {
      flag=0;
      shift(img,x0,y0,&x,&y,1);
      ix=(int)(x+1)-1; iy=(int)(y+1)-1;
      if (ix<rphot[img] || ix>=dataim[img].X-rphot[img] || iy<rphot[img] || iy>=dataim[img].Y-rphot[img]) flag=1;
      if (posOK(img,ix,iy)) {
	 nbad=nsat=ntot=0;
	 for (y1=iy-1;y1<=iy+1;y1++) for (x1=ix-1;x1<=ix+1;x1++) if (posOK(img,x1,y1)) {
	    if (y1==iy && x1==ix) n1=4;
	    else if (y1==iy || x1==ix) n1=2;
	    else n1=1;
	    ntot+=n1;
	    if (data[img][y1][x1]>=iDMAX[img]) nsat+=n1;
	    else if (data[img][y1][x1]<=iDMIN[img]) nbad+=n1;
	 }
	 if (nsat+nbad>=8) flag|=8;
	 if (nsat>=4 && nbad>=4) flag|=6;
	 else if (nsat>=4) flag|=4;
	 else if (nbad+nsat>=4) flag|=2;
	 else {
	    nbad=ntot=0;
	    cx=(int)((x+10-(int)(x+10))*50+0.5);
	    if (cx<0) cx=0; if (cx>50) cx=50;
	    cy=(int)((y+10-(int)(y+10))*50+0.5);
	    if (cy<0) cy=0; if (cy>50) cy=50;
	    for (y1=iy-rphot[img];y1<=iy+rphot[img];y1++) for (x1=ix+circ[cy][cx][img][y1-iy][0];x1<=ix+circ[cy][cx][img][y1-iy][1];x1++) {
	       ntot++;
	       if (!ppixOK(img,x1,y1)) nbad++;
	    }
	    if (nbad*2>=ntot) flag|=8;
	    if (nbad*3>=ntot) flag|=2;
	 }
      }
      else flag|=9;
      if (flag<best) best=flag;
   }
   return best;
}

Inline void eval1_apphot(int img,int ix,int iy,int cx,int cy,float*s,float*ss,float*ssky,float*cs) {
   float c=0,d=0;
   int y1,x1,dx,dy;

   *ssky=getsky(img,ix,iy,NULL);
   *cs=*s=*ss=0;
   for (y1=iy-rphot[img],dy=-rphot[img];y1<=iy+rphot[img];y1++,dy++) for (x1=ix+circ[cy][cx][img][dy][0],dx=circ[cy][cx][img][dy][0];dx<=circ[cy][cx][img][dy][1];x1++,dx++) if (ppixOK(img,x1,y1)) {
      d+=cwt[cy][cx][img][dy][dx];
      c+=psf[dy][dx]*cwt[cy][cx][img][dy][dx];
      if (!NOSUBTRACT) (*s)+=(res[img][y1][x1]-(*ssky))*cwt[cy][cx][img][dy][dx];
      else (*s)+=(data[img][y1][x1]-(*ssky))*cwt[cy][cx][img][dy][dx];
      (*ss)+=noise(img,x1,y1,*ssky)*cwt[cy][cx][img][dy][dx]*cwt[cy][cx][img][dy][dx];
   }
   if (*ss>0 && c>0 && d>0) {
      *cs=c/d;
      (*s)/=c;
      (*ss)=sqrt(*ss)/c;
   }
   else *ss=*s=*cs=0.;
   return;
}

Inline void eval1_psfphot(int img,int ix,int iy,int cx,int cy,float*s,float*ss,float*ssky,float*cs) {
   float c,d,w,n,s0=0.;
#ifdef USESKYSIG
   float dsky;
#endif
   int y1,x1,dx,dy,it;

#ifdef USESKYSIG
   *ssky=getsky(img,ix,iy,&dsky);
#else
   *ssky=getsky(img,ix,iy,NULL);
#endif
   for (it=0;it<PSFPhotIt;it++) {
      if (it) s0=*s;
      *cs=*s=*ss=c=d=0;
      for (y1=iy-rphot[img],dy=-rphot[img];y1<=iy+rphot[img];y1++,dy++) for (x1=ix+circ[cy][cx][img][dy][0],dx=circ[cy][cx][img][dy][0];dx<=circ[cy][cx][img][dy][1];x1++,dx++) if (ppixOK(img,x1,y1)) {
	 if (!it) n=noise(img,x1,y1,*ssky);
	 else n=modelnoise(img,x1,y1,dx,dy,s0,*ssky);
	 if (PSFPhot==2) w=cwt[cy][cx][img][dy][dx]*psf[dy][dx]*psf[dy][dx]/n/n;
	 else w=cwt[cy][cx][img][dy][dx]*psf[dy][dx]/n;
	 (*cs)+=psf[dy][dx]*w;
	 d+=w;
	 if (psf[dy][dx]>0) {
	    if (!NOSUBTRACT) (*s)+=(res[img][y1][x1]-(*ssky))*w;
	    else (*s)+=(data[img][y1][x1]-(*ssky))*w;
	    c+=psf[dy][dx]*w;
	    (*ss)+=w*w*n;
	 }
      }
      if (*ss>0 && c!=0 && d!=0) {
	 (*cs)/=d;
	 (*s)/=c;
#ifdef USESKYSIG
	 (*ss)=sqrt(*ss+d*d*dsky*dsky)/c;
#else
	 (*ss)=sqrt(*ss)/c;
#endif
      }
      else *ss=*s=*cs=0.;
   }
   return;
}

// FitSky = 3 routines
#define INCREMENT_PHOT \
   I+=a;								\
   P+=psf[dy][dx]*a;							\
   D+=d*a;								\
   DP+=d*psf[dy][dx]*a;							\
   PP+=psf[dy][dx]*psf[dy][dx]*a;					\
   a*=a*n;								\
   W+=a;								\
   WP+=psf[dy][dx]*a;							\
   WPP+=psf[dy][dx]*psf[dy][dx]*a;

#define INCREMENT_PHOT_NOD \
   I+=a;								\
   P+=psf[dy][dx]*a;							\
   PP+=psf[dy][dx]*psf[dy][dx]*a;					\
   a*=a*n;								\
   W+=a;								\
   WP+=psf[dy][dx]*a;							\
   WPP+=psf[dy][dx]*psf[dy][dx]*a;

#define INCREMENT_PHOT_SKY \
   I+=a;								\
   P+=psf[dy][dx]*a;							\
   D+=d*a;								\
   DP+=d*psf[dy][dx]*a;							\
   PP+=psf[dy][dx]*psf[dy][dx]*a;

Inline float eval1_sky3(int img,int ix,int iy,int cx,int cy) {
   float I=0,P=0,D=0,DP=0,PP=0,a,skywt=0.,b0=0.,X,d;
   int y1,x1,dx,dy,npix=0;

   if (fsky[img].lastoffset>=0) {
      b0=skyval(img,ix,iy);
      if (b0<0.) a=iRN[img];
      else a=iRN[img]+b0/iGAIN[img];
      skywt=25./a;
      //skywt=M_PI*(SQR(RSky1)-SQR(RSky0))/SQR(SkipSky)/a;
   }
   I = skywt;
   D = b0*skywt;
   for (y1=iy-rphot[img],dy=-rphot[img];y1<=iy+rphot[img];y1++,dy++) for (x1=ix+circ[cy][cx][img][dy][0];x1<=ix+circ[cy][cx][img][dy][1];x1++) if (ppixOK(img,x1,y1)) {
      dx=x1-ix;
      a=cwt[cy][cx][img][dy][dx]/noise(img,x1,y1,skyval(img,x1,y1));
      if (PSFPhot==2) a*=psf[dy][dx]/noise(img,x1,y1,skyval(img,x1,y1));
      d=res[img][y1][x1];
      npix++;
      INCREMENT_PHOT_SKY
   }
   X=I*PP-P*P;
   if (X!=0 && P>0 && npix>1) return (D*PP-P*DP)/X;
   return 0.;
}

Inline float eval1_snr_sky3(int img,int ix,int iy,int cx,int cy,float ct0) {
   float I=0,P=0,PP=0,W=0.,WP=0.,WPP=0.,a,skywt=0.,b0=0,n;
   int y1,x1,dx,dy,npix=0;

   //set to 1/sigma(sky)^2 for 9 pixels;
   if (fsky[img].lastoffset>=0) {
      b0=skyval(img,ix,iy);
      if (b0<0.) a=iRN[img];
      else a=iRN[img]+b0/iGAIN[img];
      skywt=25./a;
      //skywt=M_PI*(SQR(RSky1)-SQR(RSky0))/SQR(SkipSky)/a;
   }
   I = skywt;
   for (y1=iy-rphot[img],dy=-rphot[img];y1<=iy+rphot[img];y1++,dy++) for (x1=ix+circ[cy][cx][img][dy][0],dx=circ[cy][cx][img][dy][0];x1<=ix+circ[cy][cx][img][dy][1];x1++,dx++) if (ppixOK(img,x1,y1) && psf[dy][dx]>0) {
      n=noise(img,x1,y1,skyval(img,x1,y1));
      if (ct0>0.) n+=ct0*psf[dy][dx]/iGAIN[img];
      if (PSFPhot==2) a=cwt[cy][cx][img][dy][dx]*psf[dy][dx]/n/n;
      else a=cwt[cy][cx][img][dy][dx]/n;
      npix++;
      INCREMENT_PHOT_NOD
   }
   a=I*PP-P*P;
   n=I*I*WPP-2*I*P*WP+P*P*W;
   if (a>0 && n>0 && npix>1) return ct0*a/sqrt(n);
   return 0.;
}

Inline void eval1_apphot_sky3(int img,int ix,int iy,int cx,int cy,float*s,float*ss,float*ssky,float*cs) {
   float c=0.,d=0.;
   int y1,x1,dx,dy;

   if (NOSKY34) {
      eval1_apphot(img,ix,iy,cx,cy,s,ss,ssky,cs);
      return;
   }
   *ssky=eval1_sky3(img,ix,iy,cx,cy);
   *cs=*s=*ss=0;
   for (y1=iy-rphot[img],dy=-rphot[img];y1<=iy+rphot[img];y1++,dy++) for (x1=ix+circ[cy][cx][img][dy][0];x1<=ix+circ[cy][cx][img][dy][1];x1++) if (ppixOK(img,x1,y1)) {
      dx=x1-ix;
      d+=cwt[cy][cx][img][dy][dx];
      c+=psf[dy][dx]*cwt[cy][cx][img][dy][dx];
      if (!NOSUBTRACT) (*s)+=(res[img][y1][x1]-(*ssky))*cwt[cy][cx][img][dy][dx];
      else (*s)+=(data[img][y1][x1]-(*ssky))*cwt[cy][cx][img][dy][dx];
      (*ss)+=noise(img,x1,y1,*ssky)*cwt[cy][cx][img][dy][dx]*cwt[cy][cx][img][dy][dx];
   }
   if (*ss>0 && c!=0 && d!=0) {
      *cs=c/d;
      (*s)/=c;
      (*ss)=sqrt(*ss)/c;
   }
   else *ss=*s=*cs=0.;
   return;
}

Inline void eval1_psfphot_sky3(int img,int ix,int iy,int cx,int cy,float*s,float*ss,float*ssky,float*cs) {
   float I,P,D,DP,PP,W,WP,WPP,a,n,skywt=0.,b0=0.,d,s0=0.,ssky0=0.,p0=0.,Ppos;
   int y1,x1,dx,dy,it,npix;

   if (NOSKY34) {
      eval1_psfphot(img,ix,iy,cx,cy,s,ss,ssky,cs);
      return;
   }
   if (fsky[img].lastoffset>=0) {
      b0=skyval(img,ix,iy);
      if (b0<0.) a=iRN[img];
      else a=iRN[img]+b0/iGAIN[img];
      skywt=25./a;
      //skywt=M_PI*(SQR(RSky1)-SQR(RSky0))/SQR(SkipSky)/a;
   }
   for (it=0;it<PSFPhotIt;it++) {
      if (it) {s0=*s; ssky0=*ssky;}
      I=P=D=DP=PP=W=WP=WPP=0.;
      I = skywt;
      D = b0*skywt;
      npix=0;
      for (y1=iy-rphot[img],dy=-rphot[img];y1<=iy+rphot[img];y1++,dy++) for (x1=ix+circ[cy][cx][img][dy][0],dx=circ[cy][cx][img][dy][0];dx<=circ[cy][cx][img][dy][1];x1++,dx++) if (ppixOK(img,x1,y1)) {
	 d=res[img][y1][x1];
	 if (!it) n=noise(img,x1,y1,skyval(img,x1,y1));
	 else n=modelnoise(img,x1,y1,dx,dy,s0,ssky0);
	 if (PSFPhot==2) a=cwt[cy][cx][img][dy][dx]*psf[dy][dx]/n/n;
	 else a=cwt[cy][cx][img][dy][dx]/n;	
	 npix++;
	 INCREMENT_PHOT
      }
      a=I*PP-P*P;
      if (a>0 && P>0 && npix>1) {
	 *ssky=(D*PP-P*DP)/a;
	 if (!NOSUBTRACT) {
	    *s=(I*DP-D*P)/a;
	    *ss=sqrt(I*I*WPP-2*I*P*WP+P*P*W)/a;
	    *cs=PP/P;
	    //if (isnan(*ss)) printf("Caught NaN: npix=%d, %e %e\n",npix,I*I*WPP-2*I*P*WP+P*P*W,a);
	    if (isnan(*ss)) *s=*ss=*cs=*ssky=0;
	 }
	 else {
	    p0=P/I;
	    I=P=D=DP=PP=W=WP=WPP=Ppos=0.;
	    I = skywt;
	    D = b0*skywt;
	    for (y1=iy-rphot[img],dy=-rphot[img];y1<=iy+rphot[img];y1++,dy++) for (x1=ix+circ[cy][cx][img][dy][0],dx=circ[cy][cx][img][dy][0];dx<=circ[cy][cx][img][dy][1];x1++,dx++) if (ppixOK(img,x1,y1) && psf[dy][dx]>=p0) {
	       if (psf[dy][dx]>=p0) d=data[img][y1][x1]-(*ssky);
	       else d=res[img][y1][x1];
	       if (!it) n=noise(img,x1,y1,skyval(img,x1,y1));
	       else n=modelnoise(img,x1,y1,dx,dy,s0,ssky0);
	       if (PSFPhot==2) a=cwt[cy][cx][img][dy][dx]*psf[dy][dx]/n/n;
	       else a=cwt[cy][cx][img][dy][dx]/n;
	       if (psf[dy][dx]>=p0) Ppos+=psf[dy][dx]*a;
	       INCREMENT_PHOT
	    }
	    a=I*PP-P*P;
	    if (a>0 && P>0) {
	       D -= (*ssky)*skywt*Ppos/P;
	       *s=(I*DP-D*P)/a;
	       *ss=sqrt(I*I*WPP-2*I*P*WP+P*P*W)/a;
	       *cs=PP/P;
	    }
	    else *s=*ss=*cs=0.;
	 }
      }
      else *ss=*s=*cs=*ssky=0.;
   }
#ifdef NAN_PRINT
   if (isnan(*s)) printf("Uncaught s=nan\n");
   if (isnan(*ss)) printf("Uncaught ss=nan\n");
   if (isnan(*cs)) printf("Uncaught cs=nan\n");
   if (isnan(*ssky)) printf("Uncaught ssky=nan\n");
#endif
   if (isnan(*s) || isnan(*ss) || isnan(*cs) || isnan(*ssky)) *ss=*s=*cs=*ssky=0.;
   return;
}
#undef INCREMENT_PHOT
#undef INCREMENT_PHOT_NOD
#undef INCREMENT_PHOT_SKY

#define INCREMENT_PHOT \
   I+=a;								\
   P+=psf[dy][dx]*a;							\
   D+=d*a;								\
   DP+=d*psf[dy][dx]*a;							\
   DX+=d*dx*a;								\
   DY+=d*dy*a;								\
   PP+=psf[dy][dx]*psf[dy][dx]*a;					\
   X+=dx*a;								\
   Y+=dy*a;								\
   XX+=dx*dx*a;								\
   XY+=dx*dy*a;								\
   YY+=dy*dy*a;								\
   XP+=dx*psf[dy][dx]*a;						\
   YP+=dy*psf[dy][dx]*a;						\
   a*=a*n;								\
   W+=a;								\
   WP+=psf[dy][dx]*a;							\
   WPP+=psf[dy][dx]*psf[dy][dx]*a;					\
   WXP+=dx*psf[dy][dx]*a;						\
   WYP+=dy*psf[dy][dx]*a;						\
   WX+=dx*a;								\
   WY+=dy*a;								\
   WXX+=dx*dx*a;							\
   WXY+=dx*dy*a;							\
   WYY+=dy*dy*a;

#define INCREMENT_PHOT_SKY \
   I+=a;								\
   P+=psf[dy][dx]*a;							\
   D+=d*a;								\
   DP+=d*psf[dy][dx]*a;							\
   DX+=d*dx*a;								\
   DY+=d*dy*a;								\
   PP+=psf[dy][dx]*psf[dy][dx]*a;					\
   X+=dx*a;								\
   Y+=dy*a;								\
   XX+=dx*dx*a;								\
   XY+=dx*dy*a;								\
   YY+=dy*dy*a;								\
   XP+=dx*psf[dy][dx]*a;						\
   YP+=dy*psf[dy][dx]*a;						\

#define INCREMENT_PHOT_NOD \
   I+=a;								\
   P+=psf[dy][dx]*a;							\
   PP+=psf[dy][dx]*psf[dy][dx]*a;					\
   X+=dx*a;								\
   Y+=dy*a;								\
   XX+=dx*dx*a;								\
   XY+=dx*dy*a;								\
   YY+=dy*dy*a;								\
   XP+=dx*psf[dy][dx]*a;						\
   YP+=dy*psf[dy][dx]*a;						\
   a*=a*n;								\
   W+=a;								\
   WP+=psf[dy][dx]*a;							\
   WPP+=psf[dy][dx]*psf[dy][dx]*a;					\
   WXP+=dx*psf[dy][dx]*a;						\
   WYP+=dy*psf[dy][dx]*a;						\
   WX+=dx*a;								\
   WY+=dy*a;								\
   WXX+=dx*dx*a;							\
   WXY+=dx*dy*a;							\
   WYY+=dy*dy*a;

Inline float eval1_sky4(int img,int ix,int iy,int cx,int cy,float*mx,float*my) {
   float I=0,P=0,D=0,DP=0,PP=0,X=0,Y=0,XX=0,XY=0,YY=0,XP=0,YP=0,DX=0,DY=0,a,skywt=0,b0=0,d;
   int y1,x1,dx,dy,npix=0;

   if (fsky[img].lastoffset>=0) {
      b0=skyval(img,ix,iy);
      if (b0<0.) a=iRN[img];
      else a=iRN[img]+b0/iGAIN[img];
      skywt=25./a;
      //skywt=M_PI*(SQR(RSky1)-SQR(RSky0))/SQR(SkipSky)/a;
   }
   D = b0*skywt;
   for (y1=iy-rphot[img],dy=-rphot[img];y1<=iy+rphot[img];y1++,dy++) for (x1=ix+circ[cy][cx][img][dy][0];x1<=ix+circ[cy][cx][img][dy][1];x1++) if (ppixOK(img,x1,y1)) {
      dx=x1-ix;
      a=cwt[cy][cx][img][dy][dx]/noise(img,x1,y1,skyval(img,x1,y1));
      d=res[img][y1][x1];
      npix++;
      INCREMENT_PHOT_SKY
   }
   d = (XX*YY*P*P - P*P*XY*XY - 2*YY*P*X*XP + 2*P*X*XY*YP + 2*P*XP*XY*Y - 2*XX*P*Y*YP - X*X*YP*YP + PP*YY*X*X + 2*X*XP*Y*YP - 2*PP*X*XY*Y - XP*XP*Y*Y + I*YY*XP*XP - 2*I*XP*XY*YP + I*PP*XY*XY + PP*XX*Y*Y + I*XX*YP*YP - I*PP*XX*YY);
   if (d!=0 && P>0 && npix>3) {
      if (mx) *mx = (DX*(YY*P*P - 2*P*Y*YP + PP*Y*Y + I*YP*YP - I*PP*YY)-DY*(P*P*XY - I*PP*XY + I*XP*YP - P*X*YP - P*XP*Y + PP*X*Y)-DP*(XP*Y*Y - I*XP*YY + I*XY*YP + P*X*YY - P*XY*Y - X*Y*YP)-D*(X*YP*YP + P*XP*YY - P*XY*YP - PP*X*YY + PP*XY*Y - XP*Y*YP))/d;
      if (my) *my = (DY*(XX*P*P - 2*P*X*XP + PP*X*X + I*XP*XP - I*PP*XX)-DX*(P*P*XY - I*PP*XY + I*XP*YP - P*X*YP - P*XP*Y + PP*X*Y)-DP*(X*X*YP + I*XP*XY - I*XX*YP - P*X*XY + P*XX*Y - X*XP*Y)-D*(XP*XP*Y - P*XP*XY + PP*X*XY + P*XX*YP - PP*XX*Y - X*XP*YP))/d;
      return (D*(YY*XP*XP - 2*XP*XY*YP + PP*XY*XY + XX*YP*YP - PP*XX*YY)-DY*(XP*XP*Y - P*XP*XY + PP*X*XY + P*XX*YP - PP*XX*Y - X*XP*YP)-DX*(X*YP*YP + P*XP*YY - P*XY*YP - PP*X*YY + PP*XY*Y - XP*Y*YP)-DP*(P*XY*XY - P*XX*YY + X*XP*YY - X*XY*YP - XP*XY*Y + XX*Y*YP))/d;
   }
   *mx=*my=0;
   return 0.;
}

Inline float eval1_snr_sky4(int img,int ix,int iy,int cx,int cy,float ct0) {
   float I=0,P=0,PP=0,X=0,Y=0,XX=0,XY=0,YY=0,XP=0,YP=0,W=0,WP=0,WPP=0,WXP=0,WYP=0,WX=0,WY=0,WXX=0,WXY=0,WYY=0;
   float a,skywt=0.,b0=0,n,n1,n2,n3,n4,d;
   int y1,x1,dx,dy,npix=0;

   //set to 1/sigma(sky)^2 for 9 pixels;
   if (fsky[img].lastoffset>=0) {
      b0=skyval(img,ix,iy);
      if (b0<0.) a=iRN[img];
      else a=iRN[img]+b0/iGAIN[img];
      skywt=25./a;
      //skywt=M_PI*(SQR(RSky1)-SQR(RSky0))/SQR(SkipSky)/a;
   }
   I = skywt;
   for (y1=iy-rphot[img],dy=-rphot[img];y1<=iy+rphot[img];y1++,dy++) for (x1=ix+circ[cy][cx][img][dy][0],dx=circ[cy][cx][img][dy][0];x1<=ix+circ[cy][cx][img][dy][1];x1++,dx++) if (ppixOK(img,x1,y1) && psf[dy][dx]>0) {
      n=noise(img,x1,y1,skyval(img,x1,y1));
      if (ct0>0.) n+=ct0*psf[dy][dx]/iGAIN[img];
      if (PSFPhot==2) a=cwt[cy][cx][img][dy][dx]*psf[dy][dx]/n/n;
      else a=cwt[cy][cx][img][dy][dx]/n;
      npix++;
      INCREMENT_PHOT_NOD
   }
   d = (XX*YY*P*P - P*P*XY*XY - 2*YY*P*X*XP + 2*P*X*XY*YP + 2*P*XP*XY*Y - 2*XX*P*Y*YP - X*X*YP*YP + PP*YY*X*X + 2*X*XP*Y*YP - 2*PP*X*XY*Y - XP*XP*Y*Y + I*YY*XP*XP - 2*I*XP*XY*YP + I*PP*XY*XY + PP*XX*Y*Y + I*XX*YP*YP - I*PP*XX*YY);
   n1 = YY*X*X - 2*X*XY*Y + I*XY*XY + XX*Y*Y - I*XX*YY;
   n2 = X*X*YP + I*XP*XY - I*XX*YP - P*X*XY + P*XX*Y - X*XP*Y;
   n3 = XP*Y*Y - I*XP*YY + I*XY*YP + P*X*YY - P*XY*Y - X*Y*YP;
   n4 = P*XY*XY - P*XX*YY + X*XP*YY - X*XY*YP - XP*XY*Y + XX*Y*YP;
   n = WPP*n1*n1-WYP*2*n1*n2-WXP*2*n1*n3-WP*2*n1*n4+WYY*n2*n2+WXY*2*n2*n3+WY*2*n2*n4+WXX*n3*n3+WX*2*n3*n4+W*n4*n4;
   if (d!=0 && n>0 && npix>3) return ct0*fabs(d)/sqrt(n);
   return 0.;
}

Inline void eval1_apphot_sky4(int img,int ix,int iy,int cx,int cy,float*s,float*ss,float*ssky,float*mx,float*my,float*cs) {
   float c=0.,d=0.;
   int y1,x1,dx,dy;

   if (NOSKY34) {
      eval1_apphot(img,ix,iy,cx,cy,s,ss,ssky,cs);
      return;
   }
   *ssky=eval1_sky4(img,ix,iy,cx,cy,mx,my);
   *cs=*s=*ss=0;
   for (y1=iy-rphot[img],dy=-rphot[img];y1<=iy+rphot[img];y1++,dy++) for (x1=ix+circ[cy][cx][img][dy][0];x1<=ix+circ[cy][cx][img][dy][1];x1++) if (ppixOK(img,x1,y1)) {
      dx=x1-ix;
      d+=cwt[cy][cx][img][dy][dx];
      c+=psf[dy][dx]*cwt[cy][cx][img][dy][dx];
      if (!NOSUBTRACT) (*s)+=(res[img][y1][x1]-(*ssky)-(*mx)*dx-(*my)*dy)*cwt[cy][cx][img][dy][dx];
      else (*s)+=(data[img][y1][x1]-(*ssky)-(*ssky)-(*mx)*dx-(*my)*dy)*cwt[cy][cx][img][dy][dx];
      (*ss)+=noise(img,x1,y1,*ssky+(*mx)*dx+(*my)*dy)*cwt[cy][cx][img][dy][dx]*cwt[cy][cx][img][dy][dx];
   }
   if (*ss>0 && c!=0 && d!=0) {
      *cs=c/d;
      (*s)/=c;
      (*ss)=sqrt(*ss)/c;
   }
   else *ss=*s=*cs=0.;
   return;
}

Inline void eval1_psfphot_sky4(int img,int ix,int iy,int cx,int cy,float*s,float*ss,float*ssky,float*mx,float*my,float*cs) {
   float I,P,D,DP,DX,DY,PP,X,Y,XX,XY,YY,XP,YP,W,WP,WPP,WXP,WYP,WX,WY,WXX,WXY,WYY;
   float a,n,skywt=0.,b0=0.,d,s0=0.,ssky0=0.,mx0=0,my0=0,p0=0.,Ppos,n1,n2,n3,n4;
   int y1,x1,dx,dy,it,npix;

   *mx=*my=0.0;
   if (NOSKY34) {
      eval1_psfphot(img,ix,iy,cx,cy,s,ss,ssky,cs);
      return;
   }
   if (fsky[img].lastoffset>=0) {
      b0=skyval(img,ix,iy);
      if (b0<0.) a=iRN[img];
      else a=iRN[img]+b0/iGAIN[img];
      skywt=25./a;
      //skywt=M_PI*(SQR(RSky1)-SQR(RSky0))/SQR(SkipSky)/a;
   }
   for (it=0;it<PSFPhotIt;it++) {
      if (it) {s0=*s; ssky0=*ssky; mx0=*mx; my0=*my;}
      I=P=D=DP=DX=DY=PP=X=Y=XX=XY=YY=XP=YP=W=WP=WPP=WXP=WYP=WX=WY=WXX=WXY=WYY=0;
      npix=0;
      I = skywt;
      D = b0*skywt;
      for (y1=iy-rphot[img],dy=-rphot[img];y1<=iy+rphot[img];y1++,dy++) for (x1=ix+circ[cy][cx][img][dy][0],dx=circ[cy][cx][img][dy][0];dx<=circ[cy][cx][img][dy][1];x1++,dx++) if (ppixOK(img,x1,y1)) {
	 d=res[img][y1][x1];
	 if (!it) n=noise(img,x1,y1,skyval(img,x1,y1));
	 else n=modelnoise(img,x1,y1,dx,dy,s0,ssky0+mx0*dx+my0*dy);
	 if (PSFPhot==2) a=cwt[cy][cx][img][dy][dx]*psf[dy][dx]/n/n;
	 else a=cwt[cy][cx][img][dy][dx]/n;
	 npix++;
	 INCREMENT_PHOT
      }
      d = (XX*YY*P*P - P*P*XY*XY - 2*YY*P*X*XP + 2*P*X*XY*YP + 2*P*XP*XY*Y - 2*XX*P*Y*YP - X*X*YP*YP + PP*YY*X*X + 2*X*XP*Y*YP - 2*PP*X*XY*Y - XP*XP*Y*Y + I*YY*XP*XP - 2*I*XP*XY*YP + I*PP*XY*XY + PP*XX*Y*Y + I*XX*YP*YP - I*PP*XX*YY);
      if (d!=0 && P>0 && npix>3) {
	 *ssky = (D*(YY*XP*XP - 2*XP*XY*YP + PP*XY*XY + XX*YP*YP - PP*XX*YY)-DY*(XP*XP*Y - P*XP*XY + PP*X*XY + P*XX*YP - PP*XX*Y - X*XP*YP)-DX*(X*YP*YP + P*XP*YY - P*XY*YP - PP*X*YY + PP*XY*Y - XP*Y*YP)-DP*(P*XY*XY - P*XX*YY + X*XP*YY - X*XY*YP - XP*XY*Y + XX*Y*YP))/d;
	 *mx = (DX*(YY*P*P - 2*P*Y*YP + PP*Y*Y + I*YP*YP - I*PP*YY)-DY*(P*P*XY - I*PP*XY + I*XP*YP - P*X*YP - P*XP*Y + PP*X*Y)-DP*(XP*Y*Y - I*XP*YY + I*XY*YP + P*X*YY - P*XY*Y - X*Y*YP)-D*(X*YP*YP + P*XP*YY - P*XY*YP - PP*X*YY + PP*XY*Y - XP*Y*YP))/d;
	 *my = (DY*(XX*P*P - 2*P*X*XP + PP*X*X + I*XP*XP - I*PP*XX)-DX*(P*P*XY - I*PP*XY + I*XP*YP - P*X*YP - P*XP*Y + PP*X*Y)-DP*(X*X*YP + I*XP*XY - I*XX*YP - P*X*XY + P*XX*Y - X*XP*Y)-D*(XP*XP*Y - P*XP*XY + PP*X*XY + P*XX*YP - PP*XX*Y - X*XP*YP))/d;
	 if (!NOSUBTRACT) {
	    /*
	      Quick note: math assuming that X=Y=XY=XP=YP=0
	      d = (P*P-I*PP);
	      *s = (D*P-DP*I)/d;
	      *ssky = (DP*P-D*PP)/d;
	     */
	    n1 = YY*X*X - 2*X*XY*Y + I*XY*XY + XX*Y*Y - I*XX*YY;
	    n2 = X*X*YP + I*XP*XY - I*XX*YP - P*X*XY + P*XX*Y - X*XP*Y;
	    n3 = XP*Y*Y - I*XP*YY + I*XY*YP + P*X*YY - P*XY*Y - X*Y*YP;
	    n4 = P*XY*XY - P*XX*YY + X*XP*YY - X*XY*YP - XP*XY*Y + XX*Y*YP;
	    *s = (DP*n1-DY*n2-DX*n3-D*n4)/d;
	    *ss = sqrt(WPP*n1*n1+WYY*n2*n2+WXX*n3*n3+W*n4*n4-2*WYP*n1*n2-2*WXP*n1*n3-2*WP*n1*n4+2*WXY*n2*n3+2*WY*n2*n4+2*WX*n3*n4)/fabs(d);
	    *cs=PP/P;
	    //if (isnan(*ss)) printf("Caught NaN: npix=%d, %e %e\n",npix,WPP*n1*n1+WYY*n2*n2+WXX*n3*n3+W*n4*n4-2*WYP*n1*n2-2*WXP*n1*n3-2*WP*n1*n4+2*WXY*n2*n3+2*WY*n2*n4+2*WX*n3*n4,d);
	    if (isnan(*ss)) *s=*ss=*cs=*ssky=0;
	 }
	 else {
	    p0=P/(I+skywt);
	    I=P=D=DP=DX=DY=PP=X=Y=XX=XY=YY=XP=YP=W=WP=WPP=WXP=WYP=WX=WY=WXX=WXY=WYY=Ppos=0;
	    I = skywt;
	    D = b0*skywt;
	    for (y1=iy-rphot[img],dy=-rphot[img];y1<=iy+rphot[img];y1++,dy++) for (x1=ix+circ[cy][cx][img][dy][0],dx=circ[cy][cx][img][dy][0];dx<=circ[cy][cx][img][dy][1];x1++,dx++) if (ppixOK(img,x1,y1) && psf[dy][dx]>=p0) {
	       if (psf[dy][dx]>=p0) d=data[img][y1][x1]-(*ssky)-(*mx)*dx-(*my)*dy;
	       else d=res[img][y1][x1];
	       if (!it) n=noise(img,x1,y1,skyval(img,x1,y1));
	       else n=modelnoise(img,x1,y1,dx,dy,s0,ssky0+mx0*dx+my0*dy);
	       if (PSFPhot==2) a=cwt[cy][cx][img][dy][dx]*psf[dy][dx]/n/n;
	       else a=cwt[cy][cx][img][dy][dx]/n;
	       if (psf[dy][dx]>=p0) Ppos+=psf[dy][dx]*a;
	       INCREMENT_PHOT
	    }
	    d = (XX*YY*P*P - P*P*XY*XY - 2*YY*P*X*XP + 2*P*X*XY*YP + 2*P*XP*XY*Y - 2*XX*P*Y*YP - X*X*YP*YP + PP*YY*X*X + 2*X*XP*Y*YP - 2*PP*X*XY*Y - XP*XP*Y*Y + I*YY*XP*XP - 2*I*XP*XY*YP + I*PP*XY*XY + PP*XX*Y*Y + I*XX*YP*YP - I*PP*XX*YY);
	    if (d!=0 && P>0) {
	       D -= (*ssky)*skywt*Ppos/P;
	       *s = (DP*(YY*X*X - 2*X*Y*XY + XX*Y*Y + I*XY*XY - I*XX*YY) - DX*(XP*Y*Y - YP*X*Y + P*X*YY - P*Y*XY - XP*I*YY + YP*I*XY) - D*(P*XY*XY + XP*X*YY - XP*Y*XY - YP*X*XY + YP*Y*XX - P*XX*YY) - DY*(YP*X*X - XP*X*Y - P*X*XY + P*Y*XX + XP*I*XY - YP*I*XX))/d;
	       n1 = YY*X*X - 2*X*XY*Y + I*XY*XY + XX*Y*Y - I*XX*YY;
	       n2 = X*X*YP + I*XP*XY - I*XX*YP - P*X*XY + P*XX*Y - X*XP*Y;
	       n3 = XP*Y*Y - I*XP*YY + I*XY*YP + P*X*YY - P*XY*Y - X*Y*YP;
	       n4 = P*XY*XY - P*XX*YY + X*XP*YY - X*XY*YP - XP*XY*Y + XX*Y*YP;
	       *ss = sqrt(WPP*n1*n1-WYP*2*n1*n2-WXP*2*n1*n3-WP*2*n1*n4+WYY*n2*n2+WXY*2*n2*n3+WY*2*n2*n4+WXX*n3*n3+WX*2*n3*n4+W*n4*n4)/fabs(d);
	       *cs=PP/P;
	       if (isnan(*ss)) *s=*ss=*cs=*ssky=0;
	    }
	    else *s=*ss=*cs=0.;
	 }
      }
      else *ss=*s=*cs=*ssky=0.;
   }
   return;
}
#undef INCREMENT_PHOT
#undef INCREMENT_PHOT_NOD
#undef INCREMENT_PHOT_SKY

Inline float eval1_snr(int img,int ix,int iy,int cx,int cy,float ct0) {
   float c=0,w,n,ss=0;
   int y1,x1,dx,dy;

   if (FitSky==3 && PSFPhot) return eval1_snr_sky3(img,ix,iy,cx,cy,ct0);
   if (FitSky==4 && PSFPhot) return eval1_snr_sky4(img,ix,iy,cx,cy,ct0);
   for (y1=iy-rphot[img],dy=-rphot[img];y1<=iy+rphot[img];y1++,dy++) for (x1=ix+circ[cy][cx][img][dy][0],dx=circ[cy][cx][img][dy][0];dx<=circ[cy][cx][img][dy][1];x1++,dx++) if (ppixOK(img,x1,y1) && psf[dy][dx]>0) {
      n=noise(img,x1,y1,skyval(img,x1,y1));
      if (ct0>0.) n+=ct0*psf[dy][dx]/iGAIN[img];
      if (!PSFPhot) w=cwt[cy][cx][img][dy][dx];
      else if (PSFPhot==2) w=cwt[cy][cx][img][dy][dx]*psf[dy][dx]*psf[dy][dx]/n/n;
      else w=cwt[cy][cx][img][dy][dx]*psf[dy][dx]/n;
      c+=psf[dy][dx]*w;
      ss+=w*w*n;
   }
   if (c>0. && ss>0.) return ct0*c/sqrt(ss);
   return 0.;
}

#if 1
void eval1_phot(int img,int ix,int iy,int cx,int cy,int big,float*cm,float*s,float*ss,float*ssky,float*mx,float*my,float*cs)
{
   NOSKY34=big;
   *cm=apcor[img]/iEXP[img];
   *mx=*my=0.0;
   if (FitSky==3) {
      if (PSFPhot) eval1_psfphot_sky3(img,ix,iy,cx,cy,s,ss,ssky,cs);
      else eval1_apphot_sky3(img,ix,iy,cx,cy,s,ss,ssky,cs);
   }
   else if (FitSky==4) {
      if (PSFPhot) eval1_psfphot_sky4(img,ix,iy,cx,cy,s,ss,ssky,mx,my,cs);
      else eval1_apphot_sky4(img,ix,iy,cx,cy,s,ss,ssky,mx,my,cs);
   }
   else if (PSFPhot) eval1_psfphot(img,ix,iy,cx,cy,s,ss,ssky,cs);
   else eval1_apphot(img,ix,iy,cx,cy,s,ss,ssky,cs);
   (*s)*=(*cm);
   (*ss)*=(*cm);
   (*ssky)*=(*cm);
   NOSKY34=0;
}

float eval1_fom(int img,int ix,int iy,int cx,int cy,int star,float s,float ss,float ssky,float mx,float my,float cs,float*c,float*ndof,float*sh)
{
   int x1,y1,dx,dy,ct=0;
   float a,a0,b,A=0,np=0,ns,m=0,w,wchi,wmax[4],wmaxtot=0,tdn;

   wmax[0]=wmax[1]=wmax[2]=wmax[3]=0;
   for (y1=iy-rphot[img],dy=-rphot[img];y1<=iy+rphot[img];y1++,dy++) for (x1=ix+circ[cy][cx][img][dy][0];x1<=ix+circ[cy][cx][img][dy][1];x1++) if (ppixOK(img,x1,y1)) {
      dx=x1-ix;
      ct++;
      if (PSFPhot && PSFPhotIt>1) ns=modelnoise(img,x1,y1,dx,dy,s,ssky+mx*dx+my*dy);
      else ns=noise(img,x1,y1,ssky+mx*dx+my*dy);
      a0=psf[dy][dx]*s;
      if (IMDIFF && star>=0) tdn=a0+psf[dy][dx]*refcts[star]/refmult[img];
      else tdn=a0;
      if (tdn>0.) b=1./(ns+NoiseMult*tdn*tdn);
      else b=1./ns;
      if (!NOSUBTRACT) a=res[img][y1][x1];
      else a=data[img][y1][x1];
      a-=ssky+mx*dx+my*dy+a0;
      if (PSFPhot==2) {
	 if (psf[dy][dx]>0.) {
	    w=cwt[cy][cx][img][dy][dx]*psf[dy][dx]/ns;
	    wchi=chiwt[cy][cx][img][dy][dx]*psf[dy][dx]/ns;
	 }
	 else w=wchi=0.;
      }
      else {
	 w=cwt[cy][cx][img][dy][dx];
	 wchi=chiwt[cy][cx][img][dy][dx];
      }
      (*c)+=a*a*b*wchi;
      if (w>wmax[0]) {wmax[3]=wmax[2]; wmax[2]=wmax[1]; wmax[1]=wmax[0]; wmax[0]=wchi;}
      else if (w>wmax[1]) {wmax[3]=wmax[2]; wmax[2]=wmax[1]; wmax[1]=wchi;}
      else if (w>wmax[2]) {wmax[3]=wmax[2]; wmax[2]=wchi;}
      else if (w>wmax[3]) wmax[3]=wchi;
      np+=wchi;
      (*sh)+=a/sqrt(ns)*(psf[dy][dx]-cs)*w;
      if (IMDIFF && star>=0) A+=(psf[dy][dx]-cs)*(psf[dy][dx]-cs)*(s+refcts[star]/refmult[img])/sqrt(ns)*w;
      else A+=(psf[dy][dx]-cs)*(psf[dy][dx]-cs)*s/sqrt(ns)*w;
   }
   if (A>0 && ct>1) (*sh)/=A;
   else (*sh)=-9.999;
   if (*sh>9.999) *sh=9.999;
   else if (*sh<-9.999) *sh=-9.999;
   if (FitSky==4) wmaxtot=wmax[0]+wmax[1]+wmax[2]+wmax[3];
   else if (FitSky==3) wmaxtot=wmax[0]+wmax[1];
   else wmaxtot=wmax[0];
   if (ct>1 && np>wmaxtot) {
      *ndof=np-wmaxtot;
      (*c)=(*c)/(np-wmaxtot);
   }
   else {
      *c=1.;
      *ndof=0;
   }
   m=1/sqrt((*c)+0.25); // for chi minimization
   if (SearchMode==0) {
      m*=s/ss; // convert to SNR/chi maximization
   }
   if (m<0.0001) m=0.0001;
   *c = sqrt(*c);
#ifdef NAN_PRINT
   if (isnan(m)) printf("Uncaught m=nan\n");
#endif
#ifdef NAN_CRASH
   assert(!isnan(m));
#else
   if (isnan(m)) return 0;
#endif
   return m;
}

float eval1(int img,double x,double y,float pa,float pb,float pc,float*s,float*ss,float*c,float*ndof,float*sh,float*ssky,float*mx,float*my,float*cs,float*cm,int phot,int full,int star,int newpsf)
{
   int ix,iy,cx,cy,big=0;

   shift(img,x,y,&x,&y,1);
   ix=(int)(x+100)-100;
   iy=(int)(y+100)-100;
   if (ix<-rphot[img] || ix>=dataim[img].X+rphot[img] || iy<-rphot[img] || iy>=dataim[img].Y+rphot[img]) {
      *s=*ss=*c=*sh=*ssky=*ndof=0;
      return 0;
   }
   cx=(int)((x-ix)*50+0.5);
   if (cx<0) cx=0; if (cx>50) cx=50;
   cy=(int)((y-iy)*50+0.5);
   if (cy<0) cy=0; if (cy>50) cy=50;
   if (newpsf) calc1psf(img,x,y,rphot[img],pa,pb,pc,1,1);

   if (phot) {
      *mx=*my=0.0;
      if (RBig>0 && pa==RBig && pb==RBig && pc==0.) big=1;
      eval1_phot(img,ix,iy,cx,cy,big,cm,s,ss,ssky,mx,my,cs);
   }

   *c=*sh=*ndof=0;
   if (!(*ss>0)) return 0;
   if (full) return eval1_fom(img,ix,iy,cx,cy,star,(*s)/(*cm),(*ss)/(*cm),(*ssky)/(*cm),*mx,*my,*cs,c,ndof,sh);
   return ((*s)/(*cm))/((*ss)/(*cm));
}

#else
float eval1(int img,double x,double y,float pa,float pb,float pc,float*s,float*ss,float*c,float*ndof,float*sh,float*ssky,float*cm,int phot,int full,int star,int newpsf) {
   int ix,iy,cx,cy,x1,y1,dx,dy,ct=0;
   float cs,a,a0,b,A,np,ns,m=0,w,wchi,wmax[4],wmaxtot=0,tdn,mx=0,my=0;

   *cm=apcor[img]/iEXP[img];
   shift(img,x,y,&x,&y,1);
   ix=(int)(x+100)-100;
   iy=(int)(y+100)-100;
   if (ix<-rphot[img] || ix>=dataim[img].X+rphot[img] || iy<-rphot[img] || iy>=dataim[img].Y+rphot[img]) {
      *s=*ss=*c=*sh=*ssky=*ndof=0;
      return 0;
   }
   cx=(int)((x-ix)*50+0.5);
   if (cx<0) cx=0; if (cx>50) cx=50;
   cy=(int)((y-iy)*50+0.5);
   if (cy<0) cy=0; if (cy>50) cy=50;
   if (newpsf) calc1psf(img,x,y,rphot[img],pa,pb,pc,1,1);
   if (RBig>0 && pa==RBig && pb==RBig && pc==0.) NOSKY34=1;
   if (FitSky==3) {
      if (PSFPhot) eval1_psfphot_sky3(img,ix,iy,cx,cy,s,ss,ssky,&cs);
      else eval1_apphot_sky3(img,ix,iy,cx,cy,s,ss,ssky,&cs);
   }
   else if (FitSky==4) {
      if (PSFPhot) eval1_psfphot_sky4(img,ix,iy,cx,cy,s,ss,ssky,&mx,&my,&cs);
      else eval1_apphot_sky4(img,ix,iy,cx,cy,s,ss,ssky,&mx,&my,&cs);
   }
   else if (PSFPhot) eval1_psfphot(img,ix,iy,cx,cy,s,ss,ssky,&cs);
   else eval1_apphot(img,ix,iy,cx,cy,s,ss,ssky,&cs);
   NOSKY34=0;
   *c=*sh=A=np=*ndof=0;
   if (full && *ss>0) {
      ct=0;
      wmax[0]=wmax[1]=wmax[2]=wmax[3]=0;
      for (y1=iy-rphot[img],dy=-rphot[img];y1<=iy+rphot[img];y1++,dy++) for (x1=ix+circ[cy][cx][img][dy][0];x1<=ix+circ[cy][cx][img][dy][1];x1++) if (ppixOK(img,x1,y1)) {
	 dx=x1-ix;
	 ct++;
	 if (PSFPhot && PSFPhotIt>1) ns=modelnoise(img,x1,y1,dx,dy,*s,*ssky+mx*dx+my*dy);
	 else ns=noise(img,x1,y1,(*ssky)+mx*dx+my*dy);
	 a0=psf[dy][dx]*(*s);
	 if (IMDIFF && star>=0) tdn=a0+psf[dy][dx]*refcts[star]/refmult[img];
	 else tdn=a0;
	 if (tdn>0.) b=1./(ns+NoiseMult*tdn*tdn);
	 else b=1./ns;
	 if (!NOSUBTRACT) a=res[img][y1][x1];
	 else a=data[img][y1][x1];
	 a-=(*ssky)+mx*dx+my*dy+a0;
	 if (PSFPhot==2) {
	    if (psf[dy][dx]>0.) {
	       w=cwt[cy][cx][img][dy][dx]*psf[dy][dx]/ns;
	       wchi=chiwt[cy][cx][img][dy][dx]*psf[dy][dx]/ns;
	    }
	    else w=wchi=0.;
	 }
	 else {
	    w=cwt[cy][cx][img][dy][dx];
	    wchi=chiwt[cy][cx][img][dy][dx];
	 }
	 (*c)+=a*a*b*wchi;
	 if (w>wmax[0]) {wmax[3]=wmax[2]; wmax[2]=wmax[1]; wmax[1]=wmax[0]; wmax[0]=wchi;}
	 else if (w>wmax[1]) {wmax[3]=wmax[2]; wmax[2]=wmax[1]; wmax[1]=wchi;}
	 else if (w>wmax[2]) {wmax[3]=wmax[2]; wmax[2]=wchi;}
	 else if (w>wmax[3]) wmax[3]=wchi;
	 np+=wchi;
	 (*sh)+=a/sqrt(ns)*(psf[dy][dx]-cs)*w;
	 if (IMDIFF && star>=0) A+=(psf[dy][dx]-cs)*(psf[dy][dx]-cs)*(*s+refcts[star]/refmult[img])/sqrt(ns)*w;
	 else A+=(psf[dy][dx]-cs)*(psf[dy][dx]-cs)*(*s)/sqrt(ns)*w;
      }
      if (A>0 && ct>1) (*sh)/=A;
      else (*sh)=-9.999;
      if (*sh>9.999) *sh=9.999;
      else if (*sh<-9.999) *sh=-9.999;
      if (FitSky==4) wmaxtot=wmax[0]+wmax[1]+wmax[2]+wmax[3];
      else if (FitSky==3) wmaxtot=wmax[0]+wmax[1];
      else wmaxtot=wmax[0];
      if (ct>1 && np>wmaxtot) {
	 *ndof=np-wmaxtot;
	 (*c)=(*c)/(np-wmaxtot);
      }
      else {
	 *c=1.;
	 *ndof=0;
      }
      m=1/sqrt((*c)+0.25); // for chi minimization
      if (SearchMode==0) {
	 m*=(*s)/(*ss); // convert to SNR/chi maximization
      }
      if (m<0.0001) m=0.0001;
      *c = sqrt(*c);
   }
   else if (*ss>0) m=(*s)/(*ss);
   (*s)*=(*cm);
   (*ss)*=(*cm);
   (*ssky)*=(*cm);
#ifdef NAN_PRINT
   if (isnan(m)) printf("Uncaught m=nan\n");
#endif
#ifdef NAN_CRASH
   assert(!isnan(m));
#else
   if (isnan(m)) *s=*ss=*c=*sh=*ssky=*ndof=m=0;
#endif
   return m;
}
#endif

/* AEDDEBUG
  To-do:
  - bad photometry points (saturated, etc. that wouldn't be combined in out1star) shouldn't be included in eval outputs, especially combined forced photometry
 */
float eval(int IMG,double x,double y,float*PA,float*PB,float*PC,float*S,float*S0,float*SS,float*C,float*SH,float*SSKY,float*IS,float*IS0,float*ISS,float*IC,float*ISH,float*ISSKY,float*ICM,int checkstar,int star) {
   int img;
   float M=0,wt,swt,twt=0,tswt=0,tndof=0;
   static float *m=0,*s0,*s,*ss,*c,*ndof,*sh,*ssky,*mx,*my,*cs,*cm;
   static int *comb;

   if (m==0) {
      m=(float*)calloc(Timg,FLOATSIZE);
      s0=(float*)calloc(Timg,FLOATSIZE);
      s=(float*)calloc(Timg,FLOATSIZE);
      ss=(float*)calloc(Timg,FLOATSIZE);
      c=(float*)calloc(Timg,FLOATSIZE);
      ndof=(float*)calloc(Timg,FLOATSIZE);
      sh=(float*)calloc(Timg,FLOATSIZE);
      ssky=(float*)calloc(Timg,FLOATSIZE);
      mx=(float*)calloc(Timg,FLOATSIZE);
      my=(float*)calloc(Timg,FLOATSIZE);
      cs=(float*)calloc(Timg,FLOATSIZE);
      cm=(float*)calloc(Timg,FLOATSIZE);
      comb=(int*)calloc(Timg,INTSIZE);
      if (!m || !s0 || !s || !ss || !c || !ndof || !sh || !ssky || !mx || !my || !cs || !cm || !comb) merr();
   }
   *S=*SS=*C=*SH=*SSKY=0;
   if (S0) *S0=0;
   for (img=0;img<Timg;img++) {
      s[img]=ss[img]=0.0;
      comb[img]=1;
   }
   for (img=0;img<Timg && (img<Nimg || IMG>=0);img++) if ((img==IMG || IMG==-1) && (checkstar<0 || !star_flag(stars[checkstar].x,stars[checkstar].y,img))) {
      if ((IS0 || S0) && !IMDIFF) {
	 NOSUBTRACT=1;
	 eval1(img,x,y,PA[img],PB[img],PC[img],s0+img,ss+img,c+img,ndof+img,sh+img,ssky+img,mx+img,my+img,cs+img,cm+img,1,0,star,1);
	 NOSUBTRACT=0;
	 m[img]=eval1(img,x,y,PA[img],PB[img],PC[img],s+img,ss+img,c+img,ndof+img,sh+img,ssky+img,mx+img,my+img,cs+img,cm+img,1,1,star,0);
	 if (s0[img]>0. && s[img]>0.) {
	    s0[img]=-2.5*log10(s[img]/s0[img]);
	    if (s0[img]>9.999) s0[img]=9.999;
	    else if (s0[img]<0.00001) s0[img]=0.00001;
	 }
	 else if (s0[img]>0.) s0[img]=9.999;
	 else s0[img]=0.00001;
      }
      else {
	 m[img]=eval1(img,x,y,PA[img],PB[img],PC[img],s+img,ss+img,c+img,ndof+img,sh+img,ssky+img,mx+img,my+img,cs+img,cm+img,1,1,star,1);
	 s0[img]=0.00001;
      }
      if (ss[img]>0.0 && hstmode[img].inst!=NONE) comb[img]=0;
   }
   if (ForceSameMag) {
      for (img=0;img<Timg;img++) if (!comb[img]) {
	 int j;
	 float ts=0,twt=0,wt;
	 for (j=img;j<Timg;j++) if (!comb[j] && hstmode[j].inst==hstmode[img].inst && hstmode[j].filt==hstmode[img].filt) {
	    wt = 1.0/(ss[j]*ss[j]);
	    ts += s[j]*wt;
	    twt += wt;
	 }
	 if (twt>0.0) {
	    ts /= twt;
	    for (j=img;j<Timg;j++) if (!comb[j] && hstmode[j].inst==hstmode[img].inst && hstmode[j].filt==hstmode[img].filt) {
	       s[j] = ts;
	       m[j] = eval1(j,x,y,PA[j],PB[j],PC[j],s+j,ss+j,c+j,ndof+j,sh+j,ssky+j,mx+j,my+j,cs+j,cm+j,0,1,star,1);
	       comb[j] = 1;
	    }
	 }
      }
   }
   for (img=0;img<Timg && (img<Nimg || IMG>=0);img++) if ((img==IMG || IMG==-1) && (checkstar<0 || !star_flag(stars[checkstar].x,stars[checkstar].y,img))) {
      if (m[img]<=0.) wt=swt=0.;
      else if (IMDIFF && IMG<0 && star>=0) {
	 swt=wt=(s[img]+refcts[star]/refmult[img]*cm[img])/ss[img]/ss[img];
	 if (swt<0.) {
	    swt=0.;
	    wt=-wt;
	 }
      }
      else if (IMDIFF && star>=0) {
	 wt=1./ss[img]/ss[img];
	 swt=wt*(s[img]+refcts[star]/refmult[img]*cm[img]);
	 if (swt<0.) swt=0.;
      }
      else if (IMG<0) {
	 if (s[img]>0) wt=swt=s[img]/ss[img]/ss[img];
	 else {
	    wt=-s[img]/ss[img]/ss[img];
	    swt=0.;
	 }
      }
      else {
	 wt=1./ss[img]/ss[img];
	 if (s[img]>0) swt=wt*s[img];
	 else swt=0.;
      }
      twt+=wt;
      tswt+=swt;
      tndof+=ndof[img];
      if (IS) IS[img]=s[img];
      if (IS0) IS0[img]=s0[img];
      if (ISS) ISS[img]=ss[img];
      if (ICM) ICM[img]=cm[img];
      if (IC) IC[img]=c[img];
      if (ISH) ISH[img]=sh[img];
      if (ISSKY) ISSKY[img]=ssky[img];
      (*S)+=s[img]*wt;
      if (S0) (*S0)+=s0[img]*wt;
      (*SS)+=ss[img]*ss[img]*wt*wt;
      (*C)+=c[img]*c[img]*ndof[img];
      (*SH)+=sh[img]*swt;
      (*SSKY)+=ssky[img]*wt;
   }
   if (tndof) (*C)=sqrt((*C)/tndof);
   if (tswt) (*SH)/=tswt;
   if (twt) {
      *SS=sqrt(*SS)/twt;
      if (IMG<0 && *S>0) {
	 float tss;
	 tss=sqrt(*S)/twt;
	 if (tss>*SS) *SS=tss;
      }
      (*S)/=twt;
      if (S0) (*S0)/=twt;
      (*SSKY)/=twt;
      if (*S) {
	 M=1./(sqrt((*C)*(*C)+0.25));
	 if (SearchMode==0) M*=(*S)/(*SS);
      }
   }
   assert(!isnan(M));
   return M;
}

float getnc(float m[3]) {
   float x;
   if ((x=2*m[1]-m[0]-m[2])<=0) {
      if (m[1]>=m[0] && m[1]>=m[2]) return 0.;
      else if (m[0]>=m[2]) return -1.;
      return 1.;
   }
   x=(m[2]-m[0])*0.5/x;
   if (x<-1) return -1.;
   if (x>1) return 1.;
   return x;
}

void getshape(float a,float b,float c,float*r,float*e) {
   float r1,r2,x;
   c/=a*b;
   a=1./(a*a);
   b=1./(b*b);
   *r=1./sqrt(sqrt(a*b-c*c*0.25));
   x=sqrt((a-b)*(a-b)+c*c);
   r1=1./sqrt((a+b-x)*0.5);
   r2=1./sqrt((a+b+x)*0.5);
   *e=1.-r2/r1;
   return;
}

int getclass(float a,float b,float c,int force1OK) {
   float r,e;

   if (force1OK && Force1) return 1;
   getshape(a,b,c,&r,&e);
   if (r>MaxS) return 5;
   if (r<MinS) return 4;
   if (e>MaxE) return 3;
   return 1;
}

int getsclass(int i) {
   float pa=0,pb=0,pc=0;
   int img;

   for (img=0;img<Nimg;img++) {
      pa+=stars[i].pa[img];
      pb+=stars[i].pb[img];
      pc+=stars[i].pc[img];
   }
   return getclass(pa/Nimg,pb/Nimg,pc/Nimg,0);
}

/*
void checkstars(char*estr,int verb) {
   int i,img;
   for (i=0;i<Nstars;i++) for (img=0;img<Nimg;img++) if (stars[i].is[img]>0.) {
      if (stars[i].type==4) {if (fabs(stars[i].pa[img]-0.1)>0.001 || fabs(stars[i].pb[img]-0.1)>0.001 || stars[i].pc[img]!=0.) {printf("ACK %s (%d %d %f %f %f)\n",estr,i,stars[i].type,stars[i].pa[img],stars[i].pb[img],stars[i].pc[img]); exit(-1);}}
      else if (stars[i].type==5) {if (fabs(stars[i].pa[img]-RBig)>0.001 || fabs(stars[i].pb[img]-RBig)>0.001 || stars[i].pc[img]!=0.) {printf("ACK %s (%d %d %f %f %f)\n",estr,i,stars[i].type,stars[i].pa[img],stars[i].pb[img],stars[i].pc[img]); exit(-1);}}
      else {if (fabs(stars[i].pa[img]-1.0)>0.001 || fabs(stars[i].pb[img]-1.0)>0.001 || stars[i].pc[img]!=0.) {printf("ACK %s (%d %d %f %f %f)\n",estr,i,stars[i].type,stars[i].pa[img],stars[i].pb[img],stars[i].pc[img]); exit(-1);}}
   }
   if (verb) printf("passed %s\n",estr);
   return;
}
*/

void setpsfpars(int i) {
   int img;
   if (stars[i].type<4) {
      getpsfpars(-1,stars[i].x,stars[i].y,stars[i].pa,stars[i].pb,stars[i].pc);
      return;
   }
   for (img=0;img<Timg;img++) {
      stars[i].pc[img]=0.;
      if (stars[i].type==4) stars[i].pa[img]=stars[i].pb[img]=0.1;
      else stars[i].pa[img]=stars[i].pb[img]=RBig;
   }
   return;
}

void zerostar(int IMG,float*bis) {
   int img;
   for (img=0;img<Timg;img++) if ((IMG<0 && img<Nimg) || img==IMG) bis[img]=0.;
   return;
}

float photsearch(int IMG,double*x,double*y,float*pa,float*pb,float*pc,float*s,float*ss,float*chi,float*sh,float*ssky,float*bis,float*biss,float*bicm,int *cl,int new,int id) {
   float step,stepc,m[3],q,orig,is,iss,ichi,indof,ish,isky,imx,imy,ics,icm,ma=0,mb=0,mc=0,sn[3],scale=0.75,aq;
   int i,img,it=0,CONT=1,N=0;
#define amin 0.05
#define cmax 1.90
#define PAD 1

   step=PosStep;
   if (new==1) step*=2;
   if (((!fitpsf || *s<(*ss)*SigPSF) && *cl<=2) || new==1 || PSFStep<=0.) getpsfpars(IMG,*x,*y,pa,pb,pc);
   clear_sky();
   if (!WARMSTART) while (CONT) {
      if (new==1 || new==-1) scale=0.2;
      CONT=0;
      it++;
      m[0]=eval(IMG,(*x)-step,*y,pa,pb,pc,s,NULL,ss,chi,sh,ssky,NULL,NULL,NULL,NULL,NULL,NULL,NULL,-1,id);
      m[1]=eval(IMG,*x,*y,pa,pb,pc,s,NULL,ss,chi,sh,ssky,NULL,NULL,NULL,NULL,NULL,NULL,NULL,-1,id);
      m[2]=eval(IMG,(*x)+step,*y,pa,pb,pc,s,NULL,ss,chi,sh,ssky,NULL,NULL,NULL,NULL,NULL,NULL,NULL,-1,id);
      q=getnc(m);
      if (new==0 && it>5) q*=5./it;
      if ((new==1 || new==-1) && scale<(aq=fabs(q))) scale=aq;
      if (q<-0.6 || q>0.6) CONT=1;
      *x+=q*step;
      assert(!isnan(*x));
      if (*x<-PAD || *x>=XMAX+PAD) {*s=0; if (bis) zerostar(IMG,bis); return 0;}
      m[0]=eval(IMG,*x,(*y)-step,pa,pb,pc,s,NULL,ss,chi,sh,ssky,NULL,NULL,NULL,NULL,NULL,NULL,NULL,-1,id);
      m[1]=eval(IMG,*x,*y,pa,pb,pc,s,NULL,ss,chi,sh,ssky,NULL,NULL,NULL,NULL,NULL,NULL,NULL,-1,id);
      m[2]=eval(IMG,*x,(*y)+step,pa,pb,pc,s,NULL,ss,chi,sh,ssky,NULL,NULL,NULL,NULL,NULL,NULL,NULL,-1,id);
      q=getnc(m);
      if (new==0 && it>5) q*=5./it;
      if ((new==1 || new==-1) && scale<(aq=fabs(q))) scale=aq;
      if (q<-0.6 || q>0.6) CONT=1;
      *y+=q*step;
      assert(!isnan(*y));
      if (*y<-PAD || *y>=YMAX+PAD) {*s=0; if (bis) zerostar(IMG,bis); return 0;}
      if (new==-1) {
	 if (step>PosStep*0.101) CONT=1;
	 else CONT=0;
	 if (scale>0.95) scale=0.95;
	 step*=scale;
	 if (step<PosStep*0.1) step=PosStep*0.1;
      }
      if (new==1) {
	 if (step>PosStep) CONT=1;
	 else CONT=0;
	 if (scale>0.95) scale=0.95;
	 step*=scale;
	 if (step<PosStep*0.1) step=PosStep*0.1;
      }
   }
   m[1]=eval(IMG,*x,*y,pa,pb,pc,s,NULL,ss,chi,sh,ssky,NULL,NULL,NULL,NULL,NULL,NULL,NULL,-1,id);
   if (*s<=0) {if (bis) zerostar(IMG,bis); return 0;}
   sn[1]=(*s)/(*ss);
   if (((fitpsf || new==1) && sn[1]>=SigPSF) || *cl>2) {
      if (PSFStep>0.) {
	 for (img=0;img<Timg;img++) if ((IMG<0 && img<Nimg) || img==IMG) {
	    CONT=1;
	    it=0;
	    if (new==1) {
	       step=PSFStep*2;
	       if (pa[img]>step*3) step=pa[img]*0.333;
	       if (pb[img]>step*3) step=pb[img]*0.333;
	       stepc=1.2;
	    }
	    else {
	       step=PSFStep;
	       if (pa[img]>step*5) step=pa[img]*0.2;
	       if (pb[img]>step*5) step=pb[img]*0.2;
	       stepc=0.6;
	    }
	    if (EPSF) {while (CONT) {
	       CONT=0;
	       it++;
	       if (pa[img]<step+amin) {
		  orig=(pa[img]-amin)/step-1;
		  pa[img]=step+amin;
	       }
	       else orig=0;
	       for (i=-1;i<=1;i++) m[i+1]=eval1(img,*x,*y,pa[img]+i*step,pb[img],pc[img],&is,&iss,&ichi,&indof,&ish,&isky,&imx,&imy,&ics,&icm,1,1,id,1);
	       q=getnc(m);
	       if (new==0 && it>3) q=3./it*(q-orig)+orig;
	       if (q<-0.6+orig || q>0.6+orig) CONT=1;
	       pa[img]+=q*step;
	       assert(!isnan(pa[img]));
	       if (pb[img]<step+amin) {
		  orig=(pb[img]-amin)/step-1;
		  pb[img]=step+amin;
	       }
	       else orig=0;
	       for (i=-1;i<=1;i++) m[i+1]=eval1(img,*x,*y,pa[img],pb[img]+i*step,pc[img],&is,&iss,&ichi,&indof,&ish,&isky,&imx,&imy,&ics,&icm,1,1,id,1);
	       q=getnc(m);
	       if (new==0 && it>3) q=3./it*(q-orig)+orig;
	       if (q<-0.6+orig || q>0.6+orig) CONT=1;
	       pb[img]+=q*step;
	       assert(!isnan(pb[img]));
	       if (pc[img]>cmax-stepc) {
		  orig=(pc[img]-cmax)/stepc+1;
		  pc[img]=cmax-stepc;
	       }
	       else if (pc[img]<stepc-cmax) {
		  orig=(pc[img]+cmax)/stepc-1;
		  pc[img]=stepc-cmax;
	       }
	       else orig=0;
	       for (i=-1;i<=1;i++) m[i+1]=eval1(img,*x,*y,pa[img],pb[img],pc[img]+i*stepc,&is,&iss,&ichi,&indof,&ish,&isky,&imx,&imy,&ics,&icm,1,1,id,1);
	       q=getnc(m);
	       if (new==0 && it>3) q=3./it*(q-orig)+orig;
	       if (q<-0.6+orig || q>0.6+orig) CONT=1;
	       pc[img]+=q*stepc;
	       assert(!isnan(pc[img]));
	       if (new==-1) {
		  if (CONT) {step*=0.95; stepc*=0.95;}
		  else {step*=0.6; stepc*=0.6;}
		  if (step>=0.1*PSFStep || stepc>=0.1) CONT=1;
		  else CONT=0;
	       }
	       if (new==1) {
		  if (CONT) {step*=0.95; stepc*=0.95;}
		  else {step*=0.6; stepc*=0.6;}
		  if (step>=PSFStep || stepc>=0.6) CONT=1;
		  else CONT=0;
	       }
	    }}
	    else {while (CONT) {
	       CONT=0;
	       it++;
	       if (pa[img]>step*5) step=(pa[img])*0.2;
	       if (pa[img]<step+amin) {
		  orig=(pa[img]-amin)/step-1;
		  pa[img]=step+amin;
	       }
	       else orig=0;
	       for (i=-1;i<=1;i++) m[i+1]=eval1(img,*x,*y,pa[img]+i*step,pa[img]+i*step,pc[img],&is,&iss,&ichi,&indof,&ish,&isky,&imx,&imy,&ics,&icm,1,1,id,1);
	       q=getnc(m);
	       if (new==0 && it>3) q=3./it*(q-orig)+orig;
	       if (q<-0.6+orig || q>0.6+orig) CONT=1;
	       pa[img]+=q*step;
	       assert(!isnan(pa[img]));
	       pb[img]=pa[img];
	       pc[img]=0;
	       if (new==-1) {
		  if (CONT) step*=0.95;
		  else step*=0.6;
		  if (step>=0.1*PSFStep) CONT=1;
		  else CONT=0;
	       }
	       if (new==1) {
		  if (CONT) step*=0.95;
		  else step*=0.6;
		  if (step>=PSFStep) CONT=1;
		  else CONT=0;
	       }
	    }}
	    ma+=pa[img];
	    mb+=pb[img];
	    mc+=pc[img];
	    N++;
	 }
      }
      else if (!Force1) {
	 for (img=0;img<Timg;img++) if ((IMG<0 && img<Nimg) || img==IMG) {
	    pa[img]=pb[img]=0.1;
	    pc[img]=0.;
	 }
	 m[0]=eval(IMG,*x,*y,pa,pb,pc,s,NULL,ss,chi,sh,ssky,NULL,NULL,NULL,NULL,NULL,NULL,NULL,-1,id);
	 sn[0]=(*s)/(*ss);
	 if (FitSky!=3 && FitSky!=4) {
	    for (img=0;img<Timg;img++) if ((IMG<0 && img<Nimg) || img==IMG) pa[img]=pb[img]=RBig;
	    m[2]=eval(IMG,*x,*y,pa,pb,pc,s,NULL,ss,chi,sh,ssky,NULL,NULL,NULL,NULL,NULL,NULL,NULL,-1,id);
	    sn[2]=(*s)/(*ss);
	 }
	 else {
	    m[2] = m[1]-100.0;
	    sn[2] = 0.;
	 }
	 if (m[0]>m[1] && m[0]>m[2]) {
	    sn[1]=sn[0];
	    for (img=0;img<Timg;img++) if ((IMG<0 && img<Nimg) || img==IMG) {
	       pa[img]=pb[img]=0.1;
	       pc[img]=0.;
	    }
	 }
	 else if (m[2]>m[1]) {
	    sn[1]=sn[2];
	    for (img=0;img<Timg;img++) if ((IMG<0 && img<Nimg) || img==IMG) {
	       pa[img]=pb[img]=RBig;
	       pc[img]=0.;
	    }
	 }
	 else getpsfpars(IMG,*x,*y,pa,pb,pc);
	 for (img=0;img<Timg;img++) if ((IMG<0 && img<Nimg) || img==IMG) {
	    ma+=pa[img];
	    mb+=pb[img];
	    mc+=pc[img];
	    N++;
	 }
      }
      else {
	 getpsfpars(IMG,*x,*y,pa,pb,pc);
	 for (img=0;img<Timg;img++) if ((IMG<0 && img<Nimg) || img==IMG) {
	    ma+=pa[img];
	    mb+=pb[img];
	    mc+=pc[img];
	    N++;
	 }
      }
      if (!WARMSTART) {
	 *cl=getclass(ma/N,mb/N,mc/N,1);
	 if (sn[1]<SigPSF) *cl=2;
	 if (PSFStep<=0. && !Force1) {
	    for (img=0;img<Timg;img++) if ((IMG<0 && img<Nimg) || img==IMG) {
	       if (*cl==4) {
		  pa[img]=pb[img]=0.1;
		  pc[img]=0.;
	       }
	       else if (*cl==5) {
		  pa[img]=pb[img]=RBig;
		  pc[img]=0.;
	       }
	       else getpsfpars(img,*x,*y,pa,pb,pc);
	    }
	 }
	 else if (!fitpsf && (*cl==1 || *cl==2)) getpsfpars(IMG,*x,*y,pa,pb,pc);
      }
   }
   else {
      getpsfpars(IMG,*x,*y,pa,pb,pc);
      if (!WARMSTART) {
	 if (sn[1]<SigPSF) *cl=2;
	 else *cl=1;
      }
   }
   m[0]=eval(IMG,*x,*y,pa,pb,pc,s,NULL,ss,chi,sh,ssky,bis,NULL,biss,NULL,NULL,NULL,bicm,-1,id);
   sn[0]=(*s)/(*ss);
   if (!WARMSTART && sn[0]<SigPSF) {
      if (!fitpsf && *cl>2) getpsfpars(IMG,*x,*y,pa,pb,pc);
      *cl=2;
   }
   if (WARMSTART && *cl<3) {
      if (sn[0]<SigPSF) *cl=2;
      *cl=1;
   }
   return m[0];
}

void avgstar(int i,stype *s) {
   double m;
   int img;

   if (s->type==1) stars[i].type=1;
   if (s->type==2 && stars[i].type>2) stars[i].type=2;
   if (stars[i].s>0 && s->s>0) m=stars[i].s/(stars[i].s+(s->s));
   else if (stars[i].s>0) m=1.0;
   else if (s->s>0) m=0.0;
   else m=0.5;
   stars[i].x=stars[i].x*m+(s->x)*(1-m);
   stars[i].y=stars[i].y*m+(s->y)*(1-m);
   stars[i].s=stars[i].s*m+(s->s)*(1-m);
   stars[i].ss=sqrt(stars[i].ss*stars[i].ss*m+(s->ss)*(s->ss)*(1-m));
   stars[i].chi=stars[i].chi*m+(s->chi)*(1-m);
   stars[i].sh=stars[i].sh*m+(s->sh)*(1-m);
   stars[i].sky=stars[i].sky*m+(s->sky)*(1-m);
   for (img=0;img<Nimg;img++) {
      stars[i].is[img]=stars[i].is[img]*m+(s->is[img])*(1-m);
      stars[i].iss[img]=stars[i].iss[img]*m+(s->iss[img])*(1-m);
      stars[i].icm[img]=stars[i].icm[img]*m+(s->icm[img])*(1-m);
      stars[i].pa[img]=stars[i].pa[img]*m+(s->pa[img])*(1-m);
      stars[i].pb[img]=stars[i].pb[img]*m+(s->pb[img])*(1-m);
      stars[i].pc[img]=stars[i].pc[img]*m+(s->pc[img])*(1-m);
   }
   if (PSFStep<=0) setpsfpars(i);
   return;
}

void phot(int x,int y,int mark,int nit) {
   int st,i,j,img,bcl;
   double sx,sy;
   float bs,bss,bsky,bc,bsh,m;
   stype old;
   static float *oldis,*oldiss,*oldcm,*oldpa,*oldpb,*oldpc,*bis,*biss,*bcm,*bpa,*bpb,*bpc;
   static int isfirst=1;

   if (isfirst) {
      isfirst=0;
      oldis=(float*)calloc(Timg,FLOATSIZE);
      oldiss=(float*)calloc(Timg,FLOATSIZE);
      oldcm=(float*)calloc(Timg,FLOATSIZE);
      oldpa=(float*)calloc(Timg,FLOATSIZE);
      oldpb=(float*)calloc(Timg,FLOATSIZE);
      oldpc=(float*)calloc(Timg,FLOATSIZE);
      bis=(float*)calloc(Timg,FLOATSIZE);
      biss=(float*)calloc(Timg,FLOATSIZE);
      bcm=(float*)calloc(Timg,FLOATSIZE);
      bpa=(float*)calloc(Timg,FLOATSIZE);
      bpb=(float*)calloc(Timg,FLOATSIZE);
      bpc=(float*)calloc(Timg,FLOATSIZE);
      if (!oldis || !oldiss || !oldcm || !oldpa || !oldpb || !oldpc || !bis || !biss || !bcm || !bpa || !bpb || !bpc) merr();
   }
   if (mark<0) {
      centroid(-1,x,y,&sx,&sy);
      if (sx<XMIN) sx=XMIN;
      if (sx>=XMAX) sx=XMAX-0.01;
      if (sy<YMIN) sy=YMIN;
      if (sy>=YMAX) sy=YMAX-0.01;
      if (ran[(int)(sy*SubResRef)][(int)(sx*SubResRef)]) return;
      bcl=1;
      bs=bss=0;
      m=photsearch(-1,&sx,&sy,bpa,bpb,bpc,&bs,&bss,&bc,&bsh,&bsky,bis,biss,bcm,&bcl,1,-1);
   }
   else {
      memcpy(&old,stars+mark,sizeof(stype));
      memcpy(oldis,stars[mark].is,Timg*FLOATSIZE);
      memcpy(oldiss,stars[mark].iss,Timg*FLOATSIZE);
      memcpy(oldcm,stars[mark].icm,Timg*FLOATSIZE);
      memcpy(oldpa,stars[mark].pa,Timg*FLOATSIZE);
      memcpy(oldpb,stars[mark].pb,Timg*FLOATSIZE);
      memcpy(oldpc,stars[mark].pc,Timg*FLOATSIZE);
      old.is=oldis;
      old.iss=oldiss;
      old.icm=oldcm;
      old.pa=oldpa;
      old.pb=oldpb;
      old.pc=oldpc;
      sx=old.x;
      sy=old.y;
      bcl=old.type;
      bs=old.s;
      bss=old.ss;
      memcpy(bpa,old.pa,Timg*FLOATSIZE);
      memcpy(bpb,old.pb,Timg*FLOATSIZE);
      memcpy(bpc,old.pc,Timg*FLOATSIZE);
      m=photsearch(-1,&sx,&sy,bpa,bpb,bpc,&bs,&bss,&bc,&bsh,&bsky,bis,biss,bcm,&bcl,0,mark);
   }
   if (sx<0 || sx>=X || sy<0 || sy>=Y || indx[(int)(sy*SubResRef)][(int)(sx*SubResRef)]>=0) m=0;
   else if (mark<0 && ran[(int)(sy*SubResRef)][(int)(sx*SubResRef)]) m=0;
   if (m==0) bs=bss=0;
   if (bs<=0) bs=0;
   if ((bs>=bss*SigFind && bs>0 && m>0) || mark>=0) {
      if (mark<0) {
	 if (Nstars>=MAXNSTARS) {
	    printf("****Too many stars\n");
	    exit(-1);
	 }
	 st=Nstars;
	 Nstars++;
	 stars[st].type=bcl;
      }
      else {
	 st=mark;
	 //if (mark>=0 && (fabs(bx-stars[mark].x0)>dPosMax || fabs(by-stars[mark].y0)>dPosMax)) bs=0;
	 if (m==0 || bs<=0) {
	    if (old.s!=0) markstars((int)old.x,(int)old.y);
	 }
	 else if (fabs(sx-old.x)>1.5*PosStep || fabs(sy-old.y)>1.5*PosStep) {
	    markstars((int)sx,(int)sy);
	 }
	 else for (img=0,j=1;img<Nimg && j;img++) if (fabs(bis[img]-oldis[img])>1 && fabs(bis[img]-oldis[img])>(bis[img]+oldis[img])*0.0005) {
	    markstars((int)sx,(int)sy);
	    j=0;
	 }
      }
      stars[st].flag=0;
      if (m==0 || bs<=0 || indx[(int)(sy*SubResRef)][(int)(sx*SubResRef)]>=0) {
	 stars[st].x=sx;
	 stars[st].y=sy;
	 stars[st].s=0;
	 stars[st].ss=1000.;
	 stars[st].sky=0;
	 stars[st].chi=999.99;
	 stars[st].sh=0;
	 for (i=0;i<Nimg;i++) stars[st].is[i]=0;
	 for (i=0;i<Nimg;i++) {
	    stars[st].pa[i]=bpa[i];
	    stars[st].pb[i]=bpb[i];
	    stars[st].pc[i]=bpc[i];
	 }
      }
      else {
	 stars[st].x=sx;
	 stars[st].y=sy;
	 stars[st].s=bs;
	 stars[st].ss=bss;
	 stars[st].sky=bsky;
	 stars[st].chi=bc;
	 stars[st].sh=bsh;
	 stars[st].type=bcl;
	 for (i=0;i<Nimg;i++) {
	    stars[st].is[i]=bis[i];
	    stars[st].iss[i]=biss[i];
	    stars[st].icm[i]=bcm[i];
	    stars[st].pa[i]=bpa[i];
	    stars[st].pb[i]=bpb[i];
	    stars[st].pc[i]=bpc[i];
	 }
	 if (nit>MaxIT/2 && nit%3==MaxIT%3) {
	    avgstar(st,&old);
	    if (stars[st].type<3 && !fitpsf) getpsfpars(-1,stars[st].x,stars[st].y,stars[st].pa,stars[st].pb,stars[st].pc);
	 }
	 indx[(int)(stars[st].y*SubResRef)][(int)(stars[st].x*SubResRef)]=st;
	 if (mark<0) for (j=-(int)RCombine;j<=(int)RCombine;j++) for (i=-(int)RCombine;i<=(int)RCombine;i++) if ((int)(sy*SubResRef)+j>=0 && (int)(sy*SubResRef)+j<Y && (int)(sx*SubResRef)+i>=0 && (int)(sx*SubResRef)+i<X && i*i+j*j<(RCombine-0.5)*(RCombine-0.5)) ran[(int)(sy*SubResRef)+j][(int)(sx*SubResRef)+i]=1;
      }
   }
   return;
}

void imgadd1(int img,int i) {
   int x1,y1,ix,iy;
   double x,y;

   if (stars[i].is[img]<=0) return;
   shift(img,stars[i].x,stars[i].y,&x,&y,1);
   ix=(int)x; iy=(int)y;
   calc1psf(img,x,y,RPSF[img],stars[i].pa[img],stars[i].pb[img],stars[i].pc[img],stars[i].type,0);

   for (y1=iy-RPSF[img];y1<=iy+RPSF[img];y1++) for (x1=ix-RPSF[img];x1<=ix+RPSF[img];x1++) if (ppixOK(img,x1,y1)) res[img][y1][x1]+=stars[i].is[img]/stars[i].icm[img]*psf[y1-iy][x1-ix];
   return;
}

void imgadd(int i) {
   int img;
   for (img=0;img<Nimg;img++) imgadd1(img,i);
   return;
}

void imgsub1(int img,int i) {
   int x1,y1,ix,iy;
   double x,y;

   if (stars[i].is[img]<=0) return;
   shift(img,stars[i].x,stars[i].y,&x,&y,1);
   ix=(int)x; iy=(int)y;
   calc1psf(img,x,y,RPSF[img],stars[i].pa[img],stars[i].pb[img],stars[i].pc[img],stars[i].type,0);
   for (y1=iy-RPSF[img];y1<=iy+RPSF[img];y1++) for (x1=ix-RPSF[img];x1<=ix+RPSF[img];x1++) if (ppixOK(img,x1,y1)) res[img][y1][x1]-=stars[i].is[img]/stars[i].icm[img]*psf[y1-iy][x1-ix];
   return;
}

void imgsub(int i) {
   int img;
   for (img=0;img<Nimg;img++) imgsub1(img,i);
   return;
}

//Quicksort from NR; 0-index by AED
#define SWAP(a,b) {memcpy(&temp,a,ssize);memcpy(a,b,ssize);memcpy(b,&temp,ssize);}
#define M 7
#define NSTACK 50

void sortstars() {
   int i,ir,j,k,l=0;
   int jstack=0,*istack,ssize;
   stype temp,a;

   if (WARMSTART) return;
   ir=Nstars-1;
   ssize=sizeof(stype);
   if ((istack=(int*)calloc(NSTACK,INTSIZE))==NULL) merr();
   for (;;) {
      if (ir-l<M) {
	 for (j=l+1;j<=ir;j++) {
	    memcpy(&a,stars+j,ssize);
	    for (i=j-1;i>=l;i--) {
	       if (stars[i].s>=a.s) break;
	       memcpy(stars+i+1,stars+i,ssize);
	    }
	    memcpy(stars+i+1,&a,ssize);
	 }
	 if (jstack==0) break;
	 ir=istack[jstack--];
	 l=istack[jstack--];
      } else {
	 k=(l+ir)>>1;
	 SWAP(stars+k,stars+l+1)
	 if (stars[l+1].s<stars[ir].s) SWAP(stars+l+1,stars+ir)
	 if (stars[l].s<stars[ir].s) SWAP(stars+l,stars+ir)
	 if (stars[l+1].s<stars[l].s) SWAP(stars+l+1,stars+l)
	 i=l+1;
	 j=ir;
	 memcpy(&a,stars+l,ssize);
	 for (;;) {
	    do i++; while (stars[i].s>a.s);
	    do j--; while (stars[j].s<a.s);
	    if (j<i) break;
	    SWAP(stars+i,stars+j);
	 }
	 memcpy(stars+l,stars+j,ssize);
	 memcpy(stars+j,&a,ssize);
	 jstack+=2;
	 if (jstack > NSTACK) {
	    printf("****NSTACK too small in sortstars.");
	    exit(-1);
	 }
	 if (ir-i+1>=j-l) {
	    istack[jstack]=ir;
	    istack[jstack-1]=i;
	    ir=j-1;
	 } else {
	    istack[jstack]=j-1;
	    istack[jstack-1]=l;
	    l=i;
	 }
      }
   }
   free(istack);
   return;
}
#undef M
#undef NSTACK
#undef SWAP

int igoodstar(int img,int i,int redge,int badmult) {
   int x,y,ix,iy,ct=0,tpix=0;
   double tx,ty;

   shift(img,stars[i].x,stars[i].y,&tx,&ty,1);
   if (!posOK(img,tx,ty)) return 0;
   if (star_flag(stars[i].x,stars[i].y,img)) return 0;
   if (stars[i].is[img]<=0 || stars[i].is[img]<=5*stars[i].iss[img]) return 0;
   ix=(int)tx; iy=(int)ty;
   for (x=-redge;x<=redge;x++) for (y=-redge;y<=redge;y++) if (posOK(img,ix+x,iy+y)) {
      tpix++;
      if (!datafOK(img,ix+x,iy+y)) {
	 if (abs(x)<2 && abs(y)<2) return 0;
	 ct++;
	 if (ct*badmult>=tpix) return 0;
      }
   }
   return 1;
}

int goodstar(int IMG,int i,int rsep,int redge0,int strict) {
   int j,ix,iy,redge;

   if (getsclass(i)>1) return 0;
   if (strict==2) {
      if (stars[i].s<stars[i].ss*10 || stars[i].s<=0 || stars[i].sh>0.15 || stars[i].sh<-0.15 || stars[i].chi>2.5) return 0;
   }
   else if (strict==1) {
      if (stars[i].s<stars[i].ss*10 || stars[i].s<=0 || stars[i].sh>0.25 || stars[i].sh<-0.25 || stars[i].chi>15.) return 0;
   }
   else {
      if (stars[i].s<stars[i].ss*5 || stars[i].s<=0) return 0;
   }
   ix=(int)stars[i].x; iy=(int)stars[i].y;
   for (j=0;j<i;j++) if (abs(ix-(int)stars[j].x)<=rsep && abs(iy-(int)stars[j].y)<=rsep) return 0;
   redge=redge0;
   if (redge<0) redge=rpsfMax;
   for (j=i+1;j<Nstars && stars[j].s*2>=stars[i].s;j++) if (abs(ix-(int)stars[j].x)<=redge && abs(iy-(int)stars[j].y)<=redge) return 0;
   if (ix<XMIN2 || ix-redge<XMIN || ix>=XMAX2 || ix+redge>=XMAX || iy<YMIN2 || iy-redge<YMIN || iy>=YMAX2 || iy+redge>=YMAX) return 0;
   if (IMG>=0) {
      if (redge0<0) return igoodstar(IMG,i,RPSF[IMG],3);
      else return igoodstar(IMG,i,redge,3);
   }
   for (j=0;j<Nimg;j++) {
      if (redge0<0) {if (igoodstar(j,i,RPSF[j],3)) return 1;}
      else {if (igoodstar(j,i,redge,3)) return 1;}
   }
   return 0;
}

//adapted from NR by AED;
#define SWAP(a,b) {temp=(a);(a)=(b);(b)=temp;}

void gaussj(double**a,int n,double*b) {
   int *indxc,*indxr,*ipiv;
   int i,icol=0,irow=0,j,k,l,ll;
   double big,dum,pivinv,temp;

   indxc=(int*)calloc(n,INTSIZE);
   indxr=(int*)calloc(n,INTSIZE);
   ipiv=(int*)calloc(n,INTSIZE);
   if (!indxc || !indxr || !ipiv) merr();
   for (i=0;i<n;i++) {
      big=0.0;
      for (j=0;j<n;j++) if (ipiv[j]!=1) for (k=0;k<n;k++) {
	 if (!ipiv[k]) {
	    if (fabs(a[j][k])>=big) {
	       big=fabs(a[j][k]);
	       irow=j;
	       icol=k;
	    }
	 }
	 else if (ipiv[k]>1) {
	    printf("gaussj: Singular Matrix!\n");
	    exit(-1);
	 }
      }
      ipiv[icol]++;
      if (irow!=icol) {
	 for (l=0;l<n;l++) SWAP(a[irow][l],a[icol][l])
	 SWAP(b[irow],b[icol])
      }
      indxr[i]=irow;
      indxc[i]=icol;
      if (a[icol][icol]==0.) {
	 printf("gaussj: Singular Matrix!\n");
	 exit(-1);
      }
      pivinv=1./a[icol][icol];
      a[icol][icol]=1.;
      for (l=0;l<n;l++) a[icol][l]*=pivinv;
      b[icol]*=pivinv;
      for (ll=0;ll<n;ll++) if (ll!=icol) {
	 dum=a[ll][icol];
	 a[ll][icol]=0.0;
	 for (l=0;l<n;l++) a[ll][l]-=a[icol][l]*dum;
	 b[ll]-=b[icol]*dum;
      }
   }
   for (l=n-1;l>=0;l--) if (indxr[l]!=indxc[l]) for (k=0;k<n;k++) SWAP(a[k][indxr[l]],a[k][indxc[l]])
   free(ipiv);
   free(indxr);
   free(indxc);
   return;
}
#undef SWAP

typedef double ltype[5];
//adapted from NR by AED;
void lfit(ltype*data,int N,double*vec,int NV) {
   int i,j,k;
   double wt,*afunc,**covar;

   afunc=(double*)calloc(NV,DOUBLESIZE);
   covar=(double**)calloc(NV,PTRSIZE);
   if (!afunc || !covar) merr();
   for (i=0;i<NV;i++) {
      covar[i]=(double*)calloc(NV,DOUBLESIZE);
      if (!covar[i]) merr();
      for (j=0;j<NV;j++) covar[i][j]=0.;
      vec[i]=0.;
   }
   for (i=0;i<N;i++) {
      afunc[0]=1.;
      afunc[1]=data[i][2];
      afunc[2]=data[i][3];
      if (NV>3) {
	 afunc[3]=data[i][2]*data[i][2];
	 afunc[4]=data[i][3]*data[i][3];
	 afunc[5]=data[i][2]*data[i][3];
      }
      for (j=0;j<NV;j++) {
	 wt=afunc[j]*data[i][1];
	 for (k=0;k<=j;k++) covar[j][k]+=wt*afunc[k];
	 vec[j]+=data[i][0]*wt;
      }
   }
   for (i=1;i<NV;i++) for (j=0;j<i;j++) covar[j][i]=covar[i][j];
   gaussj(covar,NV,vec);
   free(afunc);
   for (i=0;i<NV;i++) free(covar[i]);
   free(covar);
   return;
}

void solvepsf(int Np,int**plist) {
   int img,i,N,xx,v;
   ltype *list;
   double sd;

   list=(ltype*)calloc(Np,sizeof(ltype));
   if (!list) merr();
   fprintf(finfo,"PSF solution\n");
   for (img=0;img<Nimg;img++) {
      printf("Image %d Solved PSF:\n",img+1);
      if (!FakeStars[0] && VerboseData>0) fprintf(fverb,"PSF image %d:",img+1);
      for (v=0;v<3;v++) {
	 N=0;
	 for (i=0;i<Np;i++) if (plist[i][img+1]) if (stars[plist[i][0]].s>0 && stars[plist[i][0]].s>stars[plist[i][0]].ss*SigPSF) {
	    list[N][1]=stars[plist[i][0]].s/stars[plist[i][0]].ss;
	    list[N][2]=stars[plist[i][0]].x-X*0.5;
	    list[N][3]=stars[plist[i][0]].y-Y*0.5;
	    if (v==0) list[N++][0]=stars[plist[i][0]].pa[img];
	    else if (v==1) list[N++][0]=stars[plist[i][0]].pb[img];
	    else list[N++][0]=stars[plist[i][0]].pc[img];
	 }
	 if (N) {
	    xx=1;
	    while (xx) {
	       double z=0;
	       xx=0;
	       for (i=0;i<6;i++) apsf[img][v][i]=0;
	       if (PSFsol==0) {
		  for (i=0;i<N;i++) {
		     apsf[img][v][0]+=list[i][0]*list[i][1];
		     z+=list[i][1];
		  }
		  if (z<=0) return;
		  apsf[img][v][0]/=z;
	       }
	       else if (PSFsol==1) lfit(list,N,apsf[img][v],3);
	       else lfit(list,N,apsf[img][v],6);
	       sd=z=0;
	       for (i=0;i<N;i++) {
		  list[i][4]=apsf[img][v][0]+list[i][2]*apsf[img][v][1]+list[i][3]*apsf[img][v][2]+list[i][2]*list[i][2]*apsf[img][v][3]+list[i][3]*list[i][3]*apsf[img][v][4]+list[i][2]*list[i][3]*apsf[img][v][5];
		  sd+=(list[i][0]-list[i][4])*(list[i][0]-list[i][4])*list[i][1]*list[i][1];
		  z+=list[i][1]*list[i][1];
	       }
	       sd=sqrt(sd/z)*2;
	       for (i=0;i<N;)
		  if (fabs(list[i][0]-list[i][4])>sd) {
		     memcpy(list[i],list[--N],sizeof(ltype));
		     xx=1;
		  }
		  else i++;
	    }
	 }
	 printf("  %c= %5.3f + %5.3fX + %5.3fY + %5.3fX^2 + %5.3fY^2 + %5.3fXY\n",97+v,apsf[img][v][0],apsf[img][v][1]*X,apsf[img][v][2]*Y,apsf[img][v][3]*X*X,apsf[img][v][4]*Y*Y,apsf[img][v][5]*X*Y);
	 fprintf(finfo," %f %f %f %f %f %f\n",apsf[img][v][0],apsf[img][v][1],apsf[img][v][2],apsf[img][v][3],apsf[img][v][4],apsf[img][v][5]);
	 if (!FakeStars[0] && VerboseData>0) fprintf(fverb," %f %f %f %f %f %f",apsf[img][v][0],apsf[img][v][1]*X,apsf[img][v][2]*Y,apsf[img][v][3]*X*X,apsf[img][v][4]*Y*Y,apsf[img][v][5]*X*Y);
      }
      printf("X = x/%d-0.5; Y = y/%d-0.5\n",X,Y);
      fflush(stdout);
      if (!FakeStars[0] && VerboseData>0) fprintf(fverb,"\n");
   }
   for (i=0;i<Np;i++) {
      imgadd(plist[i][0]);
      if (stars[plist[i][0]].s>0) indx[(int)(stars[plist[i][0]].y*SubResRef)][(int)(stars[plist[i][0]].x*SubResRef)]=-1;
      phot(stars[plist[i][0]].x,stars[plist[i][0]].y,plist[i][0],0);
      imgsub(plist[i][0]);
   }
   free(list);
   return;
}

void getimgstr(int img,char*imgstr) {
#ifdef USEWFPC2
   if (hstmode[img].inst==WFPC2) {
      sprintf(imgstr,"%s (%s, %3.1f sec)",base[img],WFPC2imagestring(img),iEXP[img]);
   }
   else
#endif
#ifdef USEACS
   if (hstmode[img].inst==ACS) {
      sprintf(imgstr,"%s (%s, %3.1f sec)",base[img],ACSimagestring(img),iEXP[img]);
   }
   else
#endif
#ifdef USEWFC3
   if (hstmode[img].inst==WFC3) {
      sprintf(imgstr,"%s (%s, %3.1f sec)",base[img],WFC3imagestring(img),iEXP[img]);
   }
   else
#endif
      sprintf(imgstr,"%s (%3.1f sec)",base[img],iEXP[img]);
}

int fix1psf(int img,int Np,int**plist,int it) {
   int x,y,xx,yy,i,j,N;
   float *list,*sig,av,sd,z,**off,**avpsf,*ssky,toff=0,toffd=0,d;
   float is,mx,my;
   double tx,ty;

   ssky=(float*)calloc(Np,FLOATSIZE);
   list=(float*)calloc(Np,FLOATSIZE);
   sig=(float*)calloc(Np,FLOATSIZE);
   off=(float**)calloc(2*RPSF[img]+1,PTRSIZE);
   avpsf=(float**)calloc(2*RPSF[img]+1,PTRSIZE);
   if (!ssky || !list || !sig || !off || !avpsf) merr();
   off+=RPSF[img];
   avpsf+=RPSF[img];
   for (y=-RPSF[img];y<=RPSF[img];y++) {
      off[y]=(float*)calloc(2*RPSF[img]+1,FLOATSIZE);
      avpsf[y]=(float*)calloc(2*RPSF[img]+1,FLOATSIZE);
      if (!off[y] || !avpsf[y]) merr();
      off[y]+=RPSF[img];
      avpsf[y]+=RPSF[img];
      for (x=-RPSF[img];x<=RPSF[img];x++) avpsf[y][x]=0;
   }
   N=0;
   for (i=0;i<Np;i++) if (plist[i][img+1]) {
      j=plist[i][0];
      shift(img,stars[j].x,stars[j].y,&tx,&ty,1);
      clear_sky();
      if (FitSky==0 && fsky[img].lastoffset>=0) ssky[i]=skyval(img,(int)tx,(int)ty);
      else if (FitSky==3 || FitSky==4) {
	 int ix,iy,cx,cy;
	 imgadd1(img,j);
	 ix=(int)(tx+100)-100;
	 iy=(int)(ty+100)-100;
	 cx=(int)((tx-ix)*50+0.5);
	 if (cx<0) cx=0; if (cx>50) cx=50;
	 cy=(int)((ty-iy)*50+0.5);
	 if (cy<0) cy=0; if (cy>50) cy=50;
	 if (FitSky==3) ssky[i]=eval1_sky3(img,ix,iy,cx,cy);
	 else ssky[i]=eval1_sky4(img,ix,iy,cx,cy,&mx,&my);
	 imgsub1(img,j);
      }
      else ssky[i]=getsky_norm(img,(int)tx,(int)ty,NULL);
      calc1psf(img,tx,ty,RPSF[img],stars[j].pa[img],stars[j].pb[img],stars[j].pc[img],stars[j].type,0);
      for (y=-RPSF[img];y<=RPSF[img];y++) for (x=-RPSF[img];x<=RPSF[img];x++) avpsf[y][x]+=psf[y][x];
      N++;
   }
   //if (img==8) printf("**IMG 8; N=%d\n",N); // AEDDEBUG
   if (N==0) return 0;
   for (y=-RPSF[img];y<=RPSF[img];y++) for (x=-RPSF[img];x<=RPSF[img];x++) avpsf[y][x]/=N;
   for (y=-RPSF[img];y<=RPSF[img];y++) for (x=-RPSF[img];x<=RPSF[img];x++) {
      N=0;
      //if (img==8) printf("  X,Y = %d,%d\n",x,y); // AEDDEBUG
      for (i=0;i<Np;i++) if (plist[i][img+1]) {
	 int ii;
	 ii=plist[i][0];
	 if ((is=stars[ii].is[img]/stars[ii].icm[img])>0 && stars[ii].s>0 && stars[ii].s>stars[ii].ss*SigFind) {
	    double tx,ty;
	    shift(img,stars[ii].x,stars[ii].y,&tx,&ty,1);
	    xx=x+(int)tx;
	    yy=y+(int)ty;
	    if (posOK(img,xx,yy)) {
	       if (datafOK(img,xx,yy)) {
		  sig[N]=is/sqrt(noise(img,xx,yy,ssky[i]));
		  if (!psfstars[0] && ssky[i]/iEXP[img]>0) sig[N]/=1+ssky[i]/iEXP[img]*100.;
		  list[N++]=(res[img][yy][xx]-ssky[i])/is;
		  //if (img==8) printf("    %g %g good (%g %g %g %g %g)\n",sig[N-1],list[N-1],is,res[img][yy][xx],noise(img,xx,yy,ssky[i]),ssky[i],iEXP[img]); // AEDDEBUG
	       }
	       else {
		  float ns=0,rs=0;
		  int ct=0,dx,dy;
		  for (dy=-1;dy<=1;dy++) for (dx=-1;dx<=1;dx++) if (ppixfOK(img,xx+dx,yy+dy)) {
		     ns+=noise(img,xx+dx,yy+dy,ssky[i]);
		     rs+=res[img][yy+dy][xx+dx];
		     ct++;
		  }
		  if (ct) {
		     sig[N]=is/sqrt(ns/ct);
		     if (!psfstars[0] && ssky[i]/iEXP[img]>0) sig[N]/=1+ssky[i]/iEXP[img]*100.;
		     list[N++]=(rs/ct-ssky[i])/is;
		     //if (img==8) printf("    %g %g bad (%g %d %g %g %g %g)\n",sig[N-1],list[N-1],is,ct,rs,ns,ssky[i],iEXP[img]); // AEDDEBUG
		  }
	       }
	    }
	 }
      }
      av=sd=0.;
      xx=1;
      while (xx) {
	 xx=0;
	 av=sd=z=0.;
	 for (i=0;i<N;i++) {
	    av+=list[i]*sig[i];
	    z+=sig[i];
	 }
	 //if (img==8) printf("  totals %g %g %g %g\n",av,z,avpsf[y][x],poff[img][y][x]); // AEDDEBUG
	 //printf("%d %d %d: %f %f %f %f -> ",img,x,y,av,z,avpsf[y][x],poff[img][y][x]);
	 //prior of S/N=100 star;
	 //av-=poff[img][y][x]*100.;
	 //z+=100.;
	 if (z>0) {
	    //av/=z;
	    //delta = [ sum(list*sig*(avpsf+delta)) - prior*poff ] / [ sum(sig*(avpsf+delta)) + prior ];
	    //delta=(av*(avpsf+delta)-prior*poff)/(z*(avpsf+delta)+prior);
	    //delta*delta*z+delta*(z*avpsf+prior-av)+prior*poff-av*avpsf
	    //delta = 0.5*(av-z*avpsf-prior+sqrt(z*z*avpsf*avpsf+prior*prior+av*av+2*prior*z*avpsf+2*z*avpsf*av-2*prior*av-4*prior*z*poff))/z;
	    d = z*z*avpsf[y][x]*avpsf[y][x]+100.+av*av+20.*z*avpsf[y][x]+2.*z*avpsf[y][x]*av-20.*av-40.*z*poff[img][y][x];
	    if (d<0.) d=0.;
	    av=0.5*(av-z*avpsf[y][x]-10.+sqrt(d))/z;
	    for (i=0;i<N;i++) sd+=(list[i]-av)*(list[i]-av)*sig[i];
	    //if (img==8) printf("  %g %g",av,sd); // AEDDEBUG
	    sd=sqrt(sd/z)*3.5;
	    //if (img==8) printf(" %g\n",sd); // AEDDEBUG
	    for (i=0;i<N;) {
	       if (fabs(list[i]-av)>sd) {
		  //if (img==8) printf("    deleting %g %g\n",list[i],sig[i]); // AEDDEBUG
		  list[i]=list[--N];
		  sig[i]=sig[N];
		  xx=1;
	       }
	       else i++;
	    }
	 }
	 else av=-poff[img][y][x];
	 //printf("%f\n",av);
      }
      off[y][x]=av;
      toff+=avpsf[y][x];
      toffd+=avpsf[y][x]+off[y][x];
      //if (img==8) printf("  final %g %g (%g %g)\n",av,avpsf[y][x],toff,toffd); // AEDDEBUG
   }
   toff/=toffd;
   //if (img==8) printf("toff = %g\n",toff); // AEDDEBUG
   for (i=0;i<Nstars;i++) imgadd1(img,i);
   for (y=-RPSF[img];y<=RPSF[img];y++) for (x=-RPSF[img];x<=RPSF[img];x++) {
      //printf("%d %d %d %d: %f %f %f ",img,x,y,it,toff,off[y][x],avpsf[y][x]);
      off[y][x]=toff*off[y][x]+(toff-1)*avpsf[y][x];
      if (it>4) off[y][x]/=1+it*0.1;
      poff[img][y][x]+=off[y][x];
      //printf("(%f)\n",poff[img][y][x]);
      //if (img==8) printf("%d,%d: off = %g (%g %g %g)\n",x,y,poff[img][y][x],off[y][x],toff,avpsf[y][x]); // AEDDEBUG
   }
   poffreset=1;
   for (i=0;i<Nstars;i++) imgsub1(img,i);
   poffreset=0;
   free(ssky);
   free(list);
   free(sig);
   i=0;
   for (y=-RPSF[img];y<=RPSF[img] && !i;y++) for (x=-RPSF[img];x<=RPSF[img] && !i;x++) if (fabs(off[y][x])>=0.001+0.0001*it) i=1;
   for (y=-RPSF[img];y<=RPSF[img];y++) {
      free(off[y]-RPSF[img]);
      free(avpsf[y]-RPSF[img]);
   }
   free(off-RPSF[img]);
   free(avpsf-RPSF[img]);
   return i;
}

int fixpsf(int Np,int**plist,int it) {
   int img,rval=0;
   for (img=0;img<Nimg;img++) if (fix1psf(img,Np,plist,it)) rval=1;
   return rval;
}

void calcpsf(int ext,int fld) {
#define minNpsf 50
#define goodNpsf 150
#define maxNpsf 300
   int **plist,*tlist;
   int Np=0,Nt=0,i,pe,pc,bi,j,it=0,img,cont=1;
   double px,py,r,br;
   FILE *f;
   char str[161],*ptr,*ptr2;

   if (PSFsol<0 && !PSFres) return;
   plist=(int**)calloc(Nstars,PTRSIZE);
   tlist=(int*)calloc(Nstars,INTSIZE);
   if (!plist || !tlist) merr();
   for (i=0;i<Nstars;i++) {
      plist[i]=(int*)calloc(Nimg+1,INTSIZE);
      if (!plist[i]) merr();
   }
   sortstars();
   setindx(-1);
   for (i=0;i<Nstars;i++) {
      indx[(int)(stars[i].y*SubResRef)][(int)(stars[i].x*SubResRef)]=i;
      stars[i].flag=0;
   }
   if (psfstars[0]) {
      if ((f=fopen(psfstars,"r"))==NULL) {
	 printf("**%s not found\n",psfstars);
	 fflush(stdout);
      }
      else {
	 while (fgets(str,161,f)) {
	    pe=strtol(str,&ptr,10);
	    pc=strtol(ptr,&ptr,10);
	    px=strtod(ptr,&ptr)-psfoff;
	    py=strtod(ptr,&ptr2)-psfoff;
	    if (pe==ext && pc==fld+1 && ptr!=ptr2) {
	       bi=-1;
	       br=RCombine*RCombine/(SubResRef*SubResRef);
	       if (br<1) br=1;
	       for (i=0;i<Nstars;i++) if ((r=(px-stars[i].x)*(px-stars[i].x)+(py-stars[i].y)*(py-stars[i].y))<br && goodstar(-1,i,rPhotPlusPSFmax,-1,0)) {
		  for (j=0;j<Np && i!=plist[j][0];j++);
		  if (j>=Np) {
		     br=r;
		     bi=i;
		  }
	       }
	       if (bi>=0) {
		  plist[Np][0]=bi;
		  fprintf(fpsfs,"%d %d %7.2f %7.2f",ext,fld+1,stars[bi].x,stars[bi].y);
		  for (j=0;j<Nimg;j++) if (igoodstar(j,bi,RPSF[j],5)) {
		     plist[Np][j+1]=1;
		     fprintf(fpsfs," 1");
		  }
		  else {
		     plist[Np][j+1]=0;
		     fprintf(fpsfs," 0");
		  }
		  fprintf(fpsfs,"\n");
		  Np++;
		  markstars((int)stars[bi].x,(int)stars[bi].y);
	       }
	    }
	 }
	 fclose(f);
      }
   }
   else {
      int NpMin = 0;
      int *NpImg = (int*)calloc(Nimg,sizeof(int)); if (NpImg==0) merr();
      for (i=0;i<Nstars;i++) if (((stars[i].s>=stars[i].ss*50 && NpMin<maxNpsf) || (stars[i].s>=stars[i].ss*10 && NpMin<goodNpsf) || NpMin<minNpsf) && goodstar(-1,i,rPhotPlusPSFmax,-1,1)) {
	 plist[Np][0]=i;
	 fprintf(fpsfs,"%d %d %6.2f %6.2f",ext,fld+1,stars[i].x,stars[i].y);
	 NpMin=0;
	 for (j=0;j<Nimg;j++) {
	    if (((stars[i].is[j]>=stars[i].iss[j]*50 && NpImg[j]<maxNpsf) || (stars[i].is[j]>=stars[i].iss[j]*10 && NpImg[j]<goodNpsf) || NpImg[j]<minNpsf) && igoodstar(j,i,RPSF[j],5)) {
	       plist[Np][j+1]=1;
	       fprintf(fpsfs," 1");
	       NpImg[j]++;
	    }
	    else {
	       plist[Np][j+1]=0;
	       fprintf(fpsfs," 0");
	    }
	    if (j==0) NpMin = NpImg[j];
	    else if (NpImg[j]<NpMin) NpMin=NpImg[j];
	 }
	 fprintf(fpsfs,"\n");
	 Np++;
      }
      free(NpImg);
      while (Np>0 && cont) {
	 double mn=0,sd=0;

	 cont=0;
	 for (i=0;i<Np;i++) mn+=stars[plist[i][0]].sh;
	 mn/=Np;
	 for (i=0;i<Np;i++) sd+=(stars[plist[i][0]].sh-mn)*(stars[plist[i][0]].sh-mn);
	 sd=2.00*sqrt(sd/Np);
	 for (i=0;i<Np;i++) if (fabs(stars[plist[i][0]].sh-mn)>sd) {
	    Np--;
	    if (i!=Np) memcpy(plist[i],plist[Np],INTSIZE*(Nimg+1));
	    i--;
	    cont=1;
	 }
      }
      for (i=0;i<Np;i++) markstars((int)stars[plist[i][0]].x,(int)stars[plist[i][0]].y);
   }
   for (i=0;i<Nstars;i++) imgadd(i);
   POSPSF=0;
   for (i=0;i<Nstars;i++) imgsub(i);
   for (i=0;i<Nstars;i++) if (stars[i].flag) tlist[Nt++]=i;
   printf("%d PSF stars; %d neighbors\n",Np,Nt-Np);
   fflush(stdout);
   if (!FakeStars[0] && VerboseData>0) fprintf(fverb,"PSF: %d %d\n",Np,Nt-Np);
   if (PSFsol>=0 && Np) solvepsf(Np,plist);
   if (PSFres && Np) while (fixpsf(Np,plist,++it)) {
      for (i=0;i<Nt;i++) {
	 imgadd(tlist[i]);
	 if (stars[tlist[i]].s>0) indx[(int)(stars[tlist[i]].y*SubResRef)][(int)(stars[tlist[i]].x*SubResRef)]=-1;
	 phot((int)stars[tlist[i]].x,(int)stars[tlist[i]].y,tlist[i],0);
	 imgsub(tlist[i]);
      }
   }
   if (it==20) printf("Poorly converging PSF\n");
   printf("Central pixel PSF adjustments:\n");
   if (Np<50) {
      fprintf(fwarn,"Only %d stars for PSF measurement\n",Np);
      fflush(fwarn);
   }
   for (img=0;img<Nimg;img++) {
      for (i=0,j=0;i<Np;i++) if (plist[i][img+1]) j++;
      printf("image %d: %d stars, %f\n",img+1,j,poff[img][0][0]);
      fflush(stdout);
      if (!FakeStars[0] && VerboseData>0) fprintf(fverb,"PSF image %d: %d %f\n",img+1,j,poff[img][0][0]);
      if (j<25 && j<Np) {
	 char imgstr[161];
	 getimgstr(img,imgstr);
	 fprintf(fwarn,"Only %d stars for PSF measurement in image %d, %s\n",j,img+1,imgstr);
	 fflush(fwarn);
      }
   }
   for (i=0;i<Nstars;i++) free(plist[i]);
   free(plist);
   free(tlist);
   for (i=0;i<Nstars;i++) imgadd(i);
   POSPSF=1;
   for (i=0;i<Nstars;i++) imgsub(i);
   fflush(stdout);
#ifdef PGPLOT
   if (DiagPlotType[0]) {
      for (img=0;img<Nimg;img++) {
	 int x,y,i=0;
	 psfDiagData[img][DIAG_EXT][DIAG_Z].N = 2*RPSF[img]+1;
	 psfDiagData[img][DIAG_EXT][DIAG_Z].mid = poff[img][0][0];
	 psfDiagData[img][DIAG_EXT][DIAG_Z].data = (float*)calloc((RPSF[img]*2+1)*(RPSF[img]*2+1),sizeof(float));
	 if (!psfDiagData[img][DIAG_EXT][DIAG_Z].data) merr();
	 for (y=-RPSF[img];y<=RPSF[img];y++) for (x=-RPSF[img];x<=RPSF[img];x++) psfDiagData[img][DIAG_EXT][DIAG_Z].data[i++]=poff[img][y][x];
      }
   }
#endif
   return;
}

int clean(int it,int show) {
   int i,ii,j,b=0,ct=0,xx,yy,rad,x3,y3;
   float x,y,r,br,RCeff,R2Ceff;

   if (WARMSTART) return 0;
   //if (it*2<=MaxIT) RCeff=RCombine;
   //else RCeff=RCombine/(1.+(it*2.-MaxIT)/MaxIT);
   RCeff=RCombine;
   R2Ceff=RCeff*RCeff;
   for (i=1;i<Nstars;i++) if (stars[i].x>=0) {
      rad=(int)(RCeff+0.999);
      br=R2Ceff;
      for (xx=-rad;xx<=rad;xx++) if ((x3=(int)(stars[i].x*SubResRef)+xx)>=XMIN && x3<XMAX) for (yy=-rad;yy<=rad;yy++) if ((xx || yy) && (y3=(int)(stars[i].y*SubResRef)+yy)>=YMIN && y3<YMAX && (j=indx[y3][x3])>=0 && j<i && (x=fabs(stars[i].x-stars[j].x)*SubResRef)<=RCeff && (y=fabs(stars[i].y-stars[j].y)*SubResRef)<=RCeff) {
	 r=x*x+y*y;
	 if (r<br) {
	    br=r;
	    b=j;
	 }
      }
      if (br<R2Ceff) {
	 imgadd(i);
	 imgadd(b);
	 if (stars[i].s>0) indx[(int)(stars[i].y*SubResRef)][(int)(stars[i].x*SubResRef)]=-1;
	 if (stars[b].s>0) indx[(int)(stars[b].y*SubResRef)][(int)(stars[b].x*SubResRef)]=-1;
	 markstars((int)stars[i].x,(int)stars[i].y);
	 markstars((int)stars[b].x,(int)stars[b].y);
	 avgstar(b,stars+i);
	 stars[b].s*=2;
	 stars[b].ss*=2;
	 phot((int)stars[b].x,(int)stars[b].y,b,0);
	 imgsub(b);
	 stars[i].x=-1;
      }
   }
   ii=0;
   for (i=0;i<Nstars;i++) if (stars[i].x>=0 && stars[i].s>0) {
      if (ii<i) {
	 stype temp;
	 memcpy(&temp,stars+ii,sizeof(stype));
	 memcpy(stars+ii,stars+i,sizeof(stype));
	 memcpy(stars+i,&temp,sizeof(stype));
	 indx[(int)(stars[ii].y*SubResRef)][(int)(stars[ii].x*SubResRef)]=ii;
      }
      ii++;
   }
   ct=Nstars-ii;
   Nstars=ii;
   if (ct && show) {printf("%d stars eliminated\n",ct); fflush(stdout);}
   return ct;
}

float imginterp(int img,float x,float y) {
   int ix,iy;
   double v=0.,tw=0.,w,mx,my;

   ix=(int)(x-0.5);
   iy=(int)(y-0.5);
   mx=x-0.5-ix;
   my=y-0.5-iy;
   if (ppixfOK(img,ix,iy)) {
      tw+=(w=(1-mx)*(1-my));
      v+=w*data[img][iy][ix];
   }
   if (ppixfOK(img,ix+1,iy)) {
      tw+=(w=mx*(1-my));
      v+=w*data[img][iy][ix+1];
   }
   if (ppixfOK(img,ix,iy+1)) {
      tw+=(w=(1-mx)*my);
      v+=w*data[img][iy+1][ix];
   }
   if (ppixfOK(img,ix+1,iy+1)) {
      tw+=(w=mx*my);
      v+=w*data[img][iy+1][ix+1];
   }
   if (tw<=0.) return iDMIN[img]-1.;
   return v/tw;
}

float**convv;
int RPSFsub;

#define DECON_MULT 0.001
float convfunc(float*c) {
   int x,y,i,j,n;
   float vv=0.,v;

   n=2*RPSFsub+1;
   // prior to reduce ringing
   for (i=1;i<=n*n;i++) vv+=c[i]*c[i]*DECON_MULT;
   // minimize difference between convv (img PSF) and convolution of c and psf (reference PSF)
   for (j=-RPSFsub;j<=RPSFsub;j++) for (i=-RPSFsub;i<=RPSFsub;i++) {
      v=-convv[j][i];
      for (y=-RPSFsub;y<=RPSFsub;y++) if (j-y>=-RPSFsub && j-y<=RPSFsub) for (x=-RPSFsub;x<=RPSFsub;x++) if (i-x>=-RPSFsub && i-x<=RPSFsub) {
	 v+=c[1+x+RPSFsub+(y+RPSFsub)*n]*psf[j-y][i-x];
      }
      vv+=v*v;
   }
   return vv;
}

void dconvfunc(float*c,float*g) {
   int x,y,i,j,n;
   float v,dv;

   n=2*RPSFsub+1;
   for (i=1;i<=n*n;i++) g[i]=2*c[i]*DECON_MULT;
   for (j=-RPSFsub;j<=RPSFsub;j++) for (i=-RPSFsub;i<=RPSFsub;i++) {
      v=-convv[j][i];
      for (y=-RPSFsub;y<=RPSFsub;y++) if (j-y>=-RPSFsub && j-y<=RPSFsub) for (x=-RPSFsub;x<=RPSFsub;x++) if (i-x>=-RPSFsub && i-x<=RPSFsub) v+=c[1+x+RPSFsub+(y+RPSFsub)*n]*psf[j-y][i-x];
      dv=2*v;
      for (y=-RPSFsub;y<=RPSFsub;y++) if (j-y>=-RPSFsub && j-y<=RPSFsub) for (x=-RPSFsub;x<=RPSFsub;x++) if (i-x>=-RPSFsub && i-x<=RPSFsub) g[1+x+RPSFsub+(y+RPSFsub)*n]+=dv*psf[j-y][i-x];
   }
   return;
}
#undef DECON_MULT

#define ITMAX 2000
#define EPS 1.0e-10
#define FREEALL free_vector(xi,1,n);free_vector(h,1,n);free_vector(g,1,n);
void frprmn(float p[], int n, float ftol, int *iter, float *fret,
	float (*func)(float []), void (*dfunc)(float [], float []))
{
	void linmin(float p[], float xi[], int n, float *fret,
		float (*func)(float []));
	int j,its;
	float gg,gam,fp,dgg;
	float *g,*h,*xi;

	g=vector(1,n);
	h=vector(1,n);
	xi=vector(1,n);
	fp=(*func)(p);
	(*dfunc)(p,xi);
	for (j=1;j<=n;j++) {
		g[j] = -xi[j];
		xi[j]=h[j]=g[j];
	}
	for (its=1;its<=ITMAX;its++) {
		*iter=its;
		linmin(p,xi,n,fret,func);
		if (2.0*fabs(*fret-fp) <= ftol*(fabs(*fret)+fabs(fp)+EPS)) {
			FREEALL
			return;
		}
		fp=(*func)(p);
		(*dfunc)(p,xi);
		dgg=gg=0.0;
		for (j=1;j<=n;j++) {
			gg += g[j]*g[j];
			dgg += (xi[j]+g[j])*xi[j];
		}
		if (gg == 0.0) {
			FREEALL
			return;
		}
		gam=dgg/gg;
		for (j=1;j<=n;j++) {
			g[j] = -xi[j];
			xi[j]=h[j]=g[j]+gam*h[j];
		}
	}
	nrerror("Too many iterations in frprmn");
	return;
}
#undef ITMAX
#undef EPS
#undef FREEALL

float medianim(chiptype im,int img,float min,float max) {
   float*list,xx[3]={0,0,0};
   int n=0,i,x,y,b,nm,np;

   list=(float*)calloc(10000,FLOATSIZE);
   b=(dataim[img].X*dataim[img].Y)/10000;
   if (b<1) b=1;
   for (y=0;y<dataim[img].Y && n<10000;y+=b) for (x=0;x<dataim[img].X && n<10000;x+=b) if (im[y][x]>min && im[y][x]<max) {
      list[n]=im[y][x];
      if (n==0 || list[n]<xx[0]) xx[0]=list[n];
      if (n==0) xx[1]=list[n];
      if (n==0 || list[n]>xx[2]) xx[2]=list[n];
      n++;
   }
   while (1) {
      nm=np=0;
      for (i=0;i<n;i++) {
	 if (list[i]<xx[1]-0.000001) nm++;
	 else if (list[i]>xx[1]+0.000001) np++;
      }
      if (nm<=n/2 && np<=n/2) {
	 free(list);
	 return xx[1];
      }
      if (nm>np) xx[2]=xx[1];
      else xx[0]=xx[1];
      xx[1]=0.5*(xx[0]+xx[2]);
   }
   free(list);
   return 0.;
}

void subreference(int img) {
   int x,y,xx,yy,i;
   double tx,ty;
   float mlt,rbg,newMAX,newMIN,fim=0.,fdif=0.;
   float s,ss,c,ndof,sh,ssky,mx,my,cs,cm,tmp;
   double E=0.,R=0.;
   float *p;
   chiptype tmpres;
   int n;

   if (Timg==Nimg) {
      printf("External reference image required for difference photometry\n");
      return;
   }
   // Set RPSFsub equal to smaller of RPSF[Nimg] and RPSF[img]
   RPSFsub=RPSF[img];
   if (RPSFsub>RPSF[Nimg]) RPSFsub=RPSF[Nimg];
   // Set res[] and tmpres[] equal to reference image, shifted to match img
   tmpres=allocchip(dataim[img].X,dataim[img].Y);
   for (y=0;y<dataim[img].Y;y++) for (x=0;x<dataim[img].X;x++) {
      shift(img,x+0.5,y+0.5,&tx,&ty,-1);
      if (tx<0.5 || ty<0.5 || tx>dataim[Nimg].X-0.5 || ty>dataim[Nimg].Y-0.5 || data[Nimg][(int)ty][(int)tx]<=iDMIN[Nimg]) res[img][y][x]=tmpres[y][x]=iDMIN[Nimg]-1.;
      else if (data[Nimg][(int)ty][(int)tx]>=FSat*iDMAX[Nimg]) res[img][y][x]=tmpres[y][x]=iDMAX[Nimg]+1.;
      else res[img][y][x]=tmpres[y][x]=imginterp(Nimg,tx,ty);
   }
   // Set convv equal to img PSF, at center of chip; fim = sum of PSF
   convv=matrix(-RPSFsub,RPSFsub,-RPSFsub,RPSFsub);
   calc1psf(img,dataim[img].X/2+0.5,dataim[img].Y/2+0.5,RPSFsub,apsf[img][0][0],apsf[img][1][0],apsf[img][2][0],1,0);
   for (y=-RPSFsub;y<=RPSFsub;y++) for (x=-RPSFsub;x<=RPSFsub;x++) {
      convv[y][x]=psf[y][x];
      fim+=psf[y][x];
   }
   // Set psf equal to reference PSF at center of chip
   calc1psf(Nimg,dataim[Nimg].X/2+0.5,dataim[Nimg].Y/2+0.5,RPSFsub,apsf[Nimg][0][0],apsf[Nimg][1][0],apsf[Nimg][2][0],1,0);
   n=2*RPSFsub+1;
   // Add PSF offset image (if defined)
   if (refpsfimg) for (y=0;y<n;y++) for (x=0;x<n;x++) psf[y-RPSFsub][x-RPSFsub]+=refpsfimg[y][x];

   // Compute kernel that convolves with reference PSF to get img PSF
   p=vector(1,n*n);
   for (y=1;y<=n*n;y++) p[y]=0.;
   frprmn(p,n*n,1.e-3,&y,&tmp,&convfunc,&dconvfunc);
   // fdif will equal sum of convolved (ref*kernel) PSF divided by sum of kernel
   for (y=-RPSFsub;y<=RPSFsub;y++) for (x=-RPSFsub;x<=RPSFsub;x++) {
      // ty = convolution of kernel and reference PSF
      // tx = sum of kernel
      tx=ty=0.;
      for (yy=-RPSFsub;yy<=RPSFsub;yy++) if (y+yy>=-RPSFsub && y+yy<RPSFsub) for (xx=-RPSFsub;xx<=RPSFsub;xx++) if (x+xx>=-RPSFsub && x+xx<RPSFsub) {
	 ty+=psf[y+yy][x+xx]*p[1+RPSFsub-xx+(RPSFsub-yy)*n];
	 tx+=p[1+RPSFsub-xx+(RPSFsub-yy)*n];
      }
      if (tx) fdif+=ty/tx;
   }
   /*
   tx=0.; for (y=RPSFsub;y>=-RPSFsub;y--) for (x=-RPSFsub;x<=RPSFsub;x++) tx+=psf[y][x]; printf("Template PSF sum: %f\n",tx);
   printf("Image PSF sum: %f\nConv. ref PSF sum: %f\n",fim,fdif);
   tx=0.;
   for (y=RPSFsub;y>=-RPSFsub;y--) {
      for (x=-RPSFsub;x<=RPSFsub;x++) tx+=p[1+x+RPSFsub+(y+RPSFsub)*n];
      for (x=-RPSFsub;x<=RPSFsub;x++) printf(" %3d",(int)(10000*p[1+x+RPSFsub+(y+RPSFsub)*n]));
      printf("\n");
   }
   printf("kernel sum: %f\n",tx);
   */
   free_matrix(convv,-RPSFsub,RPSFsub,-RPSFsub,RPSFsub);
   // tmpres equals convolved shifted reference image.  Note that "res" right now is the shifted reference image
   for (y=0;y<dataim[img].Y;y++) for (x=0;x<dataim[img].X;x++) {
      if (res[img][y][x]>iDMIN[Nimg] && res[img][y][x]<iDMAX[Nimg]) {
	 tx=ty=0.;
	 for (yy=-RPSFsub;yy<=RPSFsub;yy++) if (y+yy>=0 && y+yy<dataim[img].Y) for (xx=-RPSFsub;xx<=RPSFsub;xx++) if (x+xx>=0 && x+xx<dataim[img].X && res[img][y+yy][x+xx]>iDMIN[Nimg] && res[img][y+yy][x+xx]<iDMAX[Nimg]) {
	    ty+=res[img][y+yy][x+xx]*p[1+RPSFsub-xx+(RPSFsub-yy)*n];
	    tx+=p[1+RPSFsub-xx+(RPSFsub-yy)*n];
	 }
	 if (tx) tmpres[y][x]=ty/tx;
	 else tmpres[y][x]=iDMIN[Nimg]-1.;
      }

      // invalidate any bad points
      if (data[img][y][x]<iDMAX[img] && (tmpres[y][x]<=iDMIN[Nimg] || tmpres[y][x]>=FSat*iDMAX[Nimg])) data[img][y][x] = iDMIN[img]-1;
   }
   free_vector(p,1,n*n);
   // E = sum(imgcounts*sqrt(refcts))
   // R = sum(refcts^1.5)
   E=R=0.;
   for (i=0;i<Nstars;i++) if (stars[i].type==1 && stars[i].icm[img]>0 && refcts[i]>0) {
      shift(img,stars[i].x,stars[i].y,&tx,&ty,1);
      x=(int)(tx+1)-1; y=(int)(ty+1)-1;
      if (ppixfOK(img,x,y)) {
	 E+=stars[i].is[img]/stars[i].icm[img];
	 R+=refcts[i];
      }
   }
   // refmult = average ratio of reference/image counts
   refmult[img]=R/E;
   // mlt = average ratio of (imageCounts*imgPSFsum)/(referenceCounts*convolvedPSFsum/kernelSum) ~= imageCountsInRPSF / referenceCountsInRPSF
   mlt=E*fim/R/fdif;
   // rbg = rough background of tmpres, using median
   rbg=medianim(tmpres,img,iDMIN[Nimg],iDMAX[Nimg]);
   if (R==0. || mlt<0.) {
      printf("Math error in subreference\n");
      return;
   }
   //printf("image multiple: %f %f\n",mlt,rbg);
   IMDIFF=1;
   while (1) {
      // compute DMIN and DMAX and create subtracted reference in res[]
      newMIN=iDMIN[img]-mlt*(iDMAX[Nimg]-rbg)+mlt*mlt*(iDMIN[Nimg]+iRN[Nimg]);
      newMAX=iDMAX[img]-mlt*(iDMIN[Nimg]-rbg)+mlt*mlt*(iDMAX[Nimg]+iRN[Nimg]);
      for (y=0;y<dataim[img].Y;y++) for (x=0;x<dataim[img].X;x++) if (datafOK(img,x,y) && tmpres[y][x]>iDMIN[Nimg] && tmpres[y][x]<FSat*iDMAX[Nimg]) res[img][y][x]=data[img][y][x]-mlt*(tmpres[y][x]-rbg);
      else if (data[img][y][x]<iDMAX[img]*FSat) data[img][y][x]=res[img][y][x]=newMIN-1.;
      else data[img][y][x]=res[img][y][x]=newMAX+1.;
      // E = sum(imgcounts*sqrt(refcts)) -- unchanged from previous
      // R = sum of subtracted photometry * sqrt(refcts) -- hopefully zero
      E=R=0.;
      for (i=0;i<Nstars;i++) if (stars[i].type==1 && stars[i].icm[img]>0 && refcts[i]>0) {
	 clear_sky();
	 shift(img,stars[i].x,stars[i].y,&tx,&ty,1);
	 x=(int)(tx+1)-1; y=(int)(ty+1)-1;
	 if (ppixfOK(img,x,y)) {
	    eval1(img,stars[i].x,stars[i].y,stars[i].pa[img],stars[i].pb[img],stars[i].pc[img],&s,&ss,&c,&ndof,&sh,&ssky,&mx,&my,&cs,&cm,1,0,i,1);
	    E+=stars[i].is[img]/stars[i].icm[img];
	    R+=s/cm;
	 }
      }
      refmult[img]*=1.-R/E;
      mlt/=1.-R/E;
      if (R>=E) {
	 printf("**Fatal Math error in subreference**\n");
	 return;
      }
      //printf("%f %f\n",mlt,1./(1.-R/E));
      if (fabs(R/E)<=0.001) break;
   }
   // increase noise to sqrt( imgNoise^2 + mlt^2*refNoise^2 )
   for (y=0;y<dataim[img].Y;y++) for (x=0;x<dataim[img].X;x++) if (datafOK(img,x,y)) {
      if (data[Nimg][y][x]>0.) data[img][y][x]+=mlt*mlt*data[Nimg][y][x];
      data[img][y][x]+=mlt*mlt*iRN[Nimg];
   }
   iDMIN[img]=newMIN;
   iDMAX[img]=newMAX;
   for (x=0;x<Nstars;x++) stars[x].s=stars[x].is[img]=0.;
   freechip(tmpres,dataim[img].X,dataim[img].Y);
   printf("Subtraction made; scale factor=%6.4f\n",mlt);
   fflush(stdout);
   return;
}

int solve(int ext,int fld,int it,int show) {
   int i,reit=0,c;

   setindx(-1);
   for (i=0;i<Nstars;i++) indx[(int)(stars[i].y*SubResRef)][(int)(stars[i].x*SubResRef)]=i;
   for (i=0;i<Nstars;i++) if (stars[i].flag) {
      reit=1;
      imgadd(i);
      indx[(int)(stars[i].y*SubResRef)][(int)(stars[i].x*SubResRef)]=-1;
      phot((int)stars[i].x,(int)stars[i].y,i,it);
      imgsub(i);
   }
   c=clean(it,show);
   if (it<3 || (c && it<MaxIT/2)) return 2;
   if (fitpsf) {
      calcpsf(ext,fld);
      clean(it,show);
      fitpsf=0;
      if (WARMSTART==2) {
	 for (i=0;i<Nstars;i++) stars[i].flag=1;
	 solve(ext,fld,it,show);
	 for (i=0;i<Nstars;i++) stars[i].flag=1;
	 solve(ext,fld,it,show);
	 for (i=0;i<Nimg;i++) subreference(i);
      }
      for (i=0;i<Nstars;i++) stars[i].flag=1;
      return 2;
   }
   if (reit) return 1;
   return 0;
}

void clearsnmap(void) {
   memset(snmap[0],0,FLOATSIZE*X*Y*SubResRef*SubResRef);
   memset(snmmap[0],0,FLOATSIZE*X*Y*SubResRef*SubResRef);
   memset(snmapflag[0],0,X*Y*SubResRef*SubResRef);
   return;
}

void clearsnmapflag(void) {
   memset(snmapflag[0],0,X*Y*SubResRef*SubResRef);
   return;
}

void setsnmap(int IMG,int x0,int x1,int y0,int y1) {
   int x,y,img;
   double fx,fy;
   float s,ss,c,ndof,sh,ssky,mx,my,cs,cm,pa,pb,pc;

   for (img=0;img<Nimg || (img==Nimg && IMG==Nimg);img++) if (IMG<0 || img==IMG) {
      shift(img,0.5*(XMIN+XMAX),0.5*(YMIN+YMAX),&fx,&fy,1);
      pa=apsf[img][0][0]+(apsf[img][0][3]+apsf[img][0][4])/12.;
      pb=apsf[img][1][0]+(apsf[img][1][3]+apsf[img][1][4])/12.;
      pc=apsf[img][2][0]+(apsf[img][2][3]+apsf[img][2][4])/12.;
      calc1psf(img,0.5+(int)fx,0.5+(int)fy,rphot[img],pa,pb,pc,1,1);
      for (y=y0*SubResRef;y<=y1*SubResRef;y++) if (y>=YMIN*SubResRef && y<YMAX*SubResRef) for (x=x0*SubResRef;x<=x1*SubResRef;x++) if (x>=XMIN*SubResRef && x<XMAX*SubResRef && snmapflag[y][x]<2) {
	 clear_sky();
	 eval1(img,(x+0.5)/(double)SubResRef,(y+0.5)/(double)SubResRef,pa,pb,pc,&s,&ss,&c,&ndof,&sh,&ssky,&mx,&my,&cs,&cm,1,1,-1,0);
	 if (snmapflag[y][x]==0) {
	    snmap[y][x]=snmmap[y][x]=0.;
	    snmapflag[y][x]=1;
	 }
	 if (ss>0) {
	    snmap[y][x]+=s*fabs(s)/(ss*ss);
	    snmmap[y][x]+=s*fabs(s)/(ss*ss*(c*c+0.25));
	 }
      }
   }
   for (y=y0*SubResRef;y<=y1*SubResRef;y++) if (y>=YMIN*SubResRef && y<YMAX*SubResRef) for (x=x0*SubResRef;x<=x1*SubResRef;x++) if (x>=XMIN*SubResRef && x<XMAX*SubResRef && snmapflag[y][x]<2) {
      if (snmap[y][x]<=0) snmap[y][x]=0.;
      else snmap[y][x]=sqrt(snmap[y][x]);
      if (snmmap[y][x]>0.) snmmap[y][x]=sqrt(snmmap[y][x]);
      else snmmap[y][x]=-sqrt(-snmmap[y][x]);
      snmapflag[y][x]=2;
   }
   return;
}

Inline float peaksn(int img,int x,int y) {
   if (img<0) {
      if (x<=0 || y<=0 || x>=X*SubResRef-1 || y>=Y*SubResRef-1) return 0.;
   }
   else if (x<=0 || y<=0 || x>=dataim[img].X*SubResRef-1 || y>=dataim[img].Y*SubResRef-1) return 0.;
   if (snmmap[y][x]<snmmap[y-1][x-1] || snmmap[y][x]<snmmap[y-1][x] || snmmap[y][x]<snmmap[y-1][x+1] || snmmap[y][x]<snmmap[y][x+1] || snmmap[y][x]<snmmap[y+1][x+1] || snmmap[y][x]<snmmap[y+1][x] || snmmap[y][x]<snmmap[y+1][x-1] || snmmap[y][x]<snmmap[y][x-1]) return 0.;
   return snmap[y][x];
}

Inline float qpeak(int x,int y,int strict) {
   if (x<=0 || y<=0 || x>=X*SubResRef-1 || y>=Y*SubResRef-1) return 0.;
   if (strict && ran[y][x]) return 0;
   return peaksn(-1,x,y);
}

float alphot(int img,double x0,double y0,double*x,double*y,int finding,int strictpos) {
   float s,ss,sky,chi,sh;
   int i,j,cl;
   static float *pa,*pb,*pc;
   static int first=1;

   if (first) {
      pa=(float*)calloc(Timg,FLOATSIZE);
      pb=(float*)calloc(Timg,FLOATSIZE);
      pc=(float*)calloc(Timg,FLOATSIZE);
      if (!pa || !pb || !pc) merr();
      first=0;
   }
   if (finding!=0) {
      if (Nstars>=MAXNSTARS) return 0;
      centroid(img,(int)x0,(int)y0,x,y);
      cl=1;
      s=ss=0.0;
      getpsfpars(img,*x,*y,pa,pb,pc);
      if (photsearch(img,x,y,pa,pb,pc,&s,&ss,&chi,&sh,&sky,NULL,NULL,NULL,&cl,1,-1)<=0) return 0;
      if (finding==1) {
	 if (cl>1 || s<ss*10 || s<=0 || fabs(sh)>0.3) return 0;
	 getpsfpars(img,*x,*y,pa,pb,pc);
	 if (photsearch(img,x,y,pa,pb,pc,&s,&ss,&chi,&sh,&sky,NULL,NULL,NULL,&cl,-1,-1)<=0) return 0;
	 if (cl>1 || s<ss*10 || s<=0 || fabs(sh)>0.3) return 0;
      }
      else {
	 if (cl>1 || s<ss*5 || s<=0) return 0;
	 getpsfpars(img,*x,*y,pa,pb,pc);
	 if (photsearch(img,x,y,pa,pb,pc,&s,&ss,&chi,&sh,&sky,NULL,NULL,NULL,&cl,-1,-1)<=0) return 0;
	 if (cl>1 || s<ss*5 || s<=0) return 0;
      }
      if (*x<0 || *x>=X || *y<0 || *y>=Y || indx[(int)(*y*SubResRef)][(int)(*x*SubResRef)]>=0) return 0;
      indx[(int)(*y*SubResRef)][(int)(*x*SubResRef)]=Nstars;
      for (j=-(int)RCombine;j<=(int)RCombine;j++) for (i=-(int)RCombine;i<=(int)RCombine;i++) if ((int)(*y*SubResRef)+j>=0 && (int)(*y*SubResRef)+j<dataim[img].Y && (int)(*x*SubResRef)+i>=0 && (int)(*x*SubResRef)+i<dataim[img].X && i*i+j*j<(RCombine-0.5)*(RCombine-0.5)) ran[(int)(*y*SubResRef)+j][(int)(*x*SubResRef)+i]=1;
      stars[Nstars].x=*x;
      stars[Nstars].y=*y;
      Nstars++;
      return 1;
   }
   else {
      *x=x0;
      *y=y0;
      cl=1;
      s=ss=0.0;
      getpsfpars(img,*x,*y,pa,pb,pc);
      if (photsearch(img,x,y,pa,pb,pc,&s,&ss,&chi,&sh,&sky,NULL,NULL,NULL,&cl,-1,-1)<=0) return 0;
      if (cl>1 || s<ss*5 || s<=0 || chi>10. || fabs(sh)>0.4) return 0;
      if (strictpos) {
	 double sx0,sy0,sx1,sy1;
	 shift(img,x0,y0,&sx0,&sy0,1);
	 shift(img,*x,*y,&sx1,&sy1,1);
	 if ((sx1-sx0)*(sx1-sx0)+(sy1-sy0)*(sy1-sy0)>1+RAper[img]*RAper[img]*0.5) return 0;
      }
      return s/ss;
   }
}

typedef float dptype[5];
int NFREE=5,ALIMG=0,Nalign;
dptype *Alist;
static double Align45scale[21]={1,
			       1e3,1e3,
			       1e6,1e6,1e6,
			       1e9,1e9,1e9,1e9,
			       1e12,1e12,1e12,1e12,1e12,
			       1e15,1e15,1e15,1e15,1e15,1e15};

float aminfunc(float*c) {
   double x,y;
   float rval=0,r2;
   int i;

   if (Align<4) {
      for (i=0;i<=Align;i++) dpos[ALIMG][i]=c[i+1];
      if (Rotate) dpos[ALIMG][4]=c[Align+2];
   }
   else {
      for (i=0;i<NFREE;i++) dpos[ALIMG][i]=c[i+1]/Align45scale[i%10];
   }
   for (i=0;i<Nalign;i++) {
      shift(ALIMG,Alist[i][0],Alist[i][1],&x,&y,1);
      r2=((x-Alist[i][2])*(x-Alist[i][2])+(y-Alist[i][3])*(y-Alist[i][3]))*25;
      //rval+=r2;
      rval+=log(1+0.5*r2);
   }
   return rval/Nalign;
}

float aminfunc4(float *c) {
   int i,x0,y0,N=0,xstep,ystep;
   double x,y;
   float rval=0;

   xstep = dataim[ALIMG].X/16;
   ystep = dataim[ALIMG].Y/16;
   for (i=0;i<20;i++) dpos[ALIMG][i+20] = c[i+1]/Align45scale[i%10];
   for (x0=0;x0<=dataim[ALIMG].X;x0+=xstep) {
      for (y0=0;y0<=dataim[ALIMG].Y;y0+=ystep) {
	 shift(ALIMG,x0,y0,&x,&y,-1);
	 shift(ALIMG,x,y,&x,&y,1);
	 rval += ((x0-x)*(x0-x)+(y0-y)*(y0-y))*25;
	 N++;
      }
   }
   return rval/N;
}

float aminfuncR2I(float *c) {
   int i,x0,y0,N=0,xstep,ystep;
   double x,y;
   float rval=0;

   xstep = dataim[ALIMG].X/16;
   ystep = dataim[ALIMG].Y/16;
   for (i=0;i<20;i++) ref2img[ALIMG][i+20] = c[i+1]/Align45scale[i%10];
   for (x0=0;x0<=dataim[ALIMG].X;x0+=xstep) {
      for (y0=0;y0<=dataim[ALIMG].Y;y0+=ystep) {
	 shift(ALIMG,x0,y0,&x,&y,-2);
	 shift(ALIMG,x,y,&x,&y,2);
	 rval += ((x0-x)*(x0-x)+(y0-y)*(y0-y))*25;
	 N++;
      }
   }
   return rval/N;
}

//powell() and affiliated code from NR;
int ncom;
float *pcom,*xicom,(*nrfunc)(float []);

float f1dim(float x)
{
	int j;
	float f,*xt;

	xt=vector(1,ncom);
	for (j=1;j<=ncom;j++) xt[j]=pcom[j]+x*xicom[j];
	f=(*nrfunc)(xt);
	free_vector(xt,1,ncom);
	return f;
}

#define GOLD 1.618034
#define GLIMIT 100.0
#define TINY 1.0e-20
#define SHFT(a,b,c,d) (a)=(b);(b)=(c);(c)=(d);
void mnbrak(float *ax, float *bx, float *cx, float *fa, float *fb, float *fc,
	float (*func)(float))
{
	float ulim,u,r,q,fu,dum;

	*fa=(*func)(*ax);
	*fb=(*func)(*bx);
	if (*fb > *fa) {
		SHFT(dum,*ax,*bx,dum)
		SHFT(dum,*fb,*fa,dum)
	}
	*cx=(*bx)+GOLD*(*bx-*ax);
	*fc=(*func)(*cx);
	while (*fb > *fc) {
		r=(*bx-*ax)*(*fb-*fc);
		q=(*bx-*cx)*(*fb-*fa);
		u=(*bx)-((*bx-*cx)*q-(*bx-*ax)*r)/
			(2.0*SIGN(FMAX(fabs(q-r),TINY),q-r));
		ulim=(*bx)+GLIMIT*(*cx-*bx);
		if ((*bx-u)*(u-*cx) > 0.0) {
			fu=(*func)(u);
			if (fu < *fc) {
				*ax=(*bx);
				*bx=u;
				*fa=(*fb);
				*fb=fu;
				return;
			} else if (fu > *fb) {
				*cx=u;
				*fc=fu;
				return;
			}
			u=(*cx)+GOLD*(*cx-*bx);
			fu=(*func)(u);
		} else if ((*cx-u)*(u-ulim) > 0.0) {
			fu=(*func)(u);
			if (fu < *fc) {
				SHFT(*bx,*cx,u,*cx+GOLD*(*cx-*bx))
				SHFT(*fb,*fc,fu,(*func)(u))
			}
		} else if ((u-ulim)*(ulim-*cx) >= 0.0) {
			u=ulim;
			fu=(*func)(u);
		} else {
			u=(*cx)+GOLD*(*cx-*bx);
			fu=(*func)(u);
		}
		SHFT(*ax,*bx,*cx,u)
		SHFT(*fa,*fb,*fc,fu)
	}
	return;
}
#undef GOLD
#undef GLIMIT
#undef TINY
#undef SHFT

#define ITMAX 100
#define CGOLD 0.3819660
#define ZEPS 1.0e-10
#define SHFT(a,b,c,d) (a)=(b);(b)=(c);(c)=(d);
float brent(float ax, float bx, float cx, float (*f)(float), float tol,
	float *xmin)
{
	int iter;
	float a,b,d=0,etemp,fu,fv,fw,fx,p,q,r,tol1,tol2,u,v,w,x,xm;
	float e=0.0;

	a=(ax < cx ? ax : cx);
	b=(ax > cx ? ax : cx);
	x=w=v=bx;
	fw=fv=fx=(*f)(x);
	for (iter=1;iter<=ITMAX;iter++) {
		xm=0.5*(a+b);
		tol2=2.0*(tol1=tol*fabs(x)+ZEPS);
		if (fabs(x-xm) <= (tol2-0.5*(b-a))) {
			*xmin=x;
			return fx;
		}
		if (fabs(e) > tol1) {
			r=(x-w)*(fx-fv);
			q=(x-v)*(fx-fw);
			p=(x-v)*q-(x-w)*r;
			q=2.0*(q-r);
			if (q > 0.0) p = -p;
			q=fabs(q);
			etemp=e;
			e=d;
			if (fabs(p) >= fabs(0.5*q*etemp) || p <= q*(a-x) || p >= q*(b-x))
				d=CGOLD*(e=(x >= xm ? a-x : b-x));
			else {
				d=p/q;
				u=x+d;
				if (u-a < tol2 || b-u < tol2)
					d=SIGN(tol1,xm-x);
			}
		} else {
			d=CGOLD*(e=(x >= xm ? a-x : b-x));
		}
		u=(fabs(d) >= tol1 ? x+d : x+SIGN(tol1,d));
		fu=(*f)(u);
		if (fu <= fx) {
			if (u >= x) a=x; else b=x;
			SHFT(v,w,x,u)
			SHFT(fv,fw,fx,fu)
		} else {
			if (u < x) a=u; else b=u;
			if (fu <= fw || w == x) {
				v=w;
				w=u;
				fv=fw;
				fw=fu;
			} else if (fu <= fv || v == x || v == w) {
				v=u;
				fv=fu;
			}
		}
	}
	nrerror("Too many iterations in brent");
	*xmin=x;
	return fx;
}
#undef ITMAX
#undef CGOLD
#undef ZEPS
#undef SHFT

#define TOL 2.0e-4
void linmin(float p[], float xi[], int n, float *fret, float (*func)(float []))
{
	int j;
	float xx,xmin,fx,fb,fa,bx,ax;

	ncom=n;
	pcom=vector(1,n);
	xicom=vector(1,n);
	nrfunc=func;
	for (j=1;j<=n;j++) {
		pcom[j]=p[j];
		xicom[j]=xi[j];
	}
	ax=0.0;
	xx=1.0;
	mnbrak(&ax,&xx,&bx,&fa,&fx,&fb,f1dim);
	*fret=brent(ax,xx,bx,f1dim,TOL,&xmin);
	for (j=1;j<=n;j++) {
		xi[j] *= xmin;
		p[j] += xi[j];
	}
	free_vector(xicom,1,n);
	free_vector(pcom,1,n);
	return;
}
#undef TOL

#define ITMAX 10000
void powell(float p[], float **xi, int n, float ftol, int *iter, float *fret,
	float (*func)(float []))
{
	int i,ibig,j;
	float del,fp,fptt,t,*pt,*ptt,*xit;

	pt=vector(1,n);
	ptt=vector(1,n);
	xit=vector(1,n);
	*fret=(*func)(p);
	for (j=1;j<=n;j++) pt[j]=p[j];
	for (*iter=1;;++(*iter)) {
		fp=(*fret);
		ibig=0;
		del=0.0;
		for (i=1;i<=n;i++) {
			for (j=1;j<=n;j++) xit[j]=xi[j][i];
			fptt=(*fret);
			linmin(p,xit,n,fret,func);
			if (fabs(fptt-(*fret)) > del) {
				del=fabs(fptt-(*fret));
				ibig=i;
			}
		}
		if (2.0*fabs(fp-(*fret)) <= ftol*(fabs(fp)+fabs(*fret))) {
			free_vector(xit,1,n);
			free_vector(ptt,1,n);
			free_vector(pt,1,n);
			return;
		}
		if (*iter == ITMAX) nrerror("powell exceeding maximum iterations.");
		for (j=1;j<=n;j++) {
			ptt[j]=2.0*p[j]-pt[j];
			xit[j]=p[j]-pt[j];
			pt[j]=p[j];
		}
		fptt=(*func)(ptt);
		if (fptt < fp) {
		   t=2.0*(fp-2.0*(*fret)+fptt)*(fp-(*fret)-del)*(fp-(*fret)-del)-del*(fp-fptt)*(fp-fptt);
			if (t < 0.0) {
				linmin(p,xit,n,fret,func);
				for (j=1;j<=n;j++) {
					xi[j][ibig]=xi[j][n];
					xi[j][n]=xit[j];
				}
			}
		}
	}
	return;
}
#undef ITMAX
//end of NR stuff;

// Compute distortion correction from WCS
void headerWCS(int img,float* imgfwd) {
   char str[81];
   int sip[2],ord,i;
   for (i=0;i<40;i++) imgfwd[i]=0.0;

   // parse coordinate types
   strcpy(str,getcardval(dataim+img,"CTYPE1",1));
   if (!strcmp(str,"RA---TAN")) sip[0]=0;
   else if (!strcmp(str,"RA---TAN-SIP")) sip[0]=1;
   else {printf("Error: CTYPE1 of %s not understood\n",str); exit(-1);}
   strcpy(str,getcardval(dataim+img,"CTYPE2",1));
   if (!strcmp(str,"DEC--TAN")) sip[1]=0;
   else if (!strcmp(str,"DEC--TAN-SIP")) sip[1]=1;
   else {printf("Error: CTYPE2 of %s not understood\n",str); exit(-1);}

   // define reference pixel
   wcsref[img][0] = atof(getcardval(dataim+img,"CRPIX1",1))-0.5;
   wcsref[img][1] = atof(getcardval(dataim+img,"CRPIX2",1))-0.5;
   wcsref[img][2] = atof(getcardval(dataim+img,"CRVAL1",1));
   wcsref[img][3] = atof(getcardval(dataim+img,"CRVAL2",1));

   // parse first transformation
   imgfwd[0] = atof(getcardval(dataim+img,"CD1_1",1));
   imgfwd[1] = atof(getcardval(dataim+img,"CD1_2",1));
   if (sip[0]) {
      ord = atoi(getcardval(dataim+img,"A_ORDER",1));
      if (ord>1) {
	 imgfwd[2] = atof(getcardval(dataim+img,"A_2_0",1));
	 imgfwd[3] = atof(getcardval(dataim+img,"A_1_1",1));
	 imgfwd[4] = atof(getcardval(dataim+img,"A_0_2",1));
      }
      if (ord>2) {
	 imgfwd[5] = atof(getcardval(dataim+img,"A_3_0",1));
	 imgfwd[6] = atof(getcardval(dataim+img,"A_2_1",1));
	 imgfwd[7] = atof(getcardval(dataim+img,"A_1_2",1));
	 imgfwd[8] = atof(getcardval(dataim+img,"A_0_3",1));
      }
      if (ord>3) {
	 imgfwd[9] = atof(getcardval(dataim+img,"A_4_0",1));
	 imgfwd[10] = atof(getcardval(dataim+img,"A_3_1",1));
	 imgfwd[11] = atof(getcardval(dataim+img,"A_2_2",1));
	 imgfwd[12] = atof(getcardval(dataim+img,"A_1_3",1));
	 imgfwd[13] = atof(getcardval(dataim+img,"A_0_4",1));
      }
      if (ord>4) {
	 imgfwd[14] = atof(getcardval(dataim+img,"A_5_0",1));
	 imgfwd[15] = atof(getcardval(dataim+img,"A_4_1",1));
	 imgfwd[16] = atof(getcardval(dataim+img,"A_3_2",1));
	 imgfwd[17] = atof(getcardval(dataim+img,"A_2_3",1));
	 imgfwd[18] = atof(getcardval(dataim+img,"A_1_4",1));
	 imgfwd[19] = atof(getcardval(dataim+img,"A_0_5",1));
      }
   }

   // parse second transformation
   imgfwd[20] = atof(getcardval(dataim+img,"CD2_1",1));
   imgfwd[21] = atof(getcardval(dataim+img,"CD2_2",1));
   if (sip[1]) {
      ord = atoi(getcardval(dataim+img,"B_ORDER",1));
      if (ord>1) {
	 imgfwd[22] = atof(getcardval(dataim+img,"B_2_0",1));
	 imgfwd[23] = atof(getcardval(dataim+img,"B_1_1",1));
	 imgfwd[24] = atof(getcardval(dataim+img,"B_0_2",1));
      }
      if (ord>2) {
	 imgfwd[25] = atof(getcardval(dataim+img,"B_3_0",1));
	 imgfwd[26] = atof(getcardval(dataim+img,"B_2_1",1));
	 imgfwd[27] = atof(getcardval(dataim+img,"B_1_2",1));
	 imgfwd[28] = atof(getcardval(dataim+img,"B_0_3",1));
      }
      if (ord>3) {
	 imgfwd[29] = atof(getcardval(dataim+img,"B_4_0",1));
	 imgfwd[30] = atof(getcardval(dataim+img,"B_3_1",1));
	 imgfwd[31] = atof(getcardval(dataim+img,"B_2_2",1));
	 imgfwd[32] = atof(getcardval(dataim+img,"B_1_3",1));
	 imgfwd[33] = atof(getcardval(dataim+img,"B_0_4",1));
      }
      if (ord>4) {
	 imgfwd[34] = atof(getcardval(dataim+img,"B_5_0",1));
	 imgfwd[35] = atof(getcardval(dataim+img,"B_4_1",1));
	 imgfwd[36] = atof(getcardval(dataim+img,"B_3_2",1));
	 imgfwd[37] = atof(getcardval(dataim+img,"B_2_3",1));
	 imgfwd[38] = atof(getcardval(dataim+img,"B_1_4",1));
	 imgfwd[39] = atof(getcardval(dataim+img,"B_0_5",1));
      }
   }
}

int REFIMG;
double refpix[2]; // position of standard pixel on reference image
double refra[2];  // RA/DEC of standard pixel

float aminfuncWCSref(float *c) {
   int i,x0,y0,N=0,xstep,ystep;
   double x,y;
   float rval=0;

   xstep = dataim[REFIMG].X/16;
   ystep = dataim[REFIMG].Y/16;
   for (i=0;i<40;i++) wcs[REFIMG][i+40] = c[i+1]/Align45scale[i%20+1];
   for (x0=0;x0<=dataim[REFIMG].X;x0+=xstep) {
      for (y0=0;y0<=dataim[REFIMG].Y;y0+=ystep) {
	 shift(REFIMG,x0,y0,&x,&y,-3);  // REF pixel to REF tangent
	 shift(REFIMG,x,y,&x,&y,3);     // REF tangent to REF pixel
	 rval += ((x0-x)*(x0-x)+(y0-y)*(y0-y))*25;
	 N++;
      }
   }
   return rval/N;
}

float aminfuncWCSrev(float *c) {
   int i,x0,y0,N=0,xstep,ystep;
   double x,y;
   float rval=0;

   xstep = dataim[ALIMG].X/16;
   ystep = dataim[ALIMG].Y/16;
   for (i=0;i<40;i++) wcs[ALIMG][i+40] = c[i+1]/Align45scale[i%20+1];
   for (x0=0;x0<=dataim[ALIMG].X;x0+=xstep) {
      for (y0=0;y0<=dataim[ALIMG].Y;y0+=ystep) {
	 wcsref[ALIMG][2]=refra[0]; wcsref[ALIMG][3]=refra[1];
	 shift(ALIMG,x0,y0,&x,&y,-3); // IMG pixel to IMG tangent
	 fixRA(ALIMG,x,y,&x,-3);      // IMG tangent to RA/DEC
	 fixRA(REFIMG,x,y,&x,3);      // RA/DEC to REF tangent
	 shift(REFIMG,x,y,&x,&y,3);   // REF tangent to REF pixel
	 wcsref[ALIMG][2]=refpix[0]; wcsref[ALIMG][3]=refpix[1];
	 shift(ALIMG,x,y,&x,&y,3);    // refpix to imgpix
	 rval += ((x0-x)*(x0-x)+(y0-y)*(y0-y))*25;
	 N++;
      }
   }
   return rval/N;
}

float aminfuncWCSfwd(float *c) {
   int i,x0,y0,N=0,xstep,ystep;
   double x,y;
   float rval=0;

   xstep = dataim[ALIMG].X/16;
   ystep = dataim[ALIMG].Y/16;
   for (i=0;i<40;i++) wcs[ALIMG][i] = c[i+1]/Align45scale[i%20+1];
   for (x0=0;x0<=dataim[ALIMG].X;x0+=xstep) {
      for (y0=0;y0<=dataim[ALIMG].Y;y0+=ystep) {
	 shift(ALIMG,x0,y0,&x,&y,-3);    // imgpix to refpix
	 shift(ALIMG,x,y,&x,&y,3);       // refpix to imgpix
	 rval += ((x0-x)*(x0-x)+(y0-y)*(y0-y))*25;
	 N++;
      }
   }
   return rval/N;
}

void convertWCS2(void) {
   int img,i,j;
   double d;
   float p[40],**xi,rv;
   double m[2][2];
   xi=matrix(1,40,1,40);

   // compute reverse transform for reference image
   int refimg=0;
   if (Timg>Nimg) refimg=Nimg;
   REFIMG = refimg;

   // read header, populate imgfwd and imgpix
   headerWCS(refimg,wcs[refimg]);

   // solve for first-order reverse transformation
   d = wcs[refimg][0]*wcs[refimg][21] - wcs[refimg][1]*wcs[refimg][20];
   if (d==0.0) {
      printf("ERROR: determinant of WCS is zero\n");
      exit(-1);
   }
   memset(p,0,sizeof(p));
   p[0] = wcs[refimg][21]/d*Align45scale[1];
   p[1] = -wcs[refimg][1]/d*Align45scale[2];
   p[20] = -wcs[refimg][20]/d*Align45scale[1];
   p[21] = wcs[refimg][0]/d*Align45scale[2];
   printf("reference rmsi = %f, ",0.2*sqrt(aminfuncWCSref(p-1)));
   for (i=1;i<=40;i++) {
      for (j=1;j<=40;j++) xi[i][j]=0.;
      xi[i][i]=1.;
   }
   powell(p-1,xi,40,1.e-8,&i,&rv,&aminfuncWCSref);
   printf("rmsf = %f\n",0.2*sqrt(aminfuncWCSref(p-1)));
   for (i=0;i<40;i++) wcs[refimg][40+i] = p[i]/Align45scale[i%20+1];

   // compute transformations for non-reference images
   for (img=0;img<Nimg;img++) if (img!=refimg) {
      ALIMG = img;
      headerWCS(img,wcs[img]);
      // location of reference pixel in RA/DEC and in reference img pixels
      refra[0]=wcsref[img][2];
      refra[1]=wcsref[img][3];
      fixRA(refimg,refra[0],refra[1],&d,3); // convert RA/DEC to REF tangent
      shift(refimg,d,refra[1],refpix,refpix+1,3);

      // solve for reverse
      m[0][0] = wcs[refimg][40]*wcs[img][0] + wcs[refimg][41]*wcs[img][20];
      m[0][1] = wcs[refimg][40]*wcs[img][1] + wcs[refimg][41]*wcs[img][21];
      m[1][0] = wcs[refimg][60]*wcs[img][0] + wcs[refimg][61]*wcs[img][20];
      m[1][1] = wcs[refimg][60]*wcs[img][1] + wcs[refimg][61]*wcs[img][21];
      d = m[0][0]*m[1][1] - m[0][1]*m[1][0];
      if (d==0.0) {
	 printf("ERROR: determinant of WCS is zero\n");
	 exit(-1);
      }
      memset(p,0,sizeof(p));
      p[0] = m[1][1]/d*Align45scale[1];
      p[1] = -m[0][1]/d*Align45scale[2];
      p[20] = -m[1][0]/d*Align45scale[1];
      p[21] = m[0][0]/d*Align45scale[2];
      printf("reverse rmsi = %f, ",0.2*sqrt(aminfuncWCSrev(p-1)));
      for (i=1;i<=40;i++) {
	 for (j=1;j<=40;j++) xi[i][j]=0.;
	 xi[i][i]=1.;
      }
      powell(p-1,xi,40,1.e-8,&i,&rv,&aminfuncWCSrev);
      printf("rmsf = %f\n",0.2*sqrt(aminfuncWCSrev(p-1)));
      for (i=0;i<40;i++) wcs[img][40+i] = p[i]/Align45scale[i%20+1];

      // solve for forward
      wcsref[img][2]=refpix[0];
      wcsref[img][3]=refpix[1];
      d = wcs[img][40]*wcs[img][61] - wcs[img][41]*wcs[img][60];
      if (d==0.0) {
	 printf("ERROR: determinant of WCS is zero\n");
	 exit(-1);
      }
      memset(p,0,sizeof(p));
      p[0] = wcs[img][61]/d*Align45scale[1];
      p[1] = -wcs[img][41]/d*Align45scale[2];
      p[20] = -wcs[img][60]/d*Align45scale[1];
      p[21] = wcs[img][40]/d*Align45scale[2];
      printf("forward rmsi = %f, ",0.2*sqrt(aminfuncWCSfwd(p-1)));
      for (i=1;i<=40;i++) {
	 for (j=1;j<=40;j++) xi[i][j]=0.;
	 xi[i][i]=1.;
      }
      powell(p-1,xi,40,1.e-8,&i,&rv,&aminfuncWCSfwd);
      printf("rmsf = %f\n",0.2*sqrt(aminfuncWCSfwd(p-1)));
      for (i=0;i<40;i++) wcs[img][i] = p[i]/Align45scale[i%20+1];
   }

   free_matrix(xi,1,40,1,40);

   // set wcs[] for reference image
   memset(wcs[refimg],0,sizeof(wcstype));
   wcsref[refimg][0] = X*0.5;
   wcsref[refimg][1] = Y*0.5;
   wcsref[refimg][2] = X*0.5;
   wcsref[refimg][3] = Y*0.5;
   wcs[refimg][0] = 1.0;
   wcs[refimg][21] = 1.0;
   wcs[refimg][40] = 1.0;
   wcs[refimg][61] = 1.0;
}

// Delta to alignment info based on WCS 1st order solution
void convertWCS1(void) {
   int img;
   double x,y,xr,yr,dx,dy,t[2],s[2],tref,sref,refm[2][2],ct,st;
   float wcsdata[40];

   // read WCS data for reference image
   int refimg=0;
   if (Timg>Nimg) refimg=Nimg;
   headerWCS(refimg,wcsdata);
   shift(refimg,wcsref[refimg][0],wcsref[refimg][1],&x,&y,-3);
   shift(refimg,wcsref[refimg][0]+1.0,wcsref[refimg][1],&dx,&dy,-3);
   dx-=x; dy-=y;
   s[0] = sqrt(dx*dx+dy*dy)/sqrt(wcsdata[0]*wcsdata[0]+wcsdata[20]*wcsdata[20]);
   t[0] = atan2(dy,dx)-atan2(wcsdata[20],-wcsdata[0]);
   shift(refimg,wcsref[refimg][0],wcsref[refimg][1]+1.0,&dx,&dy,-3);
   dx-=x; dy-=y;
   s[1] = sqrt(dx*dx+dy*dy)/sqrt(wcsdata[1]*wcsdata[1]+wcsdata[21]*wcsdata[21]);
   t[1] = atan2(dy,dx)-atan2(wcsdata[21],-wcsdata[1]);
   sref = 0.5*(s[0]+s[1]);
   if (t[0]<t[1]-M_PI) t[0] += 2*M_PI;
   else if (t[0]>t[1]+M_PI) t[0] -= 2*M_PI;
   tref = 0.5*(t[0]+t[1]);

   // RA/DEC-to-X/Y
   x = wcsdata[0]*wcsdata[21] - wcsdata[1]*wcsdata[20];
   if (x==0.0) {
      printf("ERROR: determinant of reference WCS is zero\n");
      exit(-1);
   }
   refm[0][0] = wcsdata[21]/x;  // dxref/dRA
   refm[0][1] = -wcsdata[1]/x;  // dxref/dDec
   refm[1][0] = -wcsdata[20]/x; // dyref/dRA
   refm[1][1] = wcsdata[0]/x;   // dyref/dDec

   // compute transformations for non-reference images
   printf("Alignment estimates from WCS header\n");
   for (img=0;img<Nimg;img++) if (img!=refimg) {
      headerWCS(img,wcsdata);

      // compute scale and rotation
      shift(img,wcsref[img][0],wcsref[img][1],&x,&y,-3);
      shift(img,wcsref[img][0]+1.0,wcsref[img][1],&dx,&dy,-3);
      dx-=x; dy-=y;
      s[0] = sqrt(dx*dx+dy*dy)/sqrt(wcsdata[0]*wcsdata[0]+wcsdata[20]*wcsdata[20]);
      t[0] = atan2(dy,dx)-atan2(wcsdata[20],-wcsdata[0]);
      shift(img,wcsref[img][0],wcsref[img][1]+1.0,&dx,&dy,-3);
      dx-=x; dy-=y;
      s[1] = sqrt(dx*dx+dy*dy)/sqrt(wcsdata[1]*wcsdata[1]+wcsdata[21]*wcsdata[21]);
      t[1] = atan2(dy,dx)-atan2(wcsdata[21],-wcsdata[1]);

      s[0] = 0.5*(s[0]+s[1]);
      if (t[0]<t[1]-M_PI) t[0] += 2*M_PI;
      else if (t[0]>t[1]+M_PI) t[0] -= 2*M_PI;
      t[0] = 0.5*(t[0]+t[1]);
      if (t[0]<tref-M_PI) t[0] += 2*M_PI;
      else if (t[0]>tref+M_PI) t[0] -= 2*M_PI;
      s[0] /= sref;
      t[0] -= tref;

      dpos[img][2] *= s[0];
      dpos[img][4] += t[0]*180.0/M_PI;

      // compute shift

      // Target location of reference pixel on reference image
      fixRA(refimg,wcsref[img][2],wcsref[img][3],&x,3); // RA/DEC to REF tangent
      xr = wcsref[refimg][0] + refm[0][0]*(x-wcsref[refimg][2]) + refm[0][1]*(wcsref[img][3]-wcsref[refimg][3]);
      yr = wcsref[refimg][1] + refm[1][0]*(x-wcsref[refimg][2]) + refm[1][1]*(wcsref[img][3]-wcsref[refimg][3]);
      shift(img,wcsref[img][0],wcsref[img][1],&x,&y,-3); // Output from *shift (no rotation, scaling, or translation) plus (X,Y)*0.5
      x-=X*0.5;
      y-=Y*0.5;

      // Solve for dpos[img][0] and dpos[img][1]
      ct = cos(t[0]);
      st = sin(t[0]);
      xr=(xr-X*0.5)*s[0];
      yr=(yr-Y*0.5)*s[0];
      dpos[img][0] += x-ct*xr+st*yr;
      dpos[img][1] += y-st*xr-ct*yr;
      printf("image %d: shift %4.2f %4.2f, scale %8.6f, rotate %5.3f\n",img+1,x-ct*xr+st*yr,y-st*xr-ct*yr,s[0],t[0]*180.0/M_PI);
      if (!FakeStars[0] && VerboseData>0) fprintf(fverb,"WCS image %d: %f %f %f %f\n",img+1,x-ct*xr+st*yr,y-st*xr-ct*yr,s[0],t[0]*180.0/M_PI);

      // Verify
      //shift(img,wcsref[img][0],wcsref[img][1],&x,&y,-1);
   }
   fflush(stdout);
}

// Convert from parameterized alignment to plain cubic coefficients
double fitFwdAlign4(int img) {
   float p[20],**xi,rv,rval;
   int i,j;

   ALIMG=img;
   xi=matrix(1,20,1,20);
   for (i=0;i<20;i++) p[i] = dpos[img][i+20]*Align45scale[i%10];
   //printf("rmsi = %f, ",0.2*sqrt(aminfunc4(p-1)));
   for (i=1;i<=20;i++) {
      for (j=1;j<=20;j++) xi[i][j]=0.;
      xi[i][i]=1.;
   }
   powell(p-1,xi,20,1.e-8,&i,&rv,&aminfunc4);
   rval = aminfunc4(p-1);
   //printf("rmsf = %f\n",0.2*sqrt(aminfunc4(p-1)));
   for (i=0;i<20;i++) dpos[img][i+20] = p[i]/Align45scale[i%10];
   free_matrix(xi,1,20,1,20);

   return 0.2*sqrt(rval);
}

void convertAlign4(void) {
   dpostype tmpdpos;
   double c,s,d2,m,u,v,w;
   int img;

   for (img=0;img<Timg;img++) {
      memcpy(tmpdpos,dpos[img],sizeof(dpostype));
      d2 = 4.*dpos[img][3]/3./(X*X+Y*Y);
      c=cos(M_PI/180.*dpos[img][4]);
      s=sin(M_PI/180.*dpos[img][4]);
	 
      // reference -> frame (reverse), X calculation
      dpos[img][0] = tmpdpos[0]; // constant
      dpos[img][1] = c*tmpdpos[2]; // x term
      dpos[img][2] = -s*tmpdpos[2]; // y term
      dpos[img][3] = 0; // xx term
      dpos[img][4] = 0; // xy term
      dpos[img][5] = 0; // yy term
      dpos[img][6] = c*tmpdpos[2]*d2; // xxx term
      dpos[img][7] = -s*tmpdpos[2]*d2; // xxy term
      dpos[img][8] = c*tmpdpos[2]*d2; // xyy term
      dpos[img][9] = -s*tmpdpos[2]*d2; // yyy term
      // reference -> frame (reverse), Y calculation
      dpos[img][10] = tmpdpos[1]; // constant
      dpos[img][11] = s*tmpdpos[2]; // x term
      dpos[img][12] = c*tmpdpos[2]; // y term
      dpos[img][13] = 0; // xx term
      dpos[img][14] = 0; // xy term
      dpos[img][15] = 0; // yy term
      dpos[img][16] = s*tmpdpos[2]*d2; // xxx term
      dpos[img][17] = c*tmpdpos[2]*d2; // xxy term
      dpos[img][18] = s*tmpdpos[2]*d2; // xyy term
      dpos[img][19] = c*tmpdpos[2]*d2; // yyy term

      m = d2/tmpdpos[2]/tmpdpos[2]/tmpdpos[2];
      u=tmpdpos[0]*tmpdpos[0]+tmpdpos[1]*tmpdpos[1];
      v=c*tmpdpos[0]+s*tmpdpos[1];
      w=-s*tmpdpos[0]+c*tmpdpos[1];
      dpos[img][20] = -v*(1/tmpdpos[2]-u*m); // const
      dpos[img][21] = c/tmpdpos[2]-(2*v*tmpdpos[0]+c*u)*m; // *x
      dpos[img][22] = s/tmpdpos[2]-(2*v*tmpdpos[1]+s*u)*m; // *y
      dpos[img][23] = (s*tmpdpos[1]+3*c*tmpdpos[0])*m; // *xx
      dpos[img][24] = (2*s*tmpdpos[0]+2*c*tmpdpos[1])*m; // *xy
      dpos[img][25] = (c*tmpdpos[0]+3*s*tmpdpos[1])*m; // *yy
      dpos[img][26] = -c*m; // *xxx
      dpos[img][27] = -s*m; // *xxy
      dpos[img][28] = -c*m; // *xyy
      dpos[img][29] = -s*m; // *yyy
      dpos[img][30] = -w*(1/tmpdpos[2]-m*u); // const
      dpos[img][31] = -s/tmpdpos[2]-(2*w*tmpdpos[0]-s*u)*m; // *x
      dpos[img][32] = c/tmpdpos[2]-(2*w*tmpdpos[1]+c*u)*m; // *y
      dpos[img][33] = (w-2*s*tmpdpos[0])*m; // *xx
      dpos[img][34] = (2*c*tmpdpos[0]-2*s*tmpdpos[1])*m; // *xy
      dpos[img][35] = (w+2*c*tmpdpos[1])*m; // *yy
      dpos[img][36] = s*m; // *xxx
      dpos[img][37] = -c*m; // *xxy
      dpos[img][38] = s*m; // *xyy
      dpos[img][39] = -c*m; // *yyy

      // Refine solution
      fitFwdAlign4(img);
   }
}

void setRef2Img(void) {
   int img,i,j;
   float p[20],**xi,rv,d;

   xi=matrix(1,20,1,20);
   for (img=0;img<Timg;img++) {
      ALIMG = img;
      d = ref2img[img][1]*ref2img[img][12] - ref2img[img][2]*ref2img[img][11];
      memset(p,0,sizeof(p));
      p[1] = ref2img[img][12]/d*Align45scale[1];
      p[2] = -ref2img[img][2]/d*Align45scale[2];
      p[11] = -ref2img[img][11]/d*Align45scale[1];
      p[12] = ref2img[img][1]/d*Align45scale[2];
      //printf("rmsi = %f, ",0.2*sqrt(aminfuncR2I(p-1)));
      for (i=1;i<=20;i++) {
	 for (j=1;j<=20;j++) xi[i][j]=0.;
	 xi[i][i]=1.;
      }
      powell(p-1,xi,20,1.e-8,&i,&rv,&aminfuncR2I);
      //printf("rmsf = %f\n",0.2*sqrt(aminfuncR2I(p-1)));
      for (i=0;i<20;i++) ref2img[img][i+20] = p[i]/Align45scale[i%10];
   }
   free_matrix(xi,1,20,1,20);
}

#define ALIGN_SEARCH_MODE 1
void align(int ext,int fld) {
   int img,i,j,x,y,WARMSAVE,FORCE1SAVE,SEARCHMODESAVE,e,c,t,nStep,dx0,dy0,Nalign0;
   float p,sn,sig=0,**xi,*v,sigmin,sigmax;
   double dx,dy,rsig=0.0;
   FILE *f;
   char str[1024];

   WARMSAVE=WARMSTART;
   WARMSTART=0;
   FORCE1SAVE=Force1;
   Force1=0;
   SEARCHMODESAVE=SearchMode;
   SearchMode=0;
   if (Timg>Nimg) img=Nimg;
   else img=0;
   Nstars=0;
   /*
   if (alignstars[0]) {
      if ((f=fopen(alignstars,"r"))!=NULL) {
	 while (fscanf(f,"%d %d %lf %lf",&e,&c,&dx,&dy)!=EOF && Nstars<MAXNSTARS) {
	    if (e==ext && c==fld+1) {
	       alphot(img,dx,dy,&dx,&dy,2);
	    }
	 }
      }
   }
   */
   sigmax=1.e20;
   sigmin=200.;
   if (Nstars==0 && xytfile[0]) while (sigmin>5. && Nstars<1000) {
      sigmin*=0.5;
      if (sigmin<5.) sigmin=5.;
      if ((f=fopen(xytfile,"r"))!=NULL) {
	 while (fscanf(f,"%d %d %lf %lf %d %f",&e,&c,&dx,&dy,&t,&sn)!=EOF && Nstars<MAXNSTARS) {
	    if (e==ext && c==fld+1 && t==1 && dx>=XMIN && dx<XMAX && dy>=YMIN && dy<YMAX && sn>=sigmin && sn<sigmax) {
	       stars[Nstars].x=dx;
	       stars[Nstars].y=dy;
	       Nstars++;
	    }
	    fgets(str,1024,f);
	 }
	 fclose(f);
      }
      sigmax=sigmin;
   }
   NO_SKY=1;
   if (Nstars==0) {
      setindx(-1);
      memset(ran[0],0,X*Y*SubResRef*SubResRef);
      clearsnmap();
      setsnmap(img,XMIN,XMAX-1,YMIN,YMAX-1);
      sigmax=1.e20;
      sigmin=200.;
      while (sigmin>10. && Nstars<1000) {
	 sigmin*=0.7;
	 if (sigmin<10.) sigmin=10.;
	 for (y=YMIN*SubResRef+1;y<YMAX*SubResRef-1;y++) for (x=XMIN*SubResRef+1;x<XMAX*SubResRef-1;x++) if ((p=peaksn(img,x,y))>=sigmin && p<sigmax) alphot(img,x+0.5,y+0.5,&dx,&dy,1,1);
	 sigmax=sigmin;
      }
   }
   if (Nstars<500) {
      fprintf(fwarn,"Only %d stars for alignment\n",Nstars);
      fflush(fwarn);
   }
   printf("%d stars for alignment\n",Nstars);
   fflush(stdout);
   if (!FakeStars[0] && VerboseData>0) fprintf(fverb,"Align: %d\n",Nstars);
#if ALIGN_SEARCH_MODE==1
   nStep=(int)(AlignTol/AlignStep+0.5);
   Alist=(dptype*)calloc(Nstars,sizeof(dptype));
#else
   AlignStep=1.0;
   nStep=(int)(AlignTol+0.5);
   Alist=(dptype*)calloc(Nstars*(2*nStep+1)*(2*nStep+1),sizeof(dptype));
#endif
   if (!Alist) merr();
   setindx(-1);
   memset(ran[0],0,X*Y*SubResRef*SubResRef);
   if (Align==4) NFREE=20;
   else {
      NFREE=1+Align;
      if (Rotate) NFREE++;
   }
   v=vector(1,NFREE);
   xi=matrix(1,NFREE,1,NFREE);
   for (img=0;img<Nimg;img++) if (img || Timg>Nimg) {
      int it;
      float xoff=0,yoff=0;
      ALIMG=img;

      // find any large offset
      if (nStep>0) {
	 int nlist=0,nn;
	 float *xlist,*ylist;
	 int r,n,nmax=-1,xmax=0,ymax=0;
	 xlist = (float*)calloc(Nstars*(2*nStep+1)*(2*nStep+1),FLOATSIZE);
	 ylist = (float*)calloc(Nstars*(2*nStep+1)*(2*nStep+1),FLOATSIZE);
	 for (i=0;i<Nstars;i++) {
	    double sx0,sy0;
	    shift(img,stars[i].x,stars[i].y,&sx0,&sy0,1);
	    if (posOK(img,sx0,sy0)) {
#if ALIGN_SEARCH_MODE==1
	       int err;
	       float best=0,snr;
	       for (err=0;err<=nStep && best==0;err++) {
		  for (dx0=-err;dx0<=err;dx0++) {
		     dy0=err-abs(dx0);
		     while (1) {
			shift(img,sx0+dx0*AlignStep,sy0+dy0*AlignStep,&dx,&dy,-1);
			snr = alphot(img,dx,dy,&dx,&dy,0,0);
			if (snr>0.0) {
			   shift(img,dx,dy,&dx,&dy,1);
			   dx -= sx0;
			   dy -= sy0;
			   snr *= exp( -0.5*(dx*dx+dy*dy)/(AlignStep*AlignStep) );
			   if (snr>best) {
			      best = snr;
			      xlist[nlist] = dx;
			      ylist[nlist] = dy;
			   }
			}
			if (dy0<=0) break;
			dy0 = -dy0;
		     }
		  }
	       }
	       if (best>0) nlist++;
#else
	       int nlist0;
	       nlist0=nlist;
	       for (dx0=-nStep;dx0<=nStep;dx0++) for (dy0=-nStep;dy0<=nStep;dy0++) {
		  int ix,iy;
		  ix = (int)(sx0+dx0+100)-100;
		  iy = (int)(sy0+dy0+100)-100;
		  if (posOK(img,ix,iy) && peak2(img,ix,iy)) {
		     shift(img,sx0+dx0,sy0+dy0,&dx,&dy,-1);
		     if (alphot(img,dx,dy,&dx,&dy,0,0)) {
			int newStar=1;
			shift(img,dx,dy,&dx,&dy,1);
			dx -= sx0;
			dy -= sy0;
			for (j=nlist0;j<nlist && newStar==1;j++) {
			   if (fabs(xlist[j]-dx)<0.5 && fabs(ylist[j]-dy)<0.5) newStar=0;
			}
			if (newStar==1) {
			   xlist[nlist] = dx;
			   ylist[nlist] = dy;
			   nlist++;
			}
		     }
		  }
	       }
#endif
	    }
	 }
	 r = (int)(nStep*AlignStep+3.5);
	 for (dy0=-r;dy0<=r;dy0++) for (dx0=-r;dx0<=r;dx0++) {
	    n=0;
	    for (i=0;i<nlist;i++) {
	       if (fabs(xlist[i]-dx0)<=0.5 && fabs(ylist[i]-dy0)<=0.5) n+=4;
	       else if (fabs(xlist[i]-dx0)<=1.5 && fabs(ylist[i]-dy0)<=1.5) n++;
	    }
	    if (n>nmax) {
	       xmax=dx0;
	       ymax=dy0;
	       nmax=n;
	    }
	    //printf("%3d ",n); if (dx0==r) printf("\n");
	 }
	 nn=0;
	 for (i=0;i<nlist;i++) {
	    if (xlist[i]>=xmax-1.5 && xlist[i]<=xmax+1.5 && ylist[i]>=ymax-1.5 && ylist[i]<=ymax+1.5) {
	       float z;
	       z = xlist[i];
	       for (j=nn;j>0 && z<xlist[j-1];j--) xlist[j]=xlist[j-1];
	       xlist[j] = z;
	       z = ylist[i];
	       for (j=nn;j>0 && z<ylist[j-1];j--) ylist[j]=ylist[j-1];
	       ylist[j] = z;
	       nn++;
	    }
	 }
	 if (nn%2==1) {
	    xoff = xlist[nn/2];
	    yoff = ylist[nn/2];
	 }
	 else if (nn>0) {
	    xoff = 0.5*(xlist[nn/2-1]+xlist[nn/2]);
	    yoff = 0.5*(ylist[nn/2-1]+ylist[nn/2]);
	 }
	 //printf("Mode %d %d (n=%d); median %f %f (n=%d)\n",xmax,ymax,nlist,xoff,yoff,nn);
	 printf("image %d:\n",img+1);
	 printf("  coarse alignment = %4.2f %4.2f\n",xoff,yoff);
	 fflush(stdout);
	 free(xlist);
	 free(ylist);
	 if (!FakeStars[0] && VerboseData>0) fprintf(fverb,"Coarse align image %d: %f %f\n",img+1,xoff,yoff);
      }

      // fine alignment
      for (it=0;it<AlignIter;it++) {
	 Nalign=0;
	 for (i=0;i<Nstars;i++) {
	    shift(img,stars[i].x,stars[i].y,&dx,&dy,1);
	    dx+=xoff;
	    dy+=yoff;
	    if (posOK(img,dx,dy)) {
	       shift(img,dx,dy,&dx,&dy,-1);
	       if (alphot(img,dx,dy,&dx,&dy,0,1)>0) {
		  shift(img,dx,dy,&dx,&dy,1);
		  Alist[Nalign][0]=stars[i].x;
		  Alist[Nalign][1]=stars[i].y;
		  Alist[Nalign][2]=dx;
		  Alist[Nalign][3]=dy;
		  Nalign++;
	       }
	    }
	 }
	 if (it==AlignIter-1) {
	    if (nStep==0) printf("image %d: %d matched, ",img+1,Nalign);
	    else printf("  %d matched, ",Nalign);
	    fflush(stdout);
	    if (!FakeStars[0] && VerboseData>0) fprintf(fverb,"Align image %d: %d ",img+1,Nalign);
	 }
	 Nalign0 = Nalign;
	 if (Nalign) {
	    do {
	       int nit=0;
	       float rval=0.;

	       for (i=1;i<=NFREE;i++) {
		  for (j=1;j<=NFREE;j++) xi[i][j]=0.;
		  xi[i][i]=1.;
	       }
	       if (Align==4) {
		  for (i=0;i<NFREE;i++) v[i+1]=dpos[img][i]*Align45scale[i%10];
	       }
	       else {
		  for (i=0;i<NFREE;i++) v[i+1]=dpos[img][i];
		  if (Rotate) v[NFREE]=dpos[img][4];
	       }
	       powell(v,xi,NFREE,1.e-5,&nit,&rval,&aminfunc);
	       if (Align==4) {
		  for (i=0;i<NFREE;i++) dpos[img][i]=v[i+1]/Align45scale[i%10];
	       }
	       else {
		  for (i=0;i<=Align;i++) dpos[img][i]=v[i+1];
		  if (Rotate) dpos[img][4]=v[Align+2];
	       }
	       sig=0;
	       for (i=0;i<Nalign;i++) {
		  shift(img,Alist[i][0],Alist[i][1],&dx,&dy,1);
		  sig+=(Alist[i][4]=(Alist[i][2]-dx)*(Alist[i][2]-dx)+(Alist[i][3]-dy)*(Alist[i][3]-dy));
	       }
	       sig/=Nalign;
	       if (sig<0.0001) sig=0.0001;
	       x=0;
	       for (i=0;i<Nalign;i++) {
		  if (Alist[i][4]>sig*5.0625) {
		     Nalign--;
		     memcpy(Alist[i],Alist[Nalign],sizeof(dptype));
		     x=1;
		     i--;
		  }
	       }
	    } while (Nalign && x);
	 }
	 if (Align>=4) rsig=fitFwdAlign4(img);
	 if (it==AlignIter-1) {
	    if (Nalign==0) {
	       char imgstr[161];
	       getimgstr(img,imgstr);
	       fprintf(fwarn,"No alignment stars matched for image %d, %s\n",img+1,imgstr);
	       fflush(fwarn);
	    }
	    else if (Nalign==Nalign0) {
	       char imgstr[161];
	       getimgstr(img,imgstr);
	       fprintf(fwarn,"No alignment stars trimmed for image %d, %s\n",img+1,imgstr);
	       fflush(fwarn);
	    }
	    else if (Nalign<50 && Nalign<Nstars) {
	       char imgstr[161];
	       getimgstr(img,imgstr);
	       fprintf(fwarn,"Only %d alignment stars used for image %d, %s\n",Nalign,img+1,imgstr);
	       fflush(fwarn);
	    }
	    if (hstmode[img].inst!=NONE && sqrt(sig)>=0.5) {
	       char imgstr[161];
	       getimgstr(img,imgstr);
	       fprintf(fwarn,"Large alignment scatter (%5.3f pix) for image %d, %s\n",sqrt(sig),img+1,imgstr);
	       fflush(fwarn);
	    }
	    if (Align<4) {
	       printf("%d used, %4.2f %4.2f %8.6f %7.5f %5.3f, sig=%4.2f\n",Nalign,dpos[img][0],dpos[img][1],dpos[img][2],dpos[img][3],dpos[img][4],sqrt(sig));
	       if (!FakeStars[0] && VerboseData>0) fprintf(fverb,"%d %f %f %f %f %f %f\n",Nalign,dpos[img][0],dpos[img][1],dpos[img][2],dpos[img][3],dpos[img][4],sqrt(sig));
	    }
	    else {
	       printf("%d used, sig=%4.2f, rsig=%5.3f\n",Nalign,sqrt(sig),rsig);
	       if (!FakeStars[0] && VerboseData>0) fprintf(fverb,"%d %f %f\n",Nalign,sqrt(sig),rsig);
	    }
	    fflush(stdout);
#ifdef PGPLOT
	    if (DiagPlotType[0] && Nalign>0) {
	       alignDiagData[img][DIAG_EXT][DIAG_Z].N = Nalign;
	       alignDiagData[img][DIAG_EXT][DIAG_Z].XMIN = XMIN;
	       alignDiagData[img][DIAG_EXT][DIAG_Z].XMAX = XMAX;
	       alignDiagData[img][DIAG_EXT][DIAG_Z].YMIN = YMIN;
	       alignDiagData[img][DIAG_EXT][DIAG_Z].YMAX = YMAX;
	       alignDiagData[img][DIAG_EXT][DIAG_Z].sig = sqrt(sig);
	       alignDiagData[img][DIAG_EXT][DIAG_Z].data = (alignDiagDataType*)calloc(Nalign,sizeof(alignDiagDataType));
	       if (!alignDiagData[img][DIAG_EXT][DIAG_Z].data) merr();
	       for (i=0;i<Nalign;i++) {
		  shift(img,Alist[i][0],Alist[i][1],&dx,&dy,1);
		  alignDiagData[img][DIAG_EXT][DIAG_Z].data[i].x = Alist[i][0];
		  alignDiagData[img][DIAG_EXT][DIAG_Z].data[i].y = Alist[i][1];
		  alignDiagData[img][DIAG_EXT][DIAG_Z].data[i].dx = Alist[i][2]-dx;
		  alignDiagData[img][DIAG_EXT][DIAG_Z].data[i].dy = Alist[i][3]-dy;
	       }
	    }
#endif
	 }
      }
   }
   NO_SKY=0;
   free_matrix(xi,1,NFREE,1,NFREE);
   free_vector(v,1,NFREE);
   free(Alist);
   WARMSTART=WARMSAVE;
   Force1=FORCE1SAVE;
   SearchMode=SEARCHMODESAVE;
   return;
}

double getscale(int img,double X,double Y) {
   double x0,y0,x,y;
   float scale,maxscale;
   shift(img,X,Y,&x0,&y0,-1);
   shift(img,X+1,Y,&x,&y,-1);
   maxscale = sqrt((x-x0)*(x-x0)+(y-y0)*(y-y0));
   shift(img,X,Y+1,&x,&y,-1);
   scale = sqrt((x-x0)*(x-x0)+(y-y0)*(y-y0));
   if (scale>maxscale) maxscale=scale;
   return maxscale;
}

void setRMark(void) {
   int img,r,xmax,ymax;
   float scale,maxscale;
   rMark=0;
   rMarkPSF=0;
   for (img=0;img<Nimg;img++) {
      xmax=dataim[img].X-1;
      ymax=dataim[img].Y-1;
      maxscale=getscale(img,0,0);
      scale=getscale(img,xmax,0); if (scale>maxscale) maxscale=scale;
      scale=getscale(img,0,ymax); if (scale>maxscale) maxscale=scale;
      scale=getscale(img,xmax,ymax); if (scale>maxscale) maxscale=scale;
      scale=getscale(img,xmax*0.5,ymax*0.5); if (scale>maxscale) maxscale=scale;
      r=(int)((rphot[img]+RPSF[img])*maxscale+0.9999);
      if (r>rMark) rMark=r;
      r=(int)(RPSF[img]*maxscale+0.9999);
      if (r>rMarkPSF) rMarkPSF=r;
   }
}

void find(int show) {
   int xx,x,y,nf=0,nf0=0,i;
   float p,sigmax=1.e30,sigmin=200,fsigm;

   if (show) {printf("Finding stars: "); fflush(stdout);}
   setindx(-1);
   fsigm=SigFind*SigFindMult;
   memset(ran[0],0,X*Y*SubResRef*SubResRef);
   clearsnmap();
   NO_SKY=1;
   setsnmap(-1,XMIN,XMAX-1,YMIN,YMAX-1);
   NO_SKY=0;
   if (show) {printf("."); fflush(stdout);}
   while (sigmin>fsigm) {
      sigmin*=0.5;
      if (sigmin<fsigm) sigmin=fsigm;
      for (x=XMIN*SubResRef+1;x<XMAX*SubResRef-1;x++) for (y=YMIN*SubResRef+1;y<YMAX*SubResRef-1;y++) {
              if ((p=qpeak(x,y,1))>=sigmin && p<sigmax) {
                  phot(x/(double)SubResRef,y/(double)SubResRef,-1,0);
                  for (;nf<Nstars;nf++) imgsub(nf);
              }
          }
     sigmax=sigmin;
   }
   for (i=0;i<SecondPass;i++) {
      if (show) {printf("."); fflush(stdout);}
      memset(ran[0],0,X*Y*SubResRef*SubResRef);
      for (y=0;y<Y*SubResRef;y++) for (x=0;x<X*SubResRef;x++) if (indx[y][x]>=0) {
	 ran[y][x]=1;
	 if (RCombine>=1.9) {
	    if (y>0) ran[y-1][x]=1;
	    if (y<Y-1) ran[y+1][x]=1;
	    if (x>0) ran[y][x-1]=1;
	    if (x<X-1) ran[y][x+1]=1;
	 }
	 if (RCombine>=2.4) {
	    if (y>0) {
	       if (x>0) ran[y-1][x-1]=1;
	       if (x<X-1) ran[y-1][x+1]=1;
	    }
	    if (y<Y-1) {
	       if (x>0) ran[y+1][x-1]=1;
	       if (x<X-1) ran[y+1][x+1]=1;
	    }
	 }
      }
      clearsnmapflag();
      for (xx=nf0;xx<nf;xx++) {
	 NO_SKY=1;
	 setsnmap(-1,(int)stars[xx].x-rMark-1,(int)stars[xx].x+rMark+1,(int)stars[xx].y-rMark-1,(int)stars[xx].y+rMark+1);
	 NO_SKY=0;
	 for (y=(int)((stars[xx].y-rMarkPSF)*SubResRef);y<=(int)((stars[xx].y+rMarkPSF)*SubResRef);y++) if (y>=YMIN*SubResRef && y<YMAX*SubResRef) for (x=(int)((stars[xx].x-rMarkPSF)*SubResRef);x<=(int)((stars[xx].x+rMarkPSF)*SubResRef);x++) if (x>=XMIN*SubResRef && x<XMAX*SubResRef && qpeak(x,y,1)>=fsigm) {
	    phot(x/(double)SubResRef,y/(double)SubResRef,-1,0);
	    ran[y][x]=1;
	 }
      }
      //printf("%d\n",Nstars-nf);
      nf0=nf;
      for (;nf<Nstars;nf++) imgsub(nf);
   }
   if (show) {printf(".\n%d found\n",Nstars); fflush(stdout);}
   for (x=0;x<Nstars;x++) {
      y=stars[x].flag;
      markstars(stars[x].x,stars[x].y);
      stars[x].flag=y;
   }
   clean(1,show);
   return;
}

void readwarm(int ext,int fld) {
   FILE *f;
   int e,c,t;
   float x,y,z;
   char str[1024];

   setindx(-1);
   if ((f=fopen(xytfile,"r"))==NULL) {
      printf("Cannot open %s\n",xytfile);
      fflush(stdout);
      return;
   }
   printf("Reading %s ...",xytfile);
   fflush(stdout);
   while (fscanf(f,"%d %d %f %f %d %f",&e,&c,&x,&y,&t,&z)!=EOF) {
      if (e==ext && c==fld+1) {
	 if (t==2) t=1;
	 stars[Nstars].x=x;
	 stars[Nstars].y=y;
	 stars[Nstars].type=t;
	 getpsfpars(-1,x,y,stars[Nstars].pa,stars[Nstars].pb,stars[Nstars].pc);
	 if (Timg>Nimg) getpsfpars(Nimg,x,y,stars[Nstars].pa,stars[Nstars].pb,stars[Nstars].pc);
	 phot((int)x,(int)y,Nstars,0);
	 imgsub(Nstars);
	 if (WARMSTART==2) fscanf(f,"%f %f",refcts+Nstars,refmag+Nstars);
	 Nstars++;
      }
      fgets(str,1024,f);
   }
   fclose(f);
   printf(" %d stars\n",Nstars);
   fflush(stdout);
   return;
   /* should readwarm() also search for new stars (pass2)? */
}

void cleansat(int img) {
   int x,y,cont=1,x0,x1,y0,y1;
   int NPIX;
   float DDMAX;
   double tx,ty;
   char **fixed,**last;
   void *ptr;

   return;
   NPIX=dataim[img].X*dataim[img].Y;
   fixed=(char**)calloc(dataim[img].Y,PTRSIZE);
   if (!fixed) merr();
   ptr=malloc(NPIX); if (!ptr) merr();
   for (y=0;y<dataim[img].Y;y++) fixed[y]=ptr+dataim[img].X*y;
   last=(char**)calloc(dataim[img].Y,PTRSIZE);
   if (!last) merr();
   ptr=malloc(NPIX); if (!ptr) merr();
   for (y=0;y<dataim[img].Y;y++) last[y]=ptr+dataim[img].X*y;
   x0=dataim[img].X;
   x1=0;
   y0=dataim[img].Y;
   y1=0;
   for (x=XMIN;x<XMAX;x++) {
      shift(img,x,YMIN,&tx,&ty,1);
      if (tx<x0) x0=(int)tx;
      if (tx>x1) x1=(int)tx;
      if (ty<y0) y0=(int)ty;
      if (ty>y1) y1=(int)ty;
      shift(img,x,YMAX,&tx,&ty,1);
      if (tx<x0) x0=(int)tx;
      if (tx>x1) x1=(int)tx;
      if (ty<y0) y0=(int)ty;
      if (ty>y1) y1=(int)ty;
   }
   for (y=YMIN;y<YMAX;y++) {
      shift(img,XMIN,y,&tx,&ty,1);
      if (tx<x0) x0=(int)tx;
      if (tx>x1) x1=(int)tx;
      if (ty<y0) y0=(int)ty;
      if (ty>y1) y1=(int)ty;
      shift(img,XMAX,y,&tx,&ty,1);
      if (tx<x0) x0=(int)tx;
      if (tx>x1) x1=(int)tx;
      if (ty<y0) y0=(int)ty;
      if (ty>y1) y1=(int)ty;
   }
   if (x0<0) x0=0;
   if (x1>dataim[img].X) x1=dataim[img].X;
   if (y0<0) y0=0;
   if (y1>dataim[img].Y) y1=dataim[img].Y;
   DDMAX=iDMAX[img]*FSat;
   memset(fixed[0],0,NPIX);
   for (y=y0;y<y1;y++) for (x=x0;x<x1;x++) if (data[img][y][x]>=DDMAX) {
      data[img][y][x]=iDMAX[img];
      if (x>0) fixed[y][x-1]=1;
      if (x<dataim[img].X-1) fixed[y][x+1]=1;
      if (y>0) fixed[y-1][x]=1;
      if (y<dataim[img].Y-1) fixed[y+1][x]=1;
   }
   while (cont) {
      cont=0;
      memcpy(last[0],fixed[0],NPIX);
      memset(fixed[0],0,NPIX);
      for (y=y0;y<y1;y++) for (x=x0;x<x1;x++) if (dataOK(img,x,y) && last[y][x] && peak2(img,x,y)) {
	 data[img][y][x]=iDMAX[img];
	 if (x>0 && data[img][y][x-1]<iDMAX[img]) fixed[y][x-1]=1;
	 if (x<dataim[img].X-1 && data[img][y][x+1]<iDMAX[img]) fixed[y][x+1]=1;
	 if (y>0 && data[img][y-1][1]<iDMAX[img]) fixed[y-1][x]=1;
	 if (y<dataim[img].Y-1 && data[img][y+1][1]<iDMAX[img]) fixed[y+1][x]=1;
	 cont=1;
      }
   }
   free(last[0]);
   free(last);
   free(fixed[0]);
   free(fixed);
   return;
}

float minfunc(float x,int N,float3*apc,float3*sig,int*flagUse,int idx) {
   int i;
   float q=0;
   int nGood=0;
   for (i=0;i<N;i++) if (flagUse[i]) {
      q+=log(1+(x-apc[i][idx])*(x-apc[i][idx])*0.5*sig[i][idx]*sig[i][idx]);
      nGood++;
   }
   return q/nGood;
}

float dminfunc(float x,int N,float3*apc,float3*sig,int*flagUse,int idx) {
   int i;
   float q=0;
   int nGood=0;
   for (i=0;i<N;i++) if (flagUse[i]) {
      q+=(x-apc[i][idx])/(1/sig[i][idx]/sig[i][idx]+0.5*(x-apc[i][idx])*(x-apc[i][idx]));
      nGood++;
   }
   return q/nGood;
}

//adapted from Numerical Recipes;
#define ITMAX 100
#define ZEPS 1.0e-10
#define MOV3(a,b,c, d,e,f) (a)=(d);(b)=(e);(c)=(f);
#define SIGN(a,b) ((b) >= 0.0 ? fabs(a) : -fabs(a))
float dbrent(float ax, float bx, float cx, float (*f)(float,int,float3*,float3*,int*,int),
	float (*df)(float,int,float3*,float3*,int*,int), float tol, float *xmin,int fN,float3*fapc,float3*fsig,int*flagUse,int idx)
{
	int iter,ok1,ok2;
	float a,b,d=0,d1,d2,du,dv,dw,dx,e=0.0;
	float fu,fv,fw,fx,olde,tol1,tol2,u,u1,u2,v,w,x,xm;

	a=(ax < cx ? ax : cx);
	b=(ax > cx ? ax : cx);
	x=w=v=bx;
	fw=fv=fx=(*f)(x,fN,fapc,fsig,flagUse,idx);
	dw=dv=dx=(*df)(x,fN,fapc,fsig,flagUse,idx);
	for (iter=1;iter<=ITMAX;iter++) {
		xm=0.5*(a+b);
		tol1=tol*fabs(x)+ZEPS;
		tol2=2.0*tol1;
		if (fabs(x-xm) <= (tol2-0.5*(b-a))) {
			*xmin=x;
			return fx;
		}
		if (fabs(e) > tol1) {
			d1=2.0*(b-a);
			d2=d1;
			if (dw != dx) d1=(w-x)*dx/(dx-dw);
			if (dv != dx) d2=(v-x)*dx/(dx-dv);
			u1=x+d1;
			u2=x+d2;
			ok1 = (a-u1)*(u1-b) > 0.0 && dx*d1 <= 0.0;
			ok2 = (a-u2)*(u2-b) > 0.0 && dx*d2 <= 0.0;
			olde=e;
			e=d;
			if (ok1 || ok2) {
				if (ok1 && ok2)
					d=(fabs(d1) < fabs(d2) ? d1 : d2);
				else if (ok1)
					d=d1;
				else
					d=d2;
				if (fabs(d) <= fabs(0.5*olde)) {
					u=x+d;
					if (u-a < tol2 || b-u < tol2)
						d=SIGN(tol1,xm-x);
				} else {
					d=0.5*(e=(dx >= 0.0 ? a-x : b-x));
				}
			} else {
				d=0.5*(e=(dx >= 0.0 ? a-x : b-x));
			}
		} else {
			d=0.5*(e=(dx >= 0.0 ? a-x : b-x));
		}
		if (fabs(d) >= tol1) {
			u=x+d;
			fu=(*f)(u,fN,fapc,fsig,flagUse,idx);
		} else {
			u=x+SIGN(tol1,d);
			fu=(*f)(u,fN,fapc,fsig,flagUse,idx);
			if (fu > fx) {
				*xmin=x;
				return fx;
			}
		}
		du=(*df)(u,fN,fapc,fsig,flagUse,idx);
		if (fu <= fx) {
			if (u >= x) a=x; else b=x;
			MOV3(v,fv,dv, w,fw,dw)
			MOV3(w,fw,dw, x,fx,dx)
			MOV3(x,fx,dx, u,fu,du)
		} else {
			if (u < x) a=u; else b=u;
			if (fu <= fw || w == x) {
				MOV3(v,fv,dv, w,fw,dw)
				MOV3(w,fw,dw, u,fu,du)
			} else if (fu < fv || v == x || v == w) {
				MOV3(v,fv,dv, u,fu,du)
			}
		}
	}
	printf("**Aperture corrections failed to converge!\n");
	fflush(stdout);
	*xmin=u;
	return fu;
}
#undef SIGN
#undef ITMAX
#undef ZEPS
#undef MOV3

#define APMAX 200
#define SIGCAP 10
float *slist;

Inline float getapsize(int img,int i) {
#ifdef USEWFPC2
   if (hstmode[img].inst==WFPC2) {
      double tx,ty;
      shift(img,stars[i].x,stars[i].y,&tx,&ty,1);
      return wfpc2_apsize(img,tx,ty);
   }
#endif
#ifdef USEACS
   if (hstmode[img].inst==ACS) {
      double tx,ty;
      shift(img,stars[i].x,stars[i].y,&tx,&ty,1);
      return acs_apsize(img,tx,ty);
   }
#endif
#ifdef USEWFC3
   if (hstmode[img].inst==WFC3) {
      double tx,ty;
      shift(img,stars[i].x,stars[i].y,&tx,&ty,1);
      return wfc3_apsize(img,tx,ty);
   }
#endif
   return apsize[img];
}

float calcapsky(int img,float x0,float y0,float rap,float*sigsky) {
   int x,y,ix,iy,Nsky=0;
   float r;
   ix=(int)x0; iy=(int)y0;
   for (y=iy-(int)apsky[img][1]-1;y<=iy+(int)apsky[img][1]+1;y++) for (x=ix-(int)apsky[img][1]-1;x<=ix+(int)apsky[img][1]+1;x++) if (ppixOK(img,x,y)) {
      r=sqrt((x+0.5-x0)*(x+0.5-x0)+(y+0.5-y0)*(y+0.5-y0));
      if (r>rap+0.5 && r>apsky[img][0] && r<=apsky[img][1]) slist[Nsky++]=res[img][y][x];
   }
   return calcskyval(Nsky,slist,img,sigsky);
}

void calcapcor(int ext,int fld,int img,FILE *f,int i,float rap,int*N,float3*apc,float3*sig,float*sq,int*id) {
   double tx,ty;
   float xc,yc,r,ssky,sigsky,isq,spsf=0;
   float3 iapc,isig;
   int x,y,ix,iy,sOK=1;

   imgadd(i);
   shift(img,stars[i].x,stars[i].y,&tx,&ty,1);
   ssky=calcapsky(img,tx,ty,rap,&sigsky);
   ix=(int)tx; xc=tx-0.5;
   iy=(int)ty; yc=ty-0.5;
   iapc[0]=iapc[1]=isig[0]=isig[1]=0;
   if (ssky<1.e30) {
      shift(img,stars[i].x,stars[i].y,&tx,&ty,1);
      calc1psf(img,tx,ty,RPSF[img],stars[i].pa[img],stars[i].pb[img],stars[i].pc[img],1,1);
      for (x=ix-(int)(rap)-1;x<=ix+(int)(rap)+1;x++) for (y=iy-(int)(rap)-1;y<=iy+(int)(rap)+1;y++) {
	 r=0.465+rap-sqrt((x-xc)*(x-xc)+(y-yc)*(y-yc));
	 if (r>1) r=1;
	 if (r>0) {
	    if (abs(x-ix)<=RPSF[img] && abs(y-iy)<=RPSF[img]) spsf+=psf[y-iy][x-ix]*r;
	    if (ppixfOK(img,x,y)) {
	       iapc[1]+=(res[img][y][x]-ssky)*r;
	       isig[1]+=noise(img,x,y,ssky)*r;
	    }
	    else {
	       float ns=0,rs=0;
	       int ct=0,dx,dy;
	       for (dx=-1;dx<=1;dx++) for (dy=-1;dy<=1;dy++) if (ppixfOK(img,x+dx,y+dy)) {
		  ns+=noise(img,x+dx,y+dy,ssky);
		  rs+=res[img][y+dy][x+dx];
		  ct++;
	       }
	       if (ct) {
		  iapc[1]+=(rs/ct-ssky)*r;
		  isig[1]+=ns/ct*r;
	       }
	    }
	 }
      }
   }
   if (isig[1]>0 && iapc[1]>0 && stars[i].is[img]>0) {
#ifdef SINGLE_APCOR
      isig[0]=1.0;
      iapc[0]=0.0;
      isig[1]=sqrt(isig[1]/iapc[1]/iapc[1]+stars[i].iss[img]*stars[i].iss[img]/stars[i].is[img]/stars[i].is[img])*1.0857+0.01;
      iapc[1]=-2.5*log10(iapc[1]/stars[i].is[img]*stars[i].icm[img]);
      isig[2]=0.01;
      iapc[2]=0.0;
#else
      float s,ss,ssky_dum,cs;
      int ix,iy,cx,cy;
      int savePSF,saveSky;
      sky_set[img]=1;
      sky_val[img]=ssky;
      savePSF = PSFPhot; PSFPhot=0;
      saveSky = FitSky; FitSky=1;
      ix=(int)(tx+100)-100;
      iy=(int)(ty+100)-100;
      cx=(int)((tx-ix)*50+0.5);
      if (cx<0) cx=0; if (cx>50) cx=50;
      cy=(int)((ty-iy)*50+0.5);
      if (cy<0) cy=0; if (cy>50) cy=50;
      eval1_apphot(img,ix,iy,cx,cy,&s,&ss,&ssky_dum,&cs);
      PSFPhot=savePSF;
      FitSky=saveSky;
      if (s<=0) sOK=0;
      else {
	 //printf("%6.2f %6.2f: %f %f -> %f %f (%f)\n",stars[i].x,stars[i].y,stars[i].is[img]/stars[i].icm[img],stars[i].iss[img]/stars[i].icm[img],s,ss,spsf);
	 isig[0]=sqrt(ss*ss/s/s+stars[i].iss[img]*stars[i].iss[img]/stars[i].is[img]/stars[i].is[img])*1.0857+0.01;
	 iapc[0]=-2.5*log10(s/stars[i].is[img]*stars[i].icm[img]);
	 isig[1]=sqrt((isig[1]+9.8696044*rap*rap*rap*rap*sigsky*sigsky)/iapc[1]/iapc[1]+ss*ss/s/s)*1.0857+0.01;
	 iapc[1]=-2.5*log10(iapc[1]/spsf/s);
	 isig[2]=0.01;
	 iapc[2]=-2.5*log10(spsf);
      }
#endif
      if (sOK==1) {
	 isq=1/isig[1];
	 if (stars[i].sky>=0) isq/=1+stars[i].sky*0.1;
	 fprintf(f,"%2d %d %d %7.2f %7.2f %6.3f %6.3f %6.3f %6.3f %6.3f %6.3f %6.3f %6.2f\n",img+1,ext,fld+1,stars[i].x,stars[i].y,-2.5*log10(stars[i].is[img])+Zero,iapc[0],isig[0],iapc[1],isig[1],iapc[2],isig[2],isq);
	 if (fabs(iapc[0])<=0.25 && fabs(iapc[1])<=0.25) {
	    isig[0]=1./isig[0];
	    isig[1]=1./isig[1];
	    isig[2]=1./isig[2];
	    for (x=(*N);x>0 && isq>sq[x-1];x--) {
	       apc[x][0]=apc[x-1][0];
	       sig[x][0]=sig[x-1][0];
	       apc[x][1]=apc[x-1][1];
	       sig[x][1]=sig[x-1][1];
	       apc[x][2]=apc[x-1][2];
	       sig[x][2]=sig[x-1][2];
	       sq[x]=sq[x-1];
	       id[x]=id[x-1];
	    }
	    apc[x][0]=iapc[0];
	    sig[x][0]=isig[0];
	    apc[x][1]=iapc[1];
	    sig[x][1]=isig[1];
	    apc[x][2]=iapc[2];
	    sig[x][2]=isig[2];
	    sq[x]=isq;
	    id[x]=i;
	    if ((*N)<APMAX) (*N)++;
	 }
      }
   }
   imgsub(i);
   return;
}

void startapcor(void) {
   int i;

   printf("Aperture corrections:\n");
   fflush(stdout);
   sortstars();
   setindx(-1);
   for (i=0;i<Nstars;i++) {
      indx[(int)(stars[i].y*SubResRef)][(int)(stars[i].x*SubResRef)]=i;
      stars[i].flag=0;
   }
   return;
}

void getapcor(int ext,int fld,int img,FILE *f) {
   int N=0,i,j,x;
   float3 apc[APMAX+1],sig[APMAX+1];
   float sq[APMAX+1],mn[3],hisig[SIGCAP],rap,gmn,gsd,gsg,w;
   int flagUse[APMAX+1],id[APMAX+1];

   slist=(float*)calloc((int)(4*apsky[img][1]*apsky[img][1])+1,FLOATSIZE);
   if (!slist) merr();
   if (psfstars[0]) {
      FILE *fpsf;
      if ((fpsf=fopen(psfstars,"r"))==NULL) {
	 printf("**%s not found\n",psfstars);
	 fflush(stdout);
      }
      else {
	 int pe,pc,bi;
	 double px,py,r,br,brap;
	 char str[161],*ptr,*ptr2;

	 while (fgets(str,161,fpsf)) {
	    pe=strtol(str,&ptr,10);
	    pc=strtol(ptr,&ptr,10);
	    px=strtod(ptr,&ptr)-psfoff;
	    py=strtod(ptr,&ptr2)-psfoff;
	    if (pe==ext && pc==fld+1 && ptr!=ptr2) {
	       bi=-1;
	       br=RCombine*RCombine/(double)(SubResRef*SubResRef);
	       brap=10.;
	       if (br<1) br=1;
	       for (i=0;i<Nstars;i++) if (!stars[i].flag && (r=(px-stars[i].x)*(px-stars[i].x)+(py-stars[i].y)*(py-stars[i].y))<br) {
		  rap=getapsize(img,i);
		  if (goodstar(img,i,(int)rap+1,(int)rap+1,0)) {
		     br=r;
		     bi=i;
		     brap=rap;
		  }
	       }
	       if (bi>=0) {
		  calcapcor(ext,fld,img,f,bi,brap,&N,apc,sig,sq,id);
		  stars[i].flag=1;
	       }
	    }
	 }
	 fclose(fpsf);
      }
      for (i=0;i<Nstars;i++) stars[i].flag=0;
   }
   else for (i=0;i<Nstars;i++) {
      rap=getapsize(img,i);
      if (goodstar(img,i,(int)rap+1,(int)rap+1,2)) calcapcor(ext,fld,img,f,i,rap,&N,apc,sig,sq,id);
   }
   free(slist);

   if (N<50) {
      char imgstr[161];
      getimgstr(img,imgstr);
      fprintf(fwarn,"Only %d aperture stars in image %d, %s\n",N,img+1,imgstr);
      fflush(fwarn);
   }
   printf("image %d: %d total aperture stars\n",img+1,N);
   if (!FakeStars[0] && VerboseData>0) fprintf(fverb,"Apcor image %d: %d",img+1,N);
   for (j=0;j<3;j++) {
      int nGood=0;
      // truncate sigma values
      for (i=0;i<SIGCAP;i++) hisig[i]=0.0;
      for (i=0;i<N;i++) {
	 for (x=SIGCAP;x>0 && sig[i][j]>hisig[x-1];x--) if (x<SIGCAP) hisig[x]=hisig[x-1];
	 if (x<SIGCAP) hisig[x]=sig[i][j];
      }
      for (x=SIGCAP-1;x>=0 && hisig[x]==0;x--);
      if (x>=0) for (i=0;i<N;i++) if (sig[i][j]>hisig[x]) sig[i][j]=hisig[x];

      // iterate apcor solution
      for (i=0;i<N;i++) {
	 flagUse[i]=1;
	 nGood++;
      }
      i=1;
      while (nGood>0 && i==1) {
	 dbrent(-0.75,0.,0.75,&minfunc,&dminfunc,0.01,mn+j,N,apc,sig,flagUse,j);
	 i=0;
	 for (x=0;x<N;x++) if (flagUse[x] && fabs(apc[x][j]-mn[j])>2./sig[x][j]+0.03) {
	    flagUse[x]=0;
	    i=1;
	    nGood--;
	 }
      }
      if (nGood<=0) mn[j]=0.;
      gmn=gsd=gsg=0.;
      if (nGood>0) {
	 for (x=0;x<N;x++) if (flagUse[x]) {
	    w=sig[x][j]*sig[x][j];
	    gmn+=apc[x][j]*w;
	    gsg+=w;
	    gsd+=apc[x][j]*apc[x][j]*w;
	 }
	 gmn/=gsg;
	 gsd/=gsg;
	 gsg=1./sqrt(gsg);
	 if (nGood>1) gsd=gsg*sqrt((gsd-gmn*gmn)*nGood/(nGood-1.))*1.25;
	 else gsd=gsg;
      }
      printf("  %d stars used, %6.3f (%5.3f +/- %5.3f, %5.3f)\n",nGood,mn[j],gmn,gsd,gsg);
      if (!FakeStars[0] && VerboseData>0) fprintf(fverb," %d %f %f %f %f",nGood,mn[j],gmn,gsd,gsg);
   }
   if (!FakeStars[0] && VerboseData>0) fprintf(fverb,"\n");
   fflush(stdout);
   fflush(f);
   apcor[img]=pow(10,-0.4*(mn[0]+mn[1]+mn[2]));
#ifdef PGPLOT
   if (DiagPlotType[0] && N>0) {
      apcorDiagData[img][DIAG_EXT][DIAG_Z].N = N;
      apcorDiagData[img][DIAG_EXT][DIAG_Z].dmmean = mn[0]+mn[1]+mn[2];
      apcorDiagData[img][DIAG_EXT][DIAG_Z].m = (float*)calloc(N,sizeof(float));
      apcorDiagData[img][DIAG_EXT][DIAG_Z].dm = (float*)calloc(N,sizeof(float));
      getimgstr(img,apcorDiagData[img][DIAG_EXT][DIAG_Z].imgstr);
#ifdef USEWFPC2
      if (hstmode[img].inst==WFPC2) sprintf(apcorDiagData[img][DIAG_EXT][DIAG_Z].filtstr,"WFPC2 %s",WFPC2imagestring(img));
      else
#endif
#ifdef USEACS
      if (hstmode[img].inst==ACS) sprintf(apcorDiagData[img][DIAG_EXT][DIAG_Z].filtstr,"ACS %s",ACSimagestring(img));
      else
#endif
#ifdef USEWFC3
      if (hstmode[img].inst==WFC3) sprintf(apcorDiagData[img][DIAG_EXT][DIAG_Z].filtstr,"WFC3 %s",WFC3imagestring(img));
      else
#endif
	 apcorDiagData[img][DIAG_EXT][DIAG_Z].filtstr[0]=0;
      apcorDiagData[img][DIAG_EXT][DIAG_Z].inst = hstmode[img].inst;
      if (hstmode[img].inst==NONE) apcorDiagData[img][DIAG_EXT][DIAG_Z].filt=0;
      else apcorDiagData[img][DIAG_EXT][DIAG_Z].filt = hstmode[img].filt;
      if (!apcorDiagData[img][DIAG_EXT][DIAG_Z].m || !apcorDiagData[img][DIAG_EXT][DIAG_Z].dm) merr();
      for (i=0;i<N;i++) {
	 apcorDiagData[img][DIAG_EXT][DIAG_Z].m[i] = -2.5*log10(stars[id[i]].is[img])+Zero;
	 apcorDiagData[img][DIAG_EXT][DIAG_Z].dm[i] = apc[i][0]+apc[i][1]+apc[i][2];
	 if (i==0) {
	    apcorDiagData[img][DIAG_EXT][DIAG_Z].mmin = apcorDiagData[img][DIAG_EXT][DIAG_Z].mmax = apcorDiagData[img][DIAG_EXT][DIAG_Z].m[i];
	    apcorDiagData[img][DIAG_EXT][DIAG_Z].dmmin = apcorDiagData[img][DIAG_EXT][DIAG_Z].dmmax = apcorDiagData[img][DIAG_EXT][DIAG_Z].dm[i];
	 }
	 else {
	    if (apcorDiagData[img][DIAG_EXT][DIAG_Z].m[i]<apcorDiagData[img][DIAG_EXT][DIAG_Z].mmin) apcorDiagData[img][DIAG_EXT][DIAG_Z].mmin = apcorDiagData[img][DIAG_EXT][DIAG_Z].m[i];
	    else if (apcorDiagData[img][DIAG_EXT][DIAG_Z].m[i]>apcorDiagData[img][DIAG_EXT][DIAG_Z].mmax) apcorDiagData[img][DIAG_EXT][DIAG_Z].mmax = apcorDiagData[img][DIAG_EXT][DIAG_Z].m[i];
	    if (apcorDiagData[img][DIAG_EXT][DIAG_Z].dm[i]<apcorDiagData[img][DIAG_EXT][DIAG_Z].dmmin) apcorDiagData[img][DIAG_EXT][DIAG_Z].dmmin = apcorDiagData[img][DIAG_EXT][DIAG_Z].dm[i];
	    else if (apcorDiagData[img][DIAG_EXT][DIAG_Z].dm[i]>apcorDiagData[img][DIAG_EXT][DIAG_Z].dmmax) apcorDiagData[img][DIAG_EXT][DIAG_Z].dmmax = apcorDiagData[img][DIAG_EXT][DIAG_Z].dm[i];
	 }
      }
   }
#endif
   return;
}
#undef APMAX

#define PASTEP 5
#define NPA 18
typedef float rndtype[NPA];
void get1round(int i,int img,rndtype rnd,rndtype wt) {
   int ix,iy,cx,cy,x1,y1,dx,dy,pa;
   double x,y;
   float mr,mr2,trnd,ssky,s,ss,cs,fpa,s2,mx,my;

   shift(img,stars[i].x,stars[i].y,&x,&y,1);
   ix=(int)(x+100)-100;
   iy=(int)(y+100)-100;
   cx=(int)((x+10-(int)(x+10))*50+0.5);
   if (cx<0) cx=0; if (cx>50) cx=50;
   cy=(int)((y+10-(int)(y+10))*50+0.5);
   if (cy<0) cy=0; if (cy>50) cy=50;
   calc1psf(img,x,y,RPSF[img],stars[i].pa[img],stars[i].pb[img],stars[i].pc[img],1,1);
   if (RBig>0 && stars[i].pa[img]==RBig && stars[i].pb[img]==RBig && stars[i].pc[img]==0.) NOSKY34=1;
   clear_sky();
   if (FitSky==3) {
      if (PSFPhot) eval1_psfphot_sky3(img,ix,iy,cx,cy,&s,&ss,&ssky,&cs);
      else eval1_apphot_sky3(img,ix,iy,cx,cy,&s,&ss,&ssky,&cs);
   }
   else if (FitSky==4) {
      if (PSFPhot) eval1_psfphot_sky4(img,ix,iy,cx,cy,&s,&ss,&ssky,&mx,&my,&cs);
      else eval1_apphot_sky4(img,ix,iy,cx,cy,&s,&ss,&ssky,&mx,&my,&cs);
   }
   else if (PSFPhot) eval1_psfphot(img,ix,iy,cx,cy,&s,&ss,&ssky,&cs);
   else eval1_apphot(img,ix,iy,cx,cy,&s,&ss,&ssky,&cs);
   NOSKY34=0;
   if (WARMSTART==2) s2=s+refcts[i]/refmult[img];
   else s2=s;
   if (s2<=0 || ss<=0) {
      for (pa=0;pa<NPA;pa++) {
	 if (wt) wt[pa]=0.;
	 rnd[pa]=0.;
      }
      return;
   }
   for (pa=0;pa<NPA;pa++) {
      fpa=pa*0.017453293*PASTEP;
      trnd=mr=mr2=0.;
      for (y1=iy-RPSF[img],dy=-RPSF[img];y1<=iy+RPSF[img];y1++,dy++) for (x1=ix-RPSF[img],dx=-RPSF[img];x1<=ix+RPSF[img];x1++,dx++) if ((dx!=0 || dy!=0) && dx*dx+dy*dy<(RPSF[img]+0.465)*(RPSF[img]+0.465)) {
	 if (ppixOK(img,x1,y1)) {
	    if (psf[dy][dx]>0) {
	       mr+=psf[dy][dx]*s2;
	       mr2+=res[img][y1][x1]-ssky+psf[dy][dx]*(s2-s);
	       trnd+=(res[img][y1][x1]-ssky-psf[dy][dx]*s)*cos(2*(atan2(dy,dx)-fpa));
	    }
	 }
	 else if (ppixOK(img,ix-dx,iy-dy)) {
	    if (psf[-dy][-dx]>0) {
	       mr+=psf[-dy][-dx]*s2;
	       mr2+=res[img][iy-dy][ix-dx]-ssky+psf[-dy][-dx]*(s2-s);
	       trnd+=(res[img][iy-dy][ix-dx]-ssky-psf[-dy][-dx]*s)*cos(2*(atan2(-dy,-dx)-fpa));
	    }
	 }
      }
      if (mr2>mr) mr=mr2;
      if (mr>0) {
	 rnd[pa]=trnd/mr;
	 if (wt) wt[pa]=s2/ss/ss;
      }
      else {
	 rnd[pa]=0.;
	 if (wt) wt[pa]=0.;
      }
   }
   return;
}

void getround(int i,int*pa,float*rnd,float*irnd) {
   int img,ip,pm=1;
   float tr,wt;
   rndtype *rsave,*wsave;

   rsave=(rndtype*)calloc(Nimg,sizeof(rndtype));
   wsave=(rndtype*)calloc(Nimg,sizeof(rndtype));
   if (!rsave || !wsave) merr();
   *pa=-1;
   *rnd=-1.;
   for (img=0;img<Nimg;img++) irnd[img]=-1.;
   for (img=0;img<Nimg;img++) get1round(i,img,rsave[img],wsave[img]);
   for (ip=0;ip<NPA;ip++) {
      tr=wt=0.;
      for (img=0;img<Nimg;img++) {
	 tr+=rsave[img][ip]*wsave[img][ip];
	 wt+=wsave[img][ip];
      }
      if (wt>0 && (*pa<0 || fabs(tr/wt)>*rnd)) {
	 if (tr>=0) {
	    *rnd=tr/wt;
	    *pa=ip*PASTEP;
	    pm=1;
	 }
	 else {
	    *rnd=-tr/wt;
	    *pa=ip*PASTEP+90;
	    pm=-1;
	 }
	 for (img=0;img<Nimg;img++) {
	    irnd[img]=pm*rsave[img][ip];
	    if (irnd[img]<-9.999) irnd[img]=-9.999;
	    if (irnd[img]>9.999) irnd[img]=9.999;
	 }
      }
   }
   free(rsave);
   free(wsave);
   if (*rnd<-9.999) *rnd=-9.999;
   if (*rnd>9.999) *rnd=9.999;
   return;
}
#undef PASTEP
#undef NPA

void out1starinfo(void) {
   static int first=1;
   FILE *f;
   char str[161];
   int ct=0,i;

   if (!first) return;
   first=0;
   sprintf(str,"%s.columns",outfn);
   f=fopen(str,"w");
   if (f==0) return;
   fprintf(f,"%d. Extension (zero for base image)\n",++ct);
   fprintf(f,"%d. Chip (for three-dimensional FITS image)\n",++ct);
   fprintf(f,"%d. Object X position on reference image (or first image, if no reference)\n",++ct);
   fprintf(f,"%d. Object Y position on reference image (or first image, if no reference)\n",++ct);
   fprintf(f,"%d. Chi for fit\n",++ct);
   fprintf(f,"%d. Signal-to-noise\n",++ct);
   fprintf(f,"%d. Object sharpness\n",++ct);
   fprintf(f,"%d. Object roundness\n",++ct);
   fprintf(f,"%d. Direction of major axis (if not round)\n",++ct);
   fprintf(f,"%d. Crowding\n",++ct);
   fprintf(f,"%d. Object type (1=bright star, 2=faint, 3=elongated, 4=hot pixel, 5=extended)\n",++ct);
#ifdef USEWFPC2
   WFPC2outstarinfo(f,&ct);
#endif
#ifdef USEACS
   ACSoutstarinfo(f,&ct);
#endif
#ifdef USEWFC3
   WFC3outstarinfo(f,&ct);
#endif
   for (i=0;i<Nimg;i++) {
      char imgstr[161];
      getimgstr(i,imgstr);
      fprintf(f,"%d. Measured counts, %s\n",++ct,imgstr);
      fprintf(f,"%d. Measured sky level, %s\n",++ct,imgstr);
      fprintf(f,"%d. Normalized count rate, %s\n",++ct,imgstr);
      fprintf(f,"%d. Normalized count rate uncertainty, %s\n",++ct,imgstr);
#ifdef USEWFPC2
      if (hstmode[i].inst==WFPC2) {
	 fprintf(f,"%d. Instrumental VEGAMAG magnitude, %s\n",++ct,imgstr);
	 fprintf(f,"%d. Transformed UBVRI magnitude, %s\n",++ct,imgstr);
      }
      else
#endif
#ifdef USEACS
      if (hstmode[i].inst==ACS) {
	 fprintf(f,"%d. Instrumental VEGAMAG magnitude, %s\n",++ct,imgstr);
	 fprintf(f,"%d. Transformed UBVRI magnitude, %s\n",++ct,imgstr);
      }
      else
#endif
#ifdef USEWFC3
      if (hstmode[i].inst==WFC3) {
	 fprintf(f,"%d. Instrumental VEGAMAG magnitude, %s\n",++ct,imgstr);
	 fprintf(f,"%d. Transformed UBVRI magnitude, %s\n",++ct,imgstr);
      }
      else
#endif
	 fprintf(f,"%d. Instrumental magnitude, %s\n",++ct,imgstr);
      fprintf(f,"%d. Magnitude uncertainty, %s\n",++ct,imgstr);
      fprintf(f,"%d. Chi, %s\n",++ct,imgstr);
      fprintf(f,"%d. Signal-to-noise, %s\n",++ct,imgstr);
      fprintf(f,"%d. Sharpness, %s\n",++ct,imgstr);
      fprintf(f,"%d. Roundness, %s\n",++ct,imgstr);
      fprintf(f,"%d. Crowding, %s\n",++ct,imgstr);
      if (hstmode[i].inst==NONE) {
	 fprintf(f,"%d. PSF FWHM, %s\n",++ct,imgstr);
	 fprintf(f,"%d. PSF eccentricity, %s\n",++ct,imgstr);
	 fprintf(f,"%d. PSF a parameter, %s\n",++ct,imgstr);
	 fprintf(f,"%d. PSF b parameter, %s\n",++ct,imgstr);
	 fprintf(f,"%d. PSF c parameter, %s\n",++ct,imgstr);
      }
      fprintf(f,"%d. Photometry quality flag, %s\n",++ct,imgstr);
   }
   fclose(f);
}

int out1star(int ext,int fld,int i,FILE*of,int fake) {
   int j,SBAD=0,pa;
   stype t;
   double tx,ty;
   float rnd,s0;
   static float *irnd=NULL,*is,*is0,*iss,*ic,*ish,*issky,*icm;
   static photdatatype*pdata=NULL;

   if (irnd==NULL) {
      irnd=(float*)calloc(Nimg,FLOATSIZE);
      is=(float*)calloc(Nimg,FLOATSIZE);
      is0=(float*)calloc(Nimg,FLOATSIZE);
      iss=(float*)calloc(Nimg,FLOATSIZE);
      ic=(float*)calloc(Nimg,FLOATSIZE);
      ish=(float*)calloc(Nimg,FLOATSIZE);
      issky=(float*)calloc(Nimg,FLOATSIZE);
      icm=(float*)calloc(Nimg,FLOATSIZE);
      pdata=(photdatatype*)calloc(Nimg,sizeof(photdatatype));
      if (!irnd || !is || !is0 || !iss || !ic || !ish || !issky || !icm || !pdata) merr();
   }
   clear_sky();
   if (i<0) {
      t.x=t.y=0.;
      t.s=-1;
   }
   else {
      eval(-1,stars[i].x,stars[i].y,stars[i].pa,stars[i].pb,stars[i].pc,&stars[i].s,&s0,&stars[i].ss,&stars[i].chi,&stars[i].sh,&stars[i].sky,is,is0,iss,ic,ish,issky,icm,-1,i);
      getround(i,&pa,&rnd,irnd);
      memcpy(&t,stars+i,sizeof(stype));
   }
   if (i<0 || (WARMSTART!=2 && (stars[i].s<=0 || stars[i].s<SigFinal*stars[i].ss)) || stars[i].x<XMIN2 || stars[i].x>=XMAX2 || stars[i].y<YMIN2 || stars[i].y>=YMAX2) {
      if (fake || WARMSTART) SBAD=1;
      else return 0;
   }
   fprintf(of,"%d %d %7.2f %7.2f",ext,fld+1,t.x,t.y);
   if (t.s<=0 && !WARMSTART) fprintf(of,"   0.00     0.0  0.000  0.000   0 0.000 1");
   else fprintf(of," %6.2f %7.1f %6.3f %6.3f %3d %5.3f %d",t.chi,t.s/t.ss,t.sh,rnd,pa,s0,t.type);
   for (j=0;j<Nimg;j++) {
      pdata[j].ct0=pdata[j].ct=pdata[j].ctcorr=pdata[j].sky=pdata[j].dm=pdata[j].chi=pdata[j].sh=pdata[j].rnd=pdata[j].r=pdata[j].e=pdata[j].a=pdata[j].b=pdata[j].c=pdata[j].crowd=0.;
      pdata[j].dct=9999;
      pdata[j].dctcorr=9999;
      pdata[j].m=99.999;
      pdata[j].dm=9.999;
      shift(j,t.x,t.y,&tx,&ty,1);
      if (SBAD) pdata[j].flag=16;
      else pdata[j].flag=star_flag(stars[i].x,stars[i].y,j);
      if (!SBAD && posOK(j,tx,ty)) {
	 if (WARMSTART || iss[j]>0) {
	    pdata[j].ct=pdata[j].ct0=is[j]/icm[j];
	    pdata[j].dct=iss[j]/icm[j];
	    if (WARMSTART==2) {
	       pdata[j].ct0+=refcts[i]/refmult[j];
	       pdata[j].ct=is[j]=pdata[j].ct0*refmult[j];
	       pdata[j].dct*=refmult[j];
	       if (is[j]>0. && refcts[i]>0.) pdata[j].m=refmag[i]-2.5*log10(is[j]/refcts[i]);
	       else pdata[j].m=99.999;
	       pdata[j].ctcorr=pdata[j].ct/refcts[i]*pow(10,-0.4*refmag[i]);
	       pdata[j].dctcorr=pdata[j].dct/refcts[i]*pow(10,-0.4*refmag[i]);
	    }
	    else {
	       if (is[j]>0.) pdata[j].m=-2.5*log10(is[j])+Zero;
	       else pdata[j].m=99.999;
	       pdata[j].ctcorr=is[j]*pow(10,-0.4*Zero);
	       pdata[j].dctcorr=iss[j]*pow(10,-0.4*Zero);
	    }
	    pdata[j].crowd=is0[j];
	    pdata[j].sky=issky[j]/icm[j];
	    if (pdata[j].ct>0.) pdata[j].dm=pdata[j].dct/pdata[j].ct*1.0857362;
	    pdata[j].chi=ic[j];
	    pdata[j].sh=ish[j];
	    pdata[j].rnd=irnd[j];
	    pdata[j].a=t.pa[j];
	    pdata[j].b=t.pb[j];
	    pdata[j].c=t.pc[j];
	    getshape(pdata[j].a,pdata[j].b,pdata[j].c,&pdata[j].r,&pdata[j].e);
	 }
      }
   }
#ifdef USEWFPC2
   WFPC2outstar(of,t.x,t.y,pdata);
#endif
#ifdef USEACS
   ACSoutstar(of,t.x,t.y,pdata);
#endif
#ifdef USEWFC3
   WFC3outstar(of,t.x,t.y,pdata);
#endif
   for (j=0;j<Nimg;j++) {
      if (pdata[j].ct<999999.5) fprintf(of,"  %8.1f",pdata[j].ct);
      else fprintf(of,"  %8.2e",pdata[j].ct);
      if (pdata[j].sky<99999.95) fprintf(of," %8.2f",pdata[j].sky);
      else if (pdata[j].sky<999999.5) fprintf(of," %8.1f",pdata[j].sky);
      else fprintf(of," %8.1e",pdata[j].sky);

      if (pdata[j].ctcorr<0.0) fprintf(of," %8.1e",pdata[j].ctcorr);
      else if (pdata[j].ctcorr<0.99) fprintf(of," %8.2e",pdata[j].ctcorr);
      else if (pdata[j].ctcorr<9.999995) fprintf(of," %8.6f",pdata[j].ctcorr);
      else if (pdata[j].ctcorr<99.99995) fprintf(of," %8.5f",pdata[j].ctcorr);
      else if (pdata[j].ctcorr<999.9995) fprintf(of," %8.4f",pdata[j].ctcorr);
      else if (pdata[j].ctcorr<9999.995) fprintf(of," %8.3f",pdata[j].ctcorr);
      else if (pdata[j].ctcorr<99999.95) fprintf(of," %8.2f",pdata[j].ctcorr);
      else if (pdata[j].ctcorr<999999.5) fprintf(of," %8.1f",pdata[j].ctcorr);
      else fprintf(of," %8.1e",pdata[j].ctcorr);

      if (pdata[j].dctcorr<0.0) fprintf(of," %8.1e",pdata[j].dctcorr);
      else if (pdata[j].dctcorr<0.99) fprintf(of," %8.2e",pdata[j].dctcorr);
      else if (pdata[j].dctcorr<9.999995) fprintf(of," %8.6f",pdata[j].dctcorr);
      else if (pdata[j].dctcorr<99.99995) fprintf(of," %8.5f",pdata[j].dctcorr);
      else if (pdata[j].dctcorr<999.9995) fprintf(of," %8.4f",pdata[j].dctcorr);
      else if (pdata[j].dctcorr<9999.995) fprintf(of," %8.3f",pdata[j].dctcorr);
      else if (pdata[j].dctcorr<99999.95) fprintf(of," %8.2f",pdata[j].dctcorr);
      else if (pdata[j].dctcorr<999999.5) fprintf(of," %8.1f",pdata[j].dctcorr);
      else fprintf(of," %8.1e",pdata[j].dctcorr);

#ifdef USEWFPC2
      if (hstmode[j].inst==WFPC2) WFPC2outstarimg(j,of,t.x,t.y,pdata);
      else
#endif
#ifdef USEACS
      if (hstmode[j].inst==ACS) ACSoutstarimg(j,of,t.x,t.y,pdata);
      else
#endif
#ifdef USEWFC3
      if (hstmode[j].inst==WFC3) WFC3outstarimg(j,of,t.x,t.y,pdata);
      else
#endif
      {
	 fprintf(of," %6.3f",pdata[j].m);
	 if (pdata[j].dm>9.999) fprintf(of," 9.999");
	 else fprintf(of," %5.3f",pdata[j].dm);
      }

      fprintf(of," %6.2f %7.1f %6.3f %6.3f %5.3f",pdata[j].chi,pdata[j].ct/pdata[j].dct,pdata[j].sh,pdata[j].rnd,pdata[j].crowd);
      if (hstmode[j].inst==NONE) fprintf(of," %6.3f %6.3f %6.3f %6.3f %6.3f",pdata[j].r,pdata[j].e,pdata[j].a,pdata[j].b,pdata[j].c);
      fprintf(of," %2d",pdata[j].flag);
   }
   fprintf(of,"\n");
   return 1;
}

typedef struct {double x,y,sn; int t;} realstartype;

chiptype *data0;
int Nrealstar;
realstartype *realstar;

Inline int bestmatch(double x0,double y0) {
   int i=-1,j;
   float rmin,r;

   rmin=SQR(FakeMatch);
   for (j=0;j<Nstars;j++) {
      r=SQR(stars[j].x-x0)+SQR(stars[j].y-y0);
      if (r<rmin) {i=j; rmin=r;}
   }
   return i;
}

void fakestar(int ext,int z,double x0,double y0,double xw,double yw,int tw,double *ct0) {
   int img,ix,iy,x1,y1,XMIN20,XMAX20,YMIN20,YMAX20,i,j,it,cont,xmin,xmax,xsize,ymin,ymax,cx,cy,nbg,anyPosOK;
   float pa,pb,pc;
   double x,y,r,rmin,SNRIN=0.;
   static float*bg=NULL;

   // eliminate stars with illegal positions
   if (x0<XMIN2 || x0>XMAX2 || y0<YMIN2 || y0>YMAX2) return;
   anyPosOK = 0;
   for (img=0;img<Nimg;img++) {
      shift(img,x0,y0,&x,&y,1);
      if (x>FakePad && x<dataim[img].X-FakePad && y>FakePad && y<dataim[img].Y-FakePad) anyPosOK = 1;
   }
   if (anyPosOK==0) return;

   if (bg==NULL) {
      bg=(float*)calloc(Nimg,FLOATSIZE);
      if (!bg) merr();
   }
   XMIN20=XMIN2;
   XMAX20=XMAX2;
   YMIN20=YMIN2;
   YMAX20=YMAX2;
   XMIN2=(int)x0-rMarkPSF;
   XMAX2=(int)x0+rMarkPSF+1;
   YMIN2=(int)y0-rMarkPSF;
   YMAX2=(int)y0+rMarkPSF+1;
   XMIN=XMIN2-rMark; if (XMIN<0) XMIN=0;
   XMAX=XMAX2+rMark; if (XMAX>X) XMAX=X;
   YMIN=YMIN2-rMark; if (YMIN<0) YMIN=0;
   YMAX=YMAX2+rMark; if (YMAX>Y) YMAX=Y;
   for (img=0;img<Nimg;img++) {
      shift(img,x0,y0,&x,&y,1);
      ix=(int)(x+RPSF[img])-RPSF[img];
      iy=(int)(y+RPSF[img])-RPSF[img];
      getpsfpars(img,xw,yw,(&pa)-img,(&pb)-img,(&pc)-img);
      bg[img]=0.;
      nbg=0;
      for (y1=-RPSF[img];y1<=RPSF[img] && iy+y1<dataim[img].Y;y1++) if (iy+y1>=0) for (x1=-RPSF[img];x1<=RPSF[img] && ix+x1<dataim[img].X;x1++) if (ix+x1>=0 && dataOK(img,ix+x1,iy+y1) && y1>=-rphot[img] && y1<=rphot[img] && x1>=-rphot[img] && x1<=rphot[img]) {
	 bg[img]+=data[img][iy+y1][ix+x1];
	 nbg++;
      }
      if (nbg) bg[img]/=nbg;
      else bg[img]=skyval(img,ix,iy);
#ifdef USEWFPC2
      if (hstmode[img].inst==WFPC2) WFPC2fixfakemag(img,x0,y0,ct0,bg);
#endif
#ifdef USEACS
      if (hstmode[img].inst==ACS) ACSfixfakemag(img,x0,y0,ct0,bg);
#endif
#ifdef USEWFC3
      if (hstmode[img].inst==WFC3) WFC3fixfakemag(img,x0,y0,ct0,bg);
#endif
      if (ix>=-rphot[img] && ix<dataim[img].X+rphot[img] && iy>=-rphot[img] && iy<dataim[img].Y+rphot[img]) {
	 cx=(int)((x-ix)*50+0.5);
	 if (cx<0) cx=0; if (cx>50) cx=50;
	 cy=(int)((y-iy)*50+0.5);
	 if (cy<0) cy=0; if (cy>50) cy=50;
	 calc1psf(img,x,y,rphot[img],pa,pb,pc,1,1);
	 SNRIN+=SQR(eval1_snr(img,ix,iy,cx,cy,ct0[img]));
      }
      calc1psf(img,x,y,RPSF[img],pa,pb,pc,1,0);
      for (y1=-RPSF[img];y1<=RPSF[img] && iy+y1<dataim[img].Y;y1++) if (iy+y1>=0) for (x1=-RPSF[img];x1<=RPSF[img] && ix+x1<dataim[img].X;x1++) if (ix+x1>=0) {
	 data0[img][y1+RPSF[img]][x1+RPSF[img]]=data[img][iy+y1][ix+x1];
	 if (dataOK(img,ix+x1,iy+y1)) {
	    if (RandomFake) data[img][iy+y1][ix+x1]+=poiss(ct0[img]*psf[y1][x1]);
	    else data[img][iy+y1][ix+x1]+=ct0[img]*psf[y1][x1];
	 }
      }
      shift(img,XMIN,YMIN,&x,&y,1);
      xmin=xmax=(int)(x+RPSF[img])-RPSF[img];
      ymin=ymax=(int)(y+RPSF[img])-RPSF[img];
      shift(img,XMIN,YMAX,&x,&y,1);
      ix=(int)(x+RPSF[img])-RPSF[img]; if (ix<xmin) xmin=ix; if (ix>xmax) xmax=ix;
      iy=(int)(y+RPSF[img])-RPSF[img]; if (iy<ymin) ymin=iy; if (iy>ymax) ymax=iy;
      shift(img,XMAX,YMIN,&x,&y,1);
      ix=(int)(x+RPSF[img])-RPSF[img]; if (ix<xmin) xmin=ix; if (ix>xmax) xmax=ix;
      iy=(int)(y+RPSF[img])-RPSF[img]; if (iy<ymin) ymin=iy; if (iy>ymax) ymax=iy;
      shift(img,XMAX,YMAX,&x,&y,1);
      ix=(int)(x+RPSF[img])-RPSF[img]; if (ix<xmin) xmin=ix; if (ix>xmax) xmax=ix;
      iy=(int)(y+RPSF[img])-RPSF[img]; if (iy<ymin) ymin=iy; if (iy>ymax) ymax=iy;
      xmin-=RPSF[img];
      ymin-=RPSF[img];
      xmax+=RPSF[img];
      ymax+=RPSF[img];
      if (xmin<0) xmin=0;
      if (xmax>=dataim[img].X) xmax=dataim[img].X-1;
      if (ymin<0) ymin=0;
      if (ymax>=dataim[img].Y) ymax=dataim[img].Y-1;
      xsize=1+xmax-xmin;
      if (xsize>0) for (y1=ymin;y1<=ymax;y1++) memcpy(res[img][y1]+xmin,data[img][y1]+xmin,xsize*FLOATSIZE);
   }
   SNRIN=sqrt(SNRIN);
   fitpsf=1;
   Nstars=0;
   if (WARMSTART) {
      setindx(-1);
      for (i=0;i<Nrealstar;i++) if (realstar[i].x>=XMIN && realstar[i].x<XMAX && realstar[i].y>=YMIN && realstar[i].y<YMAX) {
	 stars[Nstars].x=realstar[i].x;
	 stars[Nstars].y=realstar[i].y;
	 stars[Nstars].type=realstar[i].t;
	 getpsfpars(-1,xw,yw,stars[Nstars].pa,stars[Nstars].pb,stars[Nstars].pc);
	 if (Timg>Nimg) getpsfpars(Nimg,xw,yw,stars[Nstars].pa,stars[Nstars].pb,stars[Nstars].pc);
	 phot((int)xw,(int)yw,Nstars,0);
	 imgsub(Nstars);
	 Nstars++;
      }
      stars[Nstars].x=xw;
      stars[Nstars].y=yw;
      if (tw==2) stars[Nstars].type=1;
      else stars[Nstars].type=tw;
      getpsfpars(-1,xw,yw,stars[Nstars].pa,stars[Nstars].pb,stars[Nstars].pc);
      if (Timg>Nimg) getpsfpars(Nimg,xw,yw,stars[Nstars].pa,stars[Nstars].pb,stars[Nstars].pc);
      phot((int)xw,(int)yw,Nstars,0);
      imgsub(Nstars);
      Nstars++;
   }
   else find(0);
   sortstars();
   if (FitSky==1) for (i=0;i<Nstars;i++) stars[i].flag=1;
   it=0;
   cont=1;
   fitpsf=0;
   while (cont && it<MaxIT) {
      cont=solve(ext,z,++it,0);
      if (cont>1 && it>MaxIT/2) it=MaxIT/2;
   }
   sortstars();
   if (WARMSTART) i=Nstars-1;
   else i=bestmatch(x0,y0);
   if (i>=0 && (stars[i].x<XMIN2 || stars[i].x>=XMAX2 || stars[i].y<YMIN2 || stars[i].y>=YMAX2)) i=-1;
   if (i>=0 && !WARMSTART && (stars[i].s<=0 || stars[i].s<SigFinal*stars[i].ss)) i=-1;
   if (i>=0 && !WARMSTART) {
      rmin=SNRIN*exp(-0.5*(SQR(stars[i].x-x0)+SQR(stars[i].y-y0))/FakePSF);
      for (j=0;j<Nrealstar && i>=0;j++) {
	 r=realstar[j].sn*exp(-0.5*(SQR(stars[i].x-realstar[j].x)+SQR(stars[i].y-realstar[j].y))/FakePSF);
	 if (r>rmin && bestmatch(realstar[j].x,realstar[j].y)==i) i=-1;
      }
   }
   fprintf(ffakeout,"%d %d %7.2f %7.2f ",ext,z+1,x0,y0);
   for (img=0;img<Nimg;img++) {
      if (ct0[img]<999999.5) fprintf(ffakeout,"%8.1f ",ct0[img]);
      else fprintf(ffakeout,"%8.2e ",ct0[img]);
#ifdef USEWFPC2
      if (hstmode[img].inst==WFPC2) fprintf(ffakeout,"%6.3f ",WFPC2calcmag(img,x0,y0,ct0[img],bg[img],UseCTE));
      else
#endif
#ifdef USEACS
      if (hstmode[img].inst==ACS) fprintf(ffakeout,"%6.3f ",ACScalcmag(img,x0,y0,ct0[img],bg[img],UseCTE));
      else
#endif
#ifdef USEWFC3
      if (hstmode[img].inst==WFC3) fprintf(ffakeout,"%6.3f ",WFC3calcmag(img,x0,y0,ct0[img],bg[img],UseCTE));
      else
#endif
	 fprintf(ffakeout,"%6.3f ",-2.5*log10(ct0[img]*apcor[img]/iEXP[img])+Zero);
   }
   if (i>=0) {
      imgadd(i);
      photsearch(-1,&stars[i].x,&stars[i].y,stars[i].pa,stars[i].pb,stars[i].pc,&stars[i].s,&stars[i].ss,&stars[i].chi,&stars[i].sh,&stars[i].sky,stars[i].is,stars[i].iss,stars[i].icm,&stars[i].type,-1,i);
      out1star(ext,z,i,ffakeout,1);
   }
   else out1star(ext,z,i,ffakeout,1);
   XMIN2=XMIN20;
   XMAX2=XMAX20;
   YMIN2=YMIN20;
   YMAX2=YMAX20;
   for (img=0;img<Nimg;img++) {
      shift(img,x0,y0,&x,&y,1);
      ix=(int)(x+RPSF[img])-RPSF[img];
      iy=(int)(y+RPSF[img])-RPSF[img];
      for (y1=-RPSF[img];y1<=RPSF[img] && iy+y1<dataim[img].Y;y1++) if (iy+y1>=0) for (x1=-RPSF[img];x1<=RPSF[img] && ix+x1<dataim[img].X;x1++) if (ix+x1>=0) data[img][iy+y1][ix+x1]=data0[img][y1+RPSF[img]][x1+RPSF[img]];
   }
   fflush(ffakeout);
   return;
}

#define MAX_LINE_LENGTH 30000
void fakestars(int ext,int fld) {
   FILE *f;
   int e,z,img,tw=1,i;
   double x,y,xw=0,yw=0,*ctin;
   char str[MAX_LINE_LENGTH];

   Nrealstar=i=0;
   if ((f=fopen(outfn,"r"))==NULL) {
      printf("****Error opening original photometry file %s\n",outfn);
      exit(-1);
   }
   while (fscanf(f,"%d %d",&e,&z)==2) {
      if (e==ext && z==fld+1) Nrealstar++;
      fgets(str,MAX_LINE_LENGTH,f);
   }
   printf("Read %d real stars\n",Nrealstar);
   fflush(stdout);
   realstar=(realstartype*)calloc(Nrealstar,sizeof(realstartype)); if (!realstar) merr();
   fseek(f,0,SEEK_SET);
   while (fscanf(f,"%d %d",&e,&z)==2) {
      if (e==ext && z==fld+1) {
	 fscanf(f,"%lf %lf %lf %lf %lf %lf %lf %lf %d",&(realstar[i].x),&(realstar[i].y),&x,&(realstar[i].sn),&x,&x,&x,&x,&(realstar[i].t));
	 i++;
      }
      fgets(str,MAX_LINE_LENGTH,f);
   }
   fclose(f);
   if (i!=Nrealstar) {
      printf("****Programming error; i!=Nrealstar\n");
      exit(-1);
   }
   if ((f=fopen(FakeStars,"r"))==NULL) {
      printf("****Error opening FakeStars file %s\n",FakeStars);
      exit(-1);
   }
   ctin=(double*)calloc(Nimg,DOUBLESIZE); if (!ctin) merr();
   data0=(chiptype*)calloc(Nimg,sizeof(chiptype)); if (!data0) merr();
   for (img=0;img<Nimg;img++) {
      //cleansat(img);
      data0[img]=allocchip(2*RPSF[img]+1,2*RPSF[img]+1);
   }
   while (fscanf(f,"%d %d",&e,&z)==2) {
      if (e==ext && z==fld+1) {
	 fscanf(f,"%lf %lf",&x,&y);
	 if (WARMSTART) fscanf(f,"%lf %lf %d",&xw,&yw,&tw);
#ifdef USEWFPC2
	 WFPC2readfakemag(f);
#endif
#ifdef USEACS
	 ACSreadfakemag(f);
#endif
#ifdef USEWFC3
	 WFC3readfakemag(f);
#endif
	 for (img=0;img<Nimg;img++) if (hstmode[img].inst==NONE) fscanf(f,"%lf",ctin+img);
	 fgets(str,MAX_LINE_LENGTH,f);
	 fakestar(ext,fld,x,y,xw,yw,tw,ctin);
      }
      else fgets(str,MAX_LINE_LENGTH,f);
   }
   fclose(f);
   for (img=0;img<Nimg;img++) freechip(data0[img],2*RPSF[img]+1,2*RPSF[img]+1);
   free(ctin);
   free(data0);
   free(realstar);
   return;
}
#undef MAX_LINE_LENGTH

void outpt(FILE*of,int ext,int fld) {
   int i,Nelim[2]={0,0};

   sortstars();
   for (i=0;i<Nstars;i++) if (stars[i].x>=XMIN2 && stars[i].x<XMAX2 && stars[i].y>=YMIN2 && stars[i].y<YMAX2 && (WARMSTART || (stars[i].s>0 && stars[i].s>=SigFinal*stars[i].ss))) {
      imgadd(i);
      photsearch(-1,&stars[i].x,&stars[i].y,stars[i].pa,stars[i].pb,stars[i].pc,&stars[i].s,&stars[i].ss,&stars[i].chi,&stars[i].sh,&stars[i].sky,stars[i].is,stars[i].iss,stars[i].icm,&stars[i].type,-1,i);
      imgsub(i);
   }
   for (i=0;i<Nstars;i++) {
      imgadd(i);
      Nelim[out1star(ext,fld,i,of,0)]++;
      imgsub(i);
   }
   if (WARMSTART==2) for (i=0;i<Nstars;i++) imgadd(i);
   printf("%d stars written; %d stars deleted\n",Nelim[1],Nelim[0]);
   fflush(stdout);
   fflush(of);
   return;
}

void outptinfo(void) {
   int j;

   fprintf(finfo,"%d sets of output data\n",Nimg);
   for (j=0;j<Nimg;j++) {
      fprintf(finfo,"%s\n",base[j]);
      fprintf(finfo,"%16.8f\n",iEPOCH[j]);
   }
   fprintf(finfo,"\n");
   fflush(finfo);
   return;
}

void outptalign(void) {
   int i,j;
   fprintf(finfo,"Alignment\n");
   for (i=0;i<Nimg;i++) {
      if (Align<4) {for (j=0;j<5;j++) fprintf(finfo," %f",dpos[i][j]);}
      else {for (j=0;j<40;j++) fprintf(finfo," %12.6e",dpos[i][j]);}
      if (UseWCS==2) {
	 for (j=0;j<4;j++) fprintf(finfo," %12.6e",wcsref[i][j]);
	 for (j=0;j<80;j++) fprintf(finfo," %12.6e",wcs[i][j]);
      }
      fprintf(finfo,"\n");
   }
   fflush(finfo);
   return;
}

void outptapcor(void) {
   int i;
   fprintf(finfo,"Aperture corrections\n");
   for (i=0;i<Nimg;i++) fprintf(finfo,"  %f\n",-2.5*log10(apcor[i]));
   fprintf(finfo,"\n");
   fflush(finfo);
   return;
}

void procframe(int ext) {
   int img,i,x,y,z,it,cont,skip;
   char ch[2880];
   void *ptr;
   imtype timg;

   IMDIFF=0;
   if (!isimage(dataim)) {
      for (img=0;img<Timg;img++) {
	 fopenagain(fdata+img);
	 readimage(fdata[img].f,dataim+img, fdata[img].use_mmap);
	 freclose(fdata+img);
	 if (!FakeStars[0] && img<Nimg) {
	    fopenagain(fres+img);
	    writeimage(fres[img].f,dataim+img);
	    freclose(fres+img);
	    if (!UsePhot[0]) {
	       fopenagain(fpsf+img);
	       writeimage(fpsf[img].f,dataim+img);
	       freclose(fpsf+img);
	    }
	 }
	 freeimg(dataim[img].data,dataim[img].X,dataim[img].Y,dataim[img].Z);
      }
      return;
   }
   if (Timg>Nimg) {
      X=dataim[Nimg].X;
      Y=dataim[Nimg].Y;
   }
   else {
      X=dataim[0].X;
      Y=dataim[0].Y;
   }
#ifdef PGPLOT
   allocDiagData2(ext);
#endif
   ran=(char**)calloc(Y*SubResRef,PTRSIZE);
   snmapflag=(char**)calloc(Y*SubResRef,PTRSIZE);
   snmap=(float**)calloc(Y*SubResRef,PTRSIZE);
   snmmap=(float**)calloc(Y*SubResRef,PTRSIZE);
   indx=(int**)calloc(Y*SubResRef,PTRSIZE);
   tindx=(int*)calloc(X*SubResRef,INTSIZE);
   if (!ran || !snmapflag || !snmap || !snmmap || !indx || !tindx) merr();
   ptr=malloc(X*Y*SubResRef*SubResRef); if (!ptr) merr();
   for (y=0;y<Y*SubResRef;y++) ran[y]=ptr+(X*SubResRef)*y;
   ptr=malloc(X*Y*SubResRef*SubResRef); if (!ptr) merr();
   for (y=0;y<Y*SubResRef;y++) snmapflag[y]=ptr+X*SubResRef*y;
   ptr=malloc(X*Y*SubResRef*SubResRef*FLOATSIZE); if (!ptr) merr();
   for (y=0;y<Y*SubResRef;y++) snmap[y]=ptr+X*SubResRef*y*FLOATSIZE;
   ptr=malloc(X*Y*SubResRef*SubResRef*FLOATSIZE); if (!ptr) merr();
   for (y=0;y<Y*SubResRef;y++) snmmap[y]=ptr+X*SubResRef*y*FLOATSIZE;
   ptr=malloc(X*Y*SubResRef*SubResRef*INTSIZE); if (!ptr) merr();
   for (y=0;y<Y*SubResRef;y++) indx[y]=ptr+(X*SubResRef)*y*INTSIZE;
   for (img=0;img<Timg;img++) {
       // FIXME -- full image data array allocated here
       if (fdata[img].use_mmap) {
           data[img] = (float**)calloc(dataim[img].Y, sizeof(float*));
       } else {
           data[img]=allocchip(dataim[img].X,dataim[img].Y);
       }
      if (fsky[img].lastoffset>=0) sky[img]=allocchip(dataim[img].X,dataim[img].Y);
      res[img]=allocchip(dataim[img].X,dataim[img].Y);
   }
   if (GUSE<0 && CUSE<0) {
      XMIN2=0;
      XMAX2=X;
      YMIN2=0;
      YMAX2=Y;
   }
   for (z=0;z<dataim[0].Z;z++) {
#ifdef PGPLOT
      setDiagFrame(ext,z);
#endif
      initimgpars();
      for (img=0;img<Timg;img++) {
	 if (!FakeStars[0] && !UsePhot[0]) {
	    for (i=0;i<5;i++) dpos[img][i]=dpos0[img][i];
	 }
      }
      // update params for camera
#ifdef USEWFPC2
      wfpc2initparam();
#endif
#ifdef USEACS
      acsinitparam();
#endif
#ifdef USEWFC3
      wfc3initparam();
#endif
#ifdef USEWFPC2
      wfpc2radii(FakeStars[0]==0 && UsePhot[0]==0); // accounts for plate scale vs. chip # and shifts
#endif
      initcirc();
#ifdef USEACS
      acsinitpsf();
#endif
#ifdef USEWFC3
      wfc3initpsf();
#endif
#ifdef USEWFPC2
      wfpc2initpsf();
#endif
      if (WARMSTART==2 && refpsf.Next>=0) {
	 imtype *timg;
	 if (ext==0) timg=&(refpsf.img);
	 else timg=refpsf.ext+ext-1;
	 if (timg->X!=2*RPSF[Timg-1]+1 || timg->Y!=2*RPSF[Timg-1]+1 || timg->Z!=dataim[0].Z) {
	    printf("**Provided PSF for subtraction is the wrong size; using analytic solution only\n");
	    refpsfimg=NULL;
	 }
	 else refpsfimg=timg->data[z];
      }
      if (Nimg==Timg) DRIZZLE_BASE=0;
      else if (hstmode[Nimg].inst==NONE) DRIZZLE_BASE=1;
      out1starinfo();
      for (img=0;img<Timg;img++) {
          printf("read data %i\n", img);
	 fopenagain(fdata+img);
	 readchip(fdata[img].f,data[img], dataim+img, fdata[img].use_mmap);
	 freclose(fdata+img);
	 if (fsky[img].lastoffset>=0) {
          printf("read sky %i\n", img);
	    fopenagain(fsky+img);
	    readchip(fsky[img].f,sky[img],dataim+img, fdata[img].use_mmap);
	    freclose(fsky+img);
	 }
      }
      if ((GUSE<0 || ext==GUSE) && (CUSE<0 || z==CUSE)) {
	 if (ext>0) printf("** Extension %d, Chip %d **\n",ext,z+1);
	 else printf("** Chip %d **\n",z+1);
	 fflush(stdout);
	 if (!FakeStars[0] && VerboseData>0) fprintf(fverb,"EXTENSION %d CHIP %d\n",ext,z+1);
	 if (!FakeStars[0] && !UsePhot[0]) {
	    fprintf(finfo,"EXTENSION %d CHIP %d\n",ext,z+1);
	    fprintf(finfo,"Limits\n %d %d %d %d\n",XMIN2,XMAX2,YMIN2,YMAX2);
#ifdef USEWFPC2
	    writewfpc2info();
#endif
#ifdef USEACS
	    writeacsinfo();
#endif
#ifdef USEWFC3
	    writewfc3info();
#endif
	    fflush(finfo);
	 }
	 else readinfo(finfo,ext,z);
	 for (img=0;img<Timg;img++) {
	    if (iDMIN[img]<iDMIN0[img]) {
	       float *pd,*pr,*plast;
	       plast=data[img][0]+dataim[img].Y*dataim[img].X;
	       for (pd=data[img][0],pr=res[img][0];pd<plast;pd++,pr++) {
               if (*pd<=iDMIN0[img]) {
                   *pd=*pr=iDMIN[img]-1;
               } else {
                   *pr=*pd;
               }
	       }
	    }
	    else memcpy(res[img][0],data[img][0],dataim[img].X*dataim[img].Y*FLOATSIZE);
	    if (!FakeStars[0] && !UsePhot[0]) {
	       apcor[img]=1.;
	       for (y=-RPSF[img];y<=RPSF[img];y++) for (x=-RPSF[img];x<=RPSF[img];x++) poff[img][y][x]=0;
	    }
	    else if (img<Nimg) {
	       float* pad;
	       int nwrite = RPSF[img];
	       if (imgdata[img].RPSF<nwrite) nwrite = imgdata[img].RPSF;
	       int npad = imgdata[img].RPSF - nwrite;
	       pad = (float*)calloc(imgdata[img].RPSF*2+1,sizeof(float));
	       fopenagain(fpsf+img);
	       for (y=0;y<npad;y++) ffread(pad,4,imgdata[img].RPSF*2+1,fpsf[img].f);
	       for (y=-nwrite;y<=nwrite;y++) {
		  if (npad>0) ffread(pad,4,npad,fpsf[img].f);
		  ffread(poff[img][y]-nwrite,4,nwrite*2+1,fpsf[img].f);
		  if (npad>0) ffread(pad,4,npad,fpsf[img].f);
	       }
	       for (y=0;y<npad;y++) ffread(pad,4,imgdata[img].RPSF*2+1,fpsf[img].f);
	       freclose(fpsf+img);
	       free(pad);
	       if (!FakeStarPSF) for (y=-RPSF[img];y<=RPSF[img];y++) for (x=-RPSF[img];x<=RPSF[img];x++) poff[img][y][x]=0;
	    }
	    else {
	       for (y=-RPSF[img];y<=RPSF[img];y++) for (x=-RPSF[img];x<=RPSF[img];x++) poff[img][y][x]=0;
	    }
	 }
	 if (!FakeStars[0] && !UsePhot[0]) {
	    XMIN=XMIN2;
	    if (XMIN<0) XMIN=0;
	    XMAX=XMAX2;
	    if (XMAX>X) XMAX=X;
	    YMIN=YMIN2;
	    if (YMIN<0) YMIN=0;
	    YMAX=YMAX2;
	    if (YMAX>Y) YMAX=Y;
	    fitpsf=1;
	    if (UseWCS==2) convertWCS2();
	    else if (UseWCS==1) convertWCS1();
	    setRef2Img();
	    if (Align==4) convertAlign4();
	    if (Align && Timg>1) align(ext,z);
	    outptalign();
	 }
	 if (AlignOnly==0) {
	    setRMark();
	    XMIN=XMIN2-rMark;
	    if (XMIN<0) XMIN=0;
	    XMAX=XMAX2+rMark;
	    if (XMAX>X) XMAX=X;
	    YMIN=YMIN2-rMark;
	    if (YMIN<0) YMIN=0;
	    YMAX=YMAX2+rMark;
	    if (YMAX>Y) YMAX=Y;
	    if (FakeStars[0]) fakestars(ext,z);
	    else {
	       if (UsePhot[0]) fitpsf=0;
	       Nstars=0;
	       if (WARMSTART) readwarm(ext,z);
	       else find(1);
	       sortstars();
	       if (FitSky==1) for (i=0;i<Nstars;i++) stars[i].flag=1;
	       it=0;
	       cont=1;
	       while (cont && it<MaxIT) {
		  cont=solve(ext,z,++it,1);
		  if (cont>1 && it>MaxIT/2) it=MaxIT/2;
	       }
	       if (UsePhot[0]==0) {
		  if (ApCor) startapcor();
		  for (img=0;img<Nimg;img++) {
		     //cleansat(img);
		     if (ApCor) getapcor(ext,z,img,fapcor);
		  }
		  outptapcor();
	       }
	       outpt(of,ext,z);
	    }
	 }
      }
      // still need to advance through files
      else if (FakeStars[0] || UsePhot[0]) {
	 //readinfo(finfo,ext,z); // doesn't write finfo for skipped frames
	 for (img=0;img<Nimg;img++) {
	    float* pad;
	    int nwrite = RPSF[img];
	    if (imgdata[img].RPSF<nwrite) nwrite = imgdata[img].RPSF;
	    int npad = imgdata[img].RPSF - nwrite;
	    pad = (float*)calloc(imgdata[img].RPSF*2+1,sizeof(float));
	    fopenagain(fpsf+img);
	    for (y=0;y<npad;y++) ffread(pad,4,imgdata[img].RPSF*2+1,fpsf[img].f);
	    for (y=-nwrite;y<=nwrite;y++) {
	       if (npad>0) ffread(pad,4,npad,fpsf[img].f);
	       ffread(poff[img][y]-nwrite,4,nwrite*2+1,fpsf[img].f);
	       if (npad>0) ffread(pad,4,npad,fpsf[img].f);
	    }
	    for (y=0;y<npad;y++) ffread(pad,4,imgdata[img].RPSF*2+1,fpsf[img].f);
	    free(pad);
	    freclose(fpsf+img);
	 }
      }
      if (!FakeStars[0]) for (img=0;img<Nimg;img++) {
	 if (!UsePhot[0]) {
	    float* pad;
	    int nwrite = RPSF[img];
	    if (imgdata[img].RPSF<nwrite) nwrite = imgdata[img].RPSF;
	    int npad = imgdata[img].RPSF - nwrite;
	    pad = (float*)calloc(imgdata[img].RPSF*2+1,FLOATSIZE); if (!pad) merr();
	    memset(pad,0,FLOATSIZE*(imgdata[img].RPSF*2+1));
	    fopenagain(fpsf+img);
	    for (y=0;y<npad;y++) ffwrite(pad,FLOATSIZE,imgdata[img].RPSF*2+1,fpsf[img].f);
	    for (y=-nwrite;y<=nwrite;y++) {
	       if (npad>0) ffwrite(pad,FLOATSIZE,npad,fpsf[img].f);
	       ffwrite(poff[img][y]-nwrite,FLOATSIZE,nwrite*2+1,fpsf[img].f);
	       if (npad>0) ffwrite(pad,FLOATSIZE,npad,fpsf[img].f);
	    }
	    for (y=0;y<npad;y++) ffwrite(pad,FLOATSIZE,imgdata[img].RPSF*2+1,fpsf[img].f);
	    free(pad);
	    freclose(fpsf+img);
	 }
	 memcpy(&timg,dataim+img,sizeof(imtype));
	 timg.bscale=1.;
	 timg.bzero=0.;
	 timg.bits=-32;
	 fopenagain(fres+img);
	 writechip(fres[img].f,res[img],&timg);
	 freclose(fres+img);
      }
      freecirc();
   }
   memset(ch,0,sizeof(ch));
   for (img=0;img<Timg;img++) {
      skip=fabs(dataim[img].bits)/8*dataim[img].X*dataim[img].Y*dataim[img].Z+dataim[img].pcount;
      skip=((skip+2879)/2880)*2880-skip;
      if (fsky[img].lastoffset>=0) fsky[img].lastoffset+=skip;
      fdata[img].lastoffset+=skip;
      if (img<Nimg) {
	 if (!FakeStars[0]) {
	    skip=FLOATSIZE*dataim[img].X*dataim[img].Y*dataim[img].Z+dataim[img].pcount;
	    skip=((skip+2879)/2880)*2880-skip;
	    fopenagain(fres+img);
	    fwrite(ch,1,skip,fres[img].f);
	    freclose(fres+img);
	 }
	 skip=FLOATSIZE*(2*imgdata[img].RPSF+1)*(2*imgdata[img].RPSF+1)*dataim[img].Z+dataim[img].pcount;
	 skip=((skip+2879)/2880)*2880-skip;
	 fopenagain(fpsf+img);
	 if (!FakeStars[0] && !UsePhot[0]) fwrite(ch,1,skip,fpsf[img].f);
	 else fread(ch,1,skip,fpsf[img].f);
	 freclose(fpsf+img);
	 freechip(data[img],dataim[img].X,dataim[img].Y);
	 if (fsky[img].lastoffset>=0) freechip(sky[img],dataim[img].X,dataim[img].Y);
	 freechip(res[img],dataim[img].X,dataim[img].Y);
      }
   }
#ifdef USEWFPC2
   wfpc2freepsf();
#endif
#ifdef USEACS
   acsfreepsf();
#endif
#ifdef USEWFC3
   wfc3freepsf();
#endif
   free(indx[0]);
   free(snmmap[0]);
   free(snmap[0]);
   free(snmapflag[0]);
   free(ran[0]);
   free(tindx);
   free(indx);
   free(snmmap);
   free(snmap);
   free(snmapflag);
   free(ran);
   return;
}

int main(int argc,char**argv) {
   char str[82];
   int img,Next=0,ext,i;

   if (argc<2) {
      printf("****Usage: %s <output> <<options>>\n",*argv);
      printf("  -p<name>  for parameter file\n");
      //printf("  -m        to use mmap\n");
      printf("  x=y       to set flag x to value y\n");
      return 1;
   }
   INTSIZE=sizeof(int);
   LONGSIZE=sizeof(long);
   FLOATSIZE=sizeof(float);
   DOUBLESIZE=sizeof(double);
   PTRSIZE=sizeof(void*);
   if (sizeof(char)!=1) {printf("char size != 1\n"); return -1;}
   if (INTSIZE<4) {printf("int size < 4\n"); return -1;}
   if (LONGSIZE<4) {printf("long size < 4\n"); return -1;}
   if (FLOATSIZE!=4) {printf("float size != 4\n"); return -1;}
   if (DOUBLESIZE<8) {printf("double size < 8\n"); return -1;}
   if (sizeof(char*)!=PTRSIZE) {printf("char* size != void* size\n"); return -1;}
   if (sizeof(int*)!=PTRSIZE) {printf("int* size != void* size\n"); return -1;}
   if (sizeof(float*)!=PTRSIZE) {printf("float* size != void* size\n"); return -1;}
   if (sizeof(double*)!=PTRSIZE) {printf("double* size != void* size\n"); return -1;}
   if (sizeof(FILE*)!=PTRSIZE) {printf("FILE* size != void* size\n"); return -1;}

   // read parameters
   initimgdata();
   paramfile("dolphot.param",&dolphotparam);
   for (i=2;i<argc;i++) {
      if (!strncmp(argv[i],"-p",2) || !strncmp(argv[i],"-P",2)) paramfile1(argv[i]+2,&dolphotparam);
      //if (!strncmp(argv[i],"-m",2) || !strncmp(argv[i],"-M",2)) use_mmap = 1;
      else parseparam(argv[i],&dolphotparam);
   }
   // sanity check
   if (RPSF0<RAper0) perr("RPSF must be at least as large as RAper");

   // allocate memory
   alloc_img_data();

   fdata=(reopenableFile*)calloc(Timg,sizeof(reopenableFile));
   fsky=(reopenableFile*)calloc(Timg,sizeof(reopenableFile));
   fpsf=(reopenableFile*)calloc(Timg,sizeof(reopenableFile));
   fres=(reopenableFile*)calloc(Timg,sizeof(reopenableFile));
   if (!fdata || !fsky || !fpsf || !fres) merr();
   strcpy(outfn,argv[1]);
   if (FakeStars[0]) {
      sprintf(str,"%s.fake",outfn);
      ffakeout=fopen(str,"w");
      if (!ffakeout) {
	 printf("****Error opening %s\n",str);
	 exit(-1);
      }
      sprintf(str,"%s.info",outfn);
      finfo=fopen(str,"r");
      if (!finfo) {
	 printf("****Error opening %s\n",str);
	 exit(-1);
      }
   }
   else if (UsePhot[0]) {
      of=fopen(outfn,"w");
      if (!of) {
	 printf("****Error opening %s\n",outfn);
	 exit(-1);
      }
      sprintf(str,"%s.info",UsePhot);
      finfo=fopen(str,"r");
      if (!finfo) {
	 printf("****Error opening %s\n",str);
	 exit(-1);
      }
      sprintf(str,"%s.warnings",outfn);
      fwarn=fopen(str,"w");
      if (!fwarn) {
	 printf("****Error opening %s\n",str);
	 exit(-1);
      }
   }
   else {
      of=fopen(outfn,"w");
      if (!of) {
	 printf("****Error opening %s\n",outfn);
	 exit(-1);
      }
      sprintf(str,"%s.info",outfn);
      finfo=fopen(str,"w");
      if (!finfo) {
	 printf("****Error opening %s\n",str);
	 exit(-1);
      }
      sprintf(str,"%s.apcor",outfn);
      fapcor=fopen(str,"w");
      if (!fapcor) {
	 printf("****Error opening %s\n",str);
	 exit(-1);
      }
      sprintf(str,"%s.psfs",outfn);
      fpsfs=fopen(str,"w");
      if (!fpsfs) {
	 printf("****Error opening %s\n",str);
	 exit(-1);
      }
      sprintf(str,"%s.warnings",outfn);
      fwarn=fopen(str,"w");
      if (!fwarn) {
	 printf("****Error opening %s\n",str);
	 exit(-1);
      }
      if (VerboseData>0) {
	 sprintf(str,"%s.data",outfn);
	 fverb=fopen(str,"w");
	 if (!fverb) {
	    printf("****Error opening %s\n",str);
	    exit(-1);
	 }
      }
   }
   if (WARMSTART==2) {
      if (xytpsf[0] && !access(xytpsf,F_OK)) readfits(xytpsf,&refpsf,1);
      else {
	 refpsf.Next=-1;
	 printf("Cannot open template PSF file; ignoring\n");
      }
   }

   for (img=0;img<Timg;img++) {
      ftype tfits;
      imtype skyim;

      sprintf(str,"%s.sky.fits",base[img]);
      if (!access(str,F_OK)) {
	 fopenfirst(fsky+img,str,"rb",0);
	 fsky[img].f=readfitsh(str,&tfits,0);
	 freclose(fsky+img);
	 memcpy(&skyim,&(tfits.img),sizeof(imtype));
	 if (tfits.img.Nmax) free(tfits.img.cards);
      }
      else {
	 fsky[img].lastoffset=-1;
	 skyim.Z=0;
      }
      sprintf(str,"%s.fits",base[img]);
      fopenfirst(fdata+img,str,"rb",0);
      fdata[img].f=readfitsh(str,&tfits,1);
      freclose(fdata+img);
      memcpy(dataim+img,&(tfits.img),sizeof(imtype));
      datahd[img].Nmax=datahd[img].Ncards=dataim[img].Ncards;
      datahd[img].cards=(cardtype*)calloc(datahd[img].Nmax,sizeof(cardtype));
      memcpy(datahd[img].cards,dataim[img].cards,sizeof(cardtype)*datahd[img].Ncards);
      read_cardvals(img);
      if (fsky[img].lastoffset>=0 && (skyim.X!=dataim[img].X || skyim.Y!=dataim[img].Y || skyim.Z!=dataim[img].Z)) {
	 printf("****Sky does not match the image\n");
	 exit(-1);
      }
      if (img==0) Next=tfits.Next;
      else if (Next!=tfits.Next) {
	 printf("****Number of extensions are not the same\n");
	 exit(-1);
      }
      if (imgdata[img].RPSF<=0) imgdata[img].RPSF = RPSF0;
      if (img<Nimg) {
	 if (!FakeStars[0]) {
	    tfits.img.Ncards=dataim[img].Ncards;
	    tfits.img.Nmax=dataim[img].Ncards+1;
	    tfits.img.cards=(cardtype*)calloc(tfits.img.Nmax,sizeof(cardtype));
	    memcpy(tfits.img.cards,dataim[img].cards,sizeof(cardtype)*tfits.img.Ncards);
	    sprintf(str,"%s.%d.res.fits",outfn,img+1);
	    if (isimage(&(tfits.img))) {
	       tfits.img.bscale=1.;
	       tfits.img.bzero=0.;
	       tfits.img.bits=-32;
	       if (iDMIN[img]<iDMIN0[img]) insertcards(&(tfits.img),-1.e30,-1.e30,-1.e30,iDMIN[img],-1.e30,-1.e30,-1.e30,-1.e30);
	    }
	    fopenfirst(fres+img,str,"ab",0);
	    fres[img].f=writefitsh(str,&tfits,0);
	    freclose(fres+img);
	 }
	 if (!FakeStars[0] && !UsePhot[0]) {
	    if (isimage(&(tfits.img))) {
	       tfits.img.X=2*imgdata[img].RPSF+1;
	       tfits.img.Y=2*imgdata[img].RPSF+1;
	    }
	    sprintf(str,"%s.%d.psf.fits",outfn,img+1);
	    fopenfirst(fpsf+img,str,"ab",0);
	    fpsf[img].f=writefitsh(str,&tfits,0);
	    freclose(fpsf+img);
	    free(tfits.img.cards);
	 }
	 else {
	    if (FakeStars[0]) sprintf(str,"%s.%d.psf.fits",outfn,img+1);
	    else sprintf(str,"%s.%d.psf.fits",UsePhot,img+1);
	    fopenfirst(fpsf+img,str,"rb",0);
	    fpsf[img].f=readfitsh(str,&tfits,0);
	    if (isimage(&(tfits.img))) {
	       if (tfits.img.X!=2*imgdata[img].RPSF+1 || tfits.img.Y!=2*imgdata[img].RPSF+1) {
		  printf("Mismatch in RPSF\n");
		  exit(-1);
	       }
	       if (tfits.img.Nmax) free(tfits.img.cards);
	       tfits.img.Nmax=0;
	    }
	    else {
            readimage(fpsf[img].f,&(tfits.img), 0);
	       freeimg(tfits.img.data,tfits.img.X,tfits.img.Y,tfits.img.Z);
	    }
	    freclose(fpsf+img);
	 }
      }
   }
   if (!FakeStars[0] && !UsePhot[0]) outptinfo();
#ifdef PGPLOT
   allocDiagData1(Next);
#endif

   // set up circles and photometry weights, etc.
   procframe(0);
   for (img=0;img<Timg;img++) if (dataim[img].Nmax) free(dataim[img].cards);
   for (ext=0;ext<Next;ext++) {
      for (img=0;img<Timg;img++) {
	 ftype tfits;
	 imtype skyim;

	 if (fsky[img].lastoffset>=0) {
	    fopenagain(fsky+img);
	    readexth(fsky[img].f,&skyim,0);
	    freclose(fsky+img);
	    if (skyim.Nmax) free(skyim.cards);
	 }
	 fopenagain(fdata+img);
	 readexth(fdata[img].f,&(tfits.img),1);
	 freclose(fdata+img);
	 memcpy(dataim+img,&(tfits.img),sizeof(imtype));
	 read_cardvals(img);
	 if (fsky[img].lastoffset>=0 && (skyim.X!=dataim[img].X || skyim.Y!=dataim[img].Y || skyim.Z!=dataim[img].Z)) {
	    printf("****Sky does not match the image\n");
	    exit(-1);
	 }
	 if (img!=0 && dataim[0].Z!=dataim[img].Z) {
	    printf("****Number of chips are not the same size\n");
	    exit(-1);
	 }
	 if (img<Nimg) {
	    if (!FakeStars[0]) {
	       tfits.img.Ncards=dataim[img].Ncards;
	       tfits.img.Nmax=dataim[img].Ncards+1;
	       tfits.img.cards=(cardtype*)calloc(tfits.img.Nmax,sizeof(cardtype));
	       memcpy(tfits.img.cards,dataim[img].cards,sizeof(cardtype)*tfits.img.Ncards);
	       if (isimage(&(tfits.img))) {
		  tfits.img.bscale=1.;
		  tfits.img.bzero=0.;
		  tfits.img.bits=-32;
		  if (iDMIN[img]<iDMIN0[img]) insertcards(&(tfits.img),-1.e30,-1.e30,-1.e30,iDMIN[img],-1.e30,-1.e30,-1.e30,-1.e30);
	       }
	       fopenagain(fres+img);
	       writeexth(fres[img].f,&(tfits.img),0);
	       freclose(fres+img);
	    }
	    if (!FakeStars[0] && !UsePhot[0]) {
	       if (isimage(&(tfits.img))) {
		  tfits.img.X=2*imgdata[img].RPSF+1;
		  tfits.img.Y=2*imgdata[img].RPSF+1;
	       }
	       fopenagain(fpsf+img);
	       writeexth(fpsf[img].f,&(tfits.img),0);
	       freclose(fpsf+img);
	       free(tfits.img.cards);
	    }
	    else {
	       fopenagain(fpsf+img);
	       readexth(fpsf[img].f,&(tfits.img),0);
	       if (isimage(&(tfits.img))) {
		  if (tfits.img.X!=2*imgdata[img].RPSF+1 || tfits.img.Y!=2*imgdata[img].RPSF+1) {
		     printf("Mismatch in RPSF\n");
		     exit(-1);
		  }
	       }
	       else {
               readimage(fpsf[img].f,&(tfits.img), 0);
		  freeimg(tfits.img.data,tfits.img.X,tfits.img.Y,tfits.img.Z);
	       }
	       freclose(fpsf+img);
	    }
	 }
      }
      procframe(ext+1);
      for (img=0;img<Timg;img++) if (dataim[img].Nmax) free(dataim[img].cards);
   }
   fclose(finfo);
   if (FakeStars[0]) fclose(ffakeout);
   else {
      fclose(of);
      fclose(fwarn);
      if (!UsePhot[0]) {
	 if (VerboseData>0) fclose(fverb);
	 fclose(fapcor);
	 fclose(fpsfs);
      }
   }
#ifdef PGPLOT
   plotDiagData();
#endif
   return 0;
}
