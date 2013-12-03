#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <fits.h>
#include "../dolphot_defs.h"
#include "wfc3psfdata.h"
#include "wfc3filters.h"
#include "wfc3distort.h"

static float********wfc3psflib=NULL;
static int WFC3_RPSF=34;
static int ANY_WFC3=0;

extern void shift(int img,double x0,double y0,double *x,double *y,int dir);

void wfc3initparam(void) {
   char str[161],*ptr;
   int img;
   ANY_WFC3=0;
   for (img=0;img<Timg && !ANY_WFC3;img++) {
      strcpy(str,getcardval(dataim+img,"DOL_WFC3",0));
      strtol(str,&ptr,10);
      if (ptr!=str) ANY_WFC3=1;
   }
   if (ANY_WFC3==0) return;
   FPSF=4;
   SubPixel=1;
   EPSF=0;
   PSFsol=-1;
   PSFStep=0.;
   Zero=0.;
   MinS=0.2;
   MaxS=15.;
   MaxE=0.1;
   for (img=0;img<Timg;img++) {
      if (RPSF[img]>=wfc3_rpsf[1]) {RPSF[img]=wfc3_rpsf[1]-1; printf("Lowering RPSF to %d\n",RPSF[img]);}
      if (RAper[img]>wfc3_rpsf[1]-1.) {RAper[img]=wfc3_rpsf[1]-1; printf("Lowering RAper to %d\n",wfc3_rpsf[1]-1);}
   }
   return;
}

void wfc3initpsf(void) {
   int img,i,j,y1,x1,y2,x2,y3,n2,n3,n3skip;
   char str[161],*ptr;
   FILE *f;
   float***** ptr1;
   float**** ptr2;
   float*** ptr3;
   float** ptr4;
   float* ptr5;

   if (WFC3_NFILTERS<0) WFC3initfilters();
   if (wfc3psflib==NULL) {
      wfc3psflib=(float********)calloc(sizeof(float*******),3);
      if (!wfc3psflib) merr();
      for (i=0;i<3;i++) {
	 wfc3psflib[i]=(float*******)calloc(sizeof(float******),WFC3_NFILTERS);
	 if (!wfc3psflib[i]) merr();
      }
   }
   for (i=0;i<3;i++) for (j=0;j<WFC3_NFILTERS;j++) wfc3psflib[i][j]=NULL;
   for (img=0;img<Timg;img++) {
      strcpy(str,getcardval(dataim+img,"DOL_WFC3",0));
      int cm=strtol(str,&ptr,10);
      if (ptr==str) {
	 // AEDDEBUG need warning: printf("**Image %d has not been preprocessed with wfc3mask; cannot proceed\n",img+1);
      }
      else {
	 if (cm<-2 || cm>2) {
	    printf("**Image %d's chip cannot be identified; please report bug.\n",img+1);
	    exit(-1);
	 }
	 hstmode[img].inst=WFC3;
	 hstmode[img].cm=cm;
	 if (img==Nimg) {
	    if (hstmode[img].cm==-2) { // UVIS
	       hstmode[img].cm=1;
	       DRIZZLE_BASE=1;
	    }
	    else if (hstmode[img].cm==-1) { // IR
	       hstmode[img].cm=0;
	       DRIZZLE_BASE=1;
	    }
	    else DRIZZLE_BASE=0;
	 }
	 if (hstmode[img].cm<0) {
	    printf("**Image %d has been drizzled; cannot use for photometry\n",img+1);
	    exit(-1);
	 }
	 strcpy(str,getcardval(datahd+img,"FILTER",0));
	 hstmode[img].filt=WFC3findfilt(str);
	 //printf("Image %d: cm=%d, filt=%d (%s)\n",img+1,hstmode[img].cm,hstmode[img].filt,WFC3filters[hstmode[img].filt].name);
	 if (RPSF[img]>=wfc3_rpsf[hstmode[img].cm]) {printf("ERROR: RPSF must be less than %d for WFC3/WFC data\n",wfc3_rpsf[hstmode[img].cm]); exit(-1);}
	 if (RAper[img]>wfc3_rpsf[hstmode[img].cm]-1.) {printf("ERROR: RAper must be no more than %d for WFC3/WFC data\n",wfc3_rpsf[hstmode[img].cm]-1); exit(-1);}
      }
   }
   WFC3_RPSF=0;
   for (img=0;img<Timg;img++) if (hstmode[img].inst==WFC3) {
      i=hstmode[img].cm;
      apsf[img][0][0]=1.;
      apsf[img][1][0]=1.;
      apsf[img][2][0]=0.;
      if (i==0) apsize[img]=4.;
      else apsize[img]=12.;
      for (i=1;i<5;i++) apsf[img][0][i]=apsf[img][1][i]=apsf[img][2][i]=0.;
      if (WFC3_RPSF<RPSF[img]) WFC3_RPSF=RPSF[img];
      if (WFC3_RPSF<rphot[img]) WFC3_RPSF=rphot[img];
   }
   for (img=0;img<Timg;img++) if (hstmode[img].inst==WFC3 && wfc3psflib[i=hstmode[img].cm][j=hstmode[img].filt]==NULL) {
      f=0;
      if (WFC3psfType[i]==1) {
	 sprintf(str,"%s/wfc3/data/%s_anderson.%s.psf",BASEDIR,WFC3filters[hstmode[img].filt].name,wfc3_cn[i]);
	 f=fopen(str,"rb");
      }
      if (f==0) {
	 sprintf(str,"%s/wfc3/data/%s.%s.psf",BASEDIR,WFC3filters[hstmode[img].filt].name,wfc3_cn[i]);
	 if ((f=fopen(str,"rb"))==NULL) {
	    printf("Cannot open %s\n",str);
	    exit(-1);
	 }
      }
      n2=2*wfc3_n2psf[i]+1;
      n3=2*WFC3_RPSF+1;
      n3skip=wfc3_rpsf[i]-WFC3_RPSF;
      wfc3psflib[i][j]=(float******)calloc(sizeof(float*****),wfc3_nypsfpos[i]);
      ptr1=(float*****)calloc(sizeof(float****),wfc3_nypsfpos[i]*wfc3_nxpsfpos[i]);
      ptr2=(float****)calloc(sizeof(float***),wfc3_nypsfpos[i]*wfc3_nxpsfpos[i]*n2);
      ptr3=(float***)calloc(sizeof(float**),wfc3_nypsfpos[i]*wfc3_nxpsfpos[i]*n2*n2);
      ptr4=(float**)calloc(sizeof(float*),wfc3_nypsfpos[i]*wfc3_nxpsfpos[i]*n2*n2*n3);
      ptr5=(float*)calloc(sizeof(float),wfc3_nypsfpos[i]*wfc3_nxpsfpos[i]*n2*n2*n3*n3);
      if (!wfc3psflib[i][j] || !ptr1) merr();
      for (y1=0;y1<wfc3_nypsfpos[i];y1++) {
	 wfc3psflib[i][j][y1]=ptr1;
	 ptr1+=wfc3_nxpsfpos[i];
	 for (x1=0;x1<wfc3_nxpsfpos[i];x1++) {
	    wfc3psflib[i][j][y1][x1]=ptr2+wfc3_n2psf[i];
	    ptr2+=n2;
	    for (y2=-wfc3_n2psf[i];y2<=wfc3_n2psf[i];y2++) {
	       wfc3psflib[i][j][y1][x1][y2]=ptr3+wfc3_n2psf[i];
	       ptr3+=n2;
	       for (x2=-wfc3_n2psf[i];x2<=wfc3_n2psf[i];x2++) {
		  wfc3psflib[i][j][y1][x1][y2][x2]=ptr4+WFC3_RPSF;
		  ptr4+=n3;
		  if (n3skip) fseek(f,4*n3skip*(2*wfc3_rpsf[i]+1),SEEK_CUR);
		  for (y3=-WFC3_RPSF;y3<=WFC3_RPSF;y3++) {
		     wfc3psflib[i][j][y1][x1][y2][x2][y3]=ptr5+WFC3_RPSF;
		     ptr5+=n3;
		     if (n3skip) fseek(f,4*n3skip,SEEK_CUR);
		     ffread(wfc3psflib[i][j][y1][x1][y2][x2][y3]-WFC3_RPSF,4,n3,f);
		     if (n3skip) fseek(f,4*n3skip,SEEK_CUR);
		  }
		  if (n3skip) fseek(f,4*n3skip*(2*wfc3_rpsf[i]+1),SEEK_CUR);
	       }
	    }
	 }
      }
      fclose(f);
   }
   return;
}

void wfc3freepsf(void) {
   int i,j;
   for (i=0;i<3;i++) for (j=0;j<WFC3_NFILTERS;j++) if (wfc3psflib[i][j]) {
      free(wfc3psflib[i][j][0][0][-wfc3_n2psf[i]][-wfc3_n2psf[i]][-WFC3_RPSF]-WFC3_RPSF);
      free(wfc3psflib[i][j][0][0][-wfc3_n2psf[i]][-wfc3_n2psf[i]]-WFC3_RPSF);
      free(wfc3psflib[i][j][0][0][-wfc3_n2psf[i]]-wfc3_n2psf[i]);
      free(wfc3psflib[i][j][0][0]-wfc3_n2psf[i]);
      free(wfc3psflib[i][j][0]);
      free(wfc3psflib[i][j]);
      wfc3psflib[i][j]=NULL;
   }
   return;
}

/*
void calcwfc3psf(int img,float x,float y,int r) {
   int i,j,y1,x1,y2,x2,yy,xx;
   float mx1,my1,imx1,imy1;
   float mx2,my2,imx2,imy2;
   static int first=1,lastr=0,lastimg=0;
   static float lastx=0,lasty=0;

   if (!first && lastpsftype==1 && img==lastimg && x==lastx && y==lasty && r<=lastr) return;
   first=0;
   lastpsftype=1;
   lastimg=img;
   lastx=x;
   lasty=y;
   lastr=r;
   i=hstmode[img].cm;
   j=hstmode[img].filt;
   //printf("%s/%s PSF at %f,%f; r=%d:\n",wfc3_cn[i],WFC3filters[j].name,x,y,r);
   y1=(int)(y-128)/256; if (y1<0) y1=0; if (y1>=wfc3_nypsfpos[i]-1) y1=wfc3_nypsfpos[i]-2;
   imy1=(y-128-256*y1)/256.; my1=1-imy1;
   x1=(int)(x-128)/256; if (x1<0) x1=0; if (x1>=wfc3_nxpsfpos[i]-1) x1=wfc3_nxpsfpos[i]-2;
   imx1=(x-128-256*x1)/256.; mx1=1-imx1;
   y=(y-(int)(y+50)+49.5)*wfc3_sub[i];
   y2=(int)(y+wfc3_sub[i])-wfc3_sub[i];
   imy2=y-y2; my2=1-imy2;
   x=(x-(int)(x+50)+49.5)*wfc3_sub[i];
   x2=(int)(x+wfc3_sub[i])-wfc3_sub[i];
   imx2=x-x2; mx2=1-imx2;
   //printf("y1=%d, my1=%f; x1=%d, mx1=%f\n",y1,my1,x1,mx1);
   //printf("y2=%d, my2=%f; x2=%d, mx2=%f\n",y2,my2,x2,mx2);
   for (yy=-r;yy<=r;yy++) for (xx=-r;xx<=r;xx++) {
      psf[yy][xx]=0.;
      psf[yy][xx]+=(wfc3psflib[i][j][y1][x1][y2][x2][yy][xx]*mx2*my2+wfc3psflib[i][j][y1][x1][y2][x2+1][yy][xx]*imx2*my2+wfc3psflib[i][j][y1][x1][y2+1][x2][yy][xx]*mx2*imy2+wfc3psflib[i][j][y1][x1][y2+1][x2+1][yy][xx]*imx2*imy2)*mx1*my1;
      psf[yy][xx]+=(wfc3psflib[i][j][y1][x1+1][y2][x2][yy][xx]*mx2*my2+wfc3psflib[i][j][y1][x1+1][y2][x2+1][yy][xx]*imx2*my2+wfc3psflib[i][j][y1][x1+1][y2+1][x2][yy][xx]*mx2*imy2+wfc3psflib[i][j][y1][x1+1][y2+1][x2+1][yy][xx]*imx2*imy2)*imx1*my1;
      psf[yy][xx]+=(wfc3psflib[i][j][y1+1][x1][y2][x2][yy][xx]*mx2*my2+wfc3psflib[i][j][y1+1][x1][y2][x2+1][yy][xx]*imx2*my2+wfc3psflib[i][j][y1+1][x1][y2+1][x2][yy][xx]*mx2*imy2+wfc3psflib[i][j][y1+1][x1][y2+1][x2+1][yy][xx]*imx2*imy2)*mx1*imy1;
      psf[yy][xx]+=(wfc3psflib[i][j][y1+1][x1+1][y2][x2][yy][xx]*mx2*my2+wfc3psflib[i][j][y1+1][x1+1][y2][x2+1][yy][xx]*imx2*my2+wfc3psflib[i][j][y1+1][x1+1][y2+1][x2][yy][xx]*mx2*imy2+wfc3psflib[i][j][y1+1][x1+1][y2+1][x2+1][yy][xx]*imx2*imy2)*imx1*imy1;
   }
   //for (yy=6;yy>=-6;yy--) if (yy>=-r && yy<=r) {for (xx=-6;xx<=6;xx++) if (xx>=-r && xx<=r) printf("%5d ",(int)(psf[yy][xx]*100000+0.5)); printf("\n");} fflush(stdout);
   return;
}
*/

int calcwfc3psf(int img,float x,float y,int r,int force) {
   int i,j,y1,x1,y2,x2,yy,xx;
   float mx2,my2,imx2,imy2;
   static int first=1,lastr=0,lastimg=0;
   static float lastx=0,lasty=0;

   if (hstmode[img].inst!=WFC3) {
      printf("Stupid error; called wfc3psf for non-WFC3 data\n");
      exit(-1);
   }
   if (!first && lastpsftype==1 && img==lastimg && x==lastx && y==lasty && r<=lastr && !poffreset && !force) return 0;
   first=0;lastpsftype=1;lastimg=img;lastx=x;lasty=y;lastr=r;
   i=hstmode[img].cm;
   j=hstmode[img].filt;
   //printf("%d %d %d %f %f %d",img,i,j,x,y,r); fflush(stdout);
   //printf("%s/%s PSF at %f,%f; r=%d:\n",wfc3_cn[i],WFC3filters[j].name,x,y,r);
   if (img==Nimg && DRIZZLE_BASE) {
      y1=wfc3_nypsfpos[i]/2-1;
      x1=wfc3_nxpsfpos[i]/2-1;
   }
   else {
      y1=(int)y/wfc3_psfspacing[i]; if (y1<0) y1=0; if (y1>=wfc3_nypsfpos[i]) y1=wfc3_nypsfpos[i]-1;
      x1=(int)x/wfc3_psfspacing[i]; if (x1<0) x1=0; if (x1>=wfc3_nxpsfpos[i]) x1=wfc3_nxpsfpos[i]-1;
   }
   y-=(int)y; if (y<0) y++;
   y=(y-0.5)*wfc3_sub[i];
   y2=(int)(y+wfc3_sub[i])-wfc3_sub[i];
   imy2=y-y2; my2=1-imy2;
   //printf(" (%f %d %d %f %f)",y,y1,y2,my2,imy2); fflush(stdout);
   x-=(int)x; if (x<0) x++;
   x=(x-0.5)*wfc3_sub[i];
   x2=(int)(x+wfc3_sub[i])-wfc3_sub[i];
   imx2=x-x2; mx2=1-imx2;
   //printf(" (%f %d %d %f %f)\n",x,x1,x2,mx2,imx2); fflush(stdout);
   //printf("y1=%d, my1=%f; x1=%d, mx1=%f\n",y1,my1,x1,mx1);
   //printf("y2=%d, my2=%f; x2=%d, mx2=%f\n",y2,my2,x2,mx2);
   for (yy=-r;yy<=r;yy++) for (xx=-r;xx<=r;xx++) {
      psf[yy][xx]=wfc3psflib[i][j][y1][x1][y2][x2][yy][xx]*mx2*my2+wfc3psflib[i][j][y1][x1][y2][x2+1][yy][xx]*imx2*my2+wfc3psflib[i][j][y1][x1][y2+1][x2][yy][xx]*mx2*imy2+wfc3psflib[i][j][y1][x1][y2+1][x2+1][yy][xx]*imx2*imy2;
   }
   //for (yy=6;yy>=-6;yy--) if (yy>=-r && yy<=r) {for (xx=-6;xx<=6;xx++) if (xx>=-r && xx<=r) printf("%5d ",(int)(psf[yy][xx]*100000+0.5)); printf("\n");} fflush(stdout);
   return 1;
}

double WFC3calcmag(int img,float x0,float y0,float ct0,float bg,int useCTE) {
   float cm=1.;
   double x,y,m;

   if (hstmode[img].inst!=WFC3) {
      printf("Stupid error; called wfc3calcmag for non-WFC3 data\n");
      exit(-1);
   }
   if (WFC3_NFILTERS<0) WFC3initfilters();
   m=-2.5*log10(ct0*apcor[img]/iEXP[img]*wfc3_ctmult[hstmode[img].cm])+WFC3_ZP(hstmode[img].filt,hstmode[img].cm);
   if (useCTE) {
      shift(img,x0,y0,&x,&y,1);
      if (iEXP0[img]>0.) cm=iEXP[img]/iEXP0[img];
      m-=WFC3_CTE(hstmode[img].cm,x,y,ct0,cm,iGAIN[img],bg,iEPOCH[img]);
   }
   return m;
}

void WFC3outstarinfo(FILE *f,int*ct) {
   int i,j,n;
   for (i=0;i<WFC3_NFILTERS;i++) {
      n=0;
      for (j=0;j<Nimg;j++) if (hstmode[j].inst==WFC3 && hstmode[j].filt==i) n++;
      if (n>1) {
	 fprintf(f,"%d. Total counts, %s\n",++(*ct),WFC3filters[i].name);
	 fprintf(f,"%d. Total sky level, %s\n",++(*ct),WFC3filters[i].name);
	 fprintf(f,"%d. Normalized count rate, %s\n",++(*ct),WFC3filters[i].name);
	 fprintf(f,"%d. Normalized count rate uncertainty, %s\n",++(*ct),WFC3filters[i].name);
	 fprintf(f,"%d. Instrumental VEGAMAG magnitude, %s\n",++(*ct),WFC3filters[i].name);
	 fprintf(f,"%d. Transformed UBVRI magnitude, %s\n",++(*ct),WFC3filters[i].name);
	 fprintf(f,"%d. Magnitude uncertainty, %s\n",++(*ct),WFC3filters[i].name);
	 fprintf(f,"%d. Chi, %s\n",++(*ct),WFC3filters[i].name);
	 fprintf(f,"%d. Signal-to-noise, %s\n",++(*ct),WFC3filters[i].name);
	 fprintf(f,"%d. Sharpness, %s\n",++(*ct),WFC3filters[i].name);
	 fprintf(f,"%d. Roundness, %s\n",++(*ct),WFC3filters[i].name);
	 fprintf(f,"%d. Crowding, %s\n",++(*ct),WFC3filters[i].name);
	 fprintf(f,"%d. Photometry quality flag, %s\n",++(*ct),WFC3filters[i].name);
      }
   }
}

static double*smag=NULL,*vmag=NULL,*dvmag=NULL;
void WFC3outstar(FILE *of,float x0,float y0,photdatatype*pdata) {
   int i,j;
   static int *fused=NULL;
   static photdatatype*fphot=NULL;
   double x,y,m=1.0;

   if (WFC3_NFILTERS<0) WFC3initfilters();
   if (!fphot) {
      fused=(int*)calloc(sizeof(int),WFC3_NFILTERS);
      fphot=(photdatatype*)calloc(sizeof(photdatatype),WFC3_NFILTERS);
      smag=(double*)calloc(sizeof(double),WFC3_NFILTERS);
      vmag=(double*)calloc(sizeof(double),WFC3_NFILTERS);
      dvmag=(double*)calloc(sizeof(double),WFC3_NFILTERS);
      if (!fused || !fphot || !smag || !vmag || !dvmag) merr();
   }
   for (j=0;j<Nimg;j++) if (hstmode[j].inst==WFC3) {
      double dm = -2.5*log10(wfc3_ctmult[hstmode[j].cm]) + WFC3_ZP(hstmode[j].filt,hstmode[j].cm)-Zero;
      if (pdata[j].ct>0) {
	 if (UseCTE) {
	    shift(j,x0,y0,&x,&y,1);
	    if (iEXP0[j]>0.) m=iEXP[j]/iEXP0[j];
	    else m=1.0;
	    dm -= WFC3_CTE(hstmode[j].cm,x,y,pdata[j].ct,m,iGAIN[j],pdata[j].sky,iEPOCH[j]);
	 }
	 pdata[j].m += dm;
      }
      pdata[j].ctcorr *= pow(10,-0.4*dm);
      pdata[j].dctcorr *= pow(10,-0.4*dm);
   }
   for (i=0;i<WFC3_NFILTERS;i++) {
      float wt,swt,twt=0,tswt=0,is,iss,tcm=0.;
      fphot[i].ct0=fphot[i].ct=fphot[i].chi=fphot[i].sh=fphot[i].sky=fphot[i].ctcorr=fphot[i].dctcorr=fphot[i].rnd=fphot[i].crowd=0.;
      fused[i]=0;
      fphot[i].flag=0;
      for (j=0;j<Nimg;j++) if (hstmode[j].inst==WFC3 && hstmode[j].filt==i) {
	 fused[i]++;
	 if (pdata[j].flag<8 && !(pdata[j].flag&FlagMask)) {
	    is=pdata[j].ct/iEXP[j];
	    iss=pdata[j].dct/iEXP[j];
	    wt=1./iss/iss;
	    if (is>0) swt=wt*is;
	    else swt=0.0;
	    twt+=wt;
	    tswt+=swt;
	    fphot[i].ct0+=pdata[j].ct0/iEXP[j]*wt;
	    fphot[i].ct+=is*wt;
	    fphot[i].chi+=pdata[j].chi*pdata[j].chi*wt;
	    fphot[i].sh+=pdata[j].sh*swt;
	    fphot[i].sky+=pdata[j].sky/iEXP[j]*wt;
	    fphot[i].ctcorr+=pdata[j].ctcorr*wt;
	    fphot[i].dctcorr+=pdata[j].dctcorr*pdata[j].dctcorr*wt*wt;
	    fphot[i].rnd+=pdata[j].rnd*swt;
	    fphot[i].crowd+=pdata[j].crowd*wt;
	    fphot[i].flag|=pdata[j].flag;
	 }
	 tcm+=iEXP[j];
      }
      if (twt>0.) {
	 fphot[i].ct0/=twt/tcm;
	 fphot[i].ct/=twt/tcm;
	 fphot[i].dct=tcm/sqrt(twt);
	 fphot[i].chi=sqrt(fphot[i].chi/twt);
	 fphot[i].sky/=twt/tcm;
	 fphot[i].ctcorr/=twt;
	 fphot[i].dctcorr=sqrt(fphot[i].dctcorr)/twt;
	 if (fphot[i].ctcorr>0) fphot[i].m=-2.5*log10(fphot[i].ctcorr);
	 else fphot[i].m=99.999;
	 fphot[i].crowd/=twt;
	 if (fphot[i].ct>0) fphot[i].dm=1.0857362*fphot[i].dct/fphot[i].ct;
      }
      else {
	 fphot[i].dct=9999;
	 fphot[i].dctcorr=9999;
	 fphot[i].m=99.999;
	 fphot[i].dm=9.999;
      }
      if (tswt>0.) {
	 fphot[i].sh/=tswt;
	 fphot[i].rnd/=tswt;
      }
      vmag[i]=fphot[i].m;
      dvmag[i]=fphot[i].dm;
   }
   WFC3transform(vmag,dvmag,smag);
   for (i=0;i<WFC3_NFILTERS;i++) if (fused[i]>1) {
      if (fphot[i].ct<999999.5) fprintf(of,"  %8.1f",fphot[i].ct);
      else fprintf(of,"  %8.2e",fphot[i].ct);
      if (fphot[i].sky<99999.95) fprintf(of," %8.2f",fphot[i].sky);
      else if (fphot[i].sky<999999.5) fprintf(of," %8.1f",fphot[i].sky);
      else fprintf(of," %8.1e",fphot[i].sky);

      if (fphot[i].ctcorr<0.0) fprintf(of," %8.1e",fphot[i].ctcorr);
      else if (fphot[i].ctcorr<0.99) fprintf(of," %8.2e",fphot[i].ctcorr);
      else if (fphot[i].ctcorr<9.999995) fprintf(of," %8.6f",fphot[i].ctcorr);
      else if (fphot[i].ctcorr<99.99995) fprintf(of," %8.5f",fphot[i].ctcorr);
      else if (fphot[i].ctcorr<999.9995) fprintf(of," %8.4f",fphot[i].ctcorr);
      else if (fphot[i].ctcorr<9999.995) fprintf(of," %8.3f",fphot[i].ctcorr);
      else if (fphot[i].ctcorr<99999.95) fprintf(of," %8.2f",fphot[i].ctcorr);
      else if (fphot[i].ctcorr<999999.5) fprintf(of," %8.1f",fphot[i].ctcorr);
      else fprintf(of," %8.1e",fphot[i].ctcorr);

      if (fphot[i].dctcorr<0.0) fprintf(of," %8.1e",fphot[i].dctcorr);
      else if (fphot[i].dctcorr<0.99) fprintf(of," %8.2e",fphot[i].dctcorr);
      else if (fphot[i].dctcorr<9.999995) fprintf(of," %8.6f",fphot[i].dctcorr);
      else if (fphot[i].dctcorr<99.99995) fprintf(of," %8.5f",fphot[i].dctcorr);
      else if (fphot[i].dctcorr<999.9995) fprintf(of," %8.4f",fphot[i].dctcorr);
      else if (fphot[i].dctcorr<9999.995) fprintf(of," %8.3f",fphot[i].dctcorr);
      else if (fphot[i].dctcorr<99999.95) fprintf(of," %8.2f",fphot[i].dctcorr);
      else if (fphot[i].dctcorr<999999.5) fprintf(of," %8.1f",fphot[i].dctcorr);
      else fprintf(of," %8.1e",fphot[i].dctcorr);

      fprintf(of," %6.3f %6.3f",fphot[i].m,smag[i]);
      if (fphot[i].dm>9.999) fprintf(of," 9.999");
      else fprintf(of," %5.3f",fphot[i].dm);
      fprintf(of," %6.2f %7.1f %6.3f %6.3f %5.3f %2d",fphot[i].chi,fphot[i].ct/fphot[i].dct,fphot[i].sh,fphot[i].rnd,fphot[i].crowd,fphot[i].flag);
   }
}

void WFC3outstarimg(int img,FILE *of,float x0,float y0,photdatatype*pdata) {
   float x;
   if (hstmode[img].inst!=WFC3) {
      printf("Stupid error; called outstarimg for non-WFC3 data\n");
      exit(-1);
   }
   if (pdata[img].ctcorr<=0) x=pdata[img].m;
   else if (smag[hstmode[img].filt]>99) x=smag[hstmode[img].filt];
   else x=pdata[img].m+smag[hstmode[img].filt]-vmag[hstmode[img].filt];
   fprintf(of," %6.3f %6.3f",pdata[img].m,x);
   if (pdata[img].dm>9.999) fprintf(of," 9.999");
   else fprintf(of," %5.3f",pdata[img].dm);
   return;
}

float wfc3_apsize(int img,float x,float y) {
   double area;

   if (hstmode[img].inst!=WFC3) {
      printf("Stupid error; called wfc3_apsize for non-WFC3 data\n");
      exit(-1);
   }
   int cm=hstmode[img].cm;
   if (cm==0) {
      if (x<0) x=0;
      else if (x>1014) x=1014;
      if (y<0) y=0;
      else if (y>1014) y=1014;
   }
   else {
      if (x<0) x=0;
      else if (x>4051) x=4051;
      if (y<0) y=0;
      else if (y>1024) y=1024;
   }
   area = WFC3pixsize(cm,hstmode[img].filt,x,y);
   //printf("%d %f %f %f\n",cm,x,y,area);
   if (cm==0) return 3.90/sqrt(area);
   return 12.620/sqrt(area);
}

void WFC3shift(int img,double*x,double*y) {
   double xx,yy;
   if (hstmode[img].inst!=WFC3) {
      printf("Stupid error; called wfc3shift for non-WFC3 data\n");
      exit(-1);
   }
   xx=*x; yy=*y;
   if (hstmode[img].cm>=0) WFC3fwddistort(hstmode[img].cm,hstmode[img].filt,&xx,&yy);
   *x=xx; *y=yy;
   return;
}

void WFC3unshift(int img,double*x,double*y) {
   double xx,yy;
   if (hstmode[img].inst!=WFC3) {
      printf("Stupid error; called wfc3unshift for non-WFC3 data\n");
      exit(-1);
   }
   xx=*x; yy=*y;
   if (hstmode[img].cm>=0) WFC3revdistort(hstmode[img].cm,hstmode[img].filt,&xx,&yy);
   *x=xx; *y=yy;
   return;
}

void writewfc3info(void) {
   int img;
   fprintf(finfo,"* WFC3-specific info\n");
   for (img=0;img<Nimg;img++)if (hstmode[img].inst==WFC3)  {
      fprintf(finfo,"* image %d: %s %d %f\n",img+1,WFC3filters[hstmode[img].filt].name,hstmode[img].cm,iEXP[img]);
      if (fabs(iGAIN[img]-1.)>0.001) printf("ERROR: All WFC3 data should have gain of 1\n");
   }
   return;
}

static int*ffused=NULL;
static double*fakem0;
void WFC3readfakemag(FILE*f) {
   int i,img;

   if (WFC3_NFILTERS<0) WFC3initfilters();
   if (ffused==NULL) {
      ffused=(int*)calloc(sizeof(int),WFC3_NFILTERS);
      fakem0=(double*)calloc(sizeof(double),WFC3_NFILTERS);
      if (!ffused || !fakem0) merr();
      for (img=0;img<Nimg;img++) if (hstmode[img].inst==WFC3) ffused[hstmode[img].filt]++;
   }
   for (i=0;i<WFC3_NFILTERS;i++) if (ffused[i]) fscanf(f,"%lf",fakem0+i);
   return;
}

void WFC3fixfakemag(int img,float x0,float y0,double*ct0,float*bg) {
   int i;
   double dm;

   if (hstmode[img].inst!=WFC3) {
      printf("Stupid error; called wfc3fixfakemag for non-WFC3 data\n");
      exit(-1);
   }
   dm=WFC3calcmag(img,x0,y0,1.0,bg[img],0)-fakem0[hstmode[img].filt];
   ct0[img]=pow(10,0.4*dm);
   for (i=0;i<3;i++) {
      dm=WFC3calcmag(img,x0,y0,ct0[img],bg[img],UseCTE)-fakem0[hstmode[img].filt];
      ct0[img]*=pow(10,0.4*dm);
   }
   return;
}

char *WFC3imagestring(int img) {
   if (hstmode[img].inst!=WFC3) {
      printf("Stupid error; called wfc3imagestring for non-WFC3 data\n");
      exit(-1);
   }
   return WFC3filters[hstmode[img].filt].name;
}
