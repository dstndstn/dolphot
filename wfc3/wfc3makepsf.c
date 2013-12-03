#include <fits.h>
#include "wfc3psfdata.h"

ftype fits;
double dQE=0.00;
double kernel[3][3];

void getkernel(void) {
   int i,y=0;
   char *ptr,*ptr2;
   for (i=0;i<fits.img.Ncards;i++) {
      if (!strncmp(fits.img.cards[i],"COMMENT     ",12)) {
	 if (y>=3) {printf("Too many comment lines\n"); exit(-1);}
	 kernel[y][0]=strtod(fits.img.cards[i]+12,&ptr);
	 if (ptr==fits.img.cards[i]+12) {printf("Not enough numbers in comment line\n"); exit(-1);}
	 kernel[y][1]=strtod(ptr,&ptr2);
	 if (ptr==ptr2) {printf("Not enough numbers in comment line\n"); exit(-1);}
	 kernel[y][2]=strtod(ptr2,&ptr);
	 if (ptr==ptr2) {printf("Not enough numbers in comment line\n"); exit(-1);}
	 strtod(ptr,&ptr2);
	 if (ptr!=ptr2) {printf("Too many numbers in comment line\n"); exit(-1);}
	 //printf("%f %f %f\n",kernel[y][0],kernel[y][1],kernel[y][2]);
	 y++;
      }
   }
   if (y<3) {printf("Not enough comment lines\n"); exit(-1);}
   return;
}

void getpsf(char*filt,int cm) {
   int x,y,z,cx,cy;
   int x1,y1,x2,y2,dx,dy,ct,kx,ky;
   FILE *f;
   char fn[321];
   float **psf;
   double norm,m;

   sprintf(fn,"%s/wfc3/data/%s.%s.psf",BASEDIR,filt,wfc3_cn[cm]);
   if ((f=fopen(fn,"wb"))==NULL) {
      printf("Cannot write %s\n",fn);
      exit(1);
   }
   psf=(float**)calloc(2*wfc3_rpsf[cm]+1,sizeof(float*));
   if (!psf) merr();
   psf+=wfc3_rpsf[cm];
   for (z=-wfc3_rpsf[cm];z<=wfc3_rpsf[cm];z++) {
      psf[z]=(float*)calloc(2*wfc3_rpsf[cm]+1,sizeof(float));
      if (!psf[z]) merr();
      psf[z]+=wfc3_rpsf[cm];
   }
   for (ct=0;ct<wfc3_npsfpos[cm];ct++) {
      if (cm==0) sprintf(fn,"%s/tmp/%s_WFC3/%s%d%d.fits",BASEDIR,filt,wfc3_cn[cm],ct/10,ct%10);
      else sprintf(fn,"%s/tmp/%s_WFC3/%s%d%d%d.fits",BASEDIR,filt,wfc3_cn[cm],ct/100,ct/10%10,ct%10);
      readfits(fn,&fits,0);
      getkernel();
      if (fits.Next!=0 || fits.img.Z!=1) {printf("Illegal format of PSF file\n"); exit(-1);}
      cx=cy=0;
      for (y=0;y<fits.img.Y;y++) for (x=0;x<fits.img.X;x++) if (fits.img.data[0][y][x]>fits.img.data[0][cy][cx]) {cx=x; cy=y;}
      if (wfc3_rpsf[cm]>(cx+1)/wfc3_sub[cm]-2) {printf("Insufficient room on left; max wfc3_rpsf[%d]=%d\n",cm,(cx+1)/wfc3_sub[cm]-2); exit(-1);}
      if (wfc3_rpsf[cm]>(fits.img.X-cx+1)/wfc3_sub[cm]-2) {printf("Insufficient room on right; max wfc3_rpsf[%d]=%d\n",cm,(fits.img.X-cx+1)/wfc3_sub[cm]-2); exit(-1);}
      if (wfc3_rpsf[cm]>(cy+1)/wfc3_sub[cm]-2) {printf("Insufficient room on bottom; max wfc3_rpsf[%d]=%d\n",cm,(cy+1)/wfc3_sub[cm]-2); exit(-1);}
      if (wfc3_rpsf[cm]>(fits.img.Y-cy+1)/wfc3_sub[cm]-2) {printf("Insufficient room on top; max wfc3_rpsf[%d]=%d\n",cm,(fits.img.Y-cy+1)/wfc3_sub[cm]-2); exit(-1);}
      norm=0;
      for (y=0;y<fits.img.Y;y++) for (x=0;x<fits.img.X;x++) norm+=fits.img.data[0][y][x];
      printf("%3d: Center=%d,%d; total PSF=%f\n",ct,cx,cy,norm);
      for (y=-wfc3_n2psf[cm];y<=wfc3_n2psf[cm];y++) for (x=-wfc3_n2psf[cm];x<=wfc3_n2psf[cm];x++) {
	 for (y1=-wfc3_rpsf[cm];y1<=wfc3_rpsf[cm];y1++) for (x1=-wfc3_rpsf[cm];x1<=wfc3_rpsf[cm];x1++) psf[y1][x1]=0;
	 for (y1=0;y1<fits.img.Y;y1++) for (x1=0;x1<fits.img.X;x1++) {
	    dx=x1-cx+x;
	    dy=y1-cy+y;
	    x2=(dx+wfc3_sub[cm]/2+1000*wfc3_sub[cm])/wfc3_sub[cm]-1000;
	    y2=(dy+wfc3_sub[cm]/2+1000*wfc3_sub[cm])/wfc3_sub[cm]-1000;
	    dx-=x2*wfc3_sub[cm];
	    dy-=y2*wfc3_sub[cm];
	    m=1+dQE*(0.5-3.*(float)(dx*dx+dy*dy)/(float)(wfc3_sub[cm]*wfc3_sub[cm]));
	    m*=fits.img.data[0][y1][x1];
	    //ignore normalization since TT PSFs are already normalized;
	    //m/=norm;
	    for (kx=-1;kx<=1;kx++) for (ky=-1;ky<=1;ky++) if (abs(x2+kx)<=wfc3_rpsf[cm] && abs(y2+ky)<=wfc3_rpsf[cm]) psf[y2+ky][x2+kx]+=m*kernel[ky+1][kx+1];
	 }
	 for (y1=-wfc3_rpsf[cm];y1<=wfc3_rpsf[cm];y1++) ffwrite(psf[y1]-wfc3_rpsf[cm],2*wfc3_rpsf[cm]+1,sizeof(float),f);
      }
   }
   for (z=-wfc3_rpsf[cm];z<=wfc3_rpsf[cm];z++) free(psf[z]-wfc3_rpsf[cm]);
   free(psf-wfc3_rpsf[cm]);
   fclose(f);
   return;
}

int main(int argc,char**argv) {
   int usecm[3]={0,0,0},nfilt=0,i,j;
   char filt[100][20];

   if (argc<3) {
      printf("Usage: %s <filters> <chips (0-2)>\n",*argv);
      return -1;
   }
   for (i=1;i<argc;i++) {
      if (argv[i][0]=='F') strcpy(filt[nfilt++],argv[i]);
      else if (argv[i][0]>='0' && argv[i][0]<='2' && argv[i][1]==0) usecm[argv[i][0]-'0']=1;
      else {
	 printf("Unknown option: \"%s\"\n",argv[i]);
	 return -1;
      }
   }
   for (i=0;i<nfilt;i++) for (j=0;j<3;j++) if (usecm[j]) {
      getpsf(filt[i],j);
   }
   return 0;
}
