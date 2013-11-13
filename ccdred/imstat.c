#include <unistd.h>
#include "ccdproc.h"

ftype fits;
float *list,DMIN,DMAX;
int N[3]={0,0,0};

void addallimg(imtype*img,int add) {
   int x,y,z;
   if (!isimage(img)) return;
   parsecards(&(fits.img),NULL,NULL,NULL,&DMIN,&DMAX,NULL,NULL,NULL,0,1);
   if (img!=&(fits.img)) parsecards(img,NULL,NULL,NULL,&DMIN,&DMAX,NULL,NULL,NULL,0,0);
   for (z=0;z<img->Z;z++) for (y=0;y<img->Y;y++) for (x=0;x<img->X;x++) {
      if (img->data[z][y][x]>DMIN && img->data[z][y][x]<DMAX) {
	 if (add) list[N[0]]=img->data[z][y][x];
	 N[0]++;
      }
      else if (img->data[z][y][x]<=DMIN) N[1]++;
      else N[2]++;
   }
}

void addall(void) {
   int i;
   addallimg(&(fits.img),0);
   for (i=0;i<fits.Next;i++) addallimg(fits.ext+i,0);
   list=(float*)calloc(N[0],sizeof(float));
   if (!list) merr();
   N[0]=0;
   addallimg(&(fits.img),1);
   for (i=0;i<fits.Next;i++) addallimg(fits.ext+i,1);
}

void addsecimg(imtype*img,int x1,int x2,int y1,int y2,int z1,int z2,int add) {
   int x,y,z;
   if (!isimage(img)) return;
   parsecards(&(fits.img),NULL,NULL,NULL,&DMIN,&DMAX,NULL,NULL,NULL,0,1);
   if (img!=&(fits.img)) parsecards(img,NULL,NULL,NULL,&DMIN,&DMAX,NULL,NULL,NULL,0,0);
   for (z=z1;z<z2;z++) for (y=y1;y<y2;y++) for (x=x1;x<x2;x++) if (img->data[z][y][x]>DMIN && img->data[z][y][x]<DMAX) {
      if (add) list[N[0]]=img->data[z][y][x];
      N[0]++;
   }
}

void addsec(char*str,int ext) {
   int x1,x2,y1,y2,z1,z2;
   imtype*img;

   if (!ext) img=&(fits.img);
   else if (ext>0 && ext<=fits.Next) img=fits.ext+(ext-1);
   else return;
   parsesec(img,str,&x1,&x2,&y1,&y2,&z1,&z2);
   addsecimg(img,x1,x2,y1,y2,z1,z2,0);
   list=(float*)calloc(N[0],sizeof(float));
   if (!list) merr();
   N[0]=0;
   addsecimg(img,x1,x2,y1,y2,z1,z2,1);
}

void av_sd(double*av,double*sd) {
   int i;
   //for (i=0;i<N[0];i++) printf("%f\n",list[i]);
   *av=*sd=0;
   if (N[0]<=0) return;
   for (i=0;i<N[0];i++) (*av)+=list[i];
   (*av)/=N[0];
   if (N[0]==1) return;
   for (i=0;i<N[0];i++) (*sd)+=(list[i]-(*av))*(list[i]-(*av));
   *sd=sqrt((*sd)/(N[0]-1));
}

double rav_sd(double clip) {
   int i,c=1;
   char*mask;
   double av,sd;

   if (N[0]<=1) return 0.;
   mask=(char*)calloc(N[0],1);
   if (!mask) merr();
   while (c) {
      c=0;
      av=sd=0;
      for (i=0;i<N[0];i++) if (!mask[i]) {
	 av+=list[i];
	 c++;
      }
      av/=c;
      for (i=0;i<N[0];i++) if (!mask[i]) sd+=(list[i]-av)*(list[i]-av);
      sd=sqrt(sd/(c-1))*clip;
      c=0;
      for (i=0;i<N[0];i++) if (!mask[i] && fabs(list[i]-av)>sd) {
	 mask[i]=1;
	 c=1;
      }
   }
   free(mask);
   return av;
}

void getsky(double dx,double*av,double*sd) {
   int i,c=1,n;
   char*mask;
   double smlt=3.;

   *av=*sd=0.;
   if (N[0]<=1) return;
   mask=(char*)calloc(N[0],1);
   if (!mask) merr();
   while (c) {
      n=0;
      *av=*sd=0;
      for (i=0;i<N[0];i++) if (!mask[i]) {
	 (*av)+=list[i];
	 n++;
      }
      (*av)/=n;
      for (i=0;i<N[0];i++) if (!mask[i]) (*sd)+=(list[i]-(*av))*(list[i]-(*av));
      (*sd)=sqrt(1.+(*sd)/(n-1))*smlt;
      c=0;
      for (i=0;i<N[0];i++) if (!mask[i] && fabs(list[i]-(*av))>(*sd)) {
	 mask[i]=1;
	 c=1;
      }
      if (!c && smlt>dx) {
	 smlt*=0.9;
	 if (smlt<dx) smlt=dx;
	 c=1;
      }
      (*sd)/=sqrt(n);
   }
   free(mask);
   return;
}

double median(double x0) {
   double x[3];
   int nm,np,i;

   if (N[0]<=0) return 0.;
   x[1]=x0;
   x[0]=x[2]=list[0];
   for (i=0;i<N[0];i++) {
      if (list[i]<x[0]) x[0]=list[i];
      else if (list[i]>x[2]) x[2]=list[i];
   }
   while (1) {
      nm=np=0;
      for (i=0;i<N[0];i++) {
	 if (list[i]<x[1]-0.000001) nm++;
	 else if (list[i]>x[1]+0.000001) np++;
      }
      if (nm<=N[0]/2 && np<=N[0]/2) return x[1];
      if (nm>np) x[2]=x[1];
      else x[0]=x[1];
      x[1]=0.5*(x[0]+x[2]);
   }
   return 0.;
}

double mode(void) {
   double x0,x1,t0,t1;
   int i,j,hist[10];

   x0=x1=list[0];
   if (N[0]<1) return 0.;
   if (N[0]<2) return list[0];
   for (i=0;i<N[0];i++) {
      if (list[i]<x0) x0=list[i];
      else if (list[i]>x1) x1=list[i];
   }
   x1+=1.;
   while (1) {
      for (i=0;i<10;i++) hist[i]=0;
      for (i=0;i<N[0];i++) if (list[i]>=x0 && list[i]<x1) {
	 j=(int)(10.*(list[i]-x0)/(x1-x0));
	 if (j>=0 && j<10) hist[j]++;
      }
      j=5;
      for (i=0;i<10;i++) if (hist[i]>hist[j]) j=i;
      t0=x0+0.1*j*(x1-x0);
      t1=x0+0.1*(j+1)*(x1-x0);
      x0=t0;
      x1=t1;
      if (x1-x0<1.e-5) return 0.5*(x0+x1);
   }
   return 0.;
}

int main(int argc,char**argv) {
   double x,y,z;
   if (argc!=2 && argc!=3 && argc!=4) {
      printf("Usage: %s <image1> <<section>> <<extension>>\n",*argv);
      return -1;
   }
   readfits(argv[1],&fits,1);
   if (argc==2) addall();
   else if (argc==3) addsec(argv[2],0);
   else addsec(argv[2],atoi(argv[3]));
   printf("N(good pix)=%d\n",N[0]);
   printf("N(sat pix)=%d\n",N[2]);
   printf("N(bad pix)=%d\n",N[1]);
   av_sd(&x,&y);
   printf(" Average=%f +/- %f\n",x,y/sqrt(N[0]));
   printf(" Std Dev=%f\n",y);
   /*
   printf("Clip(3)=%f\n",rav_sd(3));
   printf("Clip(2)=%f\n",rav_sd(2));
   */
   getsky(1.35,&y,&z);
   printf("Sky(1.5)=%f +/- %f\n",y,z);
   printf("  Median=%f\n",median(x));
   printf("    Mode=%f\n",mode());
   return 0;
}
