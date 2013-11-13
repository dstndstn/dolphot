#include <fits.h>

int rsky,rin2,rout2,skip;
double sighi,siglo;
ftype fits,sky;

void getsky(imtype*in,imtype*out) {
   imgtype tsky;
   double *list,sky,sig,r;
   float DMIN,DMAX;
   int i,n,xx,yy,x,y,z,YDISP;

   YDISP=((in->Z)*(in->Y)+34)/35;
   tsky=allocimg(in->X,in->Y,in->Z);
   out->data=allocimg(in->X,in->Y,in->Z);
   list=(double*)calloc((2*rsky/skip+1)*(2*rsky/skip+1),sizeof(double));
   if (!list) merr();
   parsecards(&(fits.img),NULL,NULL,NULL,&DMIN,&DMAX,NULL,NULL,NULL,1,1);
   parsecards(in,NULL,NULL,NULL,&DMIN,&DMAX,NULL,NULL,NULL,1,0);
   for (z=0;z<in->Z;z++) for (y=0;y<in->Y;y++) for (x=0;x<in->X;x++) {
      if (x==0 && (y+z*(in->Y))%YDISP==0) {printf("."); fflush(stdout);}
      n=0;
      for (yy=y-rsky;yy<=y+rsky;yy+=skip) for (xx=x-rsky;xx<=x+rsky;xx+=skip) if (xx>=0 && xx<in->X && yy>=0 && yy<in->Y && in->data[z][yy][xx]>DMIN && in->data[z][yy][xx]<DMAX) {
	 r=(xx-x)*(xx-x)+(yy-y)*(yy-y);
	 if (r>=rin2 && r<=rout2) list[n++]=in->data[z][yy][xx];
      }
      xx=1;
      while (xx) {
	 xx=0;
	 if (!n) sky=0;
	 else {
	    sky=sig=0;
	    for (i=0;i<n;i++) sky+=list[i];
	    sky/=n;
	    for (i=0;i<n;i++) sig+=(list[i]-sky)*(list[i]-sky);
	    if (n>1) sig=sqrt(1.+sig/(n-1));
	    else sig=sqrt(1.+sig);
	    for (i=0;i<n;i++) if (list[i]<sky-siglo*sig || list[i]>sky+sighi*sig) {
	       list[i--]=list[--n];
	       xx=1;
	    }
	 }
      }
      tsky[z][y][x]=sky;
   }
   free(list);
   for (z=0;z<in->Z;z++) for (y=0;y<in->Y;y++) for (x=0;x<in->X;x++) {
      if (x==0 && (y+z*(in->Y))%YDISP==0) {printf("."); fflush(stdout);}
      n=0;
      for (yy=y-skip+1;yy<=y+skip;yy++) for (xx=x-skip+1;xx<=x+skip;xx++) if (xx>=0 && xx<in->X && yy>=0 && yy<in->Y) {
	 n++;
	 out->data[z][y][x]+=tsky[z][yy][xx];
      }
      if (n) out->data[z][y][x]/=n;
      else out->data[z][y][x]=tsky[z][y][x];
   }
   freeimg(tsky,in->X,in->Y,in->Z);
   /*
   out->bscale=1.;
   out->bzero=0.;
   out->bits=-32;
   */
   printf("\n"); fflush(stdout);
   return;
}

void getsky_q(imtype*in,imtype*out) {
   chiptype tmp;
   double *list,sky=0,sig;
   float DMIN,DMAX,mx,my,v0,v1;
   int i,n,xx,yy,x,y,z,NX,NY,QUICK;

   QUICK=-skip;
   NX=(in->X+QUICK-1)/QUICK;
   NY=(in->Y+QUICK-1)/QUICK;
   tmp=allocchip(NX+1,NY+1);
   out->data=allocimg(in->X,in->Y,in->Z);
   list=(double*)calloc(QUICK*QUICK,sizeof(double));
   if (!list) merr();
   parsecards(&(fits.img),NULL,NULL,NULL,&DMIN,&DMAX,NULL,NULL,NULL,1,1);
   parsecards(in,NULL,NULL,NULL,&DMIN,&DMAX,NULL,NULL,NULL,1,0);
   for (z=0;z<in->Z;z++) {
      for (yy=0;yy<NY;yy++) for (xx=0;xx<NX;xx++) {
	 int CONT=1;
	 n=0;
	 for (y=QUICK*yy;y<QUICK*(yy+1);y++) for (x=QUICK*xx;x<QUICK*(xx+1);x++) if (x>=0 && x<in->X && y>=0 && y<in->Y && in->data[z][y][x]>DMIN && in->data[z][y][x]<DMAX) list[n++]=in->data[z][y][x];
	 while (CONT) {
	    CONT=0;
	    if (!n) sky=0;
	    else {
	       sky=sig=0;
	       for (i=0;i<n;i++) sky+=list[i];
	       sky/=n;
	       for (i=0;i<n;i++) sig+=(list[i]-sky)*(list[i]-sky);
	       if (n>1) sig=sqrt(1.+sig/(n-1));
	       else sig=sqrt(1.+sig);
	       for (i=0;i<n;i++) if (list[i]<sky-siglo*sig || list[i]>sky+sighi*sig) {
		  list[i--]=list[--n];
		  CONT=1;
	       }
	    }
	 }
	 if (n) tmp[yy][xx]=sky;
	 else tmp[yy][xx]=-2.e30;
      }
      for (y=0;y<in->Y;y++) {
	 yy=(y+QUICK/2)/QUICK-1;
	 if (yy<0) {yy=0; my=0.;}
	 else if (yy>=NY-1) {yy=NY-1; my=0.;}
	 else my=(double)y/QUICK-0.5-yy;
	 for (x=0;x<in->X;x++) {
	    xx=(x+QUICK/2)/QUICK-1;
	    if (xx<0) {xx=0; mx=0.;}
	    else if (xx>=NX-1) {xx=NX-1; mx=0.;}
	    else mx=(double)x/QUICK-0.5-xx;
	    if (mx==0) {
	       v0=tmp[yy][xx];
	       v1=tmp[yy+1][xx];
	    }
	    else {
	       if (tmp[yy][xx+1]<-1.e30) v0=tmp[yy][xx];
	       else if (tmp[yy][xx]<-1.e30) v0=tmp[yy][xx+1];
	       else v0=(1.-mx)*tmp[yy][xx]+mx*tmp[yy][xx+1];
	       if (tmp[yy+1][xx+1]<-1.e30) v1=tmp[yy+1][xx];
	       else if (tmp[yy+1][xx]<-1.e30) v1=tmp[yy+1][xx+1];
	       else v1=(1.-mx)*tmp[yy+1][xx]+mx*tmp[yy+1][xx+1];
	    }
	    if (my>0 && v1>-1.e30) {
	       if (v0<-1.e30) v0=v1;
	       else v0=(1.-my)*v0+my*v1;
	    }
	    if (v0<-1.e30) out->data[z][y][x]=0.;
	    else out->data[z][y][x]=v0;
	 }
      }
   }
   free(list);
   freechip(tmp,NX+1,NY+1);
   /*
   out->bscale=1.;
   out->bzero=0.;
   out->bits=-32;
   */
   return;
}

int main(int argc,char**argv) {
   char str[81];
   int i;

   if (argc!=7) {
      printf("Usage: %s <fits base> <inner radius> <outer radius> <skip> <lower sigma> <upper sigma>\n",*argv);
      printf("Use negative <skip> to use NxN regions and interpolation\n");
      return 1;
   }
   rin2=atoi(argv[2]);
   if (rin2<0) rin2=0;
   rin2*=rin2;
   rsky=atoi(argv[3]);
   if (rsky<1) rsky=1;
   rout2=rsky*rsky;
   skip=atoi(argv[4]);
   if (skip==0) {
      printf("<skip> must be nonzero\n");
      return -1;
   }
   sighi=atof(argv[5]);
   if (sighi<1) sighi=1;
   siglo=atof(argv[6]);
   if (siglo<1) siglo=1;
   if (rsky%skip!=0) rsky=(rsky+skip-1)/skip*skip;
   sprintf(str,"%s.fits",argv[1]);
   readfits(str,&fits,1);
   memcpy(&sky,&fits,sizeof(ftype));
   if (isimage(&(fits.img))) {
      if (skip>0) getsky(&(fits.img),&(sky.img));
      else getsky_q(&(fits.img),&(sky.img));
   }
   if (sky.Next>0) {
      sky.ext=(imtype*)calloc(sizeof(imtype),sky.Next);
      if (!sky.ext) merr();
      for (i=0;i<sky.Next;i++) {
	 memcpy(sky.ext+i,fits.ext+i,sizeof(imtype));
	 if (isimage(fits.ext+i)) {
	    if (skip>0) getsky(fits.ext+i,sky.ext+i);
	    else getsky_q(fits.ext+i,sky.ext+i);
	 }
      }
   }
   sprintf(str,"%s.sky.fits",argv[1]);
   writefits(str,&sky,1);
   return 0;
}
