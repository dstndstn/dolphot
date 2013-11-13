#include <fits.h>
#include PGHEAD

ftype *disp_fits;
int disp_xwin;
float smin=0.,smax=10000.,gmed,gsl,CONTRAST=0.25;
int BLOCK=1,IBLOCK=1;

#define NMAX 10000
void setscale(int x0,int x1,int y0,int y1,int z) {
   int x,y,i,N=0,cont=1,Ns,skip;
   float val,in,ty,sd,data[NMAX];
   char mask[NMAX];
   static float DMIN,DMAX;
   static int first=1;

   skip=(int)(sqrt((x1+1-x0)*(y1+1-y0)/NMAX)+1);
   if (skip<1) skip=1;
   if (first) {
      parsecards(&(disp_fits->img),NULL,NULL,NULL,&DMIN,&DMAX,NULL,NULL,NULL,0,0);
      first=0;
   }
   for (y=y0;y<=y1 && N<NMAX;y++) for (x=x0;x<=x1 && N<NMAX;x++) if (y>=0 && x>=0 && y<disp_fits->img.Y && x<disp_fits->img.X && y%skip==0 && x%skip==0 && (val=disp_fits->img.data[z][y][x])>DMIN && val<DMAX) {
      for (i=N;i>0 && val<data[i-1];i--) data[i]=data[i-1];
      data[i]=val;
      mask[N]=1;
      N++;
   }
   if (N==0) return;
   smin=data[0];
   smax=data[N-1];
   if (smin>=smax) return;
   if (N%2) gmed=data[N/2];
   else gmed=0.5*(data[N/2-1]+data[N/2]);
   gsl=0.;
   while (cont) {
      float I=0.,X=0.,Y=0.,XX=0.,XY=0.,x;
      for (i=0;i<N;i++) if (mask[i]) {
	 I++;
	 X+=(x=i-0.5*(N-1));
	 Y+=data[i];
	 XX+=x*x;
	 XY+=x*data[i];
      }
      in=(X*XY-Y*XX)/(X*X-I*XX);
      gsl=(X*Y-I*XY)/(X*X-I*XX)*0.5*N;
      Ns=0;
      sd=0;
      for (i=0;i<N;i++) if (mask[i]) {
	 ty=in+gsl*2./N*(i-0.5*(N-1))-data[i];
	 sd+=ty*ty;
	 Ns++;
      }
      sd=sqrt(sd/Ns);
      cont=0;
      for (i=0;i<N;i++) if (mask[i]) {
	 ty=in+gsl*2./N*(i-0.5*(N-1))-data[i];
	 if (fabs(ty)>3*sd) mask[i]=0;
      }
   }
   return;
}

void display(ftype*fits,int ext,int z,double*cx,double*cy,char*ch,int init) {
   int xc,yc,X,Y,x0,y0,fsize;
   int xb,xt,yb,yt,xs,ys,y,CONT=1;
   float tr[6]={0.,1.,0.,0.,0.,1.},*vec;
   float x1,x2,y1,y2;
   static int first=1;
   chiptype chip;
   imtype *img;
   static float fx,fy;

   if (first) {
      disp_xwin=cpgopen("/xwin");
      cpgsubp(1,1);
      cpgsvp(0,1,0,1);
      cpgask(0);
      first=0;
   }
   else cpgslct(disp_xwin);
   disp_fits=fits;
   if (ext<0 || ext>disp_fits->Next) {
      printf("Image has %d extension%s; displaying base image\n",disp_fits->Next,(disp_fits->Next!=1?"s":""));
      ext=0;
   }
   if (ext==0) img=&(disp_fits->img);
   else img=disp_fits->ext+ext;
   if (!isimage(img)) {
      printf("Specified extension is not an image\n");
      return;
   }
   if (z<=0 || z>img->X) {
      printf("Illegal chip number; displaying chip 1\n");
      z=0;
   }
   else z--;
   fsize=sizeof(float);
   chip=img->data[z];
   if (init) {
      setscale(0,img->X-1,0,img->Y-1,z);
      xc=img->X/2;
      yc=img->Y/2;
      *cx=xc+0.5;
      *cy=yc+0.5;
      cpgqvp(3,&x1,&x2,&y1,&y2);
      while ((x2-x1)*BLOCK<img->X || (y2-y1)*BLOCK<img->Y) BLOCK*=2;
   }
   else {
      *cx=fx;
      *cy=fy;
      xc=(int)*cx;
      yc=(int)*cy;
   }
   while (CONT) {
      cpgpage();
      cpgqvp(3,&x1,&x2,&y1,&y2);
      X=((int)(x2-x1+0.5)*BLOCK)/IBLOCK;
      Y=((int)(y2-y1+0.5)*BLOCK)/IBLOCK;
      x0=xc-X/2;
      y0=yc-Y/2;
      cpgwnad(0.5+x0,0.5+x0+X,0.5+y0,0.5+y0+Y);
      cpgbbuf();
      cpgbox("bc", 0.0, 0, "bc", 0.0, 0);
      xb=x0; if (xb<0) xb=0;
      xt=x0+X; if (xt>img->X) xt=img->X;
      xs=xt-xb;
      yb=y0; if (yb<0) yb=0;
      yt=y0+Y; if (yt>img->Y) yt=img->Y;
      ys=yt-yb;
      smin=gmed-gsl/CONTRAST;
      smax=gmed+gsl/CONTRAST;
      tr[0]=xb;
      tr[3]=yb;
      vec=(float*)malloc(fsize*ys*xs+fsize);
      if (!vec) merr();
      for (y=0;y<ys;y++) memcpy(vec+y*xs,chip[y+yb]+xb,fsize*xs);
      cpggray(vec,xs,ys,1,xs,1,ys,smax,smin,tr);
      free(vec);
      cpgebuf();
      fx=*cx; fy=*cy;
      cpgband(0,1,0.,0.,&fx,&fy,ch);
      *cx=fx; *cy=fy;
      if (*ch=='<') CONTRAST*=0.5;
      else if (*ch=='>') CONTRAST*=2.;
      else if (*ch=='-') {
	 if (IBLOCK==1) BLOCK*=2;
	 else IBLOCK/=2;
      }
      else if (*ch=='+') {
	 if (BLOCK==1) IBLOCK*=2;
	 else BLOCK/=2;
      }
      else {
	 xc=(int)*cx; yc=(int)*cy;
	 if (*ch!='X' && *ch!='D') CONT=0;
      }
   }
   return;
}
