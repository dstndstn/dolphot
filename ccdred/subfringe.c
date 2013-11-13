#include "ccdproc.h"

ftype fits,cal,sfits,scal;
double scale,mult;
int box;

void add_mn(imtype*img,imtype*fimg,double m,int*n,double*mn,double mn0,double sd0) {
   int dx1,dx2,dy1,dy2,dz1,dz2,x,y,z;
   double v;
   float FMIN;

   if (!isimage(img)) return;
   parsecards(fimg,NULL,NULL,NULL,&FMIN,NULL,NULL,NULL,NULL,0,1);
   if (fimg!=&(scal.img)) parsecards(fimg,NULL,NULL,NULL,&FMIN,NULL,NULL,NULL,NULL,0,1);
   getsec(img,datasec,&dx1,&dx2,&dy1,&dy2,&dz1,&dz2);
   for (z=dz1;z<dz2;z++) for (y=dy1;y<dy2;y++) for (x=dx1;x<dx2;x++) if (fimg->data[z][y][x]>FMIN && (fabs((v=img->data[z][y][x]-m*fimg->data[z][y][x])-mn0)<=sd0 || sd0<0)) {
      (*mn)+=v;
      (*n)++;
   }
}

void add_sd(imtype*img,imtype*fimg,double m,int*n,double*sd,double mn0,double sd0,double mn) {
   int dx1,dx2,dy1,dy2,dz1,dz2,x,y,z;
   double v;
   float FMIN;

   if (!isimage(img)) return;
   parsecards(fimg,NULL,NULL,NULL,&FMIN,NULL,NULL,NULL,NULL,0,1);
   if (fimg!=&(scal.img)) parsecards(fimg,NULL,NULL,NULL,&FMIN,NULL,NULL,NULL,NULL,0,1);
   getsec(img,datasec,&dx1,&dx2,&dy1,&dy2,&dz1,&dz2);
   for (z=dz1;z<dz2;z++) for (y=dy1;y<dy2;y++) for (x=dx1;x<dx2;x++) if (fimg->data[z][y][x]>FMIN && (fabs((v=img->data[z][y][x]-m*fimg->data[z][y][x])-mn0)<=sd0 || sd0<0)) (*sd)+=(v-mn)*(v-mn);
}

double eval(double m) {
   int e,cont=3,n=1;
   double mn=0,sd=-1,mn0,sd0;

   printf("."); fflush(stdout);
   while (cont && n) {
      mn0=mn;
      sd0=sd;
      mn=sd=0;
      n=0;
      add_mn(&(sfits.img),&(scal.img),m,&n,&mn,mn0,sd0);
      for (e=0;e<sfits.Next;e++) add_mn(sfits.ext+e,scal.ext+e,m,&n,&mn,mn0,sd0);
      if (n>1) {
	 mn/=n;
	 add_sd(&(sfits.img),&(scal.img),m,&n,&sd,mn0,sd0,mn);
	 for (e=0;e<sfits.Next;e++) add_sd(sfits.ext+e,scal.ext+e,m,&n,&sd,mn0,sd0,mn);
	 sd=sqrt(1.+sd/(n-1))*3.5;
      }
      cont--;
   }
   //printf("%f %f\n",m,sd/3.5);
   return sd/3.5;
}

//adapted from NR;
#define SIGN(a,b) ((b) >= 0.0 ? fabs(a) : -fabs(a))
#define FMAX(a,b) (maxarg1=(a),maxarg2=(b),(maxarg1) > (maxarg2) ? (maxarg1) : (maxarg2))
#define GOLD 1.618034
#define GLIMIT 100.0
#define TINY 1.0e-20
#define SHFT(a,b,c,d) (a)=(b);(b)=(c);(c)=(d);
void mnbrak(double *ax,double *bx,double *cx,double *fa,double *fb,double *fc,double (*func)(double)) {
   double ulim,u,r,q,fu,dum,maxarg1,maxarg2;

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
      u=(*bx)-((*bx-*cx)*q-(*bx-*ax)*r)/(2.0*SIGN(FMAX(fabs(q-r),TINY),q-r));
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
}
#undef GOLD
#undef GLIMIT
#undef TINY
#undef SHFT

#define ITMAX 100
#define CGOLD 0.3819660
#define ZEPS 1.0e-10
#define SHFT(a,b,c,d) (a)=(b);(b)=(c);(c)=(d);
double brent(double ax,double bx,double cx,double (*f)(double),double tol,double *xmin) {
   int iter;
   double a,b,d=0,etemp,fu,fv,fw,fx,p,q,r,tol1,tol2,u,v,w,x,xm;
   double e=0.0;

   a=(ax<cx?ax:cx);
   b=(ax>cx?ax:cx);
   x=w=v=bx;
   fw=fv=fx=(*f)(x);
   for (iter=1;iter<=ITMAX;iter++) {
      xm=0.5*(a+b);
      tol2=2.0*(tol1=tol*fabs(x)+ZEPS);
      if (fabs(x-xm)<=(tol2-0.5*(b-a))) {
	 *xmin=x;
	 return fx;
      }
      if (fabs(e)>tol1) {
	 r=(x-w)*(fx-fv);
	 q=(x-v)*(fx-fw);
	 p=(x-v)*q-(x-w)*r;
	 q=2.0*(q-r);
	 if (q > 0.0) p=-p;
	 q=fabs(q);
	 etemp=e;
	 e=d;
	 if (fabs(p)>=fabs(0.5*q*etemp) || p<=q*(a-x) || p>=q*(b-x)) d=CGOLD*(e=(x>=xm?a-x:b-x));
	 else {
	    d=p/q;
	    u=x+d;
	    if (u-a<tol2 || b-u<tol2) d=SIGN(tol1,xm-x);
	 }
      }
      else d=CGOLD*(e=(x>=xm?a-x:b-x));
      u=(fabs(d)>=tol1 ? x+d : x+SIGN(tol1,d));
      fu=(*f)(u);
      if (fu<=fx) {
	 if (u>=x) a=x; else b=x;
	 SHFT(v,w,x,u)
	 SHFT(fv,fw,fx,fu)
      }
      else {
	 if (u<x) a=u; else b=u;
	 if (fu<=fw || w==x) {
	    v=w;
	    w=u;
	    fv=fw;
	    fw=fu;
	 } else if (fu<=fv || v==x || v==w) {
	    v=u;
	    fv=fu;
	 }
      }
   }
   printf("Too many iterations in brent\n");
   *xmin=x;
   return fx;
}
#undef ITMAX
#undef CGOLD
#undef ZEPS
#undef SHFT
#undef SIGN

void smooth(imtype*img,imtype*simg,imtype*fimg,float FMIN) {
   int x,y,z,N,xx,yy;
   double av,sd;
   static float *list=NULL;

   if (!list) {
      list=(float*)calloc((2*box+1)*(2*box+1),sizeof(float));
      if (!list) merr();
   }
   for (z=0;z<img->Z;z++) for (y=0;y<img->Y;y++) for (x=0;x<img->X;x++) if (fimg->data[z][y][x]>FMIN) {
      N=0;
      av=0.;
      for (yy=-box;yy<=box;yy++) for (xx=-box;xx<=box;xx++) if (x+xx>=0 && x+xx<img->X && y+yy>=0 && y+yy<img->Y && fimg->data[z][y+yy][x+xx]>FMIN) av+=(list[N++]=img->data[z][y+yy][x+xx]);
      if (N==0) simg->data[z][y][x]=img->data[z][y][x];
      else if (N<3) simg->data[z][y][x]=av/N;
      else {
	 av/=N;
	 xx=1;
	 while ((xx--)>0) {
	    sd=0.;
	    for (yy=0;yy<N;yy++) sd+=(list[yy]-av)*(list[yy]-av);
	    sd=sqrt(1.+sd/(N-1))*3.0;
	    av*=N;
	    for (yy=0;yy<N;yy++) if (fabs(list[yy]-av)>sd) {
	       av-=list[yy];
	       list[yy--]=list[--N];
	    }
	    av/=N;
	 }
	 simg->data[z][y][x]=av;
      }
   }
}

void fixff(imtype*img,imtype*simg,imtype*fimg,imtype*sfimg) {
   float DMIN,DMAX,FMIN,FMAX;
   int x,y,z;
   static int first=1;

   if (!isimage(img)) return;
   if (img->X!=fimg->X || img->Y!=fimg->Y || img->Z!=fimg->Z || !isimage(fimg)) {
      printf("**fits/fringe file mismatch\n");
      return;
   }
   parsecards(&(fits.img),NULL,NULL,NULL,&DMIN,&DMAX,NULL,NULL,NULL,0,1);
   if (img!=&(fits.img)) parsecards(img,NULL,NULL,NULL,&DMIN,&DMAX,NULL,NULL,NULL,0,0);
   parsecards(fimg,NULL,NULL,NULL,&FMIN,&FMAX,NULL,NULL,NULL,0,1);
   if (fimg!=&(cal.img)) parsecards(fimg,NULL,NULL,NULL,&FMIN,&FMAX,NULL,NULL,NULL,0,1);
   for (z=0;z<img->Z;z++) for (y=0;y<img->Y;y++) for (x=0;x<img->X;x++) if (img->data[z][y][x]<=DMIN || img->data[z][y][x]>=DMAX || fimg->data[z][y][x]<=FMIN || fimg->data[z][y][x]>=FMAX) sfimg->data[z][y][x]=fimg->data[z][y][x]=FMIN-1;
   if (box<1) return;
   if (first) {
      printf("Smoothing images\n");
      first=0;
   }
   smooth(img,simg,fimg,FMIN);
   smooth(fimg,sfimg,fimg,FMIN);
}

void fixchip(imtype*img,imtype*fimg) {
   int x,y,z;
   float FMIN;

   img->bscale=1.;
   img->bzero=0.;
   img->bits=-32;
   parsecards(fimg,NULL,NULL,NULL,&FMIN,NULL,NULL,NULL,NULL,0,1);
   if (fimg!=&(cal.img)) parsecards(fimg,NULL,NULL,NULL,&FMIN,NULL,NULL,NULL,NULL,0,1);
   for (z=0;z<img->Z;z++) for (y=0;y<img->Y;y++) for (x=0;x<img->X;x++) if (fimg->data[z][y][x]>FMIN) img->data[z][y][x]-=fimg->data[z][y][x]*mult;
}

int main(int argc,char**argv) {
   if (argc!=6) {
      printf("Usage: %s <file> <fringe> <scale> <boxcar> <output>\n",*argv);
      return 1;
   }
   paramfile("ccdproc.param",&ccdprocparam);
   scale=atof(argv[3]);
   box=atoi(argv[4]);
   readfits(argv[2],&cal,1);
   fitscopy(&cal,&scal);
   readfits(argv[1],&fits,1);
   fitscopy(&fits,&sfits);
   if (cal.Next!=fits.Next) printf("**fits/fringe file mismatch\n");
   else {
      int i;
      fixff(&(fits.img),&(sfits.img),&(cal.img),&(scal.img));
      for (i=0;i<fits.Next;i++) fixff(fits.ext+i,sfits.ext+i,cal.ext+i,scal.ext+i);
      if (scale<=0) {
	 double m[3],rv[3],rval,exp=-1;

	 if (scale<0) {
	    m[0]=-0.5*scale;
	    m[1]=-scale;
	    m[2]=-1.5*scale;
	 }
	 else {
	    parsecards(&(fits.img),NULL,NULL,&exp,NULL,NULL,NULL,NULL,NULL,0,0);
	    if (exp<1) exp=1.;
	    m[0]=exp/120.;
	    m[1]=exp/12.;
	    m[2]=exp/1.2;
	 }
	 printf("Fitting fringe pattern");
	 fflush(stdout);
	 mnbrak(m,m+1,m+2,rv,rv+1,rv+2,&eval);
	 rval=brent(m[0],m[1],m[2],&eval,0.01,&mult);
	 printf("\n  sky sigma = %e\n  scale = %e\n",rval,mult);
      }
      else mult=scale;
      fixchip(&(fits.img),&(cal.img));
      for (i=0;i<fits.Next;i++) fixchip(fits.ext+i,cal.ext+i);
   }
   writefits(argv[5],&fits,1);
   freefits(&fits);
   freefits(&cal);
   return 0;
}
