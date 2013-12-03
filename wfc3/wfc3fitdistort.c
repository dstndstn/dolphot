#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <dolphot.h>
#include "wfc3filters.h"
#include "wfc3distort.h"

#define NMAX 2200
int N=0,NP=4,cm,filt;
double x[NMAX],y[NMAX],xr[NMAX],yr[NMAX];

static inline double SQR(double x) {return x*x;}
#define SIGN(a,b) ((b) >= 0.0 ? fabs(a) : -fabs(a))
#define NR_END 1
#define FREE_ARG char*
static double maxarg1,maxarg2;
#define FMAX(a,b) (maxarg1=(a),maxarg2=(b),(maxarg1) > (maxarg2) ?\
        (maxarg1) : (maxarg2))

void nrerror(char error_text[])
/* Numerical Recipes standard error handler */
{
	fprintf(stderr,"Numerical Recipes run-time error...\n");
	fprintf(stderr,"%s\n",error_text);
	fprintf(stderr,"...now exiting to system...\n");
	exit(1);
}

double *vector(long nl, long nh)
/* allocate a double vector with subscript range v[nl..nh] */
{
	double *v;

	v=(double *)malloc((size_t) ((nh-nl+1+NR_END)*sizeof(double)));
	if (!v) nrerror("allocation failure in vector()");
	return v-nl+NR_END;
}

double **matrix(long nrl, long nrh, long ncl, long nch)
/* allocate a double matrix with subscript range m[nrl..nrh][ncl..nch] */
{
	long i, nrow=nrh-nrl+1,ncol=nch-ncl+1;
	double **m;

	/* allocate pointers to rows */
	m=(double **) malloc((size_t)((nrow+NR_END)*sizeof(double*)));
	if (!m) nrerror("allocation failure 1 in matrix()");
	m += NR_END;
	m -= nrl;

	/* allocate rows and set pointers to them */
	m[nrl]=(double *) malloc((size_t)((nrow*ncol+NR_END)*sizeof(double)));
	if (!m[nrl]) nrerror("allocation failure 2 in matrix()");
	m[nrl] += NR_END;
	m[nrl] -= ncl;

	for(i=nrl+1;i<=nrh;i++) m[i]=m[i-1]+ncol;

	/* return pointer to array of pointers to rows */
	return m;
}

void free_vector(double *v, long nl, long nh)
/* free a double vector allocated with vector() */
{
	free((FREE_ARG) (v+nl-NR_END));
}

void free_matrix(double **m, long nrl, long nrh, long ncl, long nch)
/* free a double matrix allocated by matrix() */
{
	free((FREE_ARG) (m[nrl]+ncl-NR_END));
	free((FREE_ARG) (m+nrl-NR_END));
}

void xfpoint(int i,double*r,double*xy) {
   double c=cos(r[4]);
   double s=sin(r[4]);
   xy[0]=r[1]+(xr[i]*c-yr[i]*s)*r[3];
   xy[1]=r[2]+(xr[i]*s+yr[i]*c)*r[3];
   WFC3revdistort(cm,filt,xy,xy+1);
   return;
}

double func4p(double*r) {
   int i;
   double rv=0.,xy[2];

   for (i=0;i<N;i++) {
      xfpoint(i,r,xy);
      rv+=SQR(xy[0]-x[i])+SQR(xy[1]-y[i]);
   }
   return rv/(N-2.);
}

double func2p(double*r) {
   int i;
   double rv=0.,xy[2];

   for (i=0;i<N;i++) {
      xy[0]=r[1]+xr[i];
      xy[1]=r[2]+yr[i];
      WFC3revdistort(cm,filt,xy,xy+1);
      rv+=SQR(xy[0]-x[i])+SQR(xy[1]-y[i]);
   }
   return rv/(N-1.);
}

double *func2data;
double f2scale[14]={1e3,1e3,1e6,1e6,1e6,1e9,1e9,1e9,1e9,1e12,1e12,1e12,1e12,1e12};
double func2(double *r) {
   double rv=0,x0,y0,ix,iy;
   double *ptr;
   int i;

   ptr=r-5;
   for (i=0;i<N;i++) {
      x0=(xr[i]*func2data[3]-func2data[5])/1000.0;
      y0=(yr[i]*func2data[3]-func2data[4])/1000.0;
      ix=func2data[0]+ptr[6]*y0+ptr[7]*x0+ptr[8]*y0*y0+ptr[9]*x0*y0+ptr[10]*x0*x0+ptr[11]*y0*y0*y0+ptr[12]*x0*y0*y0+ptr[13]*x0*x0*y0+ptr[14]*x0*x0*x0+ptr[15]*y0*y0*y0*y0+ptr[16]*x0*y0*y0*y0+ptr[17]*x0*x0*y0*y0+ptr[18]*x0*x0*y0*x0+ptr[19]*x0*x0*x0*x0;
      iy=func2data[1]+ptr[20]*y0+ptr[21]*x0+ptr[22]*y0*y0+ptr[23]*x0*y0+ptr[24]*x0*x0+ptr[25]*y0*y0*y0+ptr[26]*x0*y0*y0+ptr[27]*x0*x0*y0+ptr[28]*x0*x0*x0+ptr[29]*y0*y0*y0*y0+ptr[30]*x0*y0*y0*y0+ptr[31]*x0*x0*y0*y0+ptr[32]*x0*x0*y0*x0+ptr[33]*x0*x0*x0*x0;
      rv+=(x[i]-ix)*(x[i]-ix)+(y[i]-iy)*(y[i]-iy);
   }
   return rv/N*100;
}

int ncom;
double *pcom,*xicom,(*nrfunc)(double []);

double f1dim(double x)
{
	int j;
	double f,*xt;

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
void mnbrak(double *ax, double *bx, double *cx, double *fa, double *fb, double *fc,
	double (*func)(double))
{
	double ulim,u,r,q,fu,dum;

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
}
#undef GOLD
#undef GLIMIT
#undef TINY
#undef SHFT

#define ITMAX 100
#define CGOLD 0.3819660
#define ZEPS 1.0e-10
#define SHFT(a,b,c,d) (a)=(b);(b)=(c);(c)=(d);
double brent(double ax, double bx, double cx, double (*f)(double), double tol,
	double *xmin)
{
	int iter;
	double a,b,d=0.,etemp,fu,fv,fw,fx,p,q,r,tol1,tol2,u,v,w,x,xm;
	double e=0.0;

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
void linmin(double p[], double xi[], int n, double *fret, double (*func)(double []))
{
	int j;
	double xx,xmin,fx,fb,fa,bx,ax;

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
}
#undef TOL

#define ITMAX 200
void powell(double p[], double **xi, int n, double ftol, int *iter, double *fret,
	double (*func)(double []))
{
	int i,ibig,j;
	double del,fp,fptt,t,*pt,*ptt,*xit;

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
			t=2.0*(fp-2.0*(*fret)+fptt)*SQR(fp-(*fret)-del)-del*SQR(fp-fptt);
			if (t < 0.0) {
				linmin(p,xit,n,fret,func);
				for (j=1;j<=n;j++) {
					xi[j][ibig]=xi[j][n];
					xi[j][n]=xit[j];
				}
			}
		}
	}
}
#undef ITMAX

void usage(char*exe) {
   printf("Usage: %s <<CM> <filt>> <<Xmax> <Ymax>> <<flags>>\n",exe);
   printf("CM: 0=IR, 1/2=UVIS\n");
   printf("No arguments to run full test\n");
   printf("No X/Y to run test and fit inverse\n");
   printf("Flags:\n");
   printf("  -file=<filename>   to read X,Y positions from file\n");
   printf("  -XY                to compute dx,dy shift only\n");
   exit(-1);
}

int main(int argc,char**argv) {
   FILE *f;
   double fp[4]={0,0,1,0},xy[2],**xi,rv,X,Y;
   int i,j;

   /*
   // test pixsize
   WFC3initfilters();
   for (j=4096-128;j>0;j-=256) for (i=128;i<4096;i+=256) {
      if (j>2048) X=WFC3pixsize(1,0,i,j-2048);
      else X=WFC3pixsize(2,0,i,j);
      printf("%5d",(int)(X*1000));
   }
   printf("\n\n");
   for (j=1024-32;j>0;j-=64) for (i=32;i<1024;i+=64) {
      X=WFC3pixsize(0,50,i,j);
      printf("%5d",(int)(X*1000));
   }
   exit(-1);
   */
   if (argc==2 || argc==4) usage(*argv);
   WFC3initfilters();
   if (argc==1) {
      for (cm=0;cm<3;cm++) for (filt=0;filt<WFC3_NFILTERS;filt++) if (WFC3filters[filt].zp[cm]>=0) {
	 N=0;
	 rv=0.;
	 for (i=0;i<=4096 && (cm || i<=1024);i+=64) for (j=0;j<=2048 && (cm || j<=1024);j+=64) {
	    X=i; Y=j;
	    WFC3fwddistort(cm,filt,&X,&Y);
	    WFC3revdistort(cm,filt,&X,&Y);
	    rv+=(X-i)*(X-i)+(Y-j)*(Y-j);
	    N++;
	 }
	 printf("%-6s %d %f\n",WFC3filters[filt].name,cm,sqrt(rv/N));
      }
      return 0;
   }
   cm=atoi(argv[1]);
   filt=WFC3findfilt(argv[2]);
   if (argc==3) {
      double p[28];
      N=0;
      for (i=0;i<=4096 && (cm || i<=1024);i+=64) for (j=0;j<=2048 && (cm || j<=1024);j+=64) {
	 x[N]=X=i; y[N]=Y=j;
	 WFC3fwddistort(cm,filt,&X,&Y);
	 xr[N]=X; yr[N]=Y;
	 N++;
	 WFC3revdistort(cm,filt,&X,&Y);
	 //printf("%f %f\n",X,Y);
      }
      func2data=WFC3filters[filt].idc[0][cm];
      for (i=0;i<28;i++) p[i]=WFC3filters[filt].idc[1][cm][6+i]*f2scale[i%14];
      printf("rms0 = %f\n",sqrt(func2(p-1)/100));
      xi=matrix(1,28,1,28);
      for (i=1;i<=28;i++) {
	 for (j=1;j<=28;j++) xi[i][j]=0.;
	 xi[i][i]=1.;
      }
      powell(p-1,xi,28,1.e-14,&i,&rv,&func2);
      printf("rmsf = %f\n",sqrt(func2(p-1)/100));
      for (i=0;i<28;i++) printf(" %13g",p[i]/f2scale[i%14]);
      printf("\n");
      return 0;
   }
   X=atof(argv[3]);
   Y=atof(argv[4]);
   f=stdin;
   for (i=5;i<argc;i++) {
      if (!strncasecmp(argv[i],"-file=",6)) {
	 f=fopen(argv[i]+6,"r");
	 if (f==NULL) {
	    printf("Cannot read \"%s\"\n",argv[i]+6);
	    return -1;
	 }
      }
      else if (!strcasecmp(argv[i],"-XY")) NP=2;
      else usage(*argv);
   }
   if (f==stdin) printf("Enter X, Y, Xref, Yref\n");
   while (fscanf(f,"%lf %lf %lf %lf",x+N,y+N,xr+N,yr+N)==4) {
      xr[N]-=X*0.5;
      yr[N]-=Y*0.5;
      N++;
   }
   if (argc==6) fclose(f);
   printf("%d points read\n",N);
   if (NP==4) {
      xi=matrix(1,4,1,4);
      for (i=1;i<=4;i++) for (j=1;j<=4;j++) xi[i][j]=0.;
      for (i=1;i<=4;i++) xi[i][i]=1.;
      powell(fp-1,xi,4,1.e-5,&i,&rv,&func4p);
   }
   else {
      xi=matrix(1,2,1,2);
      xi[1][1]=xi[2][2]=1.;
      xi[1][2]=xi[2][1]=0.;
      powell(fp-1,xi,2,1.e-5,&i,&rv,&func2p);
      fp[2]=1.;
      fp[3]=0.;
   }
   printf("rms in fit: %f\n",sqrt(rv));
   for (i=0;i<N;i++) {
      xfpoint(i,fp-1,xy);
      printf("%4d: %8.2f %8.2f  %8.2f %8.2f  %8.2f %8.2f\n",i+1,xr[i],yr[i],x[i],y[i],xy[0],xy[1]);
   }
   printf("dx=%f\n",fp[0]);
   printf("dy=%f\n",fp[1]);
   printf("sc=%f\n",fp[2]);
   printf("th=%f\n",fp[3]*180./M_PI);
   return 0;
}
