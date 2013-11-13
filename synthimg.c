#include <fits.h>
#include <time.h>

//parameters
int X=1024,Y=1024,Round=1,Noise=1,RPSF=15,SubPixel=1,FPSF=4;
double bg=1,GN=1,RN=1,EXP=100;
float apsf[3][6]={{3,0,0,0,0,0},{3,0,0,0,0,0},{0,0,0,0,0,0}};

//various data;
ftype fits;
float **psf;

extern int fitsparam(char*,char*);

void perr(char *str) {
   printf("%s\n",str);
   exit(0);
}

int synthimgparam(char*var,char*val) {
   double x;
   int i,j,k;
   char *ptr,*ptr2;

   if (!strcasecmp(var,"psfa") || !strcasecmp(var,"psfb") || !strcasecmp(var,"psfc")) {
      if (var[3]=='a' || var[3]=='A') j=0;
      else if (var[3]=='b' || var[3]=='B') j=1;
      else j=2;
      ptr2=val;
      for (k=0;k<6;k++) {ptr=ptr2; apsf[j][k]=strtod(ptr,&ptr2);}
      if (ptr==ptr2 || *ptr2) perr("img_psf requires six parameters");
      return 1;
   }
   if (!strcasecmp(var,"FPSF")) {
      if (!strcasecmp(val,"gauss")) FPSF=1;
      else if (!strcasecmp(val,"lorentz")) FPSF=2;
      else if (!strcasecmp(val,"lorentz2")) FPSF=3;
      else if (!strcasecmp(val,"g+l")) FPSF=4;
      else perr("FPSF=gauss,lorentz,lorentz2,g+l");
      return 1;
   }
   i=strtol(val,&ptr,10);
   if (!*ptr) {
      if (!strcasecmp(var,"X")) {X=i; if (X<=0) perr("X>0"); return 1;}
      if (!strcasecmp(var,"Y")) {Y=i; if (Y<=0) perr("Y>0"); return 1;}
      if (!strcasecmp(var,"RPSF")) {RPSF=i; if (RPSF<=0) perr("RPSF>0"); return 1;}
      if (!strcasecmp(var,"SubPixel")) {SubPixel=i; if (SubPixel<=0) perr("SubPixel>0"); return 1;}
      if (!strcasecmp(var,"Round")) {Round=i; if (Round<0 || Round>1) perr("Round=0,1"); return 1;}
      if (!strcasecmp(var,"Noise")) {Noise=i; if (Noise<0 || Noise>1) perr("Noise=0,1"); return 1;}
   }
   x=strtod(val,&ptr);
   if (!*ptr) {
      if (!strcasecmp(var,"GN")) {GN=x; if (GN<=0) perr("GN>0"); return 1;}
      if (!strcasecmp(var,"RN")) {RN=x; if (RN<0) perr("RN>=0"); return 1;}
      if (!strcasecmp(var,"bg")) {bg=x; if (bg<0) perr("bg>=0"); return 1;}
      if (!strcasecmp(var,"EXP")) {EXP=x; if (EXP<=0) perr("EXP>0"); return 1;}
   }
   return 0;
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

double gauss(void) {
   return (ran2()+ran2()-ran2()-ran2())/0.57735027;
}

int round2(double x) {
   int s=1,ix;
   if (x<0) {
      s=-1;
      x=-x;
   }
   if (x-((int)x)==0.5) {
      ix=(int)x;
      if (rand()<0.5) ix++;
   }
   else ix=(int)(x+0.5);
   return ix*s;
}

int poisson(float m) {
   int n,nmax;
   double lnm;

   if (m>100) return round2(m+gauss()*sqrt(m));
   nmax=(int)(10+m+5*sqrt(m));
   lnm=log(m);
   while (1) {
      n=(int)(ran2()*nmax);
      if (ran2()<exp(-m+n*lnm-lgamma(n+1))) return n;
   }
   return 0;
}

double *psfkernel;
void setpsfkernel(void) {
   int i;
   double m;

   psfkernel=(double*)calloc(sizeof(double),100001);
   if (!psfkernel) merr();
   m=1./(SubPixel*SubPixel*M_PI);
   if (FPSF==1) m*=0.693;
   else if (FPSF==3) m*=sqrt(2.)-1.;
   for (i=0;i<=100000;i++) {
      if (FPSF==1) psfkernel[i]=m*exp(-0.01*i*0.693);
      else if (FPSF==2) psfkernel[i]=m*1./(1+0.01*i);
      else if (FPSF==3) psfkernel[i]=m*1./(1+0.01*i*0.41421356)/(1+0.01*i*0.41421356);
      else if (FPSF==4) psfkernel[i]=m*(exp(-0.01*i*0.693)+1./(1+0.01*i));
   }
}

inline double evalpsf(double x) {
   int i;
   if (x<0) x=-x;
   x*=100.;
   i=(int)x;
   if (i>=100000) {
      x/=100000.;
      return psfkernel[100000]/x/x;
   }
   return (1+i-x)*psfkernel[i]+(x-i)*psfkernel[i+1];
}

void calc1psf(float x,float y,float sa,float sb,float c) {
   int i,j,ii,jj,sp1=0;
   static int first=1;
   double a,b,m,dx,dy,sps=1.,dy0;

   if (SubPixel!=1) {
      sps=1./SubPixel;
      sp1=SubPixel-1;
   }
   x-=(int)x+0.5;
   y-=(int)y+0.5;
   if (first) {
      setpsfkernel();
      psf=(float**)calloc(sizeof(float*),RPSF*2+1);
      if (!psf) merr();
      psf+=RPSF;
      for (j=-RPSF;j<=RPSF;j++) {
	 psf[j]=(float*)calloc(sizeof(float),RPSF*2+1);
	 if (!psf[j]) merr();
	 psf[j]+=RPSF;
      }
   }
   if (sa<=0 || sb<=0) {
      printf("Error in PSF generation: negative radius\n");
      exit(-1);
   }
   if (c>2 || c<-2) {
      printf("Error in PSF generation: illegal cross term\n");
      exit(-1);
   }
   c/=0.25*sa*sb;
   a=4./(sa*sa);
   b=4./(sb*sb);
   m=1.;
   m=sqrt(a*b-c*c*0.25);
   if (FPSF==2) m=m/(log(1+RPSF*RPSF*1.5*m));
   else if (FPSF==4) m=m/(1./0.693+log(1+RPSF*RPSF*1.5*m));
   if (SubPixel==1) for (j=-RPSF;j<=RPSF;j++) for (i=-RPSF;i<=RPSF;i++) {
      dx=x-i;
      dy=y-j;
      psf[j][i]=m*evalpsf(a*dx*dx+b*dy*dy+c*dx*dy);
   }
   else if (SubPixel==2) for (j=-RPSF;j<=RPSF;j++) for (i=-RPSF;i<=RPSF;i++) {
      double dxp,dxm,dyp,dym;
      dxm=x-i-0.5*(1.-sps); dxp=dxm+sps;
      dym=y-j-0.5*(1.-sps); dyp=dym+sps;
      psf[j][i]=m*(
	 evalpsf(a*dxm*dxm+b*dym*dym+c*dxm*dym)+evalpsf(a*dxp*dxp+b*dym*dym+c*dxp*dym)
	 +evalpsf(a*dxm*dxm+b*dyp*dyp+c*dxm*dyp)+evalpsf(a*dxp*dxp+b*dyp*dyp+c*dxp*dyp)
	 );
   }
   else if (SubPixel==3) for (j=-RPSF;j<=RPSF;j++) for (i=-RPSF;i<=RPSF;i++) {
      double dxp,dxm,dyp,dym;
      dx=x-i; dxp=dx+sps; dxm=dx-sps;
      dy=y-j; dyp=dy+sps; dym=dy-sps;
      psf[j][i]=m*(
	 evalpsf(a*dxm*dxm+b*dym*dym+c*dxm*dym)+evalpsf(a*dx*dx+b*dym*dym+c*dx*dym)+evalpsf(a*dxp*dxp+b*dym*dym+c*dxp*dym)
	 +evalpsf(a*dxm*dxm+b*dy*dy+c*dxm*dy)+evalpsf(a*dx*dx+b*dy*dy+c*dx*dy)+evalpsf(a*dxp*dxp+b*dy*dy+c*dxp*dy)
	 +evalpsf(a*dxm*dxm+b*dyp*dyp+c*dxm*dyp)+evalpsf(a*dx*dx+b*dyp*dyp+c*dx*dyp)+evalpsf(a*dxp*dxp+b*dyp*dyp+c*dxp*dyp)
	 );
   }
   else for (j=-RPSF;j<=RPSF;j++) for (i=-RPSF;i<=RPSF;i++) {
      psf[j][i]=0;
      dx=x-i-0.5*(1.-sps);
      dy=dy0=y-j-0.5*(1.-sps);
      psf[j][i]+=evalpsf(a*dx*dx+b*dy*dy+c*dx*dy);
      for (jj=0;jj<sp1;jj++) {
	 dy+=sps;
	 psf[j][i]+=evalpsf(a*dx*dx+b*dy*dy+c*dx*dy);
      }
      for (ii=0;ii<sp1;ii++) {
	 dx+=sps;
	 dy=dy0;
	 psf[j][i]+=evalpsf(a*dx*dx+b*dy*dy+c*dx*dy);
	 for (jj=0;jj<sp1;jj++) {
	    dy+=sps;
	    psf[j][i]+=evalpsf(a*dx*dx+b*dy*dy+c*dx*dy);
	 }
      }
      psf[j][i]*=m;
   }
   first=0;
   return;
}

void getpsfpars(float x,float y,float*a,float*b,float*c) {
   x-=X*0.5;
   y-=Y*0.5;
   *a=apsf[0][0]+x*apsf[0][1]+y*apsf[0][2]+x*x*apsf[0][3]+y*y*apsf[0][4]+x*y*apsf[0][5];
   *b=apsf[1][0]+x*apsf[1][1]+y*apsf[1][2]+x*x*apsf[1][3]+y*y*apsf[1][4]+x*y*apsf[1][5];
   *c=apsf[2][0]+x*apsf[2][1]+y*apsf[2][2]+x*x*apsf[2][3]+y*y*apsf[2][4]+x*y*apsf[2][5];
}

void addstar(double x,double y,double c) {
   int ix,iy,y1,x1;
   float pa,pb,pc;

   if (x<0) x=X*ran2();
   if (y<0) y=Y*ran2();
   printf("%7.2f %7.2f %g\n",x,y,c);
   ix=(int)x; iy=(int)y;
   getpsfpars(x,y,&pa,&pb,&pc);
   calc1psf(x,y,pa,pb,pc);
   for (y1=iy-RPSF;y1<=iy+RPSF;y1++) for (x1=ix-RPSF;x1<=ix+RPSF;x1++) if (x1>=0 && y1>=0 && x1<X && y1<Y) fits.img.data[0][y1][x1]+=c*psf[y1-iy][x1-ix];
   return;
}

int main(int argc,char**argv) {
   int i,x,y;
   FILE *f;
   float DMIN,DMAX;
   double sx,sy,c;

   if (argc<3) {
      printf("Usage: %s <stars> <output> <<pars>>\n",*argv);
      printf("Stars file: X, Y, electrons\n");
      return 1;
   }
   paramfile("fits.param",&fitsparam);
   paramfile("synthimg.param",&synthimgparam);
   for (i=3;i<argc;i++) {
      if (!strncmp(argv[i],"-p",2) || !strncmp(argv[i],"-P",2)) paramfile1(argv[i]+2,&synthimgparam);
      else parseparam(argv[i],&synthimgparam);
   }
   fits.Next=0;
   fits.img.Ncards=0;
   fits.img.Nmax=0;
   fits.img.X=X;
   fits.img.Y=Y;
   fits.img.Z=1;
   fits.img.bits=-32;
   fits.img.pcount=0;
   fits.img.bzero=0.;
   fits.img.bscale=1.;
   strcpy(fits.img.xtension,"IMAGE");
   fits.img.data=allocimg(fits.img.X,fits.img.Y,fits.img.Z);
   for (y=0;y<Y;y++) for (x=0;x<X;x++) fits.img.data[0][y][x]=bg;
   if ((f=fopen(argv[1],"r"))==NULL) {
      printf("Cannot read %s\n",argv[1]);
      exit(-1);
   }
   while (fscanf(f,"%lf %lf %lf\n",&sx,&sy,&c)==3) addstar(sx,sy,c);
   fclose(f);
   for (y=0;y<Y;y++) for (x=0;x<X;x++) {
      if (Noise) fits.img.data[0][y][x]=poisson(fits.img.data[0][y][x]);
      if (RN>0) fits.img.data[0][y][x]+=gauss()*RN;
      if (Round) fits.img.data[0][y][x]=round2(fits.img.data[0][y][x]/GN);
   }
   DMIN=DMAX=fits.img.data[0][0][0];
   for (y=0;y<Y;y++) for (x=0;x<X;x++) {
      if (fits.img.data[0][y][x]<DMIN) DMIN=fits.img.data[0][y][x];
      else if (fits.img.data[0][y][x]>DMAX) DMAX=fits.img.data[0][y][x];
   }
   DMIN-=1.e5;
   DMAX+=1.e5;
   insertcards(&(fits.img),GN,RN,EXP,DMIN,DMAX,52000,1,EXP);
   writefits(argv[2],&fits,1);
   return 0;
}
