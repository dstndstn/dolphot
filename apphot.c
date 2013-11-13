#include <fits.h>
#include <unistd.h>

#ifdef PGPLOT
#include PGHEAD
extern void display(ftype*,int,int,double*,double*,char*,int);
#endif

ftype fits,resfits,skyfits;
imtype *data,*res,*sky;
int FitSky=1,Napsize=3,RCentroid=0,Interact=0;
double GAIN,EXP,RN,apsize[100]={20,50,100},apsky=20.0,FSat=0.999,Zero=25.,SkySig=2.25;
float DMAX,DMIN;
char base[81],infn[81]="",outfn[81]="";
FILE *fin,*fout;

void perr(char *str) {
   printf("%s\n",str);
   exit(0);
}

void skyerr(void) {
   printf("****Sky does not match the image\n");
   exit(-1);
}

int apphotparam(char*var,char*val) {
   int i;
   double x;
   char *ptr,*ptr2;

   if (!strcasecmp(var,"apsize")) {
      Napsize=0;
      ptr=val;
      while (1) {
	 ptr2=ptr;
	 x=strtod(ptr2,&ptr);
	 if (ptr==ptr2) perr("Unreadable aperture size");
	 if (x<=0) perr("Aperture size must be positive");
	 apsize[Napsize++]=x;
	 if (!*ptr) break;
      }
      if (Napsize<1) perr("No aperture sizes read");
      return 1;
   }
   i=strtol(val,&ptr,10);
   if (!*ptr) {
      if (!strcasecmp(var,"FitSky")) {FitSky=i; if (FitSky<0 || FitSky>2) perr("Fit sky must be 0-2"); return 1;}
#if PGPLOT
      if (!strcasecmp(var,"Interact")) {Interact=i; if (Interact<0 || Interact>2) perr("Interactive sky fit must be 0, 1, or 2"); return 1;}
#else
      if (!strcasecmp(var,"Interact")) {Interact=i; if (Interact<0 || Interact>1) perr("Interactive sky fit must be 0 or 1"); return 1;}
#endif
   }
   x=strtod(val,&ptr);
   if (!*ptr) {
      if (!strcasecmp(var,"apsky")) {apsky=x; if (apsky<=0.5) perr("Sky annulus must be at least 0.5"); return 1;}
      if (!strcasecmp(var,"SkySig")) {SkySig=x; if (SkySig<1.) perr("Sky sigma must be 1 or more"); return 1;}
      if (!strcasecmp(var,"RCentroid")) {RCentroid=x; if (RCentroid<0) perr("Centroid radius must be positive"); return 1;}
      if (!strcasecmp(var,"FSat")) {FSat=x; if (FSat<=0) perr("Saturation fraction must be positive"); if (FSat>1) FSat=1; return 1;}
      if (!strcasecmp(var,"Zero")) {Zero=x; return 1;}
   }
   if (!strcasecmp(var,"input")) {strcpy(infn,val); return 1;}
   if (!strcasecmp(var,"output")) {strcpy(outfn,val); return 1;}
   return 0;
}

/*
Photometry info
apsize-radius of aperture (float)
apsky-sky annulus (float)
FitSky-fit sky during photometry? (int 0=no,1=yes,2=include sky image)
FSat-fraction of saturate limit for final throw-out (float)
Zero-zeropoint for 1 DN/s (float)
*/

inline int pixsat(int x,int y,int z) {
   return (x>=0 && x<data->X && y>=0 && y<data->Y && data->data[z][y][x]>=DMAX);
}

inline int ppixOK(int x,int y,int z) {
   return (x>=0 && x<data->X && y>=0 && y<data->Y && data->data[z][y][x]>DMIN && data->data[z][y][x]<DMAX);
}

inline int ppixfOK(int x,int y,int z) {
   return (x>=0 && x<data->X && y>=0 && y<data->Y && data->data[z][y][x]>DMIN && data->data[z][y][x]<DMAX*FSat);
}

inline float skyval(int x,int y,int z) {
   if (x<0) x=0; else if (x>=data->X) x=data->X-1;
   if (y<0) y=0; else if (y>=data->Y) y=data->Y-1;
   return sky->data[z][y][x];
}

inline float noise(int x,int y,int z) {
   float s=skyval(x,y,z);
   if (x<0) x=0; else if (x>=data->X) x=data->X-1;
   if (y<0) y=0; else if (y>=data->Y) y=data->Y-1;
   return RN*RN+((s>0)?s:0)+((data->data[z][y][x]>s)?data->data[z][y][x]-s:0);
}

float calcsky(double fx,double fy,int z) {
   float *slist,r,sd,rsky,ssky;
   int x,y,ix,iy,Nsky=0,irsky;

   ix=(int)(fx+100.)-100;
   iy=(int)(fy+100.)-100;
   fx-=0.5;
   fy-=0.5;
   ssky=0.;
   rsky=apsize[Napsize-1]+apsky;
   irsky=(int)(rsky+1);
   if (FitSky) {
      r=M_PI*(rsky+2)*(rsky+2)-M_PI*apsize[Napsize-1]*apsize[Napsize-1];
      slist=(float*)calloc(sizeof(float),(int)r+1);
      if (!slist) merr();
      for (y=iy-irsky;y<=iy+irsky;y++) for (x=ix-irsky;x<=ix+irsky;x++) if (ppixOK(x,y,z)) {
	 r=sqrt((x-fx)*(x-fx)+(y-fy)*(y-fy));
	 if (r>apsize[Napsize-1]+0.5 && r<=rsky+0.5) {
	    slist[Nsky]=res->data[z][y][x];
	    if (FitSky==2) slist[Nsky]+=skyval(x,y,z);
	    Nsky++;
	 }
      }
      y=1;
      while (y && Nsky) {
	 ssky=sd=0;
	 for (x=0;x<Nsky;x++) ssky+=slist[x];
	 ssky/=Nsky;
	 for (x=0;x<Nsky;x++) sd+=(slist[x]-ssky)*(slist[x]-ssky);
	 sd=sqrt(1+sd/Nsky)*SkySig;
	 y=x=0;
	 while (x<Nsky) if (fabs(slist[x]-ssky)>sd) {
	    Nsky--;
	    slist[x]=slist[Nsky];
	    y=1;
	 }
	 else x++;
      }
      if (FitSky==2) ssky-=skyval(ix,iy,z);
      free(slist);
   }
   return ssky;
}

void run1phot(double rap,double fx,double fy,int z,double*s,double*ss,double ssky,int*sat) {
   float r;
   int x,y,ix,iy;

   if (sat) *sat=0;
   ix=(int)(fx+100.)-100;
   iy=(int)(fy+100.)-100;
   fx-=0.5;
   fy-=0.5;
   *s=*ss=0.;
   for (x=ix-(int)rap-1;x<=ix+(int)rap+1;x++) for (y=iy-(int)rap-1;y<=iy+(int)rap+1;y++) {
      r=rap+0.465-sqrt((x-fx)*(x-fx)+(y-fy)*(y-fy));
      if (r>1) r=1;
      if (r>0) {
	 if (sat && pixsat(x,y,z)) *sat=1;
	 if (ppixfOK(x,y,z)) {
	    (*s)+=(res->data[z][y][x]-ssky)*r;
	    (*ss)+=noise(x,y,z)*r;
	 }
	 else {
	    float ns=0,rs=0;
	    int ct=0,dx,dy;
	    for (dx=-1;dx<=1;dx++) for (dy=-1;dy<=1;dy++) if (ppixfOK(x+dx,y+dy,z)) {
	       ns+=noise(x+dx,y+dy,z);
	       rs+=res->data[z][y+dy][x+dx];
	       ct++;
	    }
	    if (ct) {
	       (*s)+=(rs/ct-ssky)*r;
	       (*ss)+=ns/ct*r;
	    }
	 }
      }
   }
   *ss=sqrt(*ss);
}

void centroid(int x0,int y0,int z,double *X,double *Y) {
   int xx,yy;
   double min=0.,xc=0.,yc=0.,wt=0.;

   for (yy=-RCentroid;yy<=RCentroid;yy++) for (xx=-RCentroid;xx<=RCentroid;xx++) if (abs(xx)+abs(yy)<RCentroid*3/2 && ppixfOK(x0+xx,y0+yy,z) && res->data[z][y0+yy][x0+xx]<min) min=res->data[z][y0+yy][x0+xx];
   for (yy=-RCentroid;yy<=RCentroid;yy++) for (xx=-RCentroid;xx<=RCentroid;xx++) if (abs(xx)+abs(yy)<RCentroid*3/2) {
      if (ppixOK(x0+xx,y0+yy,z)) {
	 if (res->data[z][y0+yy][x0+xx]>min) {
	    xc+=xx*(res->data[z][y0+yy][x0+xx]-min);
	    yc+=yy*(res->data[z][y0+yy][x0+xx]-min);
	    wt+=(res->data[z][y0+yy][x0+xx]-min);
	 }
      }
      else if (ppixOK(x0-xx,y0-yy,z)) {
	 xc+=xx*(res->data[z][y0-yy][x0-xx]-min);
	 yc+=yy*(res->data[z][y0-yy][x0-xx]-min);
	 wt+=(res->data[z][y0-yy][x0-xx]-min);
      }
   }
   if (wt>0) {
      *X=xc/wt+x0+0.5;
      *Y=yc/wt+y0+0.5;
   }
   else {
      *X=x0+0.5;
      *Y=y0+0.5;
   }
}

#ifdef PGPLOT
int pginit=0,xwin;

void plotit(double fx,double fy,int z,double ssky) {
   int i,f=1;
   static float *r,*m;
   double s,ss,mmin=0,mmax=0;
   char str[161];

   //return;
   if (!pginit) {
      xwin=cpgopen("/xwin");
      cpgask(0);
      cpgsch(1.2);
      pginit=1;
      r=(float*)calloc(2+(int)apsize[Napsize-1],sizeof(float));
      m=(float*)calloc(2+(int)apsize[Napsize-1],sizeof(float));
      if (!r || !m) merr();
   }
   else cpgslct(xwin);
   for (i=(int)(apsize[Napsize-1]+0.5);i>=1;i--) {
      r[i-1]=i;
      run1phot(i,fx,fy,z,&s,&ss,ssky,NULL);
      m[i-1]=Zero-2.5*log10(s/EXP/GAIN);
      if (f || m[i-1]<mmin) mmin=m[i-1];
      if (f || (m[i-1]>mmax && (i>0.15*apsize[Napsize-1] || mmin>mmax-0.5))) mmax=m[i-1];
      f=0;
   }
   cpgpage();
   cpgenv(0,apsize[Napsize-1],mmax+0.1,mmin-0.1,0,0);
   cpgpt((int)(apsize[Napsize-1]+0.5),r,m,16);
   cpgmove(0,mmin);
   cpgdraw(apsize[Napsize-1],mmin);
   cpgmove(0,mmin+0.05);
   cpgdraw(apsize[Napsize-1],mmin+0.05);
   sprintf(str,"%d %7.2f %7.2f",z,fx,fy);
   cpglab("r(pix)","magnitude",str);
}
#endif

int getnewpos(int*ext,int*z,double*x,double*y) {
   char str[321];
#ifdef PGPLOT
   if (Interact==2) {
      char ch;
      static int first=1;
      *ext=0;
      *z=1;
      printf("left-click on star position or press Q to quit\n");
      do {
	 display(&fits,*ext,*z,x,y,&ch,first);
	 first=0;
	 if (ch=='Q' || ch=='q') return 0;
      } while (ch!='A');
      return 1;
   }
#endif
   if (fin==stdin) printf("Enter extension, Z, X, Y\n");
   if (fscanf(fin,"%d %d %lf %lf",ext,z,x,y)==EOF) return 0;
   fgets(str,321,fin);
   return 1;
}

void photstar(FILE *f,int ext,double fx,double fy,int z) {
   double s,ss,x0,y0,ssky;
   int i,CONT=1,sat;
   char str[161];

   if (RCentroid>0) {
      do {
	 x0=fx; y0=fy;
	 centroid((int)fx,(int)fy,z,&fx,&fy);
	 fx=0.39*x0+0.61*fx;
	 fy=0.39*y0+0.61*fy;
      } while (fabs(fx-x0)>0.01 || fabs(fy-y0)>0.01);
   }
   fprintf(f,"%d %d %7.2f %7.2f\n",ext,z+1,fx,fy);
   ssky=calcsky(fx,fy,z);
   while (CONT) {
      for (i=0;i<Napsize;i++) {
	 run1phot(apsize[i],fx,fy,z,&s,&ss,ssky,&sat);
	 fprintf(f,"  %7.2f",apsize[i]);
	 if (s>0 && ss>0) fprintf(f," %8.4f %6.4f",Zero-2.5*log10(s/EXP/GAIN),1.0857362*ss/s);
	 else fprintf(f," %8.4f %6.4f",-99.9999,9.9999);
	 fprintf(f," %8.4f",ssky);
	 if (sat) fprintf(f," sat\n");
	 else fprintf(f,"\n");
      }
      if (!Interact) CONT=0;
      else {
#ifdef PGPLOT
	 plotit(fx,fy,z,ssky);
#endif
	 printf("Enter new sky value or q to quit\n");
	 fgets(str,161,stdin);
	 if (str[0]=='q' || str[0]=='Q') CONT=0;
	 else ssky=atof(str);
      }
   }
}

int main(int argc,char**argv) {
   char str[321];
   int ext,i,x,y,z;
   double fx,fy;

   if (argc<2) {
      printf("****Usage: %s <fits base> <<parameters>>\n",*argv);
      return 1;
   }
   paramfile("apphot.param",&apphotparam);
   for (i=2;i<argc;i++) parseparam(argv[i],&apphotparam);
   strcpy(base,argv[1]);
   if (Interact==2 || !infn[0]) fin=stdin;
   else if ((fin=fopen(infn,"r"))==NULL) {
      printf("****Error opening %s\n",infn);
      exit(-1);
   }
   if (!outfn[0]) fout=stdout;
   else if ((fout=fopen(outfn,"w"))==NULL) {
      printf("****Error opening %s\n",outfn);
      exit(-1);
   }
   sprintf(str,"%s.fits",base);
   readfits(str,&fits,1);
   fitscopy(&fits,&resfits);
   fitscopy(&fits,&skyfits);
   sprintf(str,"%s.sky.fits",base);
   if (!access(str,F_OK)) {
      readfits(str,&skyfits,1);
      if (skyfits.Next!=fits.Next) skyerr();
      if (isimage(&(fits.img))) {
	 if (skyfits.img.X!=fits.img.X || skyfits.img.Y!=fits.img.Y || skyfits.img.Z!=fits.img.Z) skyerr();
	 for (z=0;z<fits.img.Z;z++) for (y=0;y<fits.img.Y;y++) for (x=0;x<fits.img.X;x++) resfits.img.data[z][y][x]=fits.img.data[z][y][x]-skyfits.img.data[z][y][x];
      }
      else if (isimage(&(skyfits.img))) skyerr();
      for (ext=0;ext<fits.Next;ext++) {
	 if (isimage(fits.ext+ext)) {
	    if (skyfits.ext[ext].X!=fits.ext[ext].X || skyfits.ext[ext].Y!=fits.ext[ext].Y || skyfits.ext[ext].Z!=fits.ext[ext].Z) skyerr();
	    for (z=0;z<fits.ext[ext].Z;z++) for (y=0;y<fits.ext[ext].Y;y++) for (x=0;x<fits.ext[ext].X;x++) resfits.ext[ext].data[z][y][x]=fits.ext[ext].data[z][y][x]-skyfits.ext[ext].data[z][y][x];
	 }
	 else if (isimage(skyfits.ext+ext)) skyerr();
      }
   }
   else {
      if (isimage(&(fits.img))) for (z=0;z<fits.img.Z;z++) for (y=0;y<fits.img.Y;y++) for (x=0;x<fits.img.X;x++) skyfits.img.data[z][y][x]=0.;
      for (ext=0;ext<fits.Next;ext++) if (isimage(fits.ext+ext)) for (z=0;z<fits.ext[ext].Z;z++) for (y=0;y<fits.ext[ext].Y;y++) for (x=0;x<fits.ext[ext].X;x++) skyfits.ext[ext].data[z][y][x]=0.;
   }
   parsecards(&(fits.img),&GAIN,&RN,&EXP,&DMIN,&DMAX,NULL,NULL,NULL,0,1);
   if (isimage(&(fits.img))) for (z=0;z<fits.img.Z;z++) for (y=0;y<fits.img.Y;y++) for (x=0;x<fits.img.X;x++) {
      skyfits.img.data[z][y][x]*=GAIN;
      resfits.img.data[z][y][x]*=GAIN;
   }
   for (ext=0;ext<fits.Next;ext++) if (isimage(fits.ext+ext)) {
      parsecards(fits.ext+ext,&GAIN,&RN,&EXP,&DMIN,&DMAX,NULL,NULL,NULL,0,0);
      for (z=0;z<fits.ext[ext].Z;z++) for (y=0;y<fits.ext[ext].Y;y++) for (x=0;x<fits.ext[ext].X;x++) {
	 skyfits.ext[ext].data[z][y][x]*=GAIN;
	 resfits.ext[ext].data[z][y][x]*=GAIN;
      }
   }
   while (getnewpos(&ext,&z,&fx,&fy)) {
      parsecards(&(fits.img),&GAIN,&RN,&EXP,&DMIN,&DMAX,NULL,NULL,NULL,0,1);
      if (ext==0 && isimage(&(fits.img))) {
	 data=&(fits.img);
	 res=&(resfits.img);
	 sky=&(skyfits.img);
      }
      else if (ext>0 && isimage(fits.ext+ext-1)) {
	 data=fits.ext+ext-1;
	 res=resfits.ext+ext-1;
	 sky=skyfits.ext+ext-1;
	 parsecards(data,&GAIN,&RN,&EXP,&DMIN,&DMAX,NULL,NULL,NULL,0,0);
      }
      photstar(fout,ext,fx,fy,z-1);
   }
   if (fin!=stdin) fclose(fin);
   if (fout!=stdout) fclose(fout);
#ifdef PGPLOT
   if (pginit) cpgend();
#endif
   return 0;
}
