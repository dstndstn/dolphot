#include "ccdproc.h"

typedef double double3[3];
typedef double double4[4];

double *dmin,*dmax,*noise,*exptime,c2,nsig,PGAIN,nrej[4]={0,0,0,0};
float PDMIN,PDMAX,minval;
int N,min,avg,scale;
imtype *data;
ftype *fits,cfits;

void cleanchip(void) {
   int i,j,x,y,z;
   double tmult=0,g;
   static double3 *list;
   static double *imult;
   static int first=1;

   g=atof(getcardval(data+0,fitsstr[0],0));
   if (!isimage(data)) return;
   for (i=0;i<N;i++) {
      if (data[i].X!=cfits.img.X || data[i].Y!=cfits.img.Y || data[i].Z!=cfits.img.Z || !isimage(data+i)) {
	 printf("**Size mismatch in images\n");
	 return;
      }
      if (atof(getcardval(data+i,fitsstr[0],0))!=g) {
	 printf("**Gain mismatch in images\n");
	 return;
      }
   }
   nrej[3]+=(double)N*data[0].X*data[0].Y*data[0].Z;
   if (first) {
      imult=(double*)calloc(sizeof(double),N);
      list=(double3*)calloc(sizeof(double3),N);
      if (!imult || !list) merr();
      first=0;
   }
   for (i=0;i<N;i++) {
      if (scale==0) imult[i]=PDMIN-1;
      else if (scale==1) imult[i]=exptime[i];
      else {
	 int dx1,dx2,dy1,dy2,dz1,dz2;

	 j=0;
	 imult[i]=0;
	 getsec(data+i,datasec,&dx1,&dx2,&dy1,&dy2,&dz1,&dz2);
	 for (z=dz1;z<dz2;z++) for (y=dy1;y<dy2;y++) for (x=dx1;x<dx2;x++) if (data[i].data[z][y][x]>=dmin[i] && data[i].data[z][y][x]<=dmax[i]) {
	    imult[i]+=data[i].data[z][y][x];
	    j++;
	 }
	 if (j) imult[i]/=j;
	 printf("Image %d scale factor: %f\n",i+1,imult[i]);
      }
      tmult+=imult[i];
   }
   for (z=0;z<cfits.img.Z;z++) for (y=0;y<cfits.img.Y;y++) for (x=0;x<cfits.img.X;x++) {
      int n=0;
      double tmp,mid;

      for (i=0;i<N;i++) {
	  if (data[i].data[z][y][x]>=dmin[i] && data[i].data[z][y][x]<=dmax[i]) {
	    tmp=data[i].data[z][y][x]/imult[i];
	    for (j=n;j>0 && list[j-1][0]>tmp;j--) {
	       list[j][0]=list[j-1][0];
	       list[j][1]=list[j-1][1];
	       list[j][2]=list[j-1][2];
	    }
	    list[j][0]=tmp;
	    list[j][1]=imult[i];
	    list[j][2]=noise[i];
	    n++;
	 }
	 else nrej[2]++;
      }
      if (n==0) cfits.img.data[z][y][x]=1.;
      else {
	 double c=0,s=0;
	 if (min) mid=list[0][0];
	 else mid=list[(n-1)/2][0];
	 for (i=0;i<n;i++) {
	    double dmid;

	    dmid=list[i][2]/PGAIN*nsig/list[i][1];
	    if (mid>0) dmid=sqrt(dmid*dmid+nsig*nsig*mid/PGAIN/list[i][1]+c2*c2*mid*mid);
	    if (fabs(list[i][0]-mid)<=dmid) {
	       c+=list[i][0]*list[i][1];
	       s+=list[i][1];
	    }
	    else if (list[i][0]>mid) nrej[0]++;
	    else nrej[1]++;
	 }
	 if (s>0) cfits.img.data[z][y][x]=c/s*tmult;
	 else cfits.img.data[z][y][x]=mid*tmult;
      }
   }
}

void add_avgain(imtype*img,double*av,long*n) {
   double GAIN;

   if (!isimage(img)) return;
   parsecards(img,&GAIN,NULL,NULL,NULL,NULL,NULL,NULL,NULL,0,1);
   (*av)+=GAIN;
   (*n)++;
}

void add_av(imtype*img,double*av,long*n,double av0,double sd0) {
   int dx1,dx2,dy1,dy2,dz1,dz2,x,y,z;
   float DMIN;

   if (!isimage(img)) return;
   parsecards(img,NULL,NULL,NULL,&DMIN,NULL,NULL,NULL,NULL,0,1);
   getsec(img,datasec,&dx1,&dx2,&dy1,&dy2,&dz1,&dz2);
   for (z=dz1;z<dz2;z++) for (y=dy1;y<dy2;y++) for (x=dx1;x<dx2;x++) if (img->data[z][y][x]>DMIN && (sd0<0 || fabs(img->data[z][y][x]-av0)<sd0)) {
      (*av)+=img->data[z][y][x];
      (*n)++;
   }
}

void add_sd(imtype*img,double*sd,double av0,double sd0,double av) {
   int dx1,dx2,dy1,dy2,dz1,dz2,x,y,z;
   float DMIN;

   if (!isimage(img)) return;
   parsecards(img,NULL,NULL,NULL,&DMIN,NULL,NULL,NULL,NULL,0,1);
   getsec(img,datasec,&dx1,&dx2,&dy1,&dy2,&dz1,&dz2);
   for (z=dz1;z<dz2;z++) for (y=dy1;y<dy2;y++) for (x=dx1;x<dx2;x++) if (img->data[z][y][x]>DMIN && (sd0<0 || fabs(img->data[z][y][x]-av0)<sd0)) (*sd)+=(img->data[z][y][x]-av)*(img->data[z][y][x]-av);
}

void chip_scale(imtype*img,double av,double avgain) {
   int x,y,z;
   double GAIN;
   float DMIN,DMAX;
   cardtype str;

   if (!isimage(img)) return;
   parsecards(img,&GAIN,NULL,NULL,&DMIN,&DMAX,NULL,NULL,NULL,0,1);
   for (z=0;z<img->Z;z++) for (y=0;y<img->Y;y++) for (x=0;x<img->X;x++) {
      if (img->data[z][y][x]<=DMIN) img->data[z][y][x]=0.;
      else img->data[z][y][x]/=av;
      if (img->data[z][y][x]<minval) img->data[z][y][x]=0.;
   }
   DMAX/=av;
   DMIN=minval;
   insertcards(img,GAIN*av,-1.e30,-1.e30,DMIN,DMAX,-1.e30,-1.e30,-1.e30);
   sprintf(str,"FLATGAIN= %20f /                                                ",avgain);
   addcard(img,str);
}

void imscale(void) {
   double av0=0.,av=1.,sd0=-1.,sd=0.,avgain=0.;
   long n=0;
   int ex,cont=1;

   add_avgain(&(cfits.img),&avgain,&n);
   for (ex=0;ex<cfits.Next;ex++) add_avgain(cfits.ext+ex,&avgain,&n);
   if (n) avgain/=n;
   while (n>0 && cont) {
      av=sd=0.;
      n=0;
      add_av(&(cfits.img),&av,&n,av0,sd0);
      for (ex=0;ex<cfits.Next;ex++) add_av(cfits.ext+ex,&av,&n,av0,sd0);
      if (n) av/=n;
      else av=1;
      add_sd(&(cfits.img),&sd,av0,sd0,av);
      for (ex=0;ex<cfits.Next;ex++) add_sd(cfits.ext+ex,&sd,av0,sd0,av);
      if (n) sd=sqrt(sd/n)*3.0;
      if (fabs(av-av0)<0.01 || fabs(av-av0)<0.001*fabs(av0)) cont=0;
      av0=av;
      sd0=sd;
   }
   printf("Mean gain = %f\n",avgain);
   printf("Scaling image by 1/%f\n",av);
   chip_scale(&(cfits.img),av,avgain);
   for (ex=0;ex<cfits.Next;ex++) chip_scale(cfits.ext+ex,av,avgain);
}

int main(int argc,char**argv) {
   int i,j,Ntot,Nfilt=0,nf,ext;
   FILE **fin,*fout=NULL;
   char str[161];
   cardtype*filt,*bfilt,*btype;

   if (argc<7) {
      printf("Usage: %s <fits files> <scale> <fudge factor> <sigma clip> <use min> <minimum flat value>\n",*argv);
      printf("Scale: 0=equal, 1=by exptime, 2=by mean value\n");
      printf("Recommend scale=2, ff=0., sig=3, use min=0\n");
      return 1;
   }
   Ntot=argc-6;
   dmin=(double*)calloc(Ntot,sizeof(double));
   dmax=(double*)calloc(Ntot,sizeof(double));
   noise=(double*)calloc(Ntot,sizeof(double));
   exptime=(double*)calloc(Ntot,sizeof(double));
   fits=(ftype*)calloc(Ntot,sizeof(ftype));
   data=(imtype*)calloc(Ntot,sizeof(imtype));
   fin=(FILE**)calloc(Ntot,sizeof(FILE*));
   filt=(cardtype*)calloc(Ntot,sizeof(cardtype));
   bfilt=(cardtype*)calloc(Ntot,sizeof(cardtype));
   btype=(cardtype*)calloc(Ntot,sizeof(cardtype));
   if (!dmin || !dmax || !noise || !exptime || !fits || !data || !fin || !filt || !bfilt || !btype) {
      printf("Memory allocation error\n");
      return 1;
   }
   paramfile("ccdproc.param",&ccdprocparam);
   scale=atof(argv[argc-5]);
   c2=atof(argv[argc-4]);
   nsig=atof(argv[argc-3]);
   min=atoi(argv[argc-2]);
   minval=atof(argv[argc-1]);
   if (minval<0) minval=0.;
   for (i=0;i<Ntot;i++) {
      int x;
      FILE *f;

      f=readfitsh(argv[i+1],fits+i,0);
      strcpy(bfilt[i],getcardval(&(fits[i].img),filtkw,0));
      strcpy(btype[i],getcardval(&(fits[i].img),typekw,0));
      for (j=0;(!bfilt[i][0] || !btype[i][0]) && j<fits[i].Next;j++) {
	 skipimage(f,&(fits[i].img));
	 if (fits[i].img.Nmax) free(fits[i].img.cards);
	 readexth(f,&(fits[i].img),0);
	 if (!bfilt[i][0]) strcpy(bfilt[i],getcardval(&(fits[i].img),filtkw,0));
	 if (!btype[i][0]) strcpy(btype[i],getcardval(&(fits[i].img),typekw,0));
      }
      if (fits[i].img.Nmax) free(fits[i].img.cards);
      if (strstr(btype[i],"flat")) {
	 for (x=0;x<Nfilt && strcmp(bfilt[i],filt[x]);x++);
	 if (x==Nfilt) strcpy(filt[Nfilt++],bfilt[i]);
      }
      if (!i) cfits.Next=fits[i].Next;
      else if (fits[i].Next!=cfits.Next) {
	 printf("Number of extensions do not match\n");
	 return 1;
      }
   }
   cfits.img.Nmax=0;
   for (nf=0;nf<Nfilt;nf++) {
      for (ext=0;ext<=cfits.Next;ext++) {
	 double RN=0,EXP=0,EPOCH=0,AIR=0,GAIN=0,EXP0=0;
	 float DMIN=0,DMAX=0;

	 if (ext>0) printf("Extension %d\n",ext);
	 for (i=0,N=0;i<Ntot;i++) if (strstr(btype[i],"flat")) {
	    char tfilt[81];

	    if (!ext) fin[i]=readfitsh(argv[1+i],fits+i,0);
	    else readexth(fin[i],&(fits[i].img),0);
	    strcpy(tfilt,getcardval(&(fits[i].img),filtkw,0));
	    if (!tfilt[0]) strcpy(tfilt,bfilt[i]);
	    if (!strcmp(filt[nf],tfilt)) {
	       double IEXP,IGAIN,IRN,IEPOCH,IAIR,IEXP0;
	       float IDMAX,IDMIN;

	       if (!ext) {
		  if (N==0) printf("Generating %s flat\n",filt[nf]);
		  printf("Using FITS file %s\n",argv[i+1]);
	       }
	       readimage(fin[i],&(fits[i].img));
	       memcpy(data+N,&(fits[i].img),sizeof(imtype));
	       parsecards(data+N,&IGAIN,&IRN,&IEXP,&IDMIN,&IDMAX,&IEPOCH,&IAIR,&IEXP0,0,1);
	       if (!N) {
		  if (!ext) strcpy(cfits.img.xtension,"IMAGE");
		  else strcpy(cfits.img.xtension,data[N].xtension);
		  cfits.img.Ncards=0;
		  cfits.img.X=data[N].X;
		  cfits.img.Y=data[N].Y;
		  cfits.img.Z=data[N].Z;
		  if (isimage(&(cfits.img))) {
		     cfits.img.bits=-32;
		     cfits.img.bzero=0.;
		     cfits.img.bscale=1.;
		  }
		  else {
		     cfits.img.bits=data[N].bits;
		     cfits.img.bzero=data[N].bzero;
		     cfits.img.bscale=data[N].bscale;
		  }
		  for (j=0;j<data[N].Ncards;j++) addcard(&cfits.img,data[N].cards[j]);
		  GAIN=IGAIN;
		  DMIN=IDMIN/IEXP;
		  DMAX=IDMAX/IEXP;
	       }
	       else {
		  if (IDMIN/IEXP<DMIN) DMIN=IDMIN/IEXP;
		  if (IDMAX/IEXP>DMAX) DMAX=IDMAX/IEXP;
	       }
	       dmin[N]=IDMIN;
	       dmax[N]=IDMAX;
	       noise[N]=IRN;
	       exptime[N]=IEXP;
	       RN+=IRN*IRN;
	       EPOCH+=IEXP*IEPOCH;
	       AIR+=IEXP*IAIR;
	       EXP+=IEXP;
	       EXP0+=IEXP*IEXP0;
	       N++;
	    }
	    else skipimage(fin[i],&(fits[i].img));
	 }
	 RN=sqrt(RN);
	 EPOCH/=EXP;
	 AIR/=EXP;
	 DMIN*=EXP;
	 DMAX*=EXP;
	 EXP0/=EXP;
	 PGAIN=GAIN;
	 PDMIN=DMIN;
	 PDMAX=DMAX;
	 cfits.img.data=allocimg(cfits.img.X,cfits.img.Y,cfits.img.Z);
	 cleanchip();
	 insertcards(&(cfits.img),GAIN,RN,EXP,DMIN,DMAX,EPOCH,AIR,EXP0);
	 if (!ext) {
	    sprintf(str,"%s.%s.fits",flatbase,filt[nf]);
	    fout=writefitsh(str,&cfits,0);
	 }
	 else writeexth(fout,&(cfits.img),0);
	 writeimage(fout,&(cfits.img));
	 for (i=0;i<Ntot;i++) if (strstr(btype[i],"flat")) {
	    if (fits[i].img.Nmax) free(fits[i].img.cards);
	    fits[i].img.Nmax=0;
	    if (!ext && fits[i].Next>0) free(fits[i].ext);
	 }
	 for (i=0;i<N;i++) freeimg(data[i].data,cfits.img.X,cfits.img.Y,cfits.img.Z);
	 freeimg(cfits.img.data,cfits.img.X,cfits.img.Y,cfits.img.Z);
      }
      for (i=0;i<Ntot;i++) if (strstr(btype[i],"flat")) fclose(fin[i]);
      fclose(fout);
      sprintf(str,"%s.%s.fits",flatbase,filt[nf]);
      readfits(str,&cfits,0);
      imscale();
      writefits(str,&cfits,1);
      freefits(&cfits);
   }
   printf("%7d (%7.4f%%) rejected as CRs\n",(int)(nrej[0]+0.5),nrej[0]/nrej[3]);
   printf("%7d (%7.4f%%) rejected as low\n",(int)(nrej[1]+0.5),nrej[1]/nrej[3]);
   printf("%7d (%7.4f%%) rejected as bad\n",(int)(nrej[2]+0.5),nrej[2]/nrej[3]);
   return 0;
}
