#include "ccdproc.h"

typedef double double3[3];
typedef double double4[4];

double *dmin,*dmax,*noise,*exptime,c2,nsig,PGAIN,nrej[4]={0,0,0,0};
float PDMIN,PDMAX;
int N,min,avg;
imtype *data;
ftype *fits,cfits;

void cleanchip(void) {
   int i,j,x,y,z;
   double g;
   static double3 *list;
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
      list=(double3*)calloc(sizeof(double3),N);
      if (!list) merr();
      first=0;
   }
   for (z=0;z<cfits.img.Z;z++) for (y=0;y<cfits.img.Y;y++) for (x=0;x<cfits.img.X;x++) {
      int n=0;
      double tmp,mid;

      for (i=0;i<N;i++) {
	  if (data[i].data[z][y][x]>=dmin[i] && data[i].data[z][y][x]<=dmax[i]) {
	    tmp=data[i].data[z][y][x]/exptime[i];
	    for (j=n;j>0 && list[j-1][0]>tmp;j--) {
	       list[j][0]=list[j-1][0];
	       list[j][1]=list[j-1][1];
	       list[j][2]=list[j-1][2];
	    }
	    list[j][0]=tmp;
	    list[j][1]=exptime[i];
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
	 if (s>0) cfits.img.data[z][y][x]=c/s;
	 else cfits.img.data[z][y][x]=mid;
      }
   }
}

int main(int argc,char**argv) {
   int i,j,Ntot,ext;
   FILE **fin,*fout=NULL;
   char str[161];
   cardtype *btype;

   if (argc<5) {
      printf("Usage: %s <fits files> <fudge factor> <sigma clip> <use min>\n",*argv);
      printf("Recommend ff=0., sig=3, use min=0\n");
      return 1;
   }
   Ntot=argc-4;
   dmin=(double*)calloc(Ntot,sizeof(double));
   dmax=(double*)calloc(Ntot,sizeof(double));
   noise=(double*)calloc(Ntot,sizeof(double));
   exptime=(double*)calloc(Ntot,sizeof(double));
   fits=(ftype*)calloc(Ntot,sizeof(ftype));
   data=(imtype*)calloc(Ntot,sizeof(imtype));
   fin=(FILE**)calloc(Ntot,sizeof(FILE*));
   btype=(cardtype*)calloc(Ntot,sizeof(cardtype));
   if (!dmin || !dmax || !noise || !exptime || !fits || !data || !fin || !btype) {
      printf("Memory allocation error\n");
      return 1;
   }
   paramfile("ccdproc.param",&ccdprocparam);
   c2=atof(argv[argc-3]);
   nsig=atof(argv[argc-2]);
   min=atoi(argv[argc-1]);
   for (i=0;i<Ntot;i++) {
      FILE *f;

      f=readfitsh(argv[i+1],fits+i,0);
      strcpy(btype[i],getcardval(&(fits[i].img),typekw,0));
      for (j=0;!btype[i][0] && j<fits[i].Next;j++) {
	 skipimage(f,&(fits[i].img));
	 if (fits[i].img.Nmax) free(fits[i].img.cards);
	 readexth(f,&(fits[i].img),0);
	 if (!btype[i][0]) strcpy(btype[i],getcardval(&(fits[i].img),typekw,0));
      }
      if (fits[i].img.Nmax) free(fits[i].img.cards);
      if (!i) cfits.Next=fits[i].Next;
      else if (fits[i].Next!=cfits.Next) {
	 printf("Number of extensions do not match\n");
	 return 1;
      }
   }
   cfits.img.Nmax=0;
   for (ext=0;ext<=cfits.Next;ext++) {
      double RN=0,EXP=0,EPOCH=0,AIR=0,GAIN=0,EXP0=0;
      float DMIN=0,DMAX=0;

      if (ext>0) printf("Extension %d\n",ext);
      for (i=0,N=0;i<Ntot;i++) if (strstr(btype[i],"dark")) {
	 double IEXP,IGAIN,IRN,IEPOCH,IAIR,IEXP0;
	 float IDMAX,IDMIN;

	 if (!ext) {
	    if (N==0) printf("Generating dark frame\n");
	    printf("Using FITS file %s\n",argv[i+1]);
	    fin[i]=readfitsh(argv[1+i],fits+i,0);
	 }
	 else readexth(fin[i],&(fits[i].img),0);
	 readimage(fin[i],&(fits[i].img));
	 memcpy(data+N,&(fits[i].img),sizeof(imtype));
	 parsecards(data+N,&IGAIN,&IRN,&IEXP,&IDMIN,&IDMAX,&IEPOCH,&IAIR,&IEXP0,0,1);
	 IEXP=atof(getcardval(&(fits[i].img),darktime,0));
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
      RN=sqrt(RN);
      EPOCH/=EXP;
      AIR/=EXP;
      //DMIN*=EXP;
      //DMAX*=EXP;
      EXP0/=EXP;
      PGAIN=GAIN;
      PDMIN=DMIN;
      PDMAX=DMAX;
      cfits.img.data=allocimg(cfits.img.X,cfits.img.Y,cfits.img.Z);
      cleanchip();
      insertcards(&(cfits.img),GAIN,RN,EXP,DMIN,DMAX,EPOCH,AIR,EXP0);
      if (!ext) {
	 sprintf(str,"%s.fits",darkbase);
	 fout=writefitsh(str,&cfits,0);
      }
      else writeexth(fout,&(cfits.img),0);
      writeimage(fout,&(cfits.img));
      for (i=0;i<Ntot;i++) if (strstr(btype[i],"dark")) {
	 if (fits[i].img.Nmax) free(fits[i].img.cards);
	 fits[i].img.Nmax=0;
	 if (!ext && fits[i].Next>0) free(fits[i].ext);
      }
      for (i=0;i<N;i++) freeimg(data[i].data,cfits.img.X,cfits.img.Y,cfits.img.Z);
      freeimg(cfits.img.data,cfits.img.X,cfits.img.Y,cfits.img.Z);
   }
   for (i=0;i<Ntot;i++) if (strstr(btype[i],"dark")) fclose(fin[i]);
   fclose(fout);
   printf("%7d (%7.4f%%) rejected as CRs\n",(int)(nrej[0]+0.5),nrej[0]/nrej[3]);
   printf("%7d (%7.4f%%) rejected as low\n",(int)(nrej[1]+0.5),nrej[1]/nrej[3]);
   printf("%7d (%7.4f%%) rejected as bad\n",(int)(nrej[2]+0.5),nrej[2]/nrej[3]);
   return 0;
}
