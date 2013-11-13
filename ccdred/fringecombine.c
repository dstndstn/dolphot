#include "ccdproc.h"

typedef double double4[4];

double *dmin,*dmax,*noise,*exptime,*iscale,*isky,c2,siglo,sighi,PGAIN,nrej[4]={0,0,0,0};
float PDMIN,PDMAX;
int N,min,avg,scale,SIZE;
imtype *data;
ftype *fits,cfits,*mask=NULL;
char maskfn[81]="";

void cleanchip(void) {
   int i,j,x,y,z;
   double tmult=0,g;
   static double4 *list;
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
      list=(double4*)calloc(sizeof(double4),N);
      if (!imult || !list) merr();
      first=0;
   }
   for (i=0;i<N;i++) {
      if (scale==0) imult[i]=iscale[i];
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
	    tmp=(data[i].data[z][y][x]-isky[i])/imult[i];
	    for (j=n;j>0 && list[j-1][0]>tmp;j--) {
	       list[j][0]=list[j-1][0];
	       list[j][1]=list[j-1][1];
	       list[j][2]=list[j-1][2];
	       list[j][3]=list[j-1][3];
	    }
	    list[j][0]=tmp;
	    list[j][1]=imult[i];
	    list[j][2]=isky[i];
	    list[j][3]=noise[i];
	    n++;
	  }
	 else nrej[2]++;
      }
      if (n==0) cfits.img.data[z][y][x]=PDMIN-1;
      else {
	 double c=0,s=0;

	 if (min==3) mid=list[n/3][0];
	 else if (min==4) mid=list[n/4][0];
	 else if (min) mid=list[0][0];
	 else mid=list[(n-1)/2][0];
	 for (i=0;i<n;i++) {
	    double tmid,tmidlo,tmidhi;

	    tmid=mid+list[i][2]/list[i][1];
	    tmidlo=list[i][3]/PGAIN*siglo/list[i][1];
	    if (tmid>0) tmidlo=sqrt(tmidlo*tmidlo+siglo*siglo*tmid/PGAIN/list[i][1]+c2*c2*tmid*tmid);
	    tmidhi=list[i][3]/PGAIN*sighi/list[i][1];
	    if (tmid>0) tmidhi=sqrt(tmidhi*tmidhi+sighi*sighi*tmid/PGAIN/list[i][1]+c2*c2*tmid*tmid);
	    tmid-=list[i][2]/list[i][1];
	    if (list[i][0]>=tmid-tmidlo && list[i][0]<=tmid+tmidhi) {
	       c+=list[i][0]*list[i][1];
	       s+=list[i][1];
	    }
	    else if (list[i][0]>tmid) nrej[0]++;
	    else nrej[1]++;
	 }
	 if (s>0) cfits.img.data[z][y][x]=c/s*tmult;
	 else cfits.img.data[z][y][x]=mid*tmult;
      }
   }
}

int calcavav=1;
double avav=0.;

void flattenchip(imtype*img,int ext) {
   int dx1,dx2,dy1,dy2,dz1,dz2,x,y,z,xx,yy,n,nav=0;
   float DMIN,DMAX;
   imgtype av;

   if (SIZE<1 || !isimage(img)) return;
   parsecards(img,NULL,NULL,NULL,&DMIN,&DMAX,NULL,NULL,NULL,0,1);
   getsec(img,datasec,&dx1,&dx2,&dy1,&dy2,&dz1,&dz2);
   av=allocimg(img->X,img->Y,1);
   for (z=dz1;z<dz2;z++) {
      if (calcavav) {
	 avav=0;
	 nav=0;
      }
      for (y=dy1;y<dy2;y++) for (x=dx1;x<dx2;x++) if (img->data[z][y][x]>DMIN && img->data[z][y][x]<DMAX && (mask==NULL || mask[ext].img.data[z][y][x]<0.5)) {
	 av[0][y][x]=0;
	 n=0;
	 for (yy=-SIZE*7;yy<=SIZE*7;yy+=SIZE) for (xx=-SIZE*7;xx<=SIZE*7;xx+=SIZE) if (x+xx>=dx1 && x+xx<dx2 && y+yy>=dy1 && y+yy<dy2 && img->data[z][y+yy][x+xx]>DMIN && img->data[z][y+yy][x+xx]<DMAX && (mask==NULL || mask[ext].img.data[z][y+yy][x+xx]<0.5)) {
	    av[0][y][x]+=img->data[z][y+yy][x+xx];
	    n++;
	 }
	 if (n) {
	    av[0][y][x]/=n;
	    if (calcavav) {
	       avav+=av[0][y][x];
	       nav++;
	    }
	 }
      }
      if (nav && calcavav) avav/=nav;
      for (y=dy1;y<dy2;y++) for (x=dx1;x<dx2;x++) if (img->data[z][y][x]>DMIN && img->data[z][y][x]<DMAX && (mask==NULL || mask[ext].img.data[z][y][x]<0.5)) img->data[z][y][x]-=av[0][y][x]-avav;
      calcavav=0;
   }
   freeimg(av,img->X,img->Y,1);
}

void flatten(void) {
   int e;

   calcavav=1;
   flattenchip(&(cfits.img),0);
   for (e=0;e<cfits.Next;e++) flattenchip(cfits.ext+e,e+1);
}

void add_av(imtype*img,int ext,double*av,long*n,double av0,double sd0) {
   int dx1,dx2,dy1,dy2,dz1,dz2,x,y,z;
   float DMIN,DMAX;

   if (!isimage(img)) return;
   parsecards(img,NULL,NULL,NULL,&DMIN,&DMAX,NULL,NULL,NULL,0,1);
   getsec(img,datasec,&dx1,&dx2,&dy1,&dy2,&dz1,&dz2);
   for (z=dz1;z<dz2;z++) for (y=dy1;y<dy2;y++) for (x=dx1;x<dx2;x++) if (img->data[z][y][x]>DMIN && img->data[z][y][x]<DMAX && (mask==NULL || mask[ext].img.data[z][y][x]<0.5) && (sd0<0 || fabs(img->data[z][y][x]-av0)<sd0)) {
      (*av)+=img->data[z][y][x];
      (*n)++;
   }
}

void add_sd(imtype*img,int ext,double*sd,double av0,double sd0,double av) {
   int dx1,dx2,dy1,dy2,dz1,dz2,x,y,z;
   float DMIN,DMAX;

   if (!isimage(img)) return;
   parsecards(img,NULL,NULL,NULL,&DMIN,&DMAX,NULL,NULL,NULL,0,1);
   getsec(img,datasec,&dx1,&dx2,&dy1,&dy2,&dz1,&dz2);
   for (z=dz1;z<dz2;z++) for (y=dy1;y<dy2;y++) for (x=dx1;x<dx2;x++) if (img->data[z][y][x]>DMIN && img->data[z][y][x]<DMAX && (mask==NULL || mask[ext].img.data[z][y][x]<0.5) && (sd0<0 || fabs(img->data[z][y][x]-av0)<sd0)) (*sd)+=(img->data[z][y][x]-av)*(img->data[z][y][x]-av);
}

void chip_scale(imtype*img,int ext,double av,double sd) {
   int x,y,z;
   double GAIN;
   float DMIN,DMAX;

   if (!isimage(img)) return;
   parsecards(img,&GAIN,NULL,NULL,&DMIN,&DMAX,NULL,NULL,NULL,0,1);
   for (z=0;z<img->Z;z++) for (y=0;y<img->Y;y++) for (x=0;x<img->X;x++) {
      if (img->data[z][y][x]>DMIN && img->data[z][y][x]<DMAX) img->data[z][y][x]=(img->data[z][y][x]-av)/sd;
      else img->data[z][y][x]=0.;
   }
   insertcards(img,GAIN*sd,-1.e30,-1.e30,(DMIN-av)/sd,(DMAX-av)/sd,-1.e30,-1.e30,-1.e30);
}

void imscale(void) {
   double av0=0.,av=1.,sd0=-1.,sd=0.;
   long n=1;
   int ex,cont=1;

   while (n>0 && cont) {
      av=sd=0.;
      n=0;
      add_av(&(cfits.img),0,&av,&n,av0,sd0);
      for (ex=0;ex<cfits.Next;ex++) add_av(cfits.ext+ex,ex+1,&av,&n,av0,sd0);
      if (n) av/=n;
      else av=1;
      add_sd(&(cfits.img),0,&sd,av0,sd0,av);
      for (ex=0;ex<cfits.Next;ex++) add_sd(cfits.ext+ex,ex+1,&sd,av0,sd0,av);
      if (n) sd=sqrt(sd/n)*3.5;
      if (fabs(av-av0)<0.01 || fabs(av-av0)<0.001*fabs(av0)) cont=0;
      av0=av;
      sd0=sd;
   }
   sd/=3.5;
   printf("Scaling image by av=%f, sd=%f\n",av,sd);
   chip_scale(&(cfits.img),0,av,sd);
   for (ex=0;ex<cfits.Next;ex++) chip_scale(cfits.ext+ex,ex+1,av,sd);
}

int main(int argc,char**argv) {
   int i,j,Ntot,ext;
   FILE **fin,*fout=NULL;
   char fn[161];

   if (argc<11 || (argc-9)%3) {
      printf("Usage: %s <fits files, scale, sky> <scale> <fudge factor> <lo sigma> <hi sigma> <pattern size> <use min> <imager mask> <output file>\n",*argv);
      printf("Scale: 0=use values, 1=by exptime, 2=by mean value\n");
      printf("Recommend scale=0, ff=0., sig=3, use min=0\n");
      return 1;
   }
   Ntot=(argc-9)/3;
   iscale=(double*)calloc(Ntot,sizeof(double));
   isky=(double*)calloc(Ntot,sizeof(double));
   dmin=(double*)calloc(Ntot,sizeof(double));
   dmax=(double*)calloc(Ntot,sizeof(double));
   noise=(double*)calloc(Ntot,sizeof(double));
   exptime=(double*)calloc(Ntot,sizeof(double));
   fits=(ftype*)calloc(Ntot,sizeof(ftype));
   data=(imtype*)calloc(Ntot,sizeof(imtype));
   fin=(FILE**)calloc(Ntot,sizeof(FILE*));
   if (!iscale || !isky || !dmin || !dmax || !noise || !exptime || !fits || !data || !fin) {
      printf("Memory allocation error\n");
      return 1;
   }
   paramfile("ccdproc.param",&ccdprocparam);
   scale=atof(argv[argc-8]);
   c2=atof(argv[argc-7]);
   siglo=atof(argv[argc-6]);
   sighi=atof(argv[argc-5]);
   SIZE=atoi(argv[argc-4]);
   min=atoi(argv[argc-3]);
   strcpy(maskfn,argv[argc-2]);
   for (i=0;i<Ntot;i++) {
      FILE *f;

      iscale[i]=atof(argv[i*3+2]);
      isky[i]=atof(argv[i*3+3]);
      f=readfitsh(argv[i*3+1],fits+i,0);
      fclose(f);
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
      double IEXP,IGAIN,IRN,IEPOCH,IAIR,IEXP0;
      float IDMAX,IDMIN,DMIN=0,DMAX=0;

      if (ext>0) printf("Extension %d\n",ext);
      for (i=0,N=0;i<Ntot;i++) {
	 if (!ext) fin[i]=readfitsh(argv[i*3+1],fits+i,0);
	 else readexth(fin[i],&(fits[i].img),0);
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
	 EXP0+=IEXP0*IEXP;
	 N++;
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
      if (!ext) fout=writefitsh(argv[argc-1],&cfits,0);
      else writeexth(fout,&(cfits.img),0);
      writeimage(fout,&(cfits.img));
      for (i=0;i<Ntot;i++) {
	 if (fits[i].img.Nmax) free(fits[i].img.cards);
	 fits[i].img.Nmax=0;
	 if (!ext && fits[i].Next>0) free(fits[i].ext);
      }
      for (i=0;i<N;i++) freeimg(data[i].data,cfits.img.X,cfits.img.Y,cfits.img.Z);
      freeimg(cfits.img.data,cfits.img.X,cfits.img.Y,cfits.img.Z);
   }
   for (i=0;i<Ntot;i++) fclose(fin[i]);
   fclose(fout);
   readfits(argv[argc-1],&cfits,0);
   if (maskfn[0]) {
      mask=(ftype*)calloc(sizeof(ftype),cfits.Next+1);
      if (!mask) merr();
      if (isimage(&(cfits.img))) {
	 sprintf(fn,"%s/masks/%s_%d.fits",BASEDIR,maskfn,0);
	 readfits(fn,mask,0);
      }
      for (ext=0;ext<cfits.Next;ext++) if (isimage(cfits.ext+ext)) {
	 sprintf(fn,"%s/masks/%s_%d.fits",BASEDIR,maskfn,ext+1);
	 readfits(fn,mask+ext+1,0);
      }
   }
   flatten();
   imscale();
   writefits(argv[argc-1],&cfits,1);
   if (maskfn[0]) {
      if (isimage(&(cfits.img))) freefits(mask);
      for (ext=0;ext<cfits.Next;ext++) if (isimage(cfits.ext+ext)) freefits(mask+ext+1);
      free(mask);
   }
   freefits(&cfits);
   printf("%7d (%7.4f%%) rejected as CRs\n",(int)(nrej[0]+0.5),nrej[0]/nrej[3]);
   printf("%7d (%7.4f%%) rejected as low\n",(int)(nrej[1]+0.5),nrej[1]/nrej[3]);
   printf("%7d (%7.4f%%) rejected as bad\n",(int)(nrej[2]+0.5),nrej[2]/nrej[3]);
   return 0;
}
