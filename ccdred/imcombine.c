#include <fits.h>
#include "ccdproc.h"

typedef double double4[4];

double *mult,*dmin,*dmax,*setsky,r2in,r2out,PGAIN;
float PDMIN,PDMAX;
int N,min,avg,*dx,*dy,rout,skip,*IX,*IY;
imtype *data;
ftype *fits,cfits;

double calcsky(int img,int x,int y,int z) {
   int xx,yy,r2,j;
   int ns=0,cont=1;
   double sky=0;
   static double *slist;
   static int first=1;

   if (rout<0) return setsky[img];
   if (first) {
      slist=(double*)calloc(sizeof(double),(rout/2+1)*(rout/2+1));
      if (!slist) merr();
      first=0;
   }
   if (r2out<=r2in || rout<=0) return 0;
   for (yy=-rout;yy<=rout;yy+=skip) for (xx=-rout;xx<=rout;xx+=skip) if (yy+y>=0 && yy+y<IY[img] && xx+x>=0 && xx+x<IX[img] && data[img].data[z][yy+y][xx+x]>dmin[img] && data[img].data[z][yy+y][xx+x]<dmax[img]) {
      r2=yy*yy+xx*xx;
      if (r2>=r2in && r2<=r2out) slist[ns++]=data[img].data[z][yy+y][xx+x];
   }
   if (ns) while (cont) {
      double sig=0;
      cont=0;
      sky=0;
      for (j=0;j<ns;j++) sky+=slist[j];
      sky/=ns;
      for (j=0;j<ns;j++) sig+=(slist[j]-sky)*(slist[j]-sky);
      if (ns>1) sig=sqrt(sig/(ns-1));
      else sig=sqrt(sig);
      for (j=0;j<ns;j++) if (slist[j]<sky-2.5*sig || slist[j]>sky+2.0*sig) {
	 cont=1;
	 slist[j]=slist[--ns];
	 j--;
      }
   }
   return sky;
}

void cleanchip(int z) {
   int i,x,y,xx;
   double tmult=0;
   static double4 *sedge;
   static double *sky;
   static int first=1,*EX,*EY;

   if (!isimage(data)) return;
   if (first) {
      sky=(double*)calloc(sizeof(double),N);
      sedge=(double4*)calloc(sizeof(double4),N);
      EX=(int*)calloc(sizeof(int),N);
      EY=(int*)calloc(sizeof(int),N);
      if (!sky || !sedge || !EX || !EY) merr();
      first=0;
   }
   xx=r2in;
   r2in=0;
   for (i=0;i<N;i++) {
      EX[i]=EY[i]=rout;
      if (rout<0) EX[i]=EY[0]=0;
      if (EX[i]>IX[i]/4) EX[i]=IX[i]/4;
      if (EY[i]>IY[i]/4) EY[i]=IY[i]/4;
      tmult+=mult[i];
      sedge[i][0]=(calcsky(i,IX[i]/4,EY[i],z)+calcsky(i,IX[i]/2,EY[i],z)+calcsky(i,3*IX[i]/4,EY[i],z))/3;
      sedge[i][1]=(calcsky(i,IX[i]/4,IY[i]-1-EY[i],z)+calcsky(i,IX[i]/2,IY[i]-1-EY[i],z)+calcsky(i,3*IX[i]/4,IY[i]-1-EY[i],z))/3;
      sedge[i][2]=(calcsky(i,EX[i],IY[i]/4,z)+calcsky(i,EX[i],IY[i]/2,z)+calcsky(i,EX[i],3*IY[i]/4,z))/3;
      sedge[i][3]=(calcsky(i,IX[i]-1-EX[i],IY[i]/4,z)+calcsky(i,IX[i]-1-EX[i],IY[i]/2,z)+calcsky(i,IX[i]-1-EX[i],3*IY[i]/4,z))/3;
   }
   r2in=xx;
   for (y=0;y<cfits.img.Y;y++) for (x=0;x<cfits.img.X;x++) {
      int sat=0;
      double c=0,s=0,tsky=0;

      for (i=0;i<N;i++) {
	  if (y+dy[i]>=EY[i] && y+dy[i]<IY[i]-EY[i] && x+dx[i]>=EX[i] && x+dx[i]<IX[i]-EX[i]) sky[i]=calcsky(i,x+dx[i],y+dy[i],z);
	  else if (y+dy[i]<EY[i]) sky[i]=sedge[i][0];
	  else if (y+dy[i]>=IY[i]-EY[i]) sky[i]=sedge[i][1];
	  else if (x+dx[i]<EX[i]) sky[i]=sedge[i][2];
	  else sky[i]=sedge[i][3];
	  tsky+=sky[i];

	  if (y+dy[i]>=0 && y+dy[i]<IY[i] && x+dx[i]>=0 && x+dx[i]<IX[i]) {
	     if (data[i].data[z][y+dy[i]][x+dx[i]]>dmin[i] && data[i].data[z][y+dy[i]][x+dx[i]]<dmax[i]) {
		c+=data[i].data[z][y+dy[i]][x+dx[i]]-sky[i];
		s+=mult[i];
	     }
	     else if (data[i].data[z][y+dy[i]][x+dx[i]]>=dmax[i]) sat=1;
	  }
      }
      if (s<=0) {
	 if (!sat) cfits.img.data[z][y][x]=PDMIN-1;
	 else cfits.img.data[z][y][x]=PDMAX+1;
      }
      else {
	 cfits.img.data[z][y][x]=c/s*tmult+tsky;
	 if (avg) cfits.img.data[z][y][x]/=N;
      }
   }
}

int main(int argc,char**argv) {
   int i,j,z,ext,NX,NY;
   FILE **fin,*fout=NULL;

   if (argc<13 || (argc-8)%5) {
      printf("Usage: %s <fits files, dx, dy, scale, sky> Xfin Yfin\n  <rsky-in> <rsky-out> <sky skip> <average> <fits output>\n",*argv);
      printf("  (%d arguments given)\n",argc);
      printf("Scale=0 to scale by exposure time\n");
      printf("Set rsky-out<0 to use sky values\n");
      return 1;
   }
   N=(argc-8)/5;
   dmin=(double*)calloc(N,sizeof(double));
   dmax=(double*)calloc(N,sizeof(double));
   fits=(ftype*)calloc(N,sizeof(ftype));
   data=(imtype*)calloc(N,sizeof(imtype));
   mult=(double*)calloc(N,sizeof(double));
   setsky=(double*)calloc(N,sizeof(double));
   dx=(int*)calloc(N,sizeof(int));
   dy=(int*)calloc(N,sizeof(int));
   IX=(int*)calloc(N,sizeof(int));
   IY=(int*)calloc(N,sizeof(int));
   fin=(FILE**)calloc(N,sizeof(FILE*));
   if (!dmin || !dmax || !fits || !data || !mult || !setsky || !dx || !dy || !fin || !IX || !IY) {
      printf("Memory allocation error\n");
      return 1;
   }
   NX=atoi(argv[argc-7]);
   NY=atoi(argv[argc-6]);
   r2in=atof(argv[argc-5])*atof(argv[argc-5]);
   r2out=atof(argv[argc-4])*atof(argv[argc-4]);
   if (atof(argv[argc-4])<0) rout=-1;
   else rout=(int)(atof(argv[argc-4])+0.999);
   skip=atoi(argv[argc-3]);
   if (rout>0) rout=(rout*skip+skip-1)/skip;
   avg=atoi(argv[argc-2]);
   for (i=0;i<N;i++) {
      FILE *f;

      f=readfitsh(argv[i*5+1],fits+i,0);
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
      float DMIN=0,DMAX=0;

      if (ext>0) printf("Extension %d\n",ext);
      for (i=0;i<N;i++) {
	 double IEXP,IGAIN,IRN,IEPOCH,IAIR,IEXP0;
	 float IDMAX,IDMIN;

	 if (!ext) fin[i]=readfitsh(argv[1+i*5],fits+i,0);
	 else readexth(fin[i],&(fits[i].img),0);
	 readimage(fin[i],&(fits[i].img));
	 memcpy(data+i,&(fits[i].img),sizeof(imtype));
	 parsecards(data+i,&IGAIN,&IRN,&IEXP,&IDMIN,&IDMAX,&IEPOCH,&IAIR,&IEXP0,0,1);
	 if (!i) {
	    if (!ext) strcpy(cfits.img.xtension,"IMAGE");
	    else strcpy(cfits.img.xtension,data[i].xtension);
	    cfits.img.Ncards=0;
	    cfits.img.X=data[i].X;
	    cfits.img.Y=data[i].Y;
	    cfits.img.Z=data[i].Z;
	    if (isimage(&(cfits.img))) {
	       cfits.img.X=NX;
	       cfits.img.Y=NY;
	       cfits.img.bits=-32;
	       cfits.img.bzero=0.;
	       cfits.img.bscale=1.;
	    }
	    else {
	       cfits.img.X=0;
	       cfits.img.Y=0;
	       cfits.img.bits=data[i].bits;
	       cfits.img.bzero=data[i].bzero;
	       cfits.img.bscale=data[i].bscale;
	    }
	    for (j=0;j<data[i].Ncards;j++) addcard(&cfits.img,data[i].cards[j]);
	    GAIN=IGAIN;
	    DMIN=IDMIN/IEXP;
	    DMAX=IDMAX/IEXP;
	 }
	 else {
	    if (IDMIN/IEXP<DMIN) DMIN=IDMIN/IEXP;
	    if (IDMAX/IEXP>DMAX) DMAX=IDMAX/IEXP;
	 }
	 dmin[i]=IDMIN;
	 dmax[i]=IDMAX;
	 IX[i]=data[i].X;
	 IY[i]=data[i].Y;
	 dx[i]=atoi(argv[2+5*i]);
	 dy[i]=atoi(argv[3+5*i]);
	 mult[i]=atof(argv[4+5*i]);
	 if (mult[i]<=0) mult[i]=IEXP;
	 setsky[i]=atof(argv[5+5*i]);
	 RN+=IRN*IRN;
	 EPOCH+=IEXP*IEPOCH;
	 AIR+=IEXP*IAIR;
	 EXP+=IEXP;
	 EXP0+=IEXP0*IEXP;
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
      for (z=0;z<cfits.img.Z;z++) cleanchip(z);
      if (avg) {
	 EXP/=N;
	 DMIN=PDMIN/N;
	 DMAX=PDMAX/N;
	 GAIN=PGAIN*N;
      }
      insertcards(&(cfits.img),GAIN,RN,EXP,DMIN,DMAX,EPOCH,AIR,EXP0);
      if (!ext) fout=writefitsh(argv[argc-1],&cfits,1);
      else writeexth(fout,&(cfits.img),0);
      writeimage(fout,&(cfits.img));
      for (i=0;i<N;i++) {
	 if (fits[i].img.Nmax) free(fits[i].img.cards);
	 fits[i].img.Nmax=0;
	 if (!ext && fits[i].Next>0) free(fits[i].ext);
      }
      for (i=0;i<N;i++) freeimg(data[i].data,data[i].X,data[i].Y,data[i].Z);
      freeimg(cfits.img.data,cfits.img.X,cfits.img.Y,cfits.img.Z);
   }
   for (i=0;i<N;i++) fclose(fin[i]);
   fclose(fout);
   return 0;
}
