#include <fits.h>
#include "ccdproc.h"

typedef double double4[4];

double *mult,*dmin,*dmax,*noise,*setsky,*gain,c2,nsig,r2in,r2out,nrej[4]={0,0,0,0};
int N,min,*dx,*dy,rout,skip,NX,NY;
imtype *data,*data0;
ftype *fits;

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
   for (yy=-rout;yy<=rout;yy+=skip) for (xx=-rout;xx<=rout;xx+=skip) if (yy+y>=0 && yy+y<data0[img].Y && xx+x>=0 && xx+x<data0[img].X && data0[img].data[z][yy+y][xx+x]>dmin[img] && data0[img].data[z][yy+y][xx+x]<dmax[img]) {
      r2=yy*yy+xx*xx;
      if (r2>=r2in && r2<=r2out) slist[ns++]=data0[img].data[z][yy+y][xx+x];
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
   int i,j,x,y,xx,yy;
   double tmult=0;
   static double4 *sedge;
   static double *list,*cmplist,*sky;
   static int first=1,*EX,*EY;

   if (first) {
      cmplist=(double*)calloc(sizeof(double),N);
      sky=(double*)calloc(sizeof(double),N);
      list=(double*)calloc(sizeof(double),N);
      sedge=(double4*)calloc(sizeof(double4),N);
      EX=(int*)calloc(sizeof(int),N);
      EY=(int*)calloc(sizeof(int),N);
      if (!cmplist || !sky || !list || !sedge || !EX || !EY) merr();
      first=0;
   }
   xx=r2in;
   r2in=0;
   for (i=0;i<N;i++) {
      EX[i]=EY[i]=rout;
      if (rout<0) EX[i]=EY[0]=0;
      if (EX[i]>data0[i].X/4) EX[i]=data0[i].X/4;
      if (EY[i]>data0[i].Y/4) EY[i]=data0[i].Y/4;
      tmult+=mult[i];
      sedge[i][0]=(calcsky(i,data0[i].X/4,EY[i],z)+calcsky(i,data0[i].X/2,EY[i],z)+calcsky(i,3*data0[i].X/4,EY[i],z))/3;
      sedge[i][1]=(calcsky(i,data0[i].X/4,data0[i].Y-1-EY[i],z)+calcsky(i,data0[i].X/2,data0[i].Y-1-EY[i],z)+calcsky(i,3*data0[i].X/4,data0[i].Y-1-EY[i],z))/3;
      sedge[i][2]=(calcsky(i,EX[i],data0[i].Y/4,z)+calcsky(i,EX[i],data0[i].Y/2,z)+calcsky(i,EX[i],3*data0[i].Y/4,z))/3;
      sedge[i][3]=(calcsky(i,data0[i].X-1-EX[i],data0[i].Y/4,z)+calcsky(i,data0[i].X-1-EX[i],data0[i].Y/2,z)+calcsky(i,data0[i].X-1-EX[i],3*data0[i].Y/4,z))/3;
   }
   r2in=xx;
   for (y=0;y<NY;y++) for (x=0;x<NX;x++) {
      int n=0;
      double tmp,hi,lo,mid,dhi,dlo;

      // if (y%100==0 && x==0) printf("  y=%d\n",y);
      for (i=0;i<N;i++) {
	  if (y+dy[i]>=EY[i] && y+dy[i]<data0[i].Y-EY[i] && x+dx[i]>=EX[i] && x+dx[i]<data0[i].X-EX[i]) sky[i]=calcsky(i,x+dx[i],y+dy[i],z);
	  else if (y+dy[i]<EY[i]) sky[i]=sedge[i][0];
	  else if (y+dy[i]>=data0[i].Y-EY[i]) sky[i]=sedge[i][1];
	  else if (x+dx[i]<EX[i]) sky[i]=sedge[i][2];
	  else sky[i]=sedge[i][3];
      }

      // cmplist contains minimum (data-sky)/mult values for this pixel. "mid" is either minimum or median of this list.
      // list contains (data-sky)/mult, mult, sky, and noise for each image
      for (i=0;i<N;i++) {
	 if (y+dy[i]>=0 && y+dy[i]<data0[i].Y && x+dx[i]>=0 && x+dx[i]<data0[i].X) {
	    if (data0[i].data[z][y+dy[i]][x+dx[i]]>dmin[i] && data0[i].data[z][y+dy[i]][x+dx[i]]<dmax[i]) {
	       tmp=(data0[i].data[z][y+dy[i]][x+dx[i]]-sky[i])/mult[i];
	       for (j=n;j>0 && cmplist[j-1]>tmp;j--) cmplist[j]=cmplist[j-1];
	       cmplist[j]=tmp;
	       list[i]=tmp;
	       n++;
	    }
	    else if (data0[i].data[z][y+dy[i]][x+dx[i]]<=dmin[i]) {
	       list[i]=2.e30;
	       nrej[3]++;
	    }
	    else {
	       list[i]=2.e30;
	       nrej[2]++;
	    }
	 }
	 else list[i]=2.e30;
      }
      if (n>1) { // only bother if there are multiple valid frames
	 if (min) mid=cmplist[0];
	 else mid=cmplist[(n-1)/2];

	 // cmplist contains minimum (data-sky)/mult values for 3x3 area centered on this pixel. "lo" is either minimum or median of this list.
	 int ncmp=0;
	 for (i=0;i<N;i++) { // should verify list[i]<1.e30
	    tmp=2.e30;
	    for (yy=y-1;yy<=y+1;yy++) for (xx=x-1;xx<=x+1;xx++) if (yy+dy[i]>=0 && yy+dy[i]<data0[i].Y && xx+dx[i]>=0 && xx+dx[i]<data0[i].X && data0[i].data[z][yy+dy[i]][xx+dx[i]]>dmin[i] && data0[i].data[z][yy+dy[i]][xx+dx[i]]<dmax[i] && (data0[i].data[z][yy+dy[i]][xx+dx[i]]-sky[i])/mult[i]<tmp) tmp=(data0[i].data[z][yy+dy[i]][xx+dx[i]]-sky[i])/mult[i];
	    if (tmp<1.e30) {
	       for (j=ncmp;j>0 && cmplist[j-1]>tmp;j--) cmplist[j]=cmplist[j-1];
	       cmplist[j]=tmp;
	       ncmp++;
	    }
	 }
	 if (ncmp==0) lo=mid;
	 else {
	    if (min) lo=cmplist[0];
	    else lo=cmplist[(ncmp-1)/2];
	    if (lo<mid) lo=mid;
	 }

	 // cmplist contains maximum (data-sky)/mult values for 3x3 area centered on this pixel. "hi" is either minimum or median of this list.
	 ncmp=0;
	 for (i=0;i<N;i++) { // should verify list[i]<1.e30
	    tmp=-2.e30;
	    for (yy=y-1;yy<=y+1;yy++) for (xx=x-1;xx<=x+1;xx++) if (yy+dy[i]>=0 && yy+dy[i]<data0[i].Y && xx+dx[i]>=0 && xx+dx[i]<data0[i].X && data0[i].data[z][yy+dy[i]][xx+dx[i]]>dmin[i] && data0[i].data[z][yy+dy[i]][xx+dx[i]]<dmax[i] && (data0[i].data[z][yy+dy[i]][xx+dx[i]]-sky[i])/mult[i]>tmp) tmp=(data0[i].data[z][yy+dy[i]][xx+dx[i]]-sky[i])/mult[i];
	    if (tmp>-1.e30) {
	       for (j=ncmp;j>0 && cmplist[j-1]>tmp;j--) cmplist[j]=cmplist[j-1];
	       cmplist[j]=tmp;
	       ncmp++;
	    }
	 }
	 if (ncmp==0) hi=mid;
	 else {
	    if (min) hi=cmplist[0];
	    else hi=cmplist[(ncmp-1)/2];
	    if (hi<mid) hi=mid;
	 }
	 for (i=0;i<n;i++) if (list[i]<1.e30) {
	    double dd,tlo,thi;

	    tlo=lo+sky[i]/mult[i];
	    thi=hi+sky[i]/mult[i];
	    dlo=noise[i]/gain[i]*nsig/mult[i]; // readout noise in scaled counts
	    if (lo>0) dd=nsig*nsig*tlo/gain[i]/mult[i]+c2*c2*lo*lo; // photon noise and linear noise in scaled counts
	    else if (sky[i]>0) dd=nsig*nsig*sky[i]/gain[i]/mult[i]/mult[i]; // photon noise in scaled counts
	    else dd=0;
	    dlo=sqrt(dlo*dlo+dd);
	    dhi=noise[i]/gain[i]*nsig/mult[i];
	    if (hi>0) dd=nsig*nsig*thi/gain[i]/mult[i]+c2*c2*hi*hi;
	    else if (sky[i]>0) dd=nsig*nsig*sky[i]/gain[i]/mult[i]/mult[i];
	    else dd=0;
	    dhi=sqrt(dhi*dhi+dd);
	    if (list[i]>hi+dhi) {
	       data[i].data[z][y+dy[i]][x+dx[i]]=dmin[i]-1;
	       nrej[0]++;
	    }
	    else if (list[i]<lo-dlo) {
	       data[i].data[z][y+dy[i]][x+dx[i]]=dmin[i]-1;
	       nrej[1]++;
	    }
	 }
      }
   }
}

int main(int argc,char**argv) {
   int i,j,z,ext;

   if (argc<14 || (argc-9)%5) {
      printf("Usage: %s <fits files, dx, dy, scale, sky> Xfin Yfin\n  <rsky-in> <rsky-out> <sky skip> <reg factor> <sigma clip> <use min>\n",*argv);
      printf("  (%d arguments given)\n",argc);
      printf("Scale=0 to scale by exposure time\n");
      printf("Recommend reg=1, sig=3, use min=0\n");
      printf("Set rsky-out<0 to use sky values\n");
      return 1;
   }
   N=(argc-9)/5;
   dmin=(double*)calloc(N,sizeof(double));
   dmax=(double*)calloc(N,sizeof(double));
   noise=(double*)calloc(N,sizeof(double));
   gain=(double*)calloc(N,sizeof(double));
   fits=(ftype*)calloc(N,sizeof(ftype));
   data=(imtype*)calloc(N,sizeof(imtype));
   data0=(imtype*)calloc(N,sizeof(imtype));
   mult=(double*)calloc(N,sizeof(double));
   setsky=(double*)calloc(N,sizeof(double));
   dx=(int*)calloc(N,sizeof(int));
   dy=(int*)calloc(N,sizeof(int));
   if (!dmin || !dmax || !noise || !gain || !fits || !data || !data0 || !mult || !setsky || !dx || !dy) {
      printf("Memory allocation error\n");
      return 1;
   }
   NX=atoi(argv[argc-8]);
   NY=atoi(argv[argc-7]);
   r2in=atof(argv[argc-6])*atof(argv[argc-6]);
   r2out=atof(argv[argc-5])*atof(argv[argc-5]);
   if (atof(argv[argc-5])<0) rout=-1;
   else rout=(int)(atof(argv[argc-5])+0.999);
   skip=atoi(argv[argc-4]);
   if (rout>0) rout=(rout*skip+skip-1)/skip;
   c2=atof(argv[argc-3]);
   nsig=atof(argv[argc-2]);
   min=atoi(argv[argc-1]);
   for (i=0;i<N;i++) {
      readfits(argv[i*5+1],fits+i,0);
      if (i>0) {
	 if (fits[i].Next!=fits[0].Next) {
	    printf("Number of extensions does not match\n");
	    return 1;
	 }
	 if (isimage(&fits[i].img)!=isimage(&fits[0].img)) {
	    printf("Number of valid images does not match\n");
	    return 1;
	 }
	 if (fits[i].img.Z!=fits[0].img.Z) {
	    printf("Number of chips does not match\n");
	    return 1;
	 }
	 for (j=0;j<fits[i].Next;j++) {
	    if (isimage(fits[i].ext+j)!=isimage(fits[0].ext+j)) {
	       printf("Number of valid images does not match\n");
	       return 1;
	    }
	    if (fits[i].ext[j].Z!=fits[0].ext[j].Z) {
	       printf("Number of chips does not match\n");
	       return 1;
	    }
	 }
      }
   }
   for (ext=0;ext<=fits[0].Next;ext++) {
      if ((ext==0 && isimage(&fits[0].img)) || (ext>0 && isimage(fits[0].ext+ext-1))) {
	 if (ext>0) printf("Extension %d\n",ext);
	 for (i=0;i<N;i++) {
	    double IEXP,IGAIN,IRN,IEPOCH,IAIR,IEXP0;
	    float IDMAX,IDMIN;
	    
	    if (ext==0) data[i]=fits[i].img;
	    else data[i]=fits[i].ext[ext-1];
	    imcopy(data+i,data0+i);
	    parsecards(data+i,&IGAIN,&IRN,&IEXP,&IDMIN,&IDMAX,&IEPOCH,&IAIR,&IEXP0,0,1);
	    dmin[i]=IDMIN;
	    dmax[i]=IDMAX;
	    noise[i]=IRN;
	    gain[i]=IGAIN;
	    dx[i]=atoi(argv[2+5*i]);
	    dy[i]=atoi(argv[3+5*i]);
	    mult[i]=atof(argv[4+5*i]);
	    if (mult[i]<=0) mult[i]=IEXP;
	    setsky[i]=atof(argv[5+5*i]);
	 }
	 for (z=0;z<data[0].Z;z++) cleanchip(z);
	 for (i=0;i<N;i++) freeim(data0+i);
      }
   }
   for (i=0;i<N;i++) {
      writefits(argv[i*5+1],fits+i,0);
      freefits(fits+i);
   }
   printf("%7d (%7.4f%%) rejected as CRs\n",(int)(nrej[0]+0.5),nrej[0]/nrej[3]);
   printf("%7d (%7.4f%%) rejected as low\n",(int)(nrej[1]+0.5),nrej[1]/nrej[3]);
   printf("%7d (%7.4f%%) rejected as bad\n",(int)(nrej[2]+0.5),nrej[2]/nrej[3]);
   return 0;
}
