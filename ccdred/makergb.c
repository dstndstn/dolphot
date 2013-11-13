#include <fits.h>
#include <unistd.h>
typedef unsigned char byte;

int SCALE;
double mm[3][2],CSCALE;
ftype fits[3];

#ifdef BYTE_SWAP
#define lfwrite fwrite
#else
size_t lfwrite(const void*ptr,size_t l,size_t n,FILE*f) {
   char t[4];
   if (n==1 && l==2) {
      memcpy(t+0,ptr+1,1);
      memcpy(t+1,ptr+0,1);
   }
   else if (n==1 && l==4) {
      memcpy(t+0,ptr+3,1);
      memcpy(t+1,ptr+2,1);
      memcpy(t+2,ptr+1,1);
      memcpy(t+3,ptr+0,1);
   }
   else return fwrite(ptr,l,n,f);
   return fwrite(t,l,n,f);
}
#endif

#define ULIM 1.0
#define LLIM -2.0
void getsky(int img) {
   int cont=1;
   double av=0,sd=-1;

   while (cont) {
      int x,y,N=0;
      double t1=0,t2=0;
      for (y=0;y<fits[img].img.Y;y++) for (x=0;x<fits[img].img.X;x++) if (fits[img].img.data[0][y][x]!=0 && (sd<0 || (fits[img].img.data[0][y][x]-av<=ULIM*sd && fits[img].img.data[0][y][x]-av>=LLIM*sd))) {
	 t1+=fits[img].img.data[0][y][x];
	 N++;
      }
      if (!N) return;
      t1/=N;
      for (y=0;y<fits[img].img.Y;y++) for (x=0;x<fits[img].img.X;x++) if (fits[img].img.data[0][y][x]!=0 && (sd<0 || (fits[img].img.data[0][y][x]-av<=ULIM*sd && fits[img].img.data[0][y][x]-av>=LLIM*sd))) t2+=(fits[img].img.data[0][y][x]-t1)*(fits[img].img.data[0][y][x]-t1);
      t2=sqrt(t1/N);
      cont=0;
      for (y=0;y<fits[img].img.Y && !cont;y++) for (x=0;x<fits[img].img.X && !cont;x++) if (fits[img].img.data[0][y][x]!=0 && (sd<0 || (fits[img].img.data[0][y][x]-av<=ULIM*sd && fits[img].img.data[0][y][x]-av>=LLIM*sd)) && (fits[img].img.data[0][y][x]-t1>t2*ULIM || fits[img].img.data[0][y][x]-t1<t2*LLIM)) cont=1;
      av=t1;
      sd=t2;
   }
   printf("Image %d: sky=%f +/- %f\n",img+1,av,sd);
}

void fakeg(double f1,double f2) {
   int x,y;
   double t;

   fitscopy(fits,fits+1);
   t=f1+f2;
   mm[1][0]=mm[0][0]*f1+mm[2][0]*f2;
   mm[1][1]=(mm[0][1]*f1+mm[2][1]*f2-mm[1][0])/t+mm[1][0];
   for (y=0;y<fits[1].img.Y;y++) for (x=0;x<fits[1].img.X;x++) {
      if (fits[0].img.data[0][y][x]<mm[0][0]) fits[1].img.data[0][y][x]=0;
      else if (fits[0].img.data[0][y][x]>mm[0][1]) fits[1].img.data[0][y][x]=f1*mm[0][1];
      else fits[1].img.data[0][y][x]=f1*fits[0].img.data[0][y][x];
      if (fits[2].img.data[0][y][x]<mm[2][0]) fits[1].img.data[0][y][x]+=0;
      else if (fits[2].img.data[0][y][x]>mm[2][1]) fits[1].img.data[0][y][x]+=f2*mm[2][1];
      else fits[1].img.data[0][y][x]+=f2*fits[2].img.data[0][y][x];
   }
}

void cval(double in[3],byte out[3]) {
   double i0,i1,ii;
   int z;
   i0=ii=(in[0]+in[1]+in[2])/3.;
   if (i0<=0 || i0>=1) i0=i1=1.;
   else switch(SCALE) {
   case 1:i1=pow(i0,0.25); break;
   case 2:i1=pow(i0,0.333); break;
   case 3:i1=pow(i0,0.5); break;
   case 4:i1=log(i0*1000.+1.)/log(1001.); break;
   default:i1=i0; break;
   }
   for (z=0;z<3;z++) {
      in[z]=i1/i0*((in[z]-ii)*CSCALE+ii);
      if (in[z]<=0) out[z]=0;
      else if (in[z]>=1) out[z]=255;
      else out[z]=(byte)(in[z]*256);
   }
}

void reduce(byte bgr[3]) {
   /*
   //256-color grayscale;
   bgr[0]=bgr[1]=bgr[2]=(bgr[0]+bgr[1]+bgr[2]+1)/3;
   */
   /*
   //32-color grayscale;
   byte x;
   x=(bgr[0]+bgr[1]+bgr[2]+1)/3;
   x=x/8*8;
   bgr[0]=bgr[1]=bgr[2]=x;
   */
   /*
   //32-color color;
   bgr[0]=bgr[0]/8*8;
   bgr[1]=bgr[1]/8*8;
   bgr[2]=bgr[2]/8*8;
   */
}

void writebmp(char *fn) {
   FILE *f;
   int x,y,z,w,h,pw;
   byte used[256][256][32];
   byte p2[8]={1,2,4,8,16,32,64,128};
   long nused=0;

   w=fits[0].img.X;
   pw=(w*3+3)/4*4;
   h=fits[0].img.Y;
   for (x=0;x<256;x++) for (y=0;y<256;y++) for (z=0;z<32;z++) used[x][y][z]=0;
   for (y=0;y<h;y++) for (x=0;x<w;x++) {
      byte bgr[3];
      double bgr0[3];
      for (z=0;z<3;z++) {
	 if (fits[z].img.data[0][y][x]<=mm[z][0]) bgr0[z]=0.;
	 else if (fits[z].img.data[0][y][x]>=mm[z][1]) bgr0[z]=1.;
	 else bgr0[z]=(fits[z].img.data[0][y][x]-mm[z][0])/(mm[z][1]-mm[z][0]);
      }
      cval(bgr0,bgr);
      reduce(bgr);
      if (!(used[bgr[0]][bgr[1]][bgr[2]/32]&p2[bgr[2]%8])) nused++;
      used[bgr[0]][bgr[1]][bgr[2]/32]|=p2[bgr[2]%8];
   }
   printf("using %ld colors\n",nused);
   if ((f=fopen(fn,"wb"))==NULL) {
      printf("Cannot write %s\n",fn);
      exit(-1);
   }
   //header;
   {
      long l; short s;
      s=19778; lfwrite(&s,2,1,f); //file type ID;
      l=54+3*pw*h; lfwrite(&l,4,1,f); //file size;
      s=0; lfwrite(&s,2,1,f); //reserved;
      s=0; lfwrite(&s,2,1,f); //reserved;
      l=54; lfwrite(&l,4,1,f); //header+info size;
   }
   //info;
   {
      long l; short s;
      l=40; lfwrite(&l,4,1,f); //size of section;
      lfwrite(&w,4,1,f); lfwrite(&h,4,1,f); //width & height;
      s=1; lfwrite(&s,2,1,f); //number of target planes;
      s=24; lfwrite(&s,2,1,f); //bits;
      l=0; lfwrite(&l,4,1,f); //compression;
      l=3*pw*h; lfwrite(&l,4,1,f); //image bytes;
      l=0; lfwrite(&l,4,1,f); //X-pixels/meter (4000?);
      l=0; lfwrite(&l,4,1,f); //Y-pixels/meter (4000?);
      lfwrite(&nused,4,1,f); //used colors;
      lfwrite(&nused,4,1,f); //important colors;
   }
   for (y=0;y<h;y++) {
      for (x=0;x<w;x++) {
	 byte bgr[3];
	 double bgr0[3];
	 for (z=0;z<3;z++) {
	    if (fits[z].img.data[0][y][x]<=mm[z][0]) bgr0[z]=0.;
	    else if (fits[z].img.data[0][y][x]>=mm[z][1]) bgr0[z]=1.;
	    else bgr0[z]=(fits[z].img.data[0][y][x]-mm[z][0])/(mm[z][1]-mm[z][0]);
	 }
	 cval(bgr0,bgr);
	 reduce(bgr);
	 fwrite(bgr,1,3,f);
      }
      for (x=w*3;x<pw;x++) {
	 byte b=0;
	 fwrite(&b,1,1,f);
      }
   }
   fclose(f);
}

int main(int argc,char **argv) {
   int img,FAKEG=0;

   if (argc!=13) {
      printf("Usage: %s <file min max> * bgr <output> <scale> <color scale>\n",*argv);
      printf("  scale: 0=linear, 1=sqrt, 2=**1/3, 3=**1/4, 4=log\n");
      printf("  color scale: 1.0 normally, 0.5 for UVI\n");
      return -1;
   }
   SCALE=atoi(argv[11]);
   CSCALE=atof(argv[12]);
   for (img=0;img<3;img++) {
      mm[img][0]=atof(argv[img*3+2]);
      mm[img][1]=atof(argv[img*3+3]);
      if (!access(argv[1+img*3],F_OK)) {
	 readfits(argv[1+img*3],fits+img,1);
	 if (fits[img].Next) {
	    printf("Extensions not supported\n");
	    return -1;
	 }
	 if (fits[img].img.Z>1) {
	    printf("Stacked images not supported\n");
	    return -1;
	 }
	 if (fits[img].img.Y<2) {
	    printf("Not an image\n");
	    return -1;
	 }
	 getsky(img);
      }
      else if (img!=1) {
	 printf("**File %s not found\n",argv[img*3+1]);
	 return -1;
      }
      else FAKEG=1;
   }
   if (fits[0].img.X!=fits[2].img.X || fits[0].img.Y!=fits[2].img.Y) {
      printf("Images are not the same size\n");
      return -1;
   }
   if (FAKEG) fakeg(mm[1][0],mm[1][1]);
   else if (fits[0].img.X!=fits[1].img.X || fits[0].img.Y!=fits[1].img.Y) {
      printf("Images are not the same size\n");
      return -1;
   }
   writebmp(argv[10]);
   return 0;
}
