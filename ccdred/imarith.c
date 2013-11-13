#include <fits.h>
#include <unistd.h>

ftype fitsin1,fitsin2,fitsout;
int im2=0;
double val2=0.;
char op;

int eqimg(imtype*im1,imtype*im2) {
   int i1,i2;

   i1=isimage(im1);
   i2=isimage(im2);
   if (!i1 && !i2) return 1;
   if (!i1 || !i2) return 0;
   if (im1->X!=im2->X || im1->Y!=im2->Y || im1->Z!=im2->Z) return 0;
   return 1;
}

float func(float f1,float f2) {
   switch(op) {
   case '+':return f1+f2;
   case '-':return f1-f2;
   case '*':return f1*f2;
   case '/':return f1/f2;
   default:
      printf("Illegal operator \"%c\"\n",op);
      exit(-1);
   }
   return 0.;
}

void img1smooth(imtype*im1) {
   float DMIN,DMAX;
   double c,ct;
   imtype im2;
   static int ksize=-1;
   static float**kernel;
   int x,y,z,xx,yy;

   if (!isimage(im1)) return;
   if (ksize<0) {
      if (op=='b' || op=='B') ksize=(int)val2;
      else ksize=(int)(val2*3+0.5);
      if (ksize<0) return;
      kernel=(float**)calloc(2*ksize+1,sizeof(float*));
      if (!kernel) merr();
      kernel+=ksize;
      for (y=-ksize;y<=ksize;y++) {
	 kernel[y]=(float*)calloc(2*ksize+1,sizeof(float));
	 if (!kernel[y]) merr();
	 kernel[y]+=ksize;
	 for (x=-ksize;x<=ksize;x++) {
	    if (op=='b' || op=='B') kernel[y][x]=1.;
	    else kernel[y][x]=exp(-0.5*(y*y+x*x)/val2/val2);
	 }
      }
   }
   imcopy(im1,&im2);
   parsecards(&(fitsout.img),NULL,NULL,NULL,&DMIN,&DMAX,NULL,NULL,NULL,0,1);
   if (im1!=&(fitsout.img)) parsecards(im1,NULL,NULL,NULL,&DMIN,&DMAX,NULL,NULL,NULL,0,0);
   for (z=0;z<im1->Z;z++) for (y=0;y<im1->Y;y++) for (x=0;x<im1->X;x++) {
      c=ct=0.;
      for (yy=y-ksize;yy<=y+ksize;yy++) if (yy>=0 && yy<im1->Y) for (xx=x-ksize;xx<=x+ksize;xx++) if (xx>=0 && xx<im1->X && im2.data[z][yy][xx]>DMIN && im2.data[z][yy][xx]<DMAX) {
	 c+=im2.data[z][yy][xx]*kernel[yy-y][xx-x];
	 ct+=kernel[yy-y][xx-x];
      }
      if (ct>0.) im1->data[z][y][x]=c/ct;
      else if (im2.data[z][y][x]>=DMAX) im1->data[z][y][x]=DMAX;
      else im1->data[z][y][x]=DMIN;
   }
   freeim(&im2);
}

void img1arith(imtype*im1) {
   float DMIN,DMAX,tDMIN,tDMAX;
   int x,y,z,flip=0;

   if (!isimage(im1)) return;
   if (val2==0 && op=='/') {
      printf("Cannot divide by zero\n");
      exit(-1);
   }
   if (val2<0 && (op=='*' || op=='/')) flip=1;
   parsecards(&(fitsout.img),NULL,NULL,NULL,&DMIN,&DMAX,NULL,NULL,NULL,0,1);
   if (im1!=&(fitsout.img)) parsecards(im1,NULL,NULL,NULL,&DMIN,&DMAX,NULL,NULL,NULL,0,0);
   if (fabs(val2)<=1 && (op=='*' || op=='/')) {
      if (!flip) {
	 tDMIN=DMIN;
	 tDMAX=DMAX;
      }
      else {
	 tDMIN=DMAX;
	 tDMAX=DMIN;
      }
   }
   else if (!flip) {
      tDMIN=func(DMIN,val2);
      tDMAX=func(DMAX,val2);
   }
   else {
      tDMIN=func(DMAX,val2);
      tDMAX=func(DMIN,val2);
   }
   for (z=0;z<im1->Z;z++) for (y=0;y<im1->Y;y++) for (x=0;x<im1->X;x++) {
      if (im1->data[z][y][x]<=DMIN) {
	 if (!flip) im1->data[z][y][x]=tDMIN-1;
	 else im1->data[z][y][x]=tDMAX+1;
      }
      else if (im1->data[z][y][x]>=DMAX) {
	 if (!flip) im1->data[z][y][x]=tDMAX+1;
	 else im1->data[z][y][x]=tDMIN-1;
      }
      else im1->data[z][y][x]=func(im1->data[z][y][x],val2);
   }
   printf("Changing minimum data value to %f\n",tDMIN);
   printf("Changing maximum data value to %f\n",tDMAX);
   insertcards(im1,-1.e30,-1.e30,-1.e30,tDMIN,tDMAX,-1.e30,-1.e30,-1.e30);
}

void img2arith(imtype*im1,imtype*im2) {
   float DMIN1,DMAX1,DMIN2,DMAX2,tDMIN,tDMAX,q;
   int x,y,z;

   if (!isimage(im1)) return;
   parsecards(&(fitsout.img),NULL,NULL,NULL,&DMIN1,&DMAX1,NULL,NULL,NULL,0,1);
   if (im1!=&(fitsout.img)) parsecards(im1,NULL,NULL,NULL,&DMIN1,&DMAX1,NULL,NULL,NULL,0,0);
   parsecards(&(fitsin2.img),NULL,NULL,NULL,&DMIN2,&DMAX2,NULL,NULL,NULL,0,1);
   if (im2!=&(fitsin2.img)) parsecards(im2,NULL,NULL,NULL,&DMIN2,&DMAX2,NULL,NULL,NULL,0,0);
   if (op=='+') {
      tDMIN=DMIN1+DMIN2;
      tDMAX=DMAX1+DMAX2;
   }
   else if (op=='-') {
      tDMIN=DMIN1-DMAX2;
      tDMAX=DMAX1-DMIN2;
   }
   else if (op=='*') {
      tDMIN=tDMAX=DMIN1*DMIN2;
      q=DMIN1*DMAX2; if (q<tDMIN) tDMIN=q; else if (q>tDMAX) tDMAX=q;
      q=DMAX1*DMIN2; if (q<tDMIN) tDMIN=q; else if (q>tDMAX) tDMAX=q;
      q=DMAX1*DMAX2; if (q<tDMIN) tDMIN=q; else if (q>tDMAX) tDMAX=q;
   }
   else if (op=='/') {
      if (DMIN1<=0 && DMAX1>=0) {tDMIN=-1.e10; tDMAX=1.e10;}
      else if (DMIN2<=0 && DMAX2>=0) {tDMIN=-1.e10; tDMAX=1.e10;}
      else {
	 tDMIN=tDMAX=DMIN1*DMIN2;
	 q=DMIN1*DMAX2; if (q<tDMIN) tDMIN=q; else if (q>tDMAX) tDMAX=q;
	 q=DMAX1*DMIN2; if (q<tDMIN) tDMIN=q; else if (q>tDMAX) tDMAX=q;
	 q=DMAX1*DMAX2; if (q<tDMIN) tDMIN=q; else if (q>tDMAX) tDMAX=q;
      }
   }
   else {
      printf("Not yet implemented\n");
      return;
   }
   for (z=0;z<im1->Z;z++) for (y=0;y<im1->Y;y++) for (x=0;x<im1->X;x++) {
      if (im1->data[z][y][x]<=DMIN1 || im2->data[z][y][x]<=DMIN2) im1->data[z][y][x]=tDMIN-1;
      else if (im1->data[z][y][x]>=DMAX1 || im2->data[z][y][x]>=DMAX2) im1->data[z][y][x]=tDMAX+1;
      else im1->data[z][y][x]=func(im1->data[z][y][x],im2->data[z][y][x]);
   }
   printf("Changing minimum data value to %f\n",tDMIN);
   printf("Changing maximum data value to %f\n",tDMAX);
   insertcards(im1,-1.e30,-1.e30,-1.e30,tDMIN,tDMAX,-1.e30,-1.e30,-1.e30);
}

int main(int argc,char**argv) {
   int i;
   char *ptr;

   if (argc!=5) {
      printf("Usage: %s <image1> <operator> <image2/val> <output>\n",*argv);
      printf("  operators: + - *\n");
      printf("  G: gaussian smooth\n");
      printf("  B: boxcar smooth by value\n");
      return -1;
   }
   readfits(argv[1],&fitsin1,1);
   if (!access(argv[3],F_OK)) {
      readfits(argv[3],&fitsin2,1);
      im2=1;
      if (fitsin1.Next!=fitsin2.Next || !eqimg(&(fitsin1.img),&(fitsin2.img))) {
	 printf("FITS images are not the same size\n");
	 return -1;
      }
      for (i=0;i<fitsin1.Next;i++) if (!eqimg(&(fitsin1.img),&(fitsin2.img))) {
      }
   }
   else {
      val2=strtod(argv[3],&ptr);
      if (ptr==argv[3]) {
	 printf("Image 2 or value unreadable\n");
	 return -1;
      }
   }
   fitscopy(&fitsin1,&fitsout);
   op=argv[2][0];
   if (im2) {
      img2arith(&(fitsout.img),&(fitsin2.img));
      for (i=0;i<fitsout.Next;i++) img2arith(fitsout.ext+i,fitsin2.ext+i);
   }
   else if (op=='G' || op=='g' || op=='B' || op=='b') {
      img1smooth(&(fitsout.img));
      for (i=0;i<fitsout.Next;i++) img1smooth(fitsout.ext+i);
   }
   else {
      img1arith(&(fitsout.img));
      for (i=0;i<fitsout.Next;i++) img1arith(fitsout.ext+i);
   }
   writefits(argv[4],&fitsout,1);
   return 0;
}
