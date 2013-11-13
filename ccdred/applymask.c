#include <fits.h>
#include <unistd.h>


int eqimg(imtype*im1,imtype*im2) {
   int i1,i2;

   i1=isimage(im1);
   i2=isimage(im2);
   if (!i1 && !i2) return 1;
   if (!i1 || !i2) return 0;
   if (im1->X!=im2->X || im1->Y!=im2->Y || im1->Z!=im2->Z) return 0;
   return 1;
}

void applymask(imtype*mask,imtype*data,float DMIN) {
   int x,y,z;
   for (z=0;z<mask->Z;z++) for (y=0;y<mask->Y;y++) for (x=0;x<mask->X;x++) {
      if (mask->data[z][y][x]==0.0) data->data[z][y][x]=DMIN-1;
   }
}

int main(int argc,char**argv) {
   int i,j;
   ftype mask,data;
   float DMIN,DMINe;

   if (argc<3) {
      printf("Usage: %s <mask> <images>\n",*argv);
      printf("  mask image defined such that zero = bad pixel\n");
      return -1;
   }
   readfits(argv[1],&mask,1);
   for (i=2;i<argc;i++) {
      readfits(argv[i],&data,1);
      if (mask.Next!=data.Next || !eqimg(&(mask.img),&(data.img))) {
	 printf("FITS images are not the same size\n");
      }
      else {
	 int sOK=1;
	 parsecards(&data.img,NULL,NULL,NULL,&DMIN,NULL,NULL,NULL,NULL,0,1);
	 if (isimage(&(mask.img))) {
	    applymask(&(mask.img),&(data.img),DMIN);
	 }
	 for (j=0;j<mask.Next;j++) {
	    DMINe=DMIN;
	    if (!eqimg(mask.ext+j,data.ext+j)) {
	       printf("FITS images are not the same size\n");
	       sOK=0;
	    }
	    else if (isimage(mask.ext+j)) {
	       parsecards(data.ext+j,NULL,NULL,NULL,&DMINe,NULL,NULL,NULL,NULL,0,0);
	       applymask(mask.ext+j,data.ext+j,DMINe);
	    }
	 }
	 if (sOK) {
	    writefits(argv[i],&data,1);
	 }
      }
   }
   return 0;
}
