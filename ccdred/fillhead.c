#include <fits.h>

ftype fits;

void imfill(imtype*img) {
   double GAIN,RN,EXP,EPOCH,AIR,EXP0;
   float DMIN,DMAX;

   if (!isimage(img)) return;
   parsecards(&(fits.img),&GAIN,&RN,&EXP,&DMIN,&DMAX,&EPOCH,&AIR,&EXP0,0,1);
   if (img!=&(fits.img)) parsecards(img,&GAIN,&RN,&EXP,&DMIN,&DMAX,&EPOCH,&AIR,&EXP0,0,0);
   insertcards(img,GAIN,RN,EXP,DMIN,DMAX,EPOCH,AIR,EXP0);
}

int main(int argc,char**argv) {
   int i;

   if (argc!=3) {
      printf("Usage: %s <input> <output>\n",*argv);
      return -1;
   }
   readfits(argv[1],&fits,1);
   imfill(&(fits.img));
   for (i=0;i<fits.Next;i++) imfill(fits.ext+i);
   writefits(argv[2],&fits,1);
   return 0;
}
