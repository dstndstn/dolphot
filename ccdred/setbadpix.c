#include <fits.h>

int main(int argc,char**argv) {
   int ext,x,y,x1,x2,y1,y2,z;
   float DMIN;
   ftype fits;
   imtype*img;

   if (argc!=8) {
      printf("Usage: %s <file> <ext> <chip> <x1> <y1> <x2> <y2>\n",*argv);
      return 1;
   }
   readfits(argv[1],&fits,1);
   ext=atoi(argv[2]);
   if (ext<0 || ext>fits.Next) {
      printf("Illegal extension number; 0-%d allowable\n",fits.Next);
      return 1;
   }
   if (!ext) img=&(fits.img);
   else img=fits.ext+ext;
   z=atoi(argv[3]);
   if (z<1 || z>img->Z) {
      printf("Illegal chip number; 1-%d allowable\n",img->Z);
      return 1;
   }
   x1=atoi(argv[4]);
   y1=atoi(argv[5]);
   x2=atoi(argv[6]);
   y2=atoi(argv[7]);
   if (x1<0 || y1<0 || x2>=img->X || y2>=img->Y) {
      printf("Illegal bounds; image is %dx%d pixels\n",img->X,img->Y);
      return 1;
   }
   if (x2<x1 || y2<y1) {
      printf("Illegal range; x2 and y2 must be greater than x1 and y1\n");
      return 1;
   }
   parsecards(&(fits.img),NULL,NULL,NULL,&DMIN,NULL,NULL,NULL,NULL,0,1);
   if (ext>0) parsecards(img,NULL,NULL,NULL,&DMIN,NULL,NULL,NULL,NULL,0,0);
   printf("Setting to %f\n",DMIN-1);
   for (y=y1;y<=y2;y++) for (x=x1;x<=x2;x++) img->data[z-1][y][x]=DMIN-1;
   writefits(argv[1],&fits,1);
   return 0;
}
