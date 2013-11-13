#include <fits.h>
extern char* fitsnum(double x);
ftype fits;

int main(int argc,char**argv) {
   int i=1,j,ver=0;
   char *str;

   if (argc>=2 && !strcasecmp(argv[1],"-full")) {
      ver=1;
      i++;
   }
   for (;i<argc;i++) {
      fclose(readfitsh(argv[i],&fits,0));
      if (!ver) {
	 if (fits.img.Z<2) printf("%s: %dx%d",argv[i],fits.img.X,fits.img.Y);
	 else printf("%s: %dx%dx%d",argv[i],fits.img.X,fits.img.Y,fits.img.Z);
	 switch(fits.img.bits) {
	 case -32:printf(" float"); break;
	 case 32:printf(" long"); break;
	 case 16:printf(" short"); break;
	 case 8:printf(" byte"); break;
	 }
	 str=getcardval(&(fits.img),"EXPTIME",0);
	 if (str[0]) printf(", %ss",str);
	 str=getcardval(&(fits.img),"FILTER",0);
	 if (str[0]) printf(", %s",str);
	 str=getcardval(&(fits.img),"OBSTYPE",0);
	 if (str[0]) printf(", type=%s",str);
	 str=getcardval(&(fits.img),"OBJECT",0);
	 if (!str[0]) str=getcardval(&(fits.img),"TARGNAME",0);
	 if (str[0]) printf(", name=%s",str);
	 printf("\n");
      }
      else {
	 printf("SIMPLE  =                    T / Standard FITS format\n");
	 printf("BITPIX  = %20d / Data type\n",fits.img.bits);
	 if (fits.img.Z>1) {
	    printf("NAXIS   =                    3 / Number of axes\n");
	    printf("NAXIS1  = %20d /\n",fits.img.X);
	    printf("NAXIS2  = %20d /\n",fits.img.Y);
	    printf("NAXIS3  = %20d /\n",fits.img.Z);
	 }
	 else if (fits.img.Y>1) {
	    printf("NAXIS   =                    2 / Number of axes\n");
	    printf("NAXIS1  = %20d /\n",fits.img.X);
	    printf("NAXIS2  = %20d /\n",fits.img.Y);
	 }
	 else if (fits.img.X>0) {
	    printf("NAXIS   =                    1 / Number of axes\n");
	    printf("NAXIS1  = %20d /\n",fits.img.X);
	 }
	 else {
	    printf("NAXIS   =                    0 / Number of axes\n");
	 }
	 if (fits.Next<=0) {printf("EXTEND  =                    F / There are no standard extensions\n");}
	 else {
	    printf("EXTEND  =                    T / There may be standard extensions\n");
	    printf("NEXTEND =       %14d /\n",fits.Next);
	 }
	 printf("BSCALE  = %20s / Image scale\n",fitsnum(fits.img.bscale));
	 printf("BZERO   = %20s / Image zero\n",fitsnum(fits.img.bzero));
	 for (j=0;j<fits.img.Ncards;j++) printf("%s\n",fits.img.cards[j]);
      }
      if (fits.img.Nmax) free(fits.img.cards);
   }
   return 0;
}
