#include <fits.h>

ftype fitsin,fitsout;

int main(int argc,char**argv) {
   int img,i,j;
   char str[161],base[161];

   if (argc<2) {
      printf("Usage: %s <input files>\n",*argv);
      return -1;
   }
   for (img=1;img<argc;img++) {
      strcpy(base,argv[img]);
      for (i=strlen(base)-1;i>0 && base[i]!='.';i--);
      if (i>0) base[i]=0;
      if (strcasecmp(base+i+1,"fits")) printf("**%s does not have a .fits extension; skipping\n",argv[img]);
      else {
	 readfits(argv[img],&fitsin,1);
	 memcpy(&fitsout,&fitsin,sizeof(ftype));
	 fitsout.Next=0;
	 if (fitsout.img.X) {
	    sprintf(str,"%s.chip0.fits",base);
	    writefits(str,&fitsout,1);
	 }
	 fitsout.img.Nmax=0;
	 for (i=0;i<fitsin.Next;i++) {
	    memcpy(&(fitsout.img),&(fitsout.ext[i]),sizeof(imtype));
	    fitsout.img.Nmax=fitsout.img.Ncards=0;
	    if (fitsin.img.X==0) for (j=0;j<fitsin.img.Ncards;j++) addcard(&(fitsout.img),fitsin.img.cards[j]);
	    for (j=0;j<fitsin.ext[i].Ncards;j++) addcard(&(fitsout.img),fitsin.ext[i].cards[j]);
	    if (fitsout.img.X) {
	       sprintf(str,"%s.chip%d.fits",base,i+1);
	       writefits(str,&fitsout,1);
	    }
	    if (fitsout.img.Nmax) free(fitsout.img.cards);
	 }
	 freefits(&fitsin);
      }
   }
   return 0;
}
