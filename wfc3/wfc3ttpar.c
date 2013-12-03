#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <unistd.h>

int main(int argc,char**argv) {
   char fn[161],pstr[321],wmag[21]=" -wmag=2";
   FILE *f;
   int x,y,c,CHAN;

   if (argc!=3) {
      printf("Usage: %s <filter> <channel> <<options>>\n",*argv);
      printf("Supported channels\n");
      printf("   UVIS  PSF is for UVIS channel\n");
      printf("   IR    PSF is for IR channel\n");
      return -1;
   }
   if (!strcasecmp(argv[2],"UVIS")) CHAN=1;
   else if (!strcasecmp(argv[2],"IR")) CHAN=0;
   else {
      printf("Unknown channel \"%s\"\n",argv[2]);
      return -1;
   }
   for (x=3;x<argc;x++) {
      //else {
	 printf("Unknown option \"%s\"\n",argv[x]);
	 return -1;
      //}
   }
   if (argv[1][strlen(argv[1])-1]=='N') wmag[0]=0;
   sprintf(fn,"%s/tmp",BASEDIR);
   mkdir(fn,00755);
   sprintf(fn,"%s/tmp/%s_WFC3",BASEDIR,argv[1]);
   mkdir(fn,00755);
   sprintf(fn,"%s/tmp/%s_WFC3/uvis.pos",BASEDIR,argv[1]);
   if ((f=fopen(fn,"w"))==NULL) {
      printf("Cannot write %s\n",fn);
      return -1;
   }
   for (y=128;y<2048;y+=256) for (x=128;x<4096;x+=256) fprintf(f,"%d %d\n",x,y);
   fclose(f);
   sprintf(fn,"%s/tmp/%s_WFC3/ir.pos",BASEDIR,argv[1]);
   if ((f=fopen(fn,"w"))==NULL) {
      printf("Cannot write %s\n",fn);
      return -1;
   }
   for (y=64;y<1024;y+=128) for (x=64;x<1024;x+=128) fprintf(f,"%d %d\n",x,y);
   fclose(f);
   for (c=0;c<3;c++) if ((c==0)==(CHAN==0)) {
      sprintf(fn,"%s/tmp/%s_WFC3/tiny1.resp",BASEDIR,argv[1]);
      if ((f=fopen(fn,"w"))==NULL) {
	 printf("Cannot write %s\n",fn);
	 return -1;
      }
      if (c) fprintf(f,"22\n%d\n@%s/tmp/%s_WFC3/uvis.pos\n",c,BASEDIR,argv[1]);
      else fprintf(f,"23\n@%s/tmp/%s_WFC3/ir.pos\n",BASEDIR,argv[1]);
      if (wmag[0]) fprintf(f,"%s\n1\n11\n",argv[1]);
      else fprintf(f,"%s\n",argv[1]);
      if (!c) fprintf(f,"5.0\nir\n");
      else fprintf(f,"3.0\nuvis%d\n",c);
      fclose(f);
      if (!c) sprintf(pstr,"cat %s/tmp/%s_WFC3/tiny1.resp | %s/tiny1 %s/tmp/%s_WFC3/ir.par%s > %s/tmp/%s_WFC3/tiny1.log",BASEDIR,argv[1],TTDIR,BASEDIR,argv[1],wmag,BASEDIR,argv[1]);
      else sprintf(pstr,"cat %s/tmp/%s_WFC3/tiny1.resp | %s/tiny1 %s/tmp/%s_WFC3/uvis%d.par%s >> %s/tmp/%s_WFC3/tiny1.log",BASEDIR,argv[1],TTDIR,BASEDIR,argv[1],c,wmag,BASEDIR,argv[1]);
      if ((f=popen(pstr,"w"))!=NULL) pclose(f);
      else printf("Cannot run script %d\n",c);
   }
   sprintf(fn,"%s/tmp/%s_WFC3/uvis.pos",BASEDIR,argv[1]);
   unlink(fn);
   sprintf(fn,"%s/tmp/%s_WFC3/ir.pos",BASEDIR,argv[1]);
   unlink(fn);
   sprintf(fn,"%s/tmp/%s_WFC3/tiny1.resp",BASEDIR,argv[1]);
   unlink(fn);
   sprintf(fn,"%s/tmp/%s_WFC3/runtiny",BASEDIR,argv[1]);
   if ((f=fopen(fn,"w"))==NULL) {
      printf("Cannot write %s\n",fn);
      return -1;
   }
   fprintf(f,"#! /bin/csh\n");
   fprintf(f,"rm tiny.log; touch tiny.log\n");
   if (CHAN==0) {
      fprintf(f,"tiny2 ir.par >> tiny.log\n");
      for (x=0;x<64;x++) fprintf(f,"tiny3 ir.par sub=9 pos=%d >> tiny.log\n",x);
   }
   else {
      fprintf(f,"tiny2 uvis1.par >> tiny.log\n");
      fprintf(f,"tiny2 uvis2.par >> tiny.log\n");
      for (x=0;x<128;x++) fprintf(f,"tiny3 uvis1.par sub=5 pos=%d >> tiny.log\n",x);
      for (x=0;x<128;x++) fprintf(f,"tiny3 uvis2.par sub=5 pos=%d >> tiny.log\n",x);
   }
   fclose(f);
   chmod(fn,00755);
   return 0;
}
