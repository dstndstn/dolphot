#include "dolphot_common.h"

// Standard photometry flags
reopenableFile *fdata;
fntype positions;

void fakestars(int ext,int fld) {
   FILE *f;
   int img0,e,z,img;
   double x0,y0,x,y;
   char str[8001];

   if ((f=fopen(positions,"r"))==NULL) {
      printf("****Error opening positions file %s\n",positions);
      exit(-1);
   }
   while (fscanf(f,"%d %d %d",&img0,&e,&z)==3) {
      if (e==ext && z==fld+1) {
	 fscanf(f,"%lf %lf",&x0,&y0);
	 printf("%3d %d %d %7.2f %7.2f",img0,e,z,x0,y0);
	 if (img0>=1 && img0<=Nimg) shift(img0-1,x0,y0,&x0,&y0,-1);
	 printf(" %d %7.2f %7.2f",0,x0,y0);
	 for (img=0;img<Nimg;img++) {
	    shift(img,x0,y0,&x,&y,1);
	    if (!posOK(img,x,y)) x=y=-1.0;
	    printf(" %d %7.2f %7.2f",img+1,x,y);
	 }
	 printf("\n");
      }
      else fgets(str,8001,f);
   }
   fclose(f);
   return;
}

void procframe(int ext) {
   int img,z,skip;

   // skip data if not an image
   if (!isimage(dataim)) {
      for (img=0;img<Timg;img++) {
	 fopenagain(fdata+img);
	 readimage(fdata[img].f,dataim+img);
	 freclose(fdata+img);
	 freeimg(dataim[img].data,dataim[img].X,dataim[img].Y,dataim[img].Z);
      }
      return;
   }
   // set up X, Y
   if (Timg>Nimg) {
      X=dataim[Nimg].X;
      Y=dataim[Nimg].Y;
   }
   else {
      X=dataim[0].X;
      Y=dataim[0].Y;
   }
   // allocate memory for images
   for (img=0;img<Timg;img++) data[img]=allocchip(dataim[img].X,dataim[img].Y);

   // loop through chips in this extension
   for (z=0;z<dataim[0].Z;z++) {
      initimgpars();
#ifdef USEWFPC2
      wfpc2radii(0); // accounts for plate scale vs. chip # and shifts
#endif
      initcirc();
#ifdef USEACS
      acsinitpsf();
#endif
#ifdef USEWFC3
      wfc3initpsf();
#endif
#ifdef USEWFPC2
      wfpc2initpsf();
#endif
      if (Nimg==Timg) DRIZZLE_BASE=0;
      else if (hstmode[Nimg].inst==NONE) DRIZZLE_BASE=1;
      for (img=0;img<Timg;img++) {
	 fopenagain(fdata+img);
	 readchip(fdata[img].f,data[img],dataim+img);
	 freclose(fdata+img);
      }
      if ((GUSE<0 || ext==GUSE) && (CUSE<0 || z==CUSE)) {
	 if (ext>0) printf("** Extension %d, Chip %d **\n",ext,z+1);
	 else printf("** Chip %d **\n",z+1);
	 fflush(stdout);
	 readinfo(finfo,ext,z);
	 fakestars(ext,z);
      }
      freecirc();
   }
   for (img=0;img<Timg;img++) {
      skip=fabs(dataim[img].bits)/8*dataim[img].X*dataim[img].Y*dataim[img].Z+dataim[img].pcount;
      skip=((skip+2879)/2880)*2880-skip;
      fdata[img].lastoffset+=skip;
      if (img<Nimg) {
	 freechip(data[img],dataim[img].X,dataim[img].Y);
      }
   }
#ifdef USEWFPC2
   wfpc2freepsf();
#endif
#ifdef USEACS
   acsfreepsf();
#endif
#ifdef USEWFC3
   wfc3freepsf();
#endif
   return;
}

int main(int argc,char**argv) {
   char str[82];
   int img,Next=0,ext,i;

   if (argc<3) {
      printf("****Usage: %s <output> <positions> <<options>>\n",*argv);
      printf("  -p<name>  for parameter file\n");
      printf("  x=y       to set flag x to value y\n");
      return 1;
   }
   INTSIZE=sizeof(int);
   LONGSIZE=sizeof(long);
   FLOATSIZE=sizeof(float);
   DOUBLESIZE=sizeof(double);
   PTRSIZE=sizeof(void*);
   if (sizeof(char)!=1) {printf("char size != 1\n"); return -1;}
   if (INTSIZE<4) {printf("int size < 4\n"); return -1;}
   if (LONGSIZE<4) {printf("long size < 4\n"); return -1;}
   if (FLOATSIZE!=4) {printf("float size != 4\n"); return -1;}
   if (DOUBLESIZE<8) {printf("double size < 8\n"); return -1;}
   if (sizeof(char*)!=PTRSIZE) {printf("char* size != void* size\n"); return -1;}
   if (sizeof(int*)!=PTRSIZE) {printf("int* size != void* size\n"); return -1;}
   if (sizeof(float*)!=PTRSIZE) {printf("float* size != void* size\n"); return -1;}
   if (sizeof(double*)!=PTRSIZE) {printf("double* size != void* size\n"); return -1;}
   if (sizeof(FILE*)!=PTRSIZE) {printf("FILE* size != void* size\n"); return -1;}

   // read parameters
   strcpy(positions,argv[2]);
   initimgdata();
   paramfile("dolphot.param",&dolphotparam);
   for (i=3;i<argc;i++) {
      if (!strncmp(argv[i],"-p",2) || !strncmp(argv[i],"-P",2)) paramfile1(argv[i]+2,&dolphotparam);
      else parseparam(argv[i],&dolphotparam);
   }

   // sanity check
   if (RPSF0<RAper0) perr("RPSF must be at least as large as RAper");

   // allocate memory
   alloc_img_data();
   // update params for camera
#ifdef USEWFPC2
   wfpc2initparam();
#endif
#ifdef USEACS
   acsinitparam();
#endif
#ifdef USEWFC3
   wfc3initparam();
#endif

   fdata=(reopenableFile*)calloc(Timg,sizeof(reopenableFile));
   if (!fdata) merr();
   sprintf(str,"%s.info",argv[1]);
   finfo=fopen(str,"r");
   if (!finfo) {
      printf("****Error opening %s\n",str);
      exit(-1);
   }
   // Process main frame
   for (img=0;img<Timg;img++) {
      ftype tfits;

      // open file and read FITS header
      sprintf(str,"%s.fits",base[img]);
      fopenfirst(fdata+img,str,"rb",0);
      fdata[img].f=readfitsh(str,&tfits,1);
      freclose(fdata+img);
      memcpy(dataim+img,&(tfits.img),sizeof(imtype));
      // copy FITS cards from dataim into datahd
      datahd[img].Nmax=datahd[img].Ncards=dataim[img].Ncards;
      datahd[img].cards=(cardtype*)calloc(datahd[img].Nmax,sizeof(cardtype));
      memcpy(datahd[img].cards,dataim[img].cards,sizeof(cardtype)*datahd[img].Ncards);
      read_cardvals(img);
      if (img==0) Next=tfits.Next;
      else if (Next!=tfits.Next) {
	 printf("****Number of extensions are not the same\n");
	 exit(-1);
      }
   }
   procframe(0);
   for (img=0;img<Timg;img++) if (dataim[img].Nmax) free(dataim[img].cards);

   // Process extensions
   for (ext=0;ext<Next;ext++) {
      for (img=0;img<Timg;img++) {
	 ftype tfits;

	 // reopen file and read FITS extension headers
	 fopenagain(fdata+img);
	 readexth(fdata[img].f,&(tfits.img),1);
	 freclose(fdata+img);
	 memcpy(dataim+img,&(tfits.img),sizeof(imtype));
	 read_cardvals(img);
	 if (img!=0 && dataim[0].Z!=dataim[img].Z) {
	    printf("****Number of chips are not the same size\n");
	    exit(-1);
	 }
      }
      procframe(ext+1);
      for (img=0;img<Timg;img++) if (dataim[img].Nmax) free(dataim[img].cards);
   }
   fclose(finfo);
   return 0;
}
