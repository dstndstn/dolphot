#include "dolphot_common.h"

// Standard photometry flags
fntype outfn;
reopenableFile *fdata;
FILE **fpsf,**fres;
long *fileRPSF;

void fakestar(int ext,int z,double x0,double y0,double xw,double yw,int tw,double *ct0) {
   int img,ix,iy,x1,y1,nbg;
   double x,y;
   float pa,pb,pc;
   static float*bg=NULL;

   if (bg==NULL) {
      bg=(float*)calloc(Nimg,FLOATSIZE);
      if (!bg) merr();
   }
   for (img=0;img<Nimg;img++) {
      shift(img,x0,y0,&x,&y,1);
      ix=(int)(x+RPSF[img])-RPSF[img];
      iy=(int)(y+RPSF[img])-RPSF[img];
      getpsfpars(img,xw,yw,(&pa)-img,(&pb)-img,(&pc)-img);
      bg[img]=0.;
      nbg=0;
      for (y1=-RPSF[img];y1<=RPSF[img] && iy+y1<dataim[img].Y;y1++) if (iy+y1>=0) for (x1=-RPSF[img];x1<=RPSF[img] && ix+x1<dataim[img].X;x1++) if (ix+x1>=0 && dataOK(img,ix+x1,iy+y1) && y1>=-rphot[img] && y1<=rphot[img] && x1>=-rphot[img] && x1<=rphot[img]) {
	 bg[img]+=data[img][iy+y1][ix+x1];
	 nbg++;
      }
      if (nbg) bg[img]/=nbg;
      else bg[img]=0.0;
#ifdef USEWFPC2
      if (hstmode[img].inst==WFPC2) WFPC2fixfakemag(img,x0,y0,ct0,bg);
#endif
#ifdef USEACS
      if (hstmode[img].inst==ACS) ACSfixfakemag(img,x0,y0,ct0,bg);
#endif
#ifdef USEWFC3
      if (hstmode[img].inst==WFC3) WFC3fixfakemag(img,x0,y0,ct0,bg);
#endif
      calc1psf(img,x,y,RPSF[img],pa,pb,pc,1,0);
      for (y1=-RPSF[img];y1<=RPSF[img] && iy+y1<dataim[img].Y;y1++) if (iy+y1>=0) for (x1=-RPSF[img];x1<=RPSF[img] && ix+x1<dataim[img].X;x1++) if (ix+x1>=0) {
	 if (dataOK(img,ix+x1,iy+y1)) {
	    if (RandomFake) data[img][iy+y1][ix+x1]+=poiss(ct0[img]*psf[y1][x1]);
	    else data[img][iy+y1][ix+x1]+=ct0[img]*psf[y1][x1];
	 }
      }
   }
   return;
}

void fakestars(int ext,int fld) {
   FILE *f;
   int e,z,img,tw=1;
   double x,y,xw=0,yw=0,*ctin;
   char str[8001];

   if ((f=fopen(FakeStars,"r"))==NULL) {
      printf("****Error opening FakeStars file %s\n",FakeStars);
      exit(-1);
   }
   ctin=(double*)calloc(Nimg,DOUBLESIZE); if (!ctin) merr();
   while (fscanf(f,"%d %d",&e,&z)==2) {
      if (e==ext && z==fld+1) {
	 fscanf(f,"%lf %lf",&x,&y);
#ifdef USEWFPC2
	 WFPC2readfakemag(f);
#endif
#ifdef USEACS
	 ACSreadfakemag(f);
#endif
#ifdef USEWFC3
	 WFC3readfakemag(f);
#endif
	 for (img=0;img<Nimg;img++) if (hstmode[img].inst==NONE) fscanf(f,"%lf",ctin+img);
	 fgets(str,8001,f);
	 fakestar(ext,fld,x,y,xw,yw,tw,ctin);
      }
      else fgets(str,8001,f);
   }
   fclose(f);
   free(ctin);
   return;
}

void procframe(int ext) {
   int img,x,y,z,skip;
   char ch[2880];
   imtype timg;
   //FILE *f;

   if (!isimage(dataim)) {
      for (img=0;img<Timg;img++) {
	 fopenagain(fdata+img);
	 readimage(fdata[img].f,dataim+img);
	 freclose(fdata+img);
	 freeimg(dataim[img].data,dataim[img].X,dataim[img].Y,dataim[img].Z);
      }
      return;
   }
   if (Timg>Nimg) {
      X=dataim[Nimg].X;
      Y=dataim[Nimg].Y;
   }
   else {
      X=dataim[0].X;
      Y=dataim[0].Y;
   }
   for (img=0;img<Timg;img++) {
      data[img]=allocchip(dataim[img].X,dataim[img].Y);
   }
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
	 for (img=0;img<Timg;img++) {
	    if (img<Nimg) {
	       float* pad;
	       int nwrite = RPSF[img];
	       if (fileRPSF[img]<nwrite) nwrite = fileRPSF[img];
	       int npad = fileRPSF[img] - nwrite;
	       pad = (float*)calloc(fileRPSF[img]*2+1,sizeof(float));
	       for (y=0;y<npad;y++) ffread(pad,fileRPSF[img]*2+1,4,fpsf[img]);
	       for (y=-nwrite;y<=nwrite;y++) {
		  if (npad>0) ffread(pad,npad,4,fpsf[img]);
		  ffread(poff[img][y]-nwrite,nwrite*2+1,4,fpsf[img]);
		  if (npad>0) ffread(pad,npad,4,fpsf[img]);
	       }
	       for (y=0;y<npad;y++) ffread(pad,fileRPSF[img]*2+1,4,fpsf[img]);
	       free(pad);
	       if (!FakeStarPSF) for (y=-RPSF[img];y<=RPSF[img];y++) for (x=-RPSF[img];x<=RPSF[img];x++) poff[img][y][x]=0;
	    }
	    else {
	       for (y=-RPSF[img];y<=RPSF[img];y++) for (x=-RPSF[img];x<=RPSF[img];x++) poff[img][y][x]=0;
	    }
	 }
	 fakestars(ext,z);
      }
      for (img=0;img<Nimg;img++) {
	 memcpy(&timg,dataim+img,sizeof(imtype));
	 timg.bscale=1.;
	 timg.bzero=0.;
	 timg.bits=-32;
	 writechip(fres[img],data[img],&timg);
      }
      freecirc();
   }
   for (img=0;img<Timg;img++) {
      skip=fabs(dataim[img].bits)/8*dataim[img].X*dataim[img].Y*dataim[img].Z+dataim[img].pcount;
      skip=((skip+2879)/2880)*2880-skip;
      fdata[img].lastoffset+=skip;
      if (img<Nimg) {
	 memset(ch,0,skip);
	 fwrite(ch,1,skip,fres[img]);
	 skip=fabs(dataim[img].bits)/8*(2*fileRPSF[img]+1)*(2*fileRPSF[img]+1)*dataim[img].Z+dataim[img].pcount;
	 skip=((skip+2879)/2880)*2880-skip;
	 fread(ch,1,skip,fpsf[img]);
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

   if (argc<2) {
      printf("****Usage: %s <output> <<options>>\n",*argv);
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
   initimgdata();
   paramfile("dolphot.param",&dolphotparam);
   for (i=2;i<argc;i++) {
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
   fpsf=(FILE**)calloc(Nimg,PTRSIZE);
   fres=(FILE**)calloc(Nimg,PTRSIZE);
   fileRPSF=(long*)calloc(Nimg,LONGSIZE);
   if (!fdata || !fpsf || !fres || !fileRPSF) merr();
   strcpy(outfn,argv[1]);
   sprintf(str,"%s.info",outfn);
   finfo=fopen(str,"r");
   if (!finfo) {
      printf("****Error opening %s\n",str);
      exit(-1);
   }
   for (img=0;img<Timg;img++) {
      ftype tfits;

      sprintf(str,"%s.fits",base[img]);
      fopenfirst(fdata+img,str,"rb",0);
      fdata[img].f=readfitsh(str,&tfits,1);
      freclose(fdata+img);
      memcpy(dataim+img,&(tfits.img),sizeof(imtype));
      datahd[img].Nmax=datahd[img].Ncards=dataim[img].Ncards;
      datahd[img].cards=(cardtype*)calloc(datahd[img].Nmax,sizeof(cardtype));
      memcpy(datahd[img].cards,dataim[img].cards,sizeof(cardtype)*datahd[img].Ncards);
      read_cardvals(img);
      if (img==0) Next=tfits.Next;
      else if (Next!=tfits.Next) {
	 printf("****Number of extensions are not the same\n");
	 exit(-1);
      }
      if (imgdata[img].RPSF<=0) imgdata[img].RPSF = RPSF0;
      if (img<Nimg) {
	 sprintf(str,"%s.%d.mod.fits",outfn,img+1);
	 if (isimage(&(tfits.img))) {
	    tfits.img.bscale=1.;
	    tfits.img.bzero=0.;
	    tfits.img.bits=-32;
	 }
	 fres[img]=writefitsh(str,&tfits,0);
	 sprintf(str,"%s.%d.psf.fits",outfn,img+1);
	 fpsf[img]=readfitsh(str,&tfits,0);
	 if (isimage(&(tfits.img))) {
	    if (tfits.img.X!=tfits.img.Y || tfits.img.X%2!=1) {
	       printf("PSF files not from DOLPHOT\n");
	       exit(-1);
	    }
	    fileRPSF[img] = (tfits.img.X-1)/2;
	    if (tfits.img.X!=2*imgdata[img].RPSF+1 || tfits.img.Y!=2*imgdata[img].RPSF+1) {
	       printf("Resizing PSF from %d x %d\n",tfits.img.X,tfits.img.Y);
	    }
	 }
	 else {
	    readimage(fpsf[img],&(tfits.img));
	    freeimg(tfits.img.data,tfits.img.X,tfits.img.Y,tfits.img.Z);
	 }
      }
   }

   // set up circles and photometry weights, etc.
   procframe(0);
   for (img=0;img<Timg;img++) if (dataim[img].Nmax) free(dataim[img].cards);
   for (ext=0;ext<Next;ext++) {
      for (img=0;img<Timg;img++) {
	 ftype tfits;

	 fopenagain(fdata+img);
	 readexth(fdata[img].f,&(tfits.img),1);
	 freclose(fdata+img);
	 memcpy(dataim+img,&(tfits.img),sizeof(imtype));
	 read_cardvals(img);
	 if (img!=0 && dataim[0].Z!=dataim[img].Z) {
	    printf("****Number of chips are not the same size\n");
	    exit(-1);
	 }
	 if (img<Nimg) {
	    if (isimage(&(tfits.img))) {
	       tfits.img.bscale=1.;
	       tfits.img.bzero=0.;
	       tfits.img.bits=-32;
	    }
	    writeexth(fres[img],&(tfits.img),0);
	    readexth(fpsf[img],&(tfits.img),0);
	    if (isimage(&(tfits.img))) {
	       if (tfits.img.X!=tfits.img.Y || tfits.img.X%2!=1) {
		  printf("PSF files not from DOLPHOT\n");
		  exit(-1);
	       }
	       fileRPSF[img] = (tfits.img.X-1)/2;
	       if (tfits.img.X!=2*imgdata[img].RPSF+1 || tfits.img.Y!=2*imgdata[img].RPSF+1) {
		  printf("Resizing PSF to %d x %d\n",tfits.img.X,tfits.img.Y);
	       }
	    }
	    else {
	       readimage(fpsf[img],&(tfits.img));
	       freeimg(tfits.img.data,tfits.img.X,tfits.img.Y,tfits.img.Z);
	    }
	 }
      }
      procframe(ext+1);
      for (img=0;img<Timg;img++) if (dataim[img].Nmax) free(dataim[img].cards);
   }
   fclose(finfo);
   for (img=0;img<Nimg;img++) {
      fclose(fpsf[img]);
      fclose(fres[img]);
   }
   return 0;
}
