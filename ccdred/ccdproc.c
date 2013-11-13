#include "ccdproc.h"

ftype fits,cal;
char obstype[81];

void convert(imtype*img,int ext) {
   if (!isimage(img)) return;
   img->bscale=1.;
   img->bzero=0.;
   img->bits=-32;
}

void special(void) {
   int i,j,x,y;

   if (!strcmp(imager,"mimo")) {
      float sat[4],bad[4];
      if (fits.Next!=4) {
	 printf("**Mini Mosaic images need to have four 1088x4096 extensions\n");
	 return;
      }
      for (i=0;i<fits.Next;i++) if (fits.ext[i].X!=1088 || fits.ext[i].Y!=4096 || fits.ext[i].Z!=1) {
	 printf("**Mini Mosaic images need to have four 1088x4096 extensions\n");
	 return;
      }
      for (i=0;i<fits.Next;i++) {
	 parsecards(&(fits.img),NULL,NULL,NULL,bad+i,sat+i,NULL,NULL,NULL,0,1);
	 parsecards(fits.ext+i,NULL,NULL,NULL,bad+i,sat+i,NULL,NULL,NULL,0,0);
      }
      for (i=0;i<fits.Next;i+=2) for (y=0;y<fits.ext[0].Y;y++) for (x=0;x<fits.ext[0].X;x++) if (fits.ext[i].data[0][y][x]>=sat[i]-0.5) {
	 for (j=0;j<fits.Next;j+=2) if (i!=j) fits.ext[j].data[0][y][x]=bad[j]-1.;
	 for (j=1;j<fits.Next;j+=2) if (i!=j) fits.ext[j].data[0][y][1087-x]=bad[j]-1.;
      }
      for (i=1;i<fits.Next;i+=2) for (y=0;y<fits.ext[0].Y;y++) for (x=0;x<fits.ext[0].X;x++) if (fits.ext[i].data[0][y][x]>=sat[i]-0.5) {
	 for (j=0;j<fits.Next;j+=2) if (i!=j) fits.ext[j].data[0][y][1087-x]=bad[j]-1.;
	 for (j=1;j<fits.Next;j+=2) if (i!=j) fits.ext[j].data[0][y][x]=bad[j]-1.;
      }
   }
}

void mask(imtype*img,int c) {
   int x,y,z;
   char fn[161];
   ftype mask;

   if (!isimage(img)) return;
   sprintf(fn,"%s/masks/%s_%d.fits",BASEDIR,imager,c);
   readfits(fn,&mask,0);
   printf("Reading bad pixel mask from %s\n",fn);
   if (mask.Next!=0 || mask.img.X!=img->X || mask.img.Y!=img->Y || mask.img.Z!=img->Z) printf("**Image size mismatch with %s\n",fn);
   else {
      float DMIN;
      parsecards(&(fits.img),NULL,NULL,NULL,&DMIN,NULL,NULL,NULL,NULL,0,1);
      if (img!=&(fits.img)) parsecards(img,NULL,NULL,NULL,&DMIN,NULL,NULL,NULL,NULL,0,0);
      for (z=0;z<img->Z;z++) for (y=0;y<img->Y;y++) for (x=0;x<img->X;x++) if (mask.img.data[z][y][x]>0.5) img->data[z][y][x]=DMIN-1;
   }
   freefits(&mask);
}

void overscan(imtype*img,int c) {
   int ox1,ox2,oy1,oy2,oz1,oz2,x,y,z,Ncorr=0;
   float DMIN,DMAX;
   double corr,avcorr=0;

   if (!isimage(img)) return;
   if (!getsec(img,oversec,&ox1,&ox2,&oy1,&oy2,&oz1,&oz2)) {
      printf("**No overscan section found\n");
      return;
   }
   else if (ox1==0 && ox1==img->X && oy1==0 && oy1==img->Y) {
      printf("**Overscan is entire chip; no correction made\n");
      return;
   }
   else if (oz1!=0 || oz2!=img->Z) printf("**Overscan selected only for certain chips\n");
   printf("Overscan section: [%d:%d,%d:%d]\n",ox1+1,ox2,oy1+1,oy2);
   parsecards(&(fits.img),NULL,NULL,NULL,&DMIN,&DMAX,NULL,NULL,NULL,0,1);
   if (img!=&(fits.img)) parsecards(img,NULL,NULL,NULL,&DMIN,&DMAX,NULL,NULL,NULL,0,0);
   for (z=oz1;z<oz2;z++) for (y=oy1;y<oy2;y++) for (x=ox1;x<ox2;x++) if (img->data[z][y][x]>DMIN && img->data[z][y][x]<DMAX) {
      avcorr+=img->data[z][y][x];
      Ncorr++;
   }
   if (Ncorr && avcorr>=0) {
      double tmp;
      avcorr/=Ncorr;
      tmp=DMIN-avcorr;
      for (z=0;z<img->Z;z++) for (y=0;y<img->Y;y++) for (x=0;x<img->X;x++) if (img->data[z][y][x]<=DMIN) img->data[z][y][x]=tmp-1;
      printf("Lowering minimum data value by %f\n",avcorr);
      DMIN=tmp;
      insertcards(img,-1.e30,-1.e30,-1.e30,DMIN,-1.e30,-1.e30,-1.e30,-1.e30);
   }
   else if (Ncorr && avcorr<0) {
      double tmp;
      avcorr/=Ncorr;
      tmp=DMAX-avcorr;
      for (z=0;z<img->Z;z++) for (y=0;y<img->Y;y++) for (x=0;x<img->X;x++) if (img->data[z][y][x]>=DMAX) img->data[z][y][x]=tmp+1;
      printf("Raising maximum data value by %f\n",-avcorr);
      DMAX=tmp;
      insertcards(img,-1.e30,-1.e30,-1.e30,-1.e30,DMAX,-1.e30,-1.e30,-1.e30);
   }
   if (OVERSCAN==1) for (z=oz1;z<oz2;z++) {
      Ncorr=0;
      corr=0.;
      for (y=oy1;y<oy2;y++) for (x=ox1;x<ox2;x++) if (img->data[z][y][x]>DMIN && img->data[z][y][x]<DMAX) {
	 corr+=img->data[z][y][x];
	 Ncorr++;
      }
      if (Ncorr) corr/=Ncorr;
      for (y=0;y<img->Y;y++) for (x=0;x<img->X;x++) if (img->data[z][y][x]>DMIN && img->data[z][y][x]<DMAX) img->data[z][y][x]-=corr;
      printf("Overscan average = %f\n",corr);
   }
   else if (oy1==0 && oy2==img->Y) for (z=oz1;z<oz2;z++) {
      if (ox2-ox1>12) {
	 ox1++;
	 ox2--;
      }
      for (y=0;y<img->Y;y++) {
	 Ncorr=0;
	 corr=0.;
	 for (x=ox1+1;x<ox2-1;x++) if (img->data[z][y][x]>DMIN && img->data[z][y][x]<DMAX) {
	    corr+=img->data[z][y][x];
	    Ncorr++;
	 }
	 if (Ncorr) corr/=Ncorr;
	 for (x=0;x<img->X;x++) if (img->data[z][y][x]>DMIN && img->data[z][y][x]<DMAX) img->data[z][y][x]-=corr;
      }
   }
   else if (ox1==0 && ox2==img->X) for (z=oz1;z<oz2;z++) {
      if (oy2-oy1>12) {
	 oy1++;
	 oy2--;
      }
      for (x=0;x<img->X;x++) {
	 Ncorr=0;
	 corr=0.;
	 for (y=oy1+1;y<oy2-1;y++) if (img->data[z][y][x]>DMIN && img->data[z][y][x]<DMAX) {
	    corr+=img->data[z][y][x];
	    Ncorr++;
	 }
	 if (Ncorr) corr/=Ncorr;
	 for (y=0;y<img->Y;y++) if (img->data[z][y][x]>DMIN && img->data[z][y][x]<DMAX) img->data[z][y][x]-=corr;
      }
   }
   else printf("**Overscan=2 requires full X or Y range in overscan region\n");
}

void zeroim(imtype*img,imtype*zimg) {
   float DMIN,DMAX,ZDMIN,ZDMAX;
   double avg=0.;
   int x,y,z,N=0;

   if (!isimage(img)) return;
   if (img->X!=zimg->X || img->Y!=zimg->Y || img->Z!=zimg->Z || !isimage(zimg)) {
      printf("**Image size mismatch in bias frame");
      return;
   }
   parsecards(&(fits.img),NULL,NULL,NULL,&DMIN,&DMAX,NULL,NULL,NULL,0,1);
   if (img!=&(fits.img)) parsecards(img,NULL,NULL,NULL,&DMIN,&DMAX,NULL,NULL,NULL,0,0);
   parsecards(zimg,NULL,NULL,NULL,&ZDMIN,&ZDMAX,NULL,NULL,NULL,0,1);
   if (zimg!=&(cal.img)) parsecards(zimg,NULL,NULL,NULL,&ZDMIN,&ZDMAX,NULL,NULL,NULL,0,1);
   for (z=0;z<img->Z;z++) for (y=0;y<img->Y;y++) for (x=0;x<img->X;x++) if (zimg->data[z][y][x]>ZDMIN && zimg->data[z][y][x]<ZDMAX) {
      avg+=zimg->data[z][y][x];
      N++;
   }
   if (N && avg>=0) {
      double tmp;
      avg/=N;
      tmp=DMIN-avg;
      for (z=0;z<img->Z;z++) for (y=0;y<img->Y;y++) for (x=0;x<img->X;x++) if (img->data[z][y][x]<=DMIN) img->data[z][y][x]=tmp-1;
      printf("Lowering minimum data value by %f\n",avg);
      DMIN=tmp;
      insertcards(img,-1.e30,-1.e30,-1.e30,DMIN,-1.e30,-1.e30,-1.e30,-1.e30);
   }
   else if (N) {
      double tmp;
      avg/=N;
      tmp=DMAX-avg;
      for (z=0;z<img->Z;z++) for (y=0;y<img->Y;y++) for (x=0;x<img->X;x++) if (img->data[z][y][x]>=DMAX) img->data[z][y][x]=tmp+1;
      printf("Raising maximum data value by %f\n",-avg);
      DMAX=tmp;
      insertcards(img,-1.e30,-1.e30,-1.e30,-1.e30,DMAX,-1.e30,-1.e30,-1.e30);
   }
   for (z=0;z<img->Z;z++) for (y=0;y<img->Y;y++) for (x=0;x<img->X;x++) if (img->data[z][y][x]>DMIN && img->data[z][y][x]<DMAX) {
      if (zimg->data[z][y][x]>ZDMIN && zimg->data[z][y][x]<ZDMAX) img->data[z][y][x]-=zimg->data[z][y][x];
      else img->data[z][y][x]=DMIN-1;
   }
}

void zero(void) {
   int i;
   char fn[161];

   sprintf(fn,"%s.fits",zerobase);
   readfits(fn,&cal,0);
   printf("Reading zero image from %s\n",fn);
   if (cal.Next!=fits.Next) printf("**Wrong number of extensions in %s\n",fn);
   else {
      zeroim(&(fits.img),&(cal.img));
      for (i=0;i<fits.Next;i++) zeroim(fits.ext+i,cal.ext+i);
   }
   freefits(&cal);
}

void flatim(imtype*img,imtype*fimg) {
   float DMIN,DMAX,FDMIN,FDMAX;
   int x,y,z;
   char str[81];

   if (!isimage(img)) return;
   if (img->X!=fimg->X || img->Y!=fimg->Y || img->Z!=fimg->Z || !isimage(fimg)) {
      printf("**Image size mismatch in bias frame");
      return;
   }
   parsecards(&(fits.img),NULL,NULL,NULL,&DMIN,&DMAX,NULL,NULL,NULL,0,1);
   if (img!=&(fits.img)) parsecards(img,NULL,NULL,NULL,&DMIN,&DMAX,NULL,NULL,NULL,0,0);
   parsecards(&(cal.img),NULL,NULL,NULL,&FDMIN,&FDMAX,NULL,NULL,NULL,0,1);
   if (fimg!=&(cal.img)) parsecards(fimg,NULL,NULL,NULL,&FDMIN,&FDMAX,NULL,NULL,NULL,0,0);
   {
      double tmp;
      if (DMIN<0) tmp=DMIN*1.5;
      else tmp=DMIN*0.75;
      for (z=0;z<img->Z;z++) for (y=0;y<img->Y;y++) for (x=0;x<img->X;x++) if (img->data[z][y][x]<=DMIN) img->data[z][y][x]=tmp-1;
      printf("Lowering minimum data value by %f\n",DMIN-tmp);
      DMIN=tmp;
      tmp=DMAX*1.5;
      for (z=0;z<img->Z;z++) for (y=0;y<img->Y;y++) for (x=0;x<img->X;x++) if (img->data[z][y][x]>=DMAX) img->data[z][y][x]=tmp+1;
      printf("Raising maximum data value by %f\n",tmp-DMAX);
      DMAX=tmp;
      insertcards(img,-1.e30,-1.e30,-1.e30,DMIN,DMAX,-1.e30,-1.e30,-1.e30);
   }
   for (z=0;z<img->Z;z++) for (y=0;y<img->Y;y++) for (x=0;x<img->X;x++) if (img->data[z][y][x]>DMIN && img->data[z][y][x]<DMAX) {
      if (fimg->data[z][y][x]>FDMIN && fimg->data[z][y][x]<FDMAX && fimg->data[z][y][x]>0) img->data[z][y][x]/=fimg->data[z][y][x];
      else img->data[z][y][x]=DMIN-1;
   }
   strcpy(str,getcardval(fimg,"FLATGAIN",0));
   if (str[0]) insertcards(img,atof(str),-1.e30,-1.e30,-1.e30,-1.e30,-1.e30,-1.e30,-1.e30);
}

void flat(void) {
   int i;
   char fn[161],filt[81];

   strcpy(filt,getcardval(&(fits.img),filtkw,0));
   for (i=0;i<fits.Next && !filt[0];i++) strcpy(filt,getcardval(fits.ext+i,filtkw,0));
   sprintf(fn,"%s.%s.fits",flatbase,filt);
   readfits(fn,&cal,0);
   printf("Reading flatfield image from %s\n",fn);
   if (cal.Next!=fits.Next) printf("**Wrong number of extensions in %s\n",fn);
   else {
      flatim(&(fits.img),&(cal.img));
      for (i=0;i<fits.Next;i++) flatim(fits.ext+i,cal.ext+i);
   }
   freefits(&cal);
}

void process(void (*func)(imtype*,int)) {
   int i;

   (*func)(&(fits.img),0);
   for (i=0;i<fits.Next;i++) (*func)(fits.ext+i,i+1);
}

void trim(imtype*img,int ext) {
   int tx1,tx2,ty1,ty2,tz1,tz2,z,y,xsize,dx,dy;
   imgtype tmp;
   cardtype str,card;

   if (!isimage(img)) return;
   if (!getsec(img,trimsec,&tx1,&tx2,&ty1,&ty2,&tz1,&tz2)) {
      printf("**No trim section found\n");
      return;
   }
   if (tz2<=tz1 || ty2<=ty1 || tx2<=tx1) {
      printf("**Illegal trim section\n");
      return;
   }
   printf("Trim section: [%d:%d,%d:%d]\n",tx1+1,tx2,ty1+1,ty2);
   tmp=allocimg(tx2-tx1,ty2-ty1,img->Z);
   xsize=sizeof(float)*(tx2-tx1);
   for (z=0;z<img->Z;z++) for (y=ty1;y<ty2;y++) memcpy(tmp[z][y-ty1],img->data[z][y]+tx1,xsize);
   freeimg(img->data,img->X,img->Y,img->Z);
   img->data=tmp;
   img->X=tx2-tx1;
   img->Y=ty2-ty1;
   sprintf(str,"\'[1:%d,1:%d]\'",tx2-tx1,ty2-ty1);
   sprintf(card,"%-8s= %-20s / Trim section                                   ",trimsec,str);
   addcard(img,card);
   dx=tx1;
   dy=ty1;
   if (getsec(img,datasec,&tx1,&tx2,&ty1,&ty2,&tz1,&tz2)) {
      sprintf(str,"\'[%d:%d,%d:%d]\'",tx1+1-dx,tx2-dx,ty1+1-dy,ty2-dy);
      sprintf(card,"%-8s= %-20s / Data section                                   ",datasec,str);
      addcard(img,card);
   }
}

int mergeamp(void) {
   int i,j,k;

   if (isimage(&(fits.img)) || fits.Next<1) {
      printf("**Can only merge amplifiers in multigroup FITS files\n");
      return 0;
   }
   for (i=0;i<fits.Next;i++) if (fits.ext[i].Z>1) {
      printf("**I don\'t know how to merge 3-diminsional images\n");
      return 0;
   }
   if (strcmp(getcardval(&(fits.img),"D_PROC6",0),"T")) {
      printf("**Must trim image before merging\n");
      return 0;
   }
   for (i=0;i<fits.Next-1;i++) for (j=i+1;j<fits.Next;j++) if (isimage(fits.ext+i) && isimage(fits.ext+j) && !strcmp(getcardval(fits.ext+i,"CCDNAME",0),getcardval(fits.ext+j,"CCDNAME",0))) {
      int cx1,cx2,cy1,cy2,cz1,cz2,sx1,sx2,sy1,sy2,sz1,sz2,x,y,dx,dy,dx1=1,dx2=1,dy1=1,dy2=1;
      imgtype tmp;
      cardtype str,card;
      double FH[3][8];
      float fFH[3][8];

      printf("Merging groups %d and %d\n",i+1,j+1);
      parsecards(&(fits.img),FH[1],FH[1]+1,FH[1]+2,fFH[1]+3,fFH[1]+4,FH[1]+5,FH[1]+6,FH[1]+7,0,1);
      parsecards(fits.ext+i,FH[1],FH[1]+1,FH[1]+2,fFH[1]+3,fFH[1]+4,FH[1]+5,FH[1]+6,FH[1]+7,0,0);
      parsecards(&(fits.img),FH[2],FH[2]+1,FH[2]+2,fFH[2]+3,fFH[2]+4,FH[2]+5,FH[2]+6,FH[1]+7,0,1);
      parsecards(fits.ext+j,FH[2],FH[2]+1,FH[2]+2,fFH[2]+3,fFH[2]+4,FH[2]+5,FH[2]+6,FH[1]+7,0,0);
      memcpy(FH[0],FH[1],sizeof(double[8]));
      memcpy(fFH[0],fFH[1],sizeof(float[8]));
      if (FH[0][0]!=FH[2][0]) {
	 printf("**Gain values do not match; using average\n");
	 FH[0][0]=(FH[0][0]+FH[2][0])*0.5;
      }
      if (FH[0][1]!=FH[2][1]) {
	 printf("**Read Noise values do not match; using average\n");
	 FH[0][1]=(FH[0][1]+FH[2][1])*0.5;
      }
      if (FH[0][2]!=FH[2][2]) {
	 printf("**Exposure time values do not match; using average\n");
	 FH[0][2]=(FH[0][2]+FH[2][2])*0.5;
      }
      if (fFH[2][3]<fFH[0][3]) fFH[0][3]=fFH[2][3];
      if (fFH[2][4]>fFH[0][4]) fFH[0][4]=fFH[2][4];
      if (FH[0][5]!=FH[2][5]) {
	 printf("**MJD values do not match; using average\n");
	 FH[0][5]=(FH[0][5]+FH[2][5])*0.5;
      }
      if (FH[0][6]!=FH[2][6]) {
	 printf("**Air mass values do not match; using average\n");
	 FH[0][6]=(FH[0][6]+FH[2][6])*0.5;
      }
      if (FH[0][7]!=FH[2][7]) {
	 printf("**Effective exposure time values do not match; using average\n");
	 FH[0][7]=(FH[0][7]+FH[2][7])*0.5;
      }
      insertcards(fits.ext+i,FH[0][0],FH[0][1],FH[0][2],fFH[0][3],fFH[0][4],FH[0][5],FH[0][6],FH[0][7]);
      if (!getsec(fits.ext+i,ccdsize,&cx1,&cx2,&cy1,&cy2,&cz1,&cz2)) {
	 printf("**No CCDSIZE found\n");
	 return 0;
      }
      if (!getsec(fits.ext+j,ccdsize,&sx1,&sx2,&sy1,&sy2,&sz1,&sz2)) {
	 printf("**No CCDSIZE found\n");
	 return 0;
      }
      if (sx1!=cx1 || sx2!=cx2 || sy1!=cy1 || sy2!=cy2 || sz1!=cz1 || sz2!=cz2) {
	 printf("**CCD sizes do not match\n");
	 return 0;
      }
      if (cz1!=0 || cz2!=1) {
	 printf("**CCD is 3-dimensional\n");
	 return 0;
      }
      if (cx1!=0 || cy1!=0) {
	 printf("**CCD size [%d:%d,%d:%d] unreadable\n",cx1+1,cx2,cy1+1,cy2);
	 return 0;
      }
      printf("  CCD size is %dx%d\n",cx2,cy2);
      if (!getsec(fits.ext+i,ccdsec,&sx1,&sx2,&sy1,&sy2,&sz1,&sz2)) {
	 printf("**CCD section not found\n");
	 return 0;
      }
      if (abs(sx2-sx1)!=fits.ext[i].X || abs(sy2-sy1)!=fits.ext[i].Y) {
	 printf("**CCD section is the wrong size\n");
	 return 0;
      }
      dx=dy=1;
      if (sx1>sx2) dx=-1;
      if (sy1>sy2) dy=-1;
      tmp=allocimg(cx2,cy2,1);
      for (y=0;y<cy2;y++) for (x=0;x<cx2;x++) tmp[0][y][x]=fFH[0][3]-1;
      for (y=0;y<fits.ext[i].Y;y++) for (x=0;x<fits.ext[i].X;x++) if (fits.ext[i].data[0][y][x]>fFH[1][3]) {
	 if (fits.ext[i].data[0][y][x]<fFH[1][4]) tmp[0][sy1+y*dy][sx1+x*dx]=fits.ext[i].data[0][y][x];
	 else tmp[0][sy1+y*dy][sx1+x*dx]=fFH[0][4];
      }
      if (!getsec(fits.ext+j,ccdsec,&sx1,&sx2,&sy1,&sy2,&sz1,&sz2)) {
	 printf("**CCD section not found\n");
	 return 0;
      }
      if (abs(sx2-sx1)!=fits.ext[j].X || abs(sy2-sy1)!=fits.ext[j].Y) {
	 printf("**CCD section is the wrong size\n");
	 return 0;
      }
      dx=dy=1;
      if (sx1>sx2) dx=-1;
      if (sy1>sy2) dy=-1;
      for (y=0;y<fits.ext[j].Y;y++) for (x=0;x<fits.ext[j].X;x++) if (fits.ext[j].data[0][y][x]>fFH[2][3]) {
	 if (fits.ext[j].data[0][y][x]<fFH[2][4]) tmp[0][sy1+y*dy][sx1+x*dx]=fits.ext[j].data[0][y][x];
	 else tmp[0][sy1+y*dy][sx1+x*dx]=fFH[0][4];
      }
      freeimg(fits.ext[i].data,fits.ext[i].X,fits.ext[i].Y,1);
      fits.ext[i].X=cx2;
      fits.ext[i].Y=cy2;
      fits.ext[i].data=tmp;
      if (getsec(fits.ext+i,detsec,&dx1,&dx2,&dy1,&dy2,&sz1,&sz2) && getsec(fits.ext+j,detsec,&sx1,&sx2,&sy1,&sy2,&sz1,&sz2)) {
	 if (dx1>sx1) dx1=sx1;
	 if (dx2<sx2) dx2=sx2;
	 if (dy1>sy1) dy1=sy1;
	 if (dy2<sy2) dy2=sy2;
	 sprintf(str,"\'[%d:%d,%d:%d]\'",dx1+1,dx2,dy1+1,dy2);
	 sprintf(card,"%-8s= %-20s / Detector section                               ",detsec,str);
	 addcard(fits.ext+i,card);
      }
      if (getsec(fits.ext+i,ccdsec,&dx1,&dx2,&dy1,&dy2,&sz1,&sz2) && getsec(fits.ext+j,ccdsec,&sx1,&sx2,&sy1,&sy2,&sz1,&sz2)) {
	 if (dx1>sx1) dx1=sx1;
	 if (dx2<sx2) dx2=sx2;
	 if (dy1>sy1) dy1=sy1;
	 if (dy2<sy2) dy2=sy2;
	 sprintf(str,"\'[%d:%d,%d:%d]\'",dx1+1,dx2,dy1+1,dy2);
	 sprintf(card,"%-8s= %-20s / CCD section                                    ",ccdsec,str);
	 addcard(fits.ext+i,card);
      }
      sprintf(str,"\'[1:%d,1:%d]\'",cx2,cy2);
      sprintf(card,"%-8s= %-20s / Data section                                   ",datasec,str);
      addcard(fits.ext+i,card);
      strcpy(card,getcardval(fits.ext+i,"WAT0_001",0));
      if (card[0] && strcasecmp(card,"system=image")) {
	 printf("Warning: WCS system does not appear to be \"image\".  I don\'t know how to\nconvert it.\n");
      }
      printf("  Deleting group %d\n",j+1);
      freeim(fits.ext+j);
      fits.Next--;
      for (k=j;k<fits.Next;k++) memcpy(fits.ext+k,fits.ext+k+1,sizeof(imtype));
      j--;
   }
   return 1;
}

int main(int argc,char**argv) {
   int i;

   if (argc<3) {
      printf("Usage: %s <type> <files>\n",*argv);
      return 1;
   }
   strcpy(class,argv[1]);
   paramfile("ccdproc.param",&ccdprocparam);
   for (i=2;i<argc;i++) {
      FILE *f;

      f=readfitsh(argv[i],&fits,1);
      strcpy(obstype,getcardval(&(fits.img),typekw,0));
      if (!class[0] || !strcasecmp(class,obstype)) {
	 readbody(f,&fits,1);
	 process(convert);
	 if (SPECIAL && !strcmp(getcardval(&(fits.img),"D_PROC1",0),"")) {
	    special();
	    addcard(&(fits.img),"D_PROC1 =                    T / Special correction applied                     ");
	 }
	 if (MASK && strcmp(obstype,"zero") && strcmp(obstype,"dark") && !strstr(obstype,"flat") && !strcmp(getcardval(&(fits.img),"D_PROC2",0),"")) {
	    process(mask);
	    addcard(&(fits.img),"D_PROC2 =                    T / Bad pixel mask applied                         ");
	 }
	 if (OVERSCAN && !strcmp(getcardval(&(fits.img),"D_PROC3",0),"")) {
	    process(overscan);
	    addcard(&(fits.img),"D_PROC3 =                    T / Overscan correction applied                    ");
	 }
	 if (ZERO && strcmp(obstype,"zero") && !strcmp(getcardval(&(fits.img),"D_PROC4",0),"")) {
	    zero();
	    addcard(&(fits.img),"D_PROC4 =                    T / Zero correction applied                        ");
	 }
	 if (DARK && strcmp(obstype,"zero") && strcmp(obstype,"dark") && !strcmp(getcardval(&(fits.img),"D_PROC8",0),"")) {
	    //zero();
	    addcard(&(fits.img),"D_PROC8 =                    T / Dark correction applied                        ");
	 }
	 if (FLAT && strcmp(obstype,"zero") && strcmp(obstype,"dark") && !strstr(obstype,"flat") && !strcmp(getcardval(&(fits.img),"D_PROC5",0),"")) {
	    flat();
	    addcard(&(fits.img),"D_PROC5 =                    T / Flatfield corretion applied                    ");
	 }
	 if (TRIM && strcmp(obstype,"zero") && strcmp(obstype,"dark") && !strstr(obstype,"flat") && !strcmp(getcardval(&(fits.img),"D_PROC6",0),"")) {
	    process(trim);
	    addcard(&(fits.img),"D_PROC6 =                    T / Trim section applied                           ");
	 }
	 if (MERGEAMP && strcmp(obstype,"zero") && strcmp(obstype,"dark") && !strstr(obstype,"flat") && !strcmp(getcardval(&(fits.img),"D_PROC7",0),"")) {
	    if (mergeamp()) addcard(&(fits.img),"D_PROC7 =                    T / Amplifiers merged                              ");
	 }
	 writefits(argv[i],&fits,1);
	 freefits(&fits);
      }
      else {
	 fclose(f);
	 if (fits.img.Nmax) free(fits.img.cards);
      }
   }
   return 0;
}
