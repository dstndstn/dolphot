#define FITS_C
#include <fits.h>
#include <stdint.h>
#include <unistd.h>
#include <sys/mman.h>
#include <errno.h>

int fitsparam(char*var,char*val) {
   double x;
   char *ptr;
   if (strlen(val)<=8) {
      if (!strcasecmp(var,"gain_kw")) {sprintf(fitsstr[0],"%-8s=",val); return 1;}
      if (!strcasecmp(var,"exptime_kw")) {sprintf(fitsstr[1],"%-8s=",val); return 1;}
      if (!strcasecmp(var,"rdnoise_kw")) {sprintf(fitsstr[2],"%-8s=",val); return 1;}
      if (!strcasecmp(var,"minval_kw")) {sprintf(fitsstr[3],"%-8s=",val); return 1;}
      if (!strcasecmp(var,"maxval_kw")) {sprintf(fitsstr[4],"%-8s=",val); return 1;}
      if (!strcasecmp(var,"mjd_kw")) {sprintf(fitsstr[5],"%-8s=",val); return 1;}
      if (!strcasecmp(var,"airmass_kw")) {sprintf(fitsstr[6],"%-8s=",val); return 1;}
      if (!strcasecmp(var,"exptime0_kw")) {sprintf(fitsstr[7],"%-8s=",val); return 1;}
   }
   x=strtod(val,&ptr);
   if (!*ptr) {
      if (!strcasecmp(var,"gain_val")) {fitsval[0]=x; return 1;}
      if (!strcasecmp(var,"exptime_val")) {fitsval[1]=x; return 1;}
      if (!strcasecmp(var,"rdnoise_val")) {fitsval[2]=x; return 1;}
      if (!strcasecmp(var,"minval_val")) {fitsval[3]=(double)((float)x); return 1;}
      if (!strcasecmp(var,"maxval_val")) {fitsval[4]=(double)((float)x); return 1;}
      if (!strcasecmp(var,"mjd_val")) {fitsval[5]=x; return 1;}
      if (!strcasecmp(var,"airmass_val")) {fitsval[6]=x; return 1;}
      if (!strcasecmp(var,"exptime0_val")) {fitsval[7]=x; return 1;}
   }
   return 0;
}

int32_t fits_round(float x) {
   return (x>=0)?(int)(x+0.5):(int)(x-0.5);
}

void addcard(imtype*img,char str[81]) {
   int i;

   if (!strncmp(str,"COMMENT  ",9)) i=img->Ncards;
   else for (i=0;i<img->Ncards && strncmp(str,img->cards[i],9);i++);
   if (i<img->Ncards) strcpy(img->cards[i],str);
   else {
      if (i>=img->Nmax) {
	 cardtype *tmp;
	 int newNmax = i+14;
	 tmp=(cardtype*)calloc(sizeof(cardtype),newNmax);
	 if (!tmp) merr();
	 if (img->Ncards) {
	    memcpy(tmp,(img->cards),sizeof(cardtype)*(img->Ncards));
	 }
	 if (img->Nmax) {
	    free(img->cards);
	 }
	 img->cards=tmp;
	 img->Nmax = newNmax;
      }
      strcpy(img->cards[i],str);
      img->Ncards=i+1;
   }
   return;
}

void skipimage(FILE *f,imtype*img) {
   int skip,z,y;
   char ch[2880],*buf;

   buf=(char*)calloc(img->X,abs(img->bits)/8);
   if (!buf) merr();
   for (z=0;z<img->Z;z++) for (y=0;y<img->Y;y++) fread(buf,abs(img->bits)/8,img->X,f);
   skip=abs(img->bits)/8*img->X*img->Y*img->Z;
   skip=((skip+2879)/2880)*2880-skip;
   fread(ch,1,skip,f);
   free(buf);
}

char* parsecardval(char*card) {
   static char str[2][81];
   static int nstr=0;
   int i;

   nstr=1-nstr;
   str[nstr][0]=0;
   if (strlen(card)<9) return str[nstr];
   i=9;
   while (card[i]==' ') i++;
   if (card[i]=='\'') {
      i++;
      strcpy(str[nstr],card+i);
      i=0;
      while (str[nstr][i] && str[nstr][i]!='\'') i++;
      while (i>0 && str[nstr][i-1]==32) i--;
      str[nstr][i]=0;
      return str[nstr];
   }
   strcpy(str[nstr],card+i);
   i=0;
   while (str[nstr][i] && str[nstr][i]!='/') i++;
   while (i>0 && str[nstr][i-1]==32) i--;
   str[nstr][i]=0;
   return str[nstr];
}

void parsecards(imtype*img,double*GAIN,double*RN,double*EXP,float*DMIN,float*DMAX,double*EPOCH,double*AIR,double*EXP0,int ver,int usedef) {
   int i;

   if (usedef) {
      if (GAIN!=NULL) *GAIN=fitsval[0];
      if (EXP!=NULL) *EXP=fitsval[1];
      if (RN!=NULL) *RN=fitsval[2];
      if (DMIN!=NULL) *DMIN=fitsval[3];
      if (DMAX!=NULL) *DMAX=fitsval[4];
      if (EPOCH!=NULL) *EPOCH=fitsval[5];
      if (AIR!=NULL) *AIR=fitsval[6];
      if (EXP0!=NULL) *EXP0=fitsval[7];
   }
   for (i=0;i<img->Ncards;i++) {
      if (GAIN!=NULL && !strncmp(img->cards[i],fitsstr[0],9)) *GAIN=atof(parsecardval(img->cards[i]));
      else if (EXP!=NULL && !strncmp(img->cards[i],fitsstr[1],9)) *EXP=atof(parsecardval(img->cards[i]));
      else if (RN!=NULL && !strncmp(img->cards[i],fitsstr[2],9)) *RN=atof(parsecardval(img->cards[i]));
      else if (DMIN!=NULL && !strncmp(img->cards[i],fitsstr[3],9)) *DMIN=atof(parsecardval(img->cards[i]));
      else if (DMAX!=NULL && !strncmp(img->cards[i],fitsstr[4],9)) *DMAX=atof(parsecardval(img->cards[i]));
      else if (EPOCH!=NULL && !strncmp(img->cards[i],fitsstr[5],9)) *EPOCH=atof(parsecardval(img->cards[i]));
      else if (AIR!=NULL && !strncmp(img->cards[i],fitsstr[6],9)) *AIR=atof(parsecardval(img->cards[i]));
      else if (EXP0!=NULL && !strncmp(img->cards[i],fitsstr[7],9)) *EXP0=atof(parsecardval(img->cards[i]));
   }
   if (ver) {
      printf(" ");
      if (GAIN!=NULL) printf(" GAIN=%4.2f",*GAIN);
      if (EXP!=NULL) printf(" EXP=%ds",(int)*EXP);
      if (RN!=NULL) printf(" NOISE=%4.2f",*RN);
      if (DMIN!=NULL) printf(" BAD=%4.2f",*DMIN);
      if (DMAX!=NULL) printf(" SAT=%4.2f",*DMAX);
      printf("\n");
   }
}

char*getcard(int Ncards,cardtype*cards,char*str) {
   char sstr[10];
   static cardtype card;
   int i;

   if (strlen(str)==9 && str[8]=='=') strcpy(sstr,str);
   else sprintf(sstr,"%-8s=",str);
   for (i=0;i<Ncards;i++) if (!strncmp(cards[i],sstr,9)) {
      strcpy(card,cards[i]);
      return card;
   }
   strcpy(card,"                                                                                ");
   return card;
}

char*getcardval(imtype*img,char*val,int err) {
   static char dumstr[1]="";
   char card[81];

   strcpy(card,getcard(img->Ncards,img->cards,val));
   if (!strcmp(card,"                                                                                ")) {
      if (err) printf("Error: FITS card \"%s\" not found\n",val);
      return dumstr;
   }
   return parsecardval(card);
}

int readcards(FILE *f,ftype*fits,imtype*img) {
   int MODE=0,line=0,tmp,EXT=0;
   cardtype str="";
   int ext_set=0;

   img->bits=-32;
   img->bzero=0.;
   img->bscale=1.;
   img->pcount=0;
   strcpy(img->xtension,"IMAGE");
   img->X=img->Y=img->Z=1;
   img->tfields=0;
   img->Ncards=0;
   img->cards=0;
   img->Nmax=0;
   if (fits) fits->Next=-1;
   while (!feof(f) && (!MODE || line%36)) {
      //fgets(str,81,f);
      fread(str,1,80,f); str[80]=0;
      line++;
      if (!MODE) {
	 if (!strncmp(str,"END     ",8)) MODE=1;
	 else if (!strncmp(str,"SIMPLE  =",9)) {
	    if (str[29]!='T' && str[10]!='T') ferr(f,"Illegal format type -- must be simple");
	    if (!fits) ferr(f,"SIMPLE should not be in extensions");
	 }
	 else if (!strncmp(str,"NAXIS   =",9)) {
	    if ((tmp=atoi(str+9))==2) img->Z=1;
	    else if (tmp==1) img->Y=img->Z=1;
	    else if (tmp==0) img->X=img->Y=img->Z=0;
	    else if (tmp!=3) ferr(f,"Illegal number of axes");
	 }
	 else if (!strncmp(str,"NAXIS",5)) {
	    tmp=atoi(str+9);
	    if (str[5]=='1') img->X=tmp;
	    else if (str[5]=='2') img->Y=tmp;
	    else if (str[5]=='3') img->Z=tmp;
	    else ferr(f,"Illegal Axis Number");
	 }
	 else if (!strncmp(str,"BITPIX  =",9)) img->bits=atoi(str+9);
	 else if (!strncmp(str,"BSCALE  =",9)) img->bscale=atof(str+9);
	 else if (!strncmp(str,"BZERO   =",9)) img->bzero=atof(str+9);
	 else if (!strncmp(str,"PCOUNT  =",9)) img->pcount=atoi(str+9);
	 else if (!strncmp(str,"EXTEND  =",9)) {
	    if (fits==NULL) ferr(f,"Extension format error");
	    else if (str[29]=='T') EXT=1;
	 }
	 else if (!strncmp(str,"NEXTEND =",9)) {
	    if (EXT) {
	       if (fits==NULL) ferr(f,"Extension format error");
	       fits->Next=atoi(str+9);
	    }
	 }
	 else if (!strncmp(str,"XTENSION=",9)) {
	    if (fits!=NULL) ferr(f,"Extension format error");
	    for (tmp=9;str[tmp] && str[tmp]!='\'';tmp++);
	    if (!str[tmp]) ferr(f,"Cannot read extension type");
	    tmp++;
	    while (str[tmp]==32) tmp++;
	    strcpy(img->xtension,str+tmp);
	    for (tmp=0;img->xtension[tmp] && img->xtension[tmp]!='\'';tmp++);
	    if (!img->xtension[tmp]) ferr(f,"Cannot read extension type");
	    tmp--;
	    while (tmp>=0 && img->xtension[tmp]==32) tmp--;
	    img->xtension[tmp+1]=0;
	    ext_set=1;
	 }
	 else if (!strncmp(str,"TFIELDS =",9)) {
	    img->tfields=atoi(str+9);
	 }
	 else addcard(img,str);
      }
   }
   // if NEXTEND not set
   if (fits && fits->Next<0) {
      fits->Next=0;
      if (EXT) {
	 long fpos;
	 imtype img;
	 fpos=ftell(f);
	 skipimage(f,&(fits->img));
	 while (readcards(f,NULL,&img)) {
	    skipimage(f,&img);
	    fits->Next++;
	 }
	 fseek(f,fpos,SEEK_SET);
      }
   }
   return ext_set;
}

FILE *readfitsh(char*fn,ftype*fits,int ver) {
   FILE *f;
   static int first=1;

   if ((f=fopen(fn,"rb"))==NULL) {
      printf("Error reading \"%s\"\n",fn);
      exit(1);
   }
   if (first) {
      paramfile("fits.param",&fitsparam);
      first=0;
   }
   fits->img.bintabdata=NULL;
   readcards(f,fits,&(fits->img));
   if (ver) {
      if (fits->img.Z>1) printf("Reading FITS file %s: %dx%dx%d\n",fn,fits->img.X,fits->img.Y,fits->img.Z);
      else if (fits->img.Y>1) printf("Reading FITS file %s: %dx%d\n",fn,fits->img.X,fits->img.Y);
      else if (fits->img.X>0) printf("Reading FITS file %s: %d\n",fn,fits->img.X);
      else printf("Reading FITS file %s\n",fn);
   }
   return f;
}

void readexth(FILE *f,imtype*img,int ver) {
   img->bintabdata=NULL;
   readcards(f,NULL,img);
   if (ver) {
      if (img->Z>1) printf("Reading %s extension: %dx%dx%d\n",img->xtension,img->X,img->Y,img->Z);
      else if (img->Y>1) printf("Reading %s extension: %dx%d\n",img->xtension,img->X,img->Y);
      else if (img->X>0) printf("Reading %s extension: %d\n",img->xtension,img->X);
      else printf("Reading %s extension\n",img->xtension);
   }
   return;
}

/*
 #ifdef USE_MMAP
 int use_mmap = 1;
 #else
 int use_mmap = 0;
 #endif
 */
#include <assert.h>

void readchip(FILE *f,chiptype chip,imtype*img, int use_mmap) {
   int i,j;
   int32_t l;
   int16_t s;
   int8_t c;
   float *fp,*fpmax;

   if (!((use_mmap == 0) || (use_mmap == 1))) {
       assert(0);
   }

   if (sizeof(float)!=4 || sizeof(double)!=8 || sizeof(int32_t)!=4 || sizeof(int16_t)!=2 || sizeof(int8_t)!=1) {
      printf("Incorrect data sizes in fits_lib.c\n");
      exit(-1);
   }
   if (img->bits==-32) {
       printf("readchip: mmap %i\n", use_mmap);
       if (use_mmap && img->bscale == 1. && img->bzero == 0.) {
           off_t nbytes = img->X * img->Y * 4;
           off_t start = ftello(f);
           off_t mapstart, mapsize;
           int delta;
           int ps = getpagesize();
           int gap = start % ps;
           void* map;
           // start must be a multiple of pagesize.
           mapstart = start  - gap;
           mapsize  = nbytes + gap;
           delta = gap;

           map = mmap(0, mapsize,
                      //PROT_READ,
                      //MAP_SHARED,
                      PROT_READ | PROT_WRITE,
                      MAP_PRIVATE,
                      fileno(f), mapstart);
           if (map == MAP_FAILED) {
               fprintf(stderr, "Failed to mmap file: %s\n", strerror(errno));
               exit(-1);
           }
           printf("Mmapped file: offset %i + size %i\n",
                  (int)mapstart, (int)mapsize);
           printf("addr: %p + %i => %p to %p\n", map, delta, map+delta,
                  map+nbytes);

           chip[0] = map + delta;
           // from allocchip
           for (i=1; i<img->Y; i++)
               chip[i] = chip[i-1] + img->X;


           printf("pixel 0: %g\n", chip[0][0]);

           // pretend we read the data
           fseeko(f, nbytes, SEEK_CUR);
           return;

       } else {
           ffread(chip[0],(img->X)*(img->Y),4,f);
           if (img->bscale==1. && img->bzero==0.) return;
           fp=chip[0]; fpmax=fp+(img->X)*(img->Y);
           while (fp<fpmax) (*fp)=(*fp)*img->bscale+img->bzero;
           return;
       }
   }
   for (i=0;i<img->Y;i++) {
      if (img->bits==32) {
	 for (j=0;j<img->X;j++) {
	    ffread(&l,1,4,f);
	    chip[i][j]=l*img->bscale+img->bzero;
	 }
      }
      else if (img->bits==16) {
	 for (j=0;j<img->X;j++) {
	    sfread(&s,1,2,f);
	    chip[i][j]=s*img->bscale+img->bzero;
	 }
      }
      else if (img->bits==8) {
	 for (j=0;j<img->X;j++) {
	    fread(&c,1,1,f);
	    chip[i][j]=c*img->bscale+img->bzero;
	 }
      }
   }
   return;
}

void endianbintable(imtype*img) {
   int i,y;
   for (i=0;i<img->tfields;i++) {
      switch(img->bintabdata[i].type) {
      case 'D': // double
	 for (y=0;y<img->Y;y++) endian8((char*)img->data[0][0]+img->X*y+img->bintabdata[i].offset,img->bintabdata[i].length);
	 break;
      case 'E': // float
      case 'J': // int
	 for (y=0;y<img->Y;y++) endian4((char*)img->data[0][0]+img->X*y+img->bintabdata[i].offset,img->bintabdata[i].length);
	 break;
      case 'I': // short
	 for (y=0;y<img->Y;y++) endian2((char*)img->data[0][0]+img->X*y+img->bintabdata[i].offset,img->bintabdata[i].length);
	 break;
      case 'A': // char
      case 'B': // byte
      default:
	 break;
      }
   }
}

void readbintable(imtype*img) {
   int i;
   int err=0;
   cardtype errstr="";
   img->bintabdata=(bintabledatatype*)calloc(img->tfields,sizeof(bintabledatatype));
   size_t offset=0;
   for (i=0;i<img->tfields && !err;i++) {
      cardtype str;
      char *fptr,*fptr2;
      img->bintabdata[i].offset=offset;
      sprintf(str,"TTYPE%d",i+1);
      strcpy(img->bintabdata[i].name,getcardval(img,str,0));
      if (img->bintabdata[i].name[0]==0) {
	 sprintf(errstr,"card %s not present",str);
	 err=1;
      }
      else {
	 sprintf(str,"TFORM%d",i+1);
	 fptr=fptr2=getcardval(img,str,0);
	 if (*fptr==0) {
	    sprintf(errstr,"card %s not present",str);
	    err=1;
	 }
	 else {
	    if (*fptr>='0' && *fptr<='9') img->bintabdata[i].length=strtol(fptr,&fptr,10);
	    else img->bintabdata[i].length=1;
	    if (*fptr=='A' || *fptr=='a') { // character
	       img->bintabdata[i].type='A';
	       offset+=(size_t)(img->bintabdata[i].length);
	    }
	    else if (*fptr=='B' || *fptr=='b') { // byte
	       img->bintabdata[i].type='B';
	       offset+=(size_t)(img->bintabdata[i].length);
	    }
	    else if (*fptr=='D' || *fptr=='d') { // double
	       img->bintabdata[i].type='D';
	       offset+=(size_t)8*(size_t)(img->bintabdata[i].length);
	    }
	    else if (*fptr=='E' || *fptr=='e') { // float
	       img->bintabdata[i].type='E';
	       offset+=(size_t)4*(size_t)(img->bintabdata[i].length);
	    }
	    else if (*fptr=='I' || *fptr=='i') { // short
	       img->bintabdata[i].type='I';
	       offset+=(size_t)2*(size_t)(img->bintabdata[i].length);
	    }
	    else if (*fptr=='J' || *fptr=='j') { // long
	       img->bintabdata[i].type='J';
	       offset+=(size_t)4*(size_t)(img->bintabdata[i].length);
	    }
	    // l=logical, x=bit, pi(13)=array of 13 i's, c=32-bit complex, m=64-bit complex
	    else {
	       sprintf(errstr,"unknown form %s",fptr2);
	       err=1;
	    }
	 }
      }
   }
   if (!err && offset!=img->X) {
      sprintf(errstr,"X axis %d does not match table length %lu",img->X,offset);
      err=1;
   }
   if (err) {
      printf("Error: unable to parse bintable (%s)\n",errstr);
      free(img->bintabdata);
      img->tfields=0;
   }
   else endianbintable(img);
   /*
   for (i=0;i<img->tfields && !err;i++) {
      int j,k;
      char*ptr;
      printf("%s (%d%c): ",img->bintabdata[i].name,img->bintabdata[i].length,img->bintabdata[i].type);
      for (j=0;j<img->Y;j++) {
	 if (j) printf(", ");
	 ptr=(char*)img->data[0][0] + j*img->X + img->bintabdata[i].offset;
	 for (k=0;k<img->bintabdata[i].length;k++) {
	 switch(img->bintabdata[i].type) {
	 case 'A':
	    printf("%c",ptr[k]);
	    break;
	 case 'D':
	    printf("%g",*((double*)(ptr+k)));
	    break;
	 case 'E':
	    printf("%g",*((float*)(ptr+k)));
	    break;
	 case 'J':
	    printf("%d",*((int32_t*)(ptr+k)));
	    break;
	 }
	 }
      }
      printf("\n");
   }
   */
}

int gettableindx(imtype*img,char*val,int err) {
   int indx;

   for (indx=0;indx<img->tfields && strcmp(img->bintabdata[indx].name,val);indx++);
   if (indx>=img->tfields) {
      if (err) printf("Error: TABLE card \"%s\" not found\n",val);
      return -1;
   }
   return indx;
}

void* gettableindxcheck(imtype*img,char*val,int row,int i,int err,char *type) {
   int indx;
   if (row<0 || row>=img->Y) {
      if (err) printf("Error: TABLE row %d out of bounds\n",row);
      return 0;
   }
   indx=gettableindx(img,val,err);
   if (indx<0) return 0;
   if (i<0 || i>=img->bintabdata[indx].length) {
      if (err) printf("Error: TABLE index %d out of bounds\n",i);
      return 0;
   }
   *type = img->bintabdata[indx].type;
   return (void*) ( (char*)img->data[0][0]+row*img->X+img->bintabdata[indx].offset );
}

int32_t gettablevalint(imtype*img,char*val,int row,int i,int err) {
   void *ptr;
   char type;
   ptr = gettableindxcheck(img,val,row,i,err,&type);
   if (ptr==0) return 0;
   switch(type) {
   case 'B':
      return (int32_t)(*((int8_t*)ptr+i));
      break;
   case 'D':
      return fits_round(*((double*)ptr+i));
      break;
   case 'E':
      return fits_round(*((float*)ptr+i));
      break;
   case 'I':
      return (int32_t)(*((int16_t*)ptr+i));
      break;
   case 'J':
      return *((int32_t*)ptr+i);
      break;
   case 'A':
   default:
      return (int32_t)(*((char*)ptr+i));
      break;
   }
}

float gettablevalfloat(imtype*img,char*val,int row,int i,int err) {
   void *ptr;
   char type;
   ptr = gettableindxcheck(img,val,row,i,err,&type);
   if (ptr==0) return 0;
   switch(type) {
   case 'B':
      return (float)(*((int8_t*)ptr+i));
      break;
   case 'D':
      return (float)(*((double*)ptr+i));
      break;
   case 'E':
      return (*((float*)ptr+i));
      break;
   case 'I':
      return (float)(*((int16_t*)ptr+i));
      break;
   case 'J':
      return (float)(*((int32_t*)ptr+i));
      break;
   case 'A':
   default:
      return (float)(*((char*)ptr+i));
      break;
   }
}

double gettablevaldouble(imtype*img,char*val,int row,int i,int err) {
   void *ptr;
   char type;
   ptr = gettableindxcheck(img,val,row,i,err,&type);
   if (ptr==0) return 0;
   switch(type) {
   case 'B':
      return (double)(*((int8_t*)ptr+i));
      break;
   case 'D':
      return (*((double*)ptr+i));
      break;
   case 'E':
      return (double)(*((float*)ptr+i));
      break;
   case 'I':
      return (double)(*((int16_t*)ptr+i));
      break;
   case 'J':
      return (double)(*((int32_t*)ptr+i));
      break;
   case 'A':
   default:
      return (double)(*((char*)ptr+i));
      break;
   }
}

#define MAXTABSTR 80
char*gettablevalstring(imtype*img,char*val,int row,int err) {
   static char dumstr[1]="";
   int indx,len;
   static char retstr[MAXTABSTR+1];
   if (row<0 || row>=img->Y) {
      if (err) printf("Error: TABLE row %d out of bounds\n",row);
      return dumstr;
   }
   indx=gettableindx(img,val,err);
   if (indx<0) return dumstr;
   if (img->bintabdata[indx].type!='A') {
      if (err) printf("Error: TABLE %s is not a string\n",val);
      return dumstr;
   }
   len = img->bintabdata[indx].length;
   if (img->bintabdata[indx].length>MAXTABSTR) {
      if (err) printf("Error: TABLE %s length %d exceeds maximum; increase MAXTABSTR\n",val,img->bintabdata[indx].length);
      len = MAXTABSTR;
   }
   memcpy(retstr,(char*)img->data[0][0]+row*img->X+img->bintabdata[indx].offset,len);
   while (len>0 && retstr[len-1]==' ') len--;
   retstr[len]=0;
   return retstr;
}
#undef MAXTABSTR

void readimage(FILE *f,imtype*img, int use_mmap) {
   int skip,z;
   char ch[2880];

   printf("readimage: use_mmap %i\n", use_mmap);

   img->bintabdata=NULL;
   if (!strcasecmp(img->xtension,"BINTABLE")) {
      img->data=allocimgchar(img->X,img->Y,img->Z);
      for (z=0;z<img->Z;z++) fread(img->data[z][0],1,img->X*img->Y,f);
      readbintable(img);
   }
   else {
      img->bintabdata=0;
      img->data=allocimg(img->X,img->Y,img->Z);
      for (z=0;z<img->Z;z++) readchip(f,img->data[z],img, use_mmap);
   }
   skip=abs(img->bits)/8*img->X*img->Y*img->Z+img->pcount;
   skip=((skip+2879)/2880)*2880-skip;
   fread(ch,1,skip,f);
}

void readbody(FILE *f,ftype*fits,int ver) {
   int i;

   printf("readbody\n");
   readimage(f,&(fits->img), 0);
   if (fits->Next>0) {
      fits->ext=(imtype*)calloc(sizeof(imtype),fits->Next);
      if (!fits->ext) merr();
      for (i=0;i<fits->Next;i++) {
	 readexth(f,fits->ext+i,ver);
	 readimage(f,&(fits->ext[i]), 0);
      }
   }
   fclose(f);
   return;
}

void readfits(char*fn,ftype*fits,int ver) {
   readbody(readfitsh(fn,fits,ver),fits,ver);
   return;
}

void readfitsinfo(char*fn,ftype*fits,int ver) {
   FILE *f;
   int i;

   f=readfitsh(fn,fits,ver);
   skipimage(f,&(fits->img));
   if (fits->Next>0) {
      fits->ext=(imtype*)calloc(sizeof(imtype),fits->Next);
      if (!fits->ext) merr();
      for (i=0;i<fits->Next;i++) {
	 readexth(f,fits->ext+i,ver);
	 skipimage(f,&(fits->ext[i]));
      }
   }
   fclose(f);
   return;
}

char* fitsnum(double x) {
   static char str[401];
   char *ptr;
   sprintf(str,"%20f",x);
   if (strlen(str)<=20) return str;
   sprintf(str,"%20.13e",x);
   if (strlen(str)>20) sprintf(str,"%20.12e",x);
   if (strlen(str)>20) sprintf(str,"%20.10e",x);
   ptr=strchr(str,'e');
   if (ptr!=NULL) *ptr='E';
   return str;
}

void insertcards(imtype*img,double GAIN,double RN,double EXP,float DMIN,float DMAX,double EPOCH,double AIR,double EXP0) {
   cardtype str;

   if (GAIN>-0.9e30) {
      sprintf(str,"%-9s %20s / Gain (electrons/ADU)                           ",fitsstr[0],fitsnum(GAIN));
      addcard(img,str);
   }
   if (EXP>-0.9e30) {
      sprintf(str,"%-9s %20s / Exposure time (s)                              ",fitsstr[1],fitsnum(EXP));
      addcard(img,str);
   }
   if (RN>-0.9e30) {
      sprintf(str,"%-9s %20s / Read noise (electrons)                         ",fitsstr[2],fitsnum(RN));
      addcard(img,str);
   }
   if (DMIN>-0.9e30) {
      sprintf(str,"%-9s %20s / Low bad data value (DN)                        ",fitsstr[3],fitsnum(DMIN));
      addcard(img,str);
   }
   if (DMAX>-0.9e30) {
      sprintf(str,"%-9s %20s / High bad data value (DN)                       ",fitsstr[4],fitsnum(DMAX));
      addcard(img,str);
   }
   if (EPOCH>-0.9e30) {
      sprintf(str,"%-9s %20s / Epoch (MJD)                                    ",fitsstr[5],fitsnum(EPOCH));
      addcard(img,str);
   }
   if (AIR>-0.9e30) {
      sprintf(str,"%-9s %20s / Air mass                                       ",fitsstr[6],fitsnum(AIR));
      addcard(img,str);
   }
   if (EXP0>-0.9e30) {
      sprintf(str,"%-9s %20s / Single image exposure time (s)                 ",fitsstr[7],fitsnum(EXP0));
      addcard(img,str);
   }
}

void writecards(FILE *f,ftype*fits,imtype*img) {
   char str[81];
   int line=1,i;

   if (fits) {
      fputs("SIMPLE  =                    T / Standard FITS format                           ",f); line++;
   }
   else {
      sprintf(str,"XTENSION= '%-18s' / Extension type                                 ",img->xtension); fputs(str,f); line++;
   }
   sprintf(str,"BITPIX  = %20d / Data type                                      ",img->bits); fputs(str,f); line++;
   if (img->Z>1) {
      fputs("NAXIS   =                    3 / Number of axes                                 ",f); line++;
      sprintf(str,"NAXIS1  = %20d /                                                ",img->X); fputs(str,f); line++;
      sprintf(str,"NAXIS2  = %20d /                                                ",img->Y); fputs(str,f); line++;
      sprintf(str,"NAXIS3  = %20d /                                                ",img->Z); fputs(str,f); line++;
   }
   else if (img->Y>1) {
      fputs("NAXIS   =                    2 / Number of axes                                 ",f); line++;
      sprintf(str,"NAXIS1  = %20d /                                                ",img->X); fputs(str,f); line++;
      sprintf(str,"NAXIS2  = %20d /                                                ",img->Y); fputs(str,f); line++;
   }
   else if (img->X>0) {
      fputs("NAXIS   =                    1 / Number of axes                                 ",f); line++;
      sprintf(str,"NAXIS1  = %20d /                                                ",img->X); fputs(str,f); line++;
   }
   else {
      fputs("NAXIS   =                    0 / Number of axes                                 ",f); line++;
   }
   if (fits) {
      if (fits->Next<=0) {fputs("EXTEND  =                    F / There are no standard extensions               ",f); line++;}
      else {
	 fputs("EXTEND  =                    T / There may be standard extensions               ",f); line++;
	 sprintf(str,"NEXTEND =       %14d /                                                ",fits->Next); fputs(str,f); line++;
      }
   }
   sprintf(str,"BSCALE  = %20s / Image scale                                    ",fitsnum(img->bscale)); fputs(str,f); line++;
   sprintf(str,"BZERO   = %20s / Image zero                                     ",fitsnum(img->bzero)); fputs(str,f); line++;
   for (i=0;i<img->Ncards;i++) {
      sprintf(str,"%-80s",img->cards[i]);
      fputs(str,f);
      line++;
   }
   for (;line%36;line++) fputs("                                                                                ",f);
   fputs("END                                                                             ",f);
}

FILE *writefitsh(char*fn,ftype*fits,int ver) {
   FILE *f;

   if ((f=fopen(fn,"wb"))==NULL) {
      printf("Error writing to \"%s\"\n",fn);
      exit(1);
   }
   writecards(f,fits,&(fits->img));
   if (ver) {
      if (fits->img.Z>1) printf("Writing FITS file %s: %dx%dx%d\n",fn,fits->img.X,fits->img.Y,fits->img.Z);
      else if (fits->img.Y>1) printf("Writing FITS file %s: %dx%d\n",fn,fits->img.X,fits->img.Y);
      else if (fits->img.X>0) printf("Writing FITS file %s: %d\n",fn,fits->img.X);
      else printf("Writing FITS file %s\n",fn);
   }
   return f;
}

void writeexth(FILE *f,imtype*img,int ver) {
   writecards(f,NULL,img);
   if (ver) {
      if (img->Z>1) printf("Writing %s extension: %dx%dx%d\n",img->xtension,img->X,img->Y,img->Z);
      else if (img->Y>1) printf("Writing %s extension: %dx%d\n",img->xtension,img->X,img->Y);
      else if (img->X>0) printf("Writing %s extension: %d\n",img->xtension,img->X);
      else printf("Writing %s extension\n",img->xtension);
   }
   return;
}

void writechip(FILE *f,chiptype chip,imtype*img) {
   int i,j;
   int32_t l;
   int16_t s;
   int8_t c;
   float fl;

   if (sizeof(float)!=4 || sizeof(int32_t)!=4 || sizeof(int16_t)!=2 || sizeof(int8_t)!=1) {
      printf("Incorrect data sizes in fits_lib.c\n");
      exit(-1);
   }
   if (img->bits==-32 && img->bscale==1. && img->bzero==0.) {
      ffwrite(chip[0],(img->X)*(img->Y),4,f);
      return;
   }
   for (i=0;i<img->Y;i++) {
      if (img->bits==-32) {
	 for (j=0;j<img->X;j++) {
	    fl=(chip[i][j]-img->bzero)/img->bscale;
	    ffwrite(&fl,1,4,f);
	 }
      }
      else if (img->bits==32) {
	 for (j=0;j<img->X;j++) {
	    l=(int32_t)fits_round((chip[i][j]-img->bzero)/img->bscale);
	    ffwrite(&l,1,4,f);
	 }
      }
      else if (img->bits==16) {
	 for (j=0;j<img->X;j++) {
	    s=(int16_t)fits_round((chip[i][j]-img->bzero)/img->bscale);
	    sfwrite(&s,1,2,f);
	 }
      }
      else if (img->bits==8) {
	 for (j=0;j<img->X;j++) {
	    c=(int8_t)fits_round((chip[i][j]-img->bzero)/img->bscale);
	    fwrite(&c,1,1,f);
	 }
      }
   }
   fflush(f);
}

void writeimage(FILE *f,imtype*img) {
   int skip,z;
   char ch[2880];

   if (img->bintabdata) endianbintable(img);
   for (z=0;z<img->Z;z++) writechip(f,img->data[z],img);
   skip=abs(img->bits)/8*img->X*img->Y*img->Z;
   skip=((skip+2879)/2880)*2880-skip;
   memset(ch,0,skip);
   fwrite(ch,1,skip,f);
   fflush(f);
}

void writefits(char*fn,ftype*fits,int ver) {
   FILE *f;
   int i;

   f=writefitsh(fn,fits,ver);
   writeimage(f,&(fits->img));
   for (i=0;i<fits->Next;i++) {
      writeexth(f,fits->ext+i,ver);
      writeimage(f,fits->ext+i);
   }
   fclose(f);
   return;
}

void freeim(imtype*img) {
   if (img->Nmax) free(img->cards);
   img->Nmax=0;
   freeimg(img->data,img->X,img->Y,img->Z);
   if (img->bintabdata) free(img->bintabdata);
}

void imcopy(imtype*oldimg,imtype*newimg) {
   int y,z;
   memcpy(newimg,oldimg,sizeof(imtype));
   if (newimg->Nmax) {
      newimg->cards=(cardtype*)calloc(sizeof(cardtype),newimg->Nmax);
      memcpy(newimg->cards,oldimg->cards,sizeof(cardtype)*(newimg->Nmax));
   }
   newimg->data=allocimg(newimg->X,newimg->Y,newimg->Z);
   for (z=0;z<newimg->Z;z++) for (y=0;y<newimg->Y;y++) memcpy(newimg->data[z][y],oldimg->data[z][y],sizeof(float)*(newimg->X));
   if (oldimg->bintabdata) {
      newimg->bintabdata=(bintabledatatype*)calloc(newimg->tfields,sizeof(bintabledatatype));
      memcpy(newimg->bintabdata,oldimg->bintabdata,newimg->tfields*sizeof(bintabledatatype));
   }
}

void fitscopy(ftype*oldfits,ftype*newfits) {
   int i;

   newfits->Next=oldfits->Next;
   imcopy(&(oldfits->img),&(newfits->img));
   if (newfits->Next) {
      newfits->ext=(imtype*)calloc(sizeof(imtype),newfits->Next);
      if (!newfits->ext) merr();
      for (i=0;i<newfits->Next;i++) imcopy(oldfits->ext+i,newfits->ext+i);
   }
}

void freefits(ftype*fits) {
   int i;
   freeim(&(fits->img));
   for (i=0;i<fits->Next;i++) freeim(fits->ext+i);
   if (fits->Next) free(fits->ext);
}

void freeinfo(ftype*fits) {
   int i;
   if (fits->img.Nmax) free(fits->img.cards);
   fits->img.Nmax=0;
   for (i=0;i<fits->Next;i++) if (fits->ext[i].Nmax) free(fits->ext[i].cards);
   if (fits->Next) free(fits->ext);
}

int isimage(imtype*img) {
   if (img->bintabdata) return 0;
   if (img->Z>1 || img->Y>1) return 1;
   if (strcasecmp(img->xtension,"image")) return 1;
   return 0;
}
