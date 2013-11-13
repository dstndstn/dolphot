#define DOLPHOT_C
#include <dolphot.h>
char camera[81]="";

int ferr(FILE *f,char *str) {
   printf("%s\n",str);
   if (f!=NULL) fclose(f);
   exit(1);
}

void merr(void) {
   printf("Memory allocation error\n");
   exit(1);
}

chiptype allocchip(int X,int Y) {
   chiptype tmp;
   int i;

   if (X<1) X=1;
   if (Y<1) Y=1;
   tmp=(float**)calloc(Y,sizeof(float*));
   if (!tmp) merr();
   tmp[0]=(float*)calloc((long)Y*(long)X,sizeof(float));
   if (!tmp[0]) merr();
   for (i=1;i<Y;i++) tmp[i]=tmp[i-1]+X;
   return tmp;
}

chiptype allocchipchar(int X,int Y) {
   chiptype tmp;
   int i;

   if (X<1) X=1;
   if (Y<1) Y=1;
   tmp=(float**)calloc(Y,sizeof(char*));
   if (!tmp) merr();
   tmp[0]=(float*)calloc((long)Y*(long)X,sizeof(char));
   if (!tmp[0]) merr();
   for (i=1;i<Y;i++) tmp[i]=tmp[i-1]+X;
   return tmp;
}

void freechip(chiptype ptr,int X,int Y) {
   free(ptr[0]);
   free(ptr);
   return;
}

imgtype allocimg(int X,int Y,int Z) {
   imgtype tmp;
   int i;

   if (Z<1) Z=1;
   tmp=(float***)calloc(Z,sizeof(float**));
   if (!tmp) merr();
   for (i=0;i<Z;i++) tmp[i]=allocchip(X,Y);
   return tmp;
}

imgtype allocimgchar(int X,int Y,int Z) {
   imgtype tmp;
   int i;

   if (Z<1) Z=1;
   tmp=(float***)calloc(Z,sizeof(char**));
   if (!tmp) merr();
   for (i=0;i<Z;i++) tmp[i]=allocchipchar(X,Y);
   return tmp;
}

void freeimg(imgtype tmp,int X,int Y,int Z) {
   int i;

   if (Z<1) Z=1;
   for (i=0;i<Z;i++) freechip(tmp[i],X,Y);
   free(tmp);
}

// endian-testing code
short int ENDIAN_TEST_WORD=1;
char *IS_LITTLE_ENDIAN=(char*)(&ENDIAN_TEST_WORD);

void endian2(char*cp,int N) {
   if (*IS_LITTLE_ENDIAN==0) return;
   int i;
   char x;
   for (i=0;i<N;i++) {
      x=cp[0]; cp[0]=cp[1]; cp[1]=x;
      cp+=2;
   }
}

void endian4(char*cp,int N) {
   if (*IS_LITTLE_ENDIAN==0) return;
   int i;
   char x;
   for (i=0;i<N;i++) {
      x=cp[0]; cp[0]=cp[3]; cp[3]=x;
      x=cp[1]; cp[1]=cp[2]; cp[2]=x;
      cp+=4;
   }
}

void endian8(char*cp,int N) {
   if (*IS_LITTLE_ENDIAN==0) return;
   int i;
   char x;
   for (i=0;i<N;i++) {
      x=cp[0]; cp[0]=cp[7]; cp[7]=x;
      x=cp[1]; cp[1]=cp[6]; cp[6]=x;
      x=cp[2]; cp[2]=cp[5]; cp[5]=x;
      x=cp[3]; cp[3]=cp[4]; cp[4]=x;
      cp+=8;
   }
}

size_t ffread(void *ptr,size_t size,size_t N,FILE *f) {
   size_t rv;
   rv=fread(ptr,size,N,f);
   endian4((char*)ptr,(N*size)/4);
   return rv;
}

size_t ffwrite(const void *ptr,size_t size,size_t N,FILE *f) {
   if (*IS_LITTLE_ENDIAN==0) {
      return fwrite(ptr,size,N,f);
   }
   char *cp;
   long tot;
   size_t rv;
   tot=N*size;
   cp=(char*)calloc(tot,1);
   if (!cp) {
      printf("Memory error in ffwrite\n");
      exit(0);
   }
   memcpy(cp,ptr,tot);
   endian4(cp,tot/4);
   rv=fwrite(cp,size,N,f);
   free(cp);
   return rv;
}

size_t sfread(void *ptr,size_t size,size_t N,FILE *f) {
   size_t rv;
   rv=fread(ptr,size,N,f);
   endian2((char*)ptr,(N*size)/2);
   return rv;
}

size_t sfwrite(const void *ptr,size_t size,size_t N,FILE *f) {
   if (*IS_LITTLE_ENDIAN==0) {
      return fwrite(ptr,size,N,f);
   }
   char *cp;
   long tot;
   size_t rv;
   tot=N*size;
   cp=(char*)calloc(tot,1);
   if (!cp) {
      printf("Memory error in sfwrite\n");
      exit(0);
   }
   memcpy(cp,ptr,tot);
   endian2(cp,tot/2);
   rv=fwrite(cp,size,N,f);
   free(cp);
   return rv;
}

int parseparam(char*str,int(*func)(char*,char*)) {
   char var[801],val[801],*ptr;

   for (ptr=str;*ptr==32 || *ptr==9;ptr++);
   if (*ptr=='#' || *ptr==0) return 0;
   strcpy(var,ptr);
   if ((ptr=strchr(var,'='))==NULL) {
      printf("Cannot parse parameter line:\n%s\n",str);
      return -2;
   }
   *ptr=0;
   ptr++;
   for (;*ptr==32 || *ptr==9;ptr++);
   strcpy(val,ptr);
   if ((ptr=strchr(val,'#'))!=NULL) *ptr=0;
   for (ptr=var+strlen(var);ptr>var && (*(ptr-1)==32 || *(ptr-1)==9);ptr--);
   *ptr=0;
   for (ptr=val+strlen(val);ptr>val && (*(ptr-1)==32 || *(ptr-1)==9 || *(ptr-1)==10 || *(ptr-1)==13);ptr--);
   *ptr=0;
   if ((*func)(var,val)) return 0;
   printf("Cannot match parameter name:\n%s\n",var);
   return -1;
}

void paramfile(char*fn,int(*func)(char*,char*)) {
   char str[801];
   FILE *f;
   int any=0;

   sprintf(str,"%s/param/%s",BASEDIR,fn);
   if ((f=fopen(str,"r"))!=NULL) {
      any=1;
      while (fgets(str,801,f)) parseparam(str,func);
      fclose(f);
   }
   if (camera[0]) {
      sprintf(str,"%s/%s/%s",BASEDIR,camera,fn);
      if ((f=fopen(str,"r"))!=NULL) {
	 any=1;
	 while (fgets(str,801,f)) parseparam(str,func);
	 fclose(f);
      }
   }
   if ((f=fopen(fn,"r"))!=NULL) {
      any=1;
      while (fgets(str,801,f)) parseparam(str,func);
      fclose(f);
   }
   if (!any) printf("No parameter file \"%s\" found; using defaults\n",fn);
   return;
}

void paramfile1(char*fn,int(*func)(char*,char*)) {
   char str[801];
   FILE *f;

   if ((f=fopen(fn,"r"))!=NULL) {
      while (fgets(str,801,f)) parseparam(str,func);
      fclose(f);
   }
   else printf("No parameter file \"%s\" found; using defaults\n",fn);
   return;
}
