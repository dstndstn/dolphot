#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <dolphot.h>
#include "wfc3filters.h"
#include "wfc3distort.h"

int main(int argc,char**argv) {
   char*ptr,*ptr2,str[8001];
   int cm,filt;
   double x,y;

   if (argc!=1) {
      printf("Usage: %s\n",*argv);
      printf("Data given via stdin, in column format:\n");
      printf("  camera filter X Y <data>\n");
      printf("  camera: 0=IR, 1=UVIS1, 2=UVIS2\n");
      printf("  filter: F555W, etc\n");
      return -1;
   }
   WFC3initfilters();
   while (fgets(str,8001,stdin)) {
      cm=strtol(str,&ptr,10);
      if (ptr==str) {
	 printf("No camera found:\n%s",str);
	 return -1;
      }
      if (cm<0 || cm>2) {
	 printf("Illegal camera (0, 1, 2 OK):\n%s",str);
	 return -1;
      }
      while (*ptr==' ') ptr++;
      for (ptr2=ptr;*ptr2 && *ptr2!=' ';ptr2++);
      if (!*ptr2) {
	 printf("No filter found:\n%s",str);
	 return -1;
      }
      *ptr2=0;
      filt=WFC3findfilt(ptr);
      *ptr2=' ';
      ptr=ptr2;
      x=strtod(ptr,&ptr2);
      if (ptr2==ptr) {
	 printf("No X value found:\n%s",str);
	 return -1;
      }
      if (x<0 || (cm==0 && x>1014) || (cm>0 && x>4096)) {
	 printf("Illegal X value (0-1014 IR, 0-4096 UVIS):\n%s",str);
	 return -1;
      }
      ptr=ptr2;
      y=strtod(ptr,&ptr2);
      if (ptr2==ptr) {
	 printf("No Y value found:\n%s",str);
	 return -1;
      }
      if (y<0 || (cm==0 && y>1014) || (cm>0 && y>2051)) {
	 printf("Illegal Y value (0-1014 IR, 0-2051 UVIS):\n%s",str);
	 return -1;
      }
      WFC3fwddistort(cm,filt,&x,&y);
      if (cm==0) {
	 x=-x*WFC3filters[filt].idc[0][cm][4]+WFC3filters[filt].idc[0][cm][2];
	 y=y*WFC3filters[filt].idc[0][cm][4]+WFC3filters[filt].idc[0][cm][3];
      }
      else {
	 x=x*WFC3filters[filt].idc[0][cm][4]+WFC3filters[filt].idc[0][cm][2];
	 y=-y*WFC3filters[filt].idc[0][cm][4]+WFC3filters[filt].idc[0][cm][3];
      }
      printf("%8.3f %8.3f\n",x,y);
   }
   return 0;
}
