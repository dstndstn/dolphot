#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <dolphot.h>
#include "wfc3filters.h"
#include "wfc3distort.h"

void usage(char*exe) {
   printf("Usage: %s <CM> <filt> <X> <Y> <fwd/rev>\n",exe);
   printf("  CM: 0=IR, 1/2=UVIS\n");
   printf("  fwd = forward (raw image to distortion corrected\n");
   printf("  rev = reverse (distortion corrected to raw image\n");
   exit(-1);
}

int main(int argc,char**argv) {
   int cm,filt;
   double X,Y;

   if (argc!=6) usage(*argv);
   WFC3initfilters();
   cm=atoi(argv[1]);
   filt=WFC3findfilt(argv[2]);
   X=atof(argv[3]);
   Y=atof(argv[4]);
   if (!strcasecmp(argv[5],"fwd")) WFC3fwddistort(cm,filt,&X,&Y);
   else if (!strcasecmp(argv[5],"rev")) WFC3revdistort(cm,filt,&X,&Y);
   else usage(*argv);
   printf("%f %f\n",X,Y);
   return 0;
}
