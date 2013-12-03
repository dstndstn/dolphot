#include <stdio.h>
#include <math.h>
#include "wfc3filters.h"

void WFC3fwddistort(int cm,int filt,double*x,double*y) {
   double x0,y0;
   x0=*x-WFC3filters[filt].idc[0][cm][0];
   y0=*y-WFC3filters[filt].idc[0][cm][1];
   *x=(WFC3filters[filt].idc[0][cm][6]*y0+WFC3filters[filt].idc[0][cm][7]*x0+WFC3filters[filt].idc[0][cm][8]*y0*y0+WFC3filters[filt].idc[0][cm][9]*x0*y0+WFC3filters[filt].idc[0][cm][10]*x0*x0+WFC3filters[filt].idc[0][cm][11]*y0*y0*y0+WFC3filters[filt].idc[0][cm][12]*x0*y0*y0+WFC3filters[filt].idc[0][cm][13]*x0*x0*y0+WFC3filters[filt].idc[0][cm][14]*x0*x0*x0+WFC3filters[filt].idc[0][cm][15]*y0*y0*y0*y0+WFC3filters[filt].idc[0][cm][16]*x0*y0*y0*y0+WFC3filters[filt].idc[0][cm][17]*x0*x0*y0*y0+WFC3filters[filt].idc[0][cm][18]*x0*x0*y0*x0+WFC3filters[filt].idc[0][cm][19]*x0*x0*x0*x0+WFC3filters[filt].idc[0][cm][5])/WFC3filters[filt].idc[0][cm][3];
   *y=(WFC3filters[filt].idc[0][cm][20]*y0+WFC3filters[filt].idc[0][cm][21]*x0+WFC3filters[filt].idc[0][cm][22]*y0*y0+WFC3filters[filt].idc[0][cm][23]*x0*y0+WFC3filters[filt].idc[0][cm][24]*x0*x0+WFC3filters[filt].idc[0][cm][25]*y0*y0*y0+WFC3filters[filt].idc[0][cm][26]*x0*y0*y0+WFC3filters[filt].idc[0][cm][27]*x0*x0*y0+WFC3filters[filt].idc[0][cm][28]*x0*x0*x0+WFC3filters[filt].idc[0][cm][29]*y0*y0*y0*y0+WFC3filters[filt].idc[0][cm][30]*x0*y0*y0*y0+WFC3filters[filt].idc[0][cm][31]*x0*x0*y0*y0+WFC3filters[filt].idc[0][cm][32]*x0*x0*y0*x0+WFC3filters[filt].idc[0][cm][33]*x0*x0*x0*x0+WFC3filters[filt].idc[0][cm][4])/WFC3filters[filt].idc[0][cm][3];
   return;
}

void WFC3revdistort(int cm,int filt,double*x,double*y) {
   double x0,y0;
   x0=(*x)*WFC3filters[filt].idc[0][cm][3]-WFC3filters[filt].idc[0][cm][5];
   y0=(*y)*WFC3filters[filt].idc[0][cm][3]-WFC3filters[filt].idc[0][cm][4];
   *x=WFC3filters[filt].idc[0][cm][0]+WFC3filters[filt].idc[1][cm][6]*y0+WFC3filters[filt].idc[1][cm][7]*x0+WFC3filters[filt].idc[1][cm][8]*y0*y0+WFC3filters[filt].idc[1][cm][9]*x0*y0+WFC3filters[filt].idc[1][cm][10]*x0*x0+WFC3filters[filt].idc[1][cm][11]*y0*y0*y0+WFC3filters[filt].idc[1][cm][12]*x0*y0*y0+WFC3filters[filt].idc[1][cm][13]*x0*x0*y0+WFC3filters[filt].idc[1][cm][14]*x0*x0*x0+WFC3filters[filt].idc[1][cm][15]*y0*y0*y0*y0+WFC3filters[filt].idc[1][cm][16]*x0*y0*y0*y0+WFC3filters[filt].idc[1][cm][17]*x0*x0*y0*y0+WFC3filters[filt].idc[1][cm][18]*x0*x0*y0*x0+WFC3filters[filt].idc[1][cm][19]*x0*x0*x0*x0;
   *y=WFC3filters[filt].idc[0][cm][1]+WFC3filters[filt].idc[1][cm][20]*y0+WFC3filters[filt].idc[1][cm][21]*x0+WFC3filters[filt].idc[1][cm][22]*y0*y0+WFC3filters[filt].idc[1][cm][23]*x0*y0+WFC3filters[filt].idc[1][cm][24]*x0*x0+WFC3filters[filt].idc[1][cm][25]*y0*y0*y0+WFC3filters[filt].idc[1][cm][26]*x0*y0*y0+WFC3filters[filt].idc[1][cm][27]*x0*x0*y0+WFC3filters[filt].idc[1][cm][28]*x0*x0*x0+WFC3filters[filt].idc[1][cm][29]*y0*y0*y0*y0+WFC3filters[filt].idc[1][cm][30]*x0*y0*y0*y0+WFC3filters[filt].idc[1][cm][31]*x0*x0*y0*y0+WFC3filters[filt].idc[1][cm][32]*x0*x0*y0*x0+WFC3filters[filt].idc[1][cm][33]*x0*x0*x0*x0;
   return;
}

// returns in standard pixel size
double WFC3pixsize(int cm,int filt,double x,double y) {
   double x1=x-0.5; double y1=y-0.5;
   double x2=x+0.5; double y2=y+0.5;
   double xa=x-0.5; double ya=y+0.5;
   double xb=x+0.5; double yb=y-0.5;
   WFC3fwddistort(cm,filt,&x1,&y1);
   WFC3fwddistort(cm,filt,&x2,&y2);
   WFC3fwddistort(cm,filt,&xa,&ya);
   WFC3fwddistort(cm,filt,&xb,&yb);
   y2-=y1;
   x2-=x1;
   return 0.5*fabs((xa-x1)*y2-(ya-y1)*x2)+0.5*fabs((xb-x1)*y2-(yb-y1)*x2);
}
