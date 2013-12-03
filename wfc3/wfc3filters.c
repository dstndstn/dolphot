#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <assert.h>
#include <math.h>

typedef struct {
   char name[11],color;
   double zp[3];
   char xorder[8];
   int xformc[8];
   double xform[8][5];
   double idc[2][3][34];
} WFC3filttype;

int WFC3_NFILTERS=-1;
WFC3filttype *WFC3filters;

/*
// New CTE Parameters (cold only);
//background "contamination";
double bgcorr[WFC3_NFILTERS][2]={
   {0.000000e-06,2.869069e-04}, //F122W
   {0.000000e-06,2.869069e-04}, //F160BW
   {0.000000e-06,2.869069e-04}, //F170W
   {0.000000e-06,2.869069e-04}, //F185W
   {0.000000e-06,2.869069e-04}, //F218W
   {0.000000e-06,2.869069e-04}, //F255W
   {8.540038e-07,1.010904e-04}, //F300W
   {8.540038e-07,1.171397e-04}, //F336W
   {8.540038e-07,1.171397e-04}, //F343N
   {8.540038e-07,9.556493e-05}, //F375N
   {8.540038e-07,9.556493e-05}, //F380W
   {8.540038e-07,9.556493e-05}, //F390N
   {1.085865e-05,3.178520e-05}, //F410M
   {1.085865e-05,3.178520e-05}, //F437N
   {1.085865e-05,3.178520e-05}, //F439W
   {1.085865e-05,3.178520e-05}, //F450W
   {1.085865e-05,3.178520e-05}, //F467M
   {1.085865e-05,3.178520e-05}, //F469N
   {1.085865e-05,3.178520e-05}, //F487N
   {1.077151e-05,3.013676e-05}, //F502N
   {1.077151e-05,3.013676e-05}, //F547M
   {1.077151e-05,3.013676e-05}, //F555W
   {1.077151e-05,3.013676e-05}, //F569W
   {1.077151e-05,3.013676e-05}, //F588N
   {1.077151e-05,3.013676e-05}, //F606W
   {7.776181e-06,2.758164e-05}, //F622W
   {7.776181e-06,2.758164e-05}, //F631N
   {7.776181e-06,2.758164e-05}, //F656N
   {7.776181e-06,2.758164e-05}, //F658N
   {7.776181e-06,2.758164e-05}, //F673N
   {7.776181e-06,2.758164e-05}, //F675W
   {7.776181e-06,2.758164e-05}, //F702W
   {7.917809e-06,2.726974e-05}, //F785LP
   {7.917809e-06,2.726974e-05}, //F791W
   {7.917809e-06,2.726974e-05}, //F814W
   {7.917809e-06,2.726974e-05}, //F850LP
   {7.917809e-06,2.726974e-05}, //F953N
   {7.917809e-06,2.726974e-05}  //F1042M
};

//zero point corrections;
double fcorr[WFC3_NFILTERS][2][4]={
   {{0.000000,0.000000,0.000000,0.000000},{0.000000,0.000000,0.000000,0.000000}}, //F122M
   {{0.000000,0.000000,0.000000,0.000000},{0.000000,0.000000,0.000000,0.000000}}, //F160BW
   {{0.000000,0.000000,0.000000,0.000000},{0.000000,0.000000,0.000000,0.000000}}, //F170W
   {{0.000000,0.000000,0.000000,0.000000},{0.000000,0.000000,0.000000,0.000000}}, //F185W
   {{0.000000,0.000000,0.000000,0.000000},{0.000000,0.000000,0.000000,0.000000}}, //F218W
   {{0.000000,0.000000,0.000000,0.000000},{0.000000,0.000000,0.000000,0.000000}}, //F255W
   {{0.000000,0.000000,0.000000,0.000000},{0.000000,0.000000,0.000000,0.000000}}, //F300W
   {{0.000000,0.000000,0.000000,0.000000},{0.000000,0.000000,0.000000,0.000000}}, //F336W
   {{0.000000,0.000000,0.000000,0.000000},{0.000000,0.000000,0.000000,0.000000}}, //F343N
   {{0.000000,0.000000,0.000000,0.000000},{0.000000,0.000000,0.000000,0.000000}}, //F375N
   {{-0.094952,-0.114463,-0.134266,-0.141137},{-0.091563,-0.124748,-0.126559,-0.129744}}, //F380W
   {{0.000000,0.000000,0.000000,0.000000},{0.000000,0.000000,0.000000,0.000000}}, //F390N
   {{-0.080175,-0.099687,-0.119490,-0.126360},{-0.076787,-0.109971,-0.111782,-0.114968}}, //F410M
   {{0.000000,0.000000,0.000000,0.000000},{0.000000,0.000000,0.000000,0.000000}}, //F437N
   {{-0.016375,-0.035886,-0.055689,-0.062560},{-0.012986,-0.046170,-0.047982,-0.051167}}, //F439W
   {{-0.010060,-0.029571,-0.049374,-0.056245},{-0.006671,-0.039855,-0.041666,-0.044852}}, //F450W
   {{0.042616,0.023105,0.003302,-0.003569},{0.046005,0.012821,0.011009,0.007824}}, //F467M
   {{0.000000,0.000000,0.000000,0.000000},{0.000000,0.000000,0.000000,0.000000}}, //F469N
   {{0.000000,0.000000,0.000000,0.000000},{0.000000,0.000000,0.000000,0.000000}}, //F487N
   {{0.000000,0.000000,0.000000,0.000000},{0.000000,0.000000,0.000000,0.000000}}, //F502N
   {{0.020200,-0.016446,-0.032966,-0.017032},{-0.000606,-0.017946,-0.032272,-0.029324}}, //F547M
   {{0.012311,-0.024335,-0.040855,-0.024921},{-0.008495,-0.025835,-0.040160,-0.037213}}, //F555W
   {{0.035053,-0.001594,-0.018113,-0.002179},{0.014246,-0.003094,-0.017419,-0.014471}}, //F569W
   {{0.000000,0.000000,0.000000,0.000000},{0.000000,0.000000,0.000000,0.000000}}, //F588N
   {{0.035776,-0.000870,-0.017390,-0.001456},{0.014970,-0.002370,-0.016696,-0.013748}}, //F606W
   {{-0.022031,-0.036042,-0.033863,-0.030176},{0.002250,-0.042243,-0.045907,-0.024163}}, //F622W
   {{0.000000,0.000000,0.000000,0.000000},{0.000000,0.000000,0.000000,0.000000}}, //F631N
   {{0.000000,0.000000,0.000000,0.000000},{0.000000,0.000000,0.000000,0.000000}}, //F656N
   {{0.000000,0.000000,0.000000,0.000000},{0.000000,0.000000,0.000000,0.000000}}, //F658N
   {{0.000000,0.000000,0.000000,0.000000},{0.000000,0.000000,0.000000,0.000000}}, //F673N
   {{-0.021277,-0.035288,-0.033109,-0.029422},{0.003003,-0.041489,-0.045153,-0.023409}}, //F675W
   {{-0.022286,-0.036297,-0.034118,-0.030431},{0.001994,-0.042498,-0.046162,-0.024418}}, //F702W
   {{0.021840,-0.013261,-0.022826,-0.010518},{-0.008700,-0.017838,-0.018059,-0.011740}}, //F785LP
   {{0.051464,0.016363,0.006798,0.019106},{0.020924,0.011786,0.011565,0.017884}}, //F791W
   {{0.023574,-0.011527,-0.021092,-0.008784},{-0.006966,-0.016104,-0.016325,-0.010006}}, //F814W
   {{0.021980,-0.013121,-0.022686,-0.010378},{-0.008560,-0.017698,-0.017919,-0.011600}}, //F850LP
   {{0.000000,0.000000,0.000000,0.000000},{0.000000,0.000000,0.000000,0.000000}}, //F953N
   {{0.023574,-0.011527,-0.021092,-0.008784},{-0.006966,-0.016104,-0.016325,-0.010006}} //F1042M
};
*/

static void readdata(void) {
   int i,j,n,c;
   FILE *f;
   char str[1025];
   double eecorr[2];

   sprintf(str,"%s/wfc3/data/filters.dat",BASEDIR);
   if ((f=fopen(str,"r"))==NULL) {
      printf("Could not find %s\n",str);
      exit(0);
   }
   fscanf(f,"%d",&WFC3_NFILTERS);
   fgets(str,81,f);
   WFC3filters=(WFC3filttype*)calloc(sizeof(WFC3filttype),WFC3_NFILTERS);
   assert(WFC3filters!=NULL);
   for (i=0;i<WFC3_NFILTERS;i++) {
      fgets(WFC3filters[i].name,11,f); WFC3filters[i].name[strlen(WFC3filters[i].name)-1]=0;
      WFC3filters[i].color=fgetc(f);
      fgets(str,81,f);
      fscanf(f,"%lf %lf",WFC3filters[i].zp,WFC3filters[i].zp+1);
      WFC3filters[i].zp[2]=WFC3filters[i].zp[1];
      fscanf(f,"%lf %lf",eecorr,eecorr+1);
      fgets(str,81,f);
      for (j=0;j<8;j++) WFC3filters[i].xorder[j]='X';
      fscanf(f,"%d",&n);
      fgets(str,81,f);
      for (j=0;j<n;j++) {
	 fscanf(f,"%d",&c);
	 do {WFC3filters[i].xorder[j]=fgetc(f);} while (WFC3filters[i].xorder[j]==' ');
	 fscanf(f,"%lf %lf %lf %lf %lf",WFC3filters[i].xform[j]+2,WFC3filters[i].xform[j],WFC3filters[i].xform[j]+1,WFC3filters[i].xform[j]+3,WFC3filters[i].xform[j]+4);
	 WFC3filters[i].xform[j][2]-=WFC3filters[i].zp[c];
	 fgets(str,81,f);
      }
      //applying encircled energy correction afterwards since .xform is differential;
      WFC3filters[i].zp[0]-=eecorr[0];
      WFC3filters[i].zp[1]-=eecorr[1];
      WFC3filters[i].zp[2]-=eecorr[1];
   }
   fclose(f);
   return;
}

int WFC3findfilt(char*str) {
   int i;
   void WFC3initfilters();
   if (WFC3_NFILTERS<0) WFC3initfilters();
   for (i=0;i<WFC3_NFILTERS;i++) if (!strcmp(WFC3filters[i].name,str)) return i;
   printf("Illegal filter %s\n",str);
   exit(0);
   return -1;
}

static void readidcfile(char*base) {
   int i,j,c,d;
   FILE *f;
   char str[1025],tmp[81];

   sprintf(str,"%s/wfc3/data/%s_idctab.dat",BASEDIR,base);
   if ((f=fopen(str,"r"))==NULL) {
      printf("Could not find %s\n",str);
      exit(0);
   }
   while (fgets(str,1025,f)) {
      c=atoi(str);
      if (c<1 || c>2) {
	 printf("Stupid parse error\n%s",str);
	 exit(-1);
      }
      if (!strcmp(base,"ir")) c=0;
      if (!strncmp(str+10,"FORWARD",7)) d=0;
      else if (!strncmp(str+6,"INVERSE",7)) d=1;
      else {
	 printf("Stupid parse error\n%s",str);
	 exit(-1);
      }
      if (d==1) { // inverse
	 memcpy(tmp,str+17,12);
	 for (i=12;i>0 && tmp[i-1]==' ';i--);
	 tmp[i]=0;
	 i=WFC3findfilt(tmp);
	 for (j=0;j<34;j++) WFC3filters[i].idc[d][c][j]=atof(str+51+j*14);
      }
      else { // forward
	 memcpy(tmp,str+20,12);
	 for (i=12;i>0 && tmp[i-1]==' ';i--);
	 tmp[i]=0;
	 i=WFC3findfilt(tmp);
	 WFC3filters[i].idc[d][c][0]=atof(str+46) - 0.5; // xref
	 WFC3filters[i].idc[d][c][1]=atof(str+57) - 0.5; // yref
	 WFC3filters[i].idc[d][c][2]=atof(str+68); // rotation
	 WFC3filters[i].idc[d][c][3]=atof(str+81); // plate scale
	 WFC3filters[i].idc[d][c][4]=atof(str+94); // v2ref
	 WFC3filters[i].idc[d][c][5]=atof(str+107); // v3ref
	 for (j=6;j<34;j++) WFC3filters[i].idc[d][c][j]=atof(str+120+(j-6)*21);
      }
   }
   fclose(f);
   return;
}

static void readidc(void) {
   int i,c;
   readidcfile("ir");
   readidcfile("uvis");
   for (i=0;i<WFC3_NFILTERS;i++) {
      WFC3filters[i].idc[0][0][4] = 0;
      WFC3filters[i].idc[0][0][5] = 0;
      double v2_1 = 0.5*(WFC3filters[i].idc[0][1][4]-WFC3filters[i].idc[0][2][4]);
      double v3_1 = 0.5*(WFC3filters[i].idc[0][1][5]-WFC3filters[i].idc[0][2][5]);
      double theta = WFC3filters[i].idc[0][1][2]*M_PI/180.0;
      WFC3filters[i].idc[0][1][4] = v2_1*cos(theta) + v3_1*sin(theta); // delta y
      WFC3filters[i].idc[0][1][5] = v3_1*cos(theta) - v2_1*sin(theta); // delta x
      WFC3filters[i].idc[0][2][4] = -WFC3filters[i].idc[0][1][4];
      WFC3filters[i].idc[0][2][5] = -WFC3filters[i].idc[0][1][5];
      for (c=0;c<3;c++) {
	 WFC3filters[i].idc[1][c][4] = WFC3filters[i].idc[0][c][4];
	 WFC3filters[i].idc[1][c][5] = WFC3filters[i].idc[0][c][5];
      }
   }
   return;
}

void WFC3initfilters(void) {
   readdata();
   readidc();
   return;
}

static double findcolor(double ci[],char color,int *n) {
   int i;
   if (WFC3_NFILTERS<0) WFC3initfilters();
   for (i=0;i<WFC3_NFILTERS;i++) if (WFC3filters[i].color==color && ci[i]<99) {
      if (n!=NULL) *n=i;
      return ci[i];
   }
   return 99.999;
}

void WFC3transform(double vmag[],double dvmag[],double smag[]) {
   int f,i,j,CONT=1,it=0;
   double c,tm,tn;
   static double*smag0=NULL;

   if (WFC3_NFILTERS<0) WFC3initfilters();
   if (smag0==NULL) {
      smag0=(double*)calloc(sizeof(double),WFC3_NFILTERS);
      assert(smag0!=NULL);
   }
   memcpy(smag,vmag,8*WFC3_NFILTERS);
   while (CONT) {
      CONT=0;
      it++;
      memcpy(smag0,smag,8*WFC3_NFILTERS);
      for (f=0;f<WFC3_NFILTERS;f++) smag[f]=99.999;
      for (f=0;f<WFC3_NFILTERS;f++) if (smag0[f]<99) {
	 tm=tn=0;
	 for (i=0;i<8;i++) if (WFC3filters[f].xorder[i]!='X' && (c=findcolor(smag0,WFC3filters[f].xorder[i],&j))<99) {
	    if (j<f) c=c-smag0[f];
	    else c=smag0[f]-c;
	    if ((WFC3filters[f].xform[i][3]==-99 || c>WFC3filters[f].xform[i][3]) && (WFC3filters[f].xform[i][4]==99 || c<=WFC3filters[f].xform[i][4])) {
	       tm+=(WFC3filters[f].xform[i][0]*c+WFC3filters[f].xform[i][1]*c*c+WFC3filters[f].xform[i][2])/(dvmag[j]+0.01);
	       tn+=1./(dvmag[j]+0.01);
	    }
	 }
	 if (tn>0) {
	    smag[f]=vmag[f]+tm/tn;
	    if (it>5) smag[f]=(smag[f]*5+smag0[f]*(it-5))/it;
	    if (fabs(smag[f]-smag0[f])>0.001) CONT=1;
	 }
      }
   }
}

double WFC3untransform(int filt,double smag[7]) {
   int i,j,f0=-1;
   double c,tm=0.,tn=0.;
   char UBVRI[8]="UBVRIJH";

   if (WFC3_NFILTERS<0) WFC3initfilters();
   for (i=0;i<7;i++) if (WFC3filters[filt].color==UBVRI[i]) f0=i;
   if (f0<0 || smag[f0]>90) return 99999.;
   for (i=0;i<8;i++) for (j=1;j<7;j++) if (j!=f0 && WFC3filters[filt].xorder[i]==UBVRI[j] && smag[j]<99) {
      if (j<f0) c=smag[j]-smag[f0];
      else c=smag[f0]-smag[j];
      if ((WFC3filters[filt].xform[i][3]==-99 || c>WFC3filters[filt].xform[i][3]) && (WFC3filters[filt].xform[i][4]==99 || c<=WFC3filters[filt].xform[i][4])) {
	 tm-=WFC3filters[filt].xform[i][0]*c+WFC3filters[filt].xform[i][1]*c*c+WFC3filters[filt].xform[i][2];
	 tn++;
      }
   }
   if (!tn) for (i=0;i<8;i++) for (j=0;j<1;j++) if (j!=f0 && WFC3filters[filt].xorder[i]==UBVRI[j] && smag[j]<99) {
      if (j<f0) c=smag[j]-smag[f0];
      else c=smag[f0]-smag[j];
      if ((WFC3filters[filt].xform[i][3]==-99 || c>WFC3filters[filt].xform[i][3]) && (WFC3filters[filt].xform[i][4]==99 || c<=WFC3filters[filt].xform[i][4])) {
	 tm-=WFC3filters[filt].xform[i][0]*c+WFC3filters[filt].xform[i][1]*c*c+WFC3filters[filt].xform[i][2];
	 tn++;
      }
   }
   if (tn>0) return tm/tn+smag[f0];
   return 99999.;
}

double WFC3_CTE(int chip,float x,float y,float cts,float ctmult,float gn,float sky,float epoch) {
   /*
   double YCTE,SKY,FLUX;

   SKY=sky*gn/ctmult;
   if (SKY<0.2) SKY=0.2;
   FLUX=cts*gn/ctmult;
   if (FLUX<0.) FLUX=50.;
   else FLUX=sqrt(FLUX*FLUX+2500.);
   if (chip==0) {
      // new correction from Chiaberge et al.
      YCTE=0.363*pow(SKY,-0.15)*pow(FLUX,-0.36)*y/1000.*(epoch-52333.)/365.;
      // old correction from ISR03-09
      //YCTE=0.87*pow(SKY,-0.27)*pow(FLUX,-0.21)*y/2048.*(epoch-52333.)/381.;
   }
   else {
      // new correction from Chiaberge et al.
      YCTE=0.708*pow(SKY,-0.25)*pow(FLUX,-0.44)*(epoch-52333.)/365.;
      if (chip==1) YCTE*=(2048.-y)/2000.;
      else YCTE*=y/2000.;
      // old correction from ISR04-06
      //YCTE=0.0037818747*pow(SKY,-0.31)*pow(FLUX,-0.64)*(epoch-52333.);
      //if (chip==1) YCTE*=1.-y/2048.;
      //else YCTE*=y/2048.;
   }
   //printf("chip=%d, mlt=%3.1f, SKY=%4.2f, FLUX=%4.2f, EPOCH=%3.1f, Y=%3.1f, YCTE=%f\n",chip,ctmult,SKY,FLUX,epoch,y,YCTE);
   return YCTE;
   */
   return 0.0;
}

double WFC3_ZP(int filt,int cm) {
   if (WFC3_NFILTERS<0) WFC3initfilters();
   return WFC3filters[filt].zp[cm];
}
