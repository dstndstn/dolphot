typedef struct {
   char name[11],color;
   double zp[3];
   char xorder[8];
   int xformc[8];
   double xform[8][5];
   double idc[2][3][34];
} WFC3filttype;

extern int WFC3_NFILTERS;
extern WFC3filttype *WFC3filters;

extern int WFC3findfilt(char*);
extern void WFC3initfilters(void);

extern void WFC3transform(double vmag[],double dvmag[],double smag[]);
extern double WFC3untransform(int filt,double smag[7]);
extern double WFC3_CTE(int chip,float x,float y,float cts,float ctmult,float gn,float sky,float epoch);
extern double WFC3_ZP(int filt,int cm);
