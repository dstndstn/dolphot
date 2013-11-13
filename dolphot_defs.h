typedef enum {
   NONE,
   WFPC2,
   ACS,
   WFC3
} instrument;

typedef struct {
   instrument inst;
   int cm;
   int filt;
   int i1; // temperature (WFPC2)
   int i2; // gain (WFPC2)
} hstmodetype;

typedef struct{float ct0,ct,dct,ctcorr,dctcorr,sky,m,dm,chi,sh,rnd,r,e,a,b,c,crowd;int flag;} photdatatype;

typedef double apsftype[3][6];
typedef float apskytype[2];
typedef float dpostype[40];
typedef float wcstype[80];
typedef double wcsreftype[4];
typedef int hstoffsettype[2];

#ifdef DOLPHOT_MAIN
#define EXTERN
#else
#define EXTERN extern
#endif

EXTERN hstmodetype *hstmode;
EXTERN apsftype *apsf;
EXTERN apskytype *apsky;
EXTERN dpostype *dpos,*dpos0,*ref2img;
EXTERN wcstype *wcs;
EXTERN wcsreftype *wcsref;
EXTERN hstoffsettype *hstoffset;
EXTERN imtype *datahd,*dataim;
EXTERN int Timg,Nimg;
EXTERN int *rphot,*RPSF,SubPixel,FPSF,EPSF;
EXTERN double *RAper,*RChi,MinS,MaxS,MaxE,PSFStep,Zero,*apsize;
EXTERN int PSFsol;
EXTERN int FlagMask,UseCTE;
EXTERN double *iGAIN,*iEXP,*iEXP0,*iEPOCH,*apcor;
EXTERN int lastpsftype;
EXTERN float **psf;
EXTERN FILE *finfo;
EXTERN int poffreset;
EXTERN double *RSky0,*RSky1;
EXTERN int DEBUG;
EXTERN int DRIZZLE_BASE;
EXTERN int ACSpsfType;
EXTERN int WFC3psfType[3];
