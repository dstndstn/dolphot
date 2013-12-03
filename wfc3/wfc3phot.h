extern int WFC3_DRIZZLE_BASE;

//cm: HRC=0, WFC1=1, WFC2=2;
//vmag: VEGA magnitudes;
//smag: UBVRI magnitudes;

//correct distortion;
extern void WFC3shift(int img,double*x,double*y);

//add distortion;
extern void WFC3unshift(int img,double*x,double*y);

//transforms, returns values in smag[];
//extern void WFC3transform(int cm,double vmag[],double dvmag[],double smag[]);

//returns VEGAMAG;
//extern double WFC3untransform(int cm,int filt,double smag[]);

//set parameters needed for WFC3;
extern void wfc3initparam(void);

//initialize PSF libraries;
extern void wfc3initpsf(void);

//dispose PSF libraries;
extern void wfc3freepsf(void);

//calculate PSF from WFC3 library;
extern int calcwfc3psf(int img,float x,float y,int r,int force);

//do CTE corrections, photometric transformations, and output mag only;
extern double WFC3calcmag(int img,float x0,float y0,float ct0,float bg,int useCTE);

//do CTE corrections, photometric transformations, and output;
extern void WFC3outstar(FILE*,float x,float y,photdatatype*);
extern void WFC3outstarimg(int img,FILE*,float x,float y,photdatatype*);

//output column headers and image descriptions
extern void WFC3outstarinfo(FILE*,int*ct);
extern char *WFC3imagestring(int img);

//aperture correction radius (0.5 arcsec WFC/0.3 arcsec HRC)
extern float wfc3_apsize(int img,float x,float y);

//put lines in .info file;
extern void writewfc3info(void);

//correct counts in fake stars;
extern void WFC3readfakemag(FILE*);
extern void WFC3fixfakemag(int img,float x0,float y0,double*ct0,float*bg);
