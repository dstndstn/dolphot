typedef double double4[4];
extern int X0,Y0,nxy[2],ncmd[2];
extern double NSTAR,xystep[2],CMDSTEP,mmin,cmin;
extern char xyfn[321],cmdfn[321];

extern void process(int ext,int chip,int*Nfake,double4**fake);
