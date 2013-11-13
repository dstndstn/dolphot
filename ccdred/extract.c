#include <fits.h>

ftype fitsin,fitsout;

int main(int argc,char**argv) {
   int C,X0,X1,Y0,Y1,x,y;
   FILE *f;

   if (argc!=4 && argc!=8 && argc!=9) {
      printf("Usage: %s <input> <output> <chip> <<X0> <X1> <Y0> <Y1>> <<mirror>>\n",*argv);
      return 1;
   }
   C=atoi(argv[3])-1;
   rreadfits(argv[1],&data,1);
   if (argc==4) {
      f=writefitsh(argv[2],X,Y,1,Ncards,cards);
      writechip(f,data[C],X,Y);
   }
   else {
      X0=atoi(argv[4]);
      X1=atoi(argv[5]);
      Y0=atoi(argv[6]);
      Y1=atoi(argv[7]);
      f=writefitsh(argv[2],X1-X0+1,Y1-Y0+1,1,Ncards,cards);
      if (argc==8) for (y=Y0;y<=Y1;y++) ffwrite(data[C][y]+X0,X1-X0+1,4,f);
      else for (y=Y0;y<=Y1;y++) for (x=X1;x>=X0;x--) ffwrite(data[C][y]+x,1,4,f);
   }
   fclose(f);
   return 0;
}
