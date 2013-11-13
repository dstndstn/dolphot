#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

int N=0,fixed[3]={0,0,0},Nfree=0;
double Xd[1000],Cd[1000],Sd[1000],Md[1000],dMd[1000],val[3]={0,0,0},dval[3]={0,0,0};
double avX=0,avC=0,chi2=0;

int main(int argc,char**argv) {
   int i;
   double I=0,S=0,X=0,C=0,XX=0,XC=0,CC=0,SX=0,SC=0,q;
   char dum[1024];

   for (i=1;i<argc;i++) {
      if (!strncmp(argv[i],"a=",2) || !strncmp(argv[i],"A=",2)) {
	 fixed[0]=1;
	 val[0]=atof(argv[i]+2);
      }
      else if (!strncmp(argv[i],"b=",2) || !strncmp(argv[i],"B=",2)) {
	 fixed[1]=1;
	 val[1]=atof(argv[i]+2);
      }
      else if (!strncmp(argv[i],"c=",2) || !strncmp(argv[i],"C=",2)) {
	 fixed[2]=1;
	 val[2]=atof(argv[i]+2);
      }
      else {
	 printf("Usage: %s <<options>>\n",*argv);
	 printf("  a=10  to set zero point\n");
	 printf("  b=10  to set airmass correction\n");
	 printf("  c=10  to set color correction\n");
	 return -1;
      }
   }
   for (i=0;i<3;i++) if (!fixed[i]) Nfree++;
   while (fscanf(stdin,"%lf %lf %lf %lf %lf",Sd+N,Xd+N,Cd+N,Md+N,dMd+N)!=EOF) {
      avC+=Cd[N];
      avX+=Xd[N];
      dMd[N]=1./dMd[N]/dMd[N];
      N++;
      if (N>=1000) {
	 printf("Too many points\n");
	 return -1;
      }
      fgets(dum,1024,stdin);
   }
   avC/=N;
   avX/=N;
   for (i=0;i<N;i++) {
      double dm;

      Xd[i]-=avX;
      Cd[i]-=avC;
      dm=Sd[i]-Md[i];
      if (fixed[0]) dm+=val[0];
      if (fixed[1]) dm+=val[1]*Xd[i];
      if (fixed[2]) dm+=val[2]*Cd[i];
      I+=dMd[i];
      S+=dm*dMd[i];
      X+=Xd[i]*dMd[i];
      C+=Cd[i]*dMd[i];
      XX+=Xd[i]*Xd[i]*dMd[i];
      CC+=Cd[i]*Cd[i]*dMd[i];
      SX+=dm*Xd[i]*dMd[i];
      SC+=dm*Cd[i]*dMd[i];
      XC+=Xd[i]*Cd[i]*dMd[i];
   }
   if (fixed[0]==0 && fixed[1]==0 && fixed[2]==0) {
      // 0=S+a*I+b*X+c*C;
      // 0=SX+a*X+b*XX+c*XC;
      // 0=SC+a*C+b*XC+c*CC;
      q=C*C*X*X-C*C*I*XX-I*CC*X*X+I*CC*I*XX+C*X*I*XC-I*XC*I*XC-C*X*X*C+I*XC*X*C;
      val[0]=(X*I*SX*CC-S*CC*I*XX+S*I*XC*XC-X*I*XC*SC+X*C*X*SC-C*SC*X*X+C*I*SC*XX-C*I*SX*XC)/q;
      val[1]=(I*SX*C*C-I*I*SX*CC-S*X*C*C+S*X*I*CC+I*XC*I*SC-I*XC*S*C-C*X*I*SC+C*X*S*C)/q;
      val[2]=(I*SC*X*X-I*I*SC*XX-S*C*X*X+S*C*I*XX+I*I*SX*XC-I*SX*X*C+S*X*X*C-S*X*I*XC)/q;
      for (i=0;i<N;i++) {
	 dval[0]+=(X*I*Xd[i]*CC-CC*I*XX+I*XC*XC-X*I*XC*Cd[i]+X*C*X*Cd[i]-C*Cd[i]*X*X+C*I*Cd[i]*XX-C*I*Xd[i]*XC)*(X*I*Xd[i]*CC-CC*I*XX+I*XC*XC-X*I*XC*Cd[i]+X*C*X*Cd[i]-C*Cd[i]*X*X+C*I*Cd[i]*XX-C*I*Xd[i]*XC)*dMd[i];
	 dval[1]+=(I*Xd[i]*C*C-I*I*Xd[i]*CC-X*C*C+X*I*CC+I*XC*I*Cd[i]-I*XC*C-C*X*I*Cd[i]+C*X*C)*(I*Xd[i]*C*C-I*I*Xd[i]*CC-X*C*C+X*I*CC+I*XC*I*Cd[i]-I*XC*C-C*X*I*Cd[i]+C*X*C)*dMd[i];
	 dval[2]+=(I*Cd[i]*X*X-I*I*Cd[i]*XX-C*X*X+C*I*XX+I*I*Xd[i]*XC-I*Xd[i]*X*C+X*X*C-X*I*XC)*(I*Cd[i]*X*X-I*I*Cd[i]*XX-C*X*X+C*I*XX+I*I*Xd[i]*XC-I*Xd[i]*X*C+X*X*C-X*I*XC)*dMd[i];
      }
      dval[0]=sqrt(dval[0])/q;
      dval[1]=sqrt(dval[1])/q;
      dval[2]=sqrt(dval[2])/q;
   }
   else if (!fixed[0] && !fixed[1]) {
      val[0]=(X*SX-S*XX)/(I*XX-X*X);
      val[1]=(S*X-I*SX)/(I*XX-X*X);
      dval[0]=sqrt(I*XX*XX-X*X*XX)/(I*XX-X*X);
      dval[1]=sqrt(I*I*XX-I*X*X)/(I*XX-X*X);
   }
   else if (!fixed[0] && !fixed[2]) {
      val[0]=(C*SC-S*CC)/(I*CC-C*C);
      val[2]=(S*C-I*SC)/(I*CC-C*C);
      dval[0]=sqrt(I*CC*CC-C*C*CC)/(I*CC-C*C);
      dval[2]=sqrt(I*I*CC-I*C*C)/(I*CC-C*C);
   }
   else if (!fixed[1] && !fixed[2]) {
      val[1]=(SC*XC-SX*CC)/(XX*CC-XC*XC);
      val[2]=(SX*XC-SC*XX)/(XX*CC-XC*XC);
      dval[1]=sqrt(XX*CC*CC-CC*XC*XC)/(XX*CC-XC*XC);
      dval[2]=sqrt(XX*XX*CC-XX*XC*XC)/(XX*CC-XC*XC);
   }
   else if (!fixed[0]) {
      val[0]=-S/I;
      dval[0]=1/sqrt(I);
   }
   else if (!fixed[1]) {
      val[1]=-SX/XX;
      dval[1]=1/sqrt(XX);
   }
   else if (!fixed[2]) {
      val[2]=-SC/CC;
      dval[1]=1/sqrt(CC);
   }
   printf("Fit: OBSMAG = STDMAG + (%f +/- %f)\n",val[0],dval[0]);
   printf("            + (%f +/- %f) (%f+X)\n",val[1],dval[1],-avX);
   printf("            + (%f +/- %f) (%f+STDCOL)\n",val[2],dval[2],-avC);
   printf("errors:\n");
   for (i=0;i<N;i++) {
      printf("%f %f\n",Md[i]-(Sd[i]+val[0]+val[1]*Xd[i]+val[2]*Cd[i]),1/sqrt(dMd[i]));
      chi2+=(Sd[i]+val[0]+val[1]*Xd[i]+val[2]*Cd[i]-Md[i])*(Sd[i]+val[0]+val[1]*Xd[i]+val[2]*Cd[i]-Md[i])*dMd[i];
   }
   printf("chi = %f\n",sqrt(chi2/(N-Nfree)));
   return 0;
}
