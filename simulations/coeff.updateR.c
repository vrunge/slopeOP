/*attempt to code up coeff.update function*/

#include <math.h>
#include <stdio.h>
#include <R.h>

/*void coeffupdate(double *coeffs[40], int *taustar, double *sigsquared, double *beta);

int main() 
{
  double coeffs[40]={1,1,1,2,3,
		       1,2,1,2,4,
		       2,2,1,2,3,
		       3,2,1,2,4,
		       4,3,1,2,5,
		       4,1,2,3,4,
		       4,2,3,4,1,
		       5,2,1,3,4};
  coeffupdate(coeffs,9,1,3);
  return (0);
  }*/

  void coeffupdate(double *coeffs, double *S, double *SJ, double *SS, int *taustar, double *sigsquared, double *beta, int *nrows, double *coeffnew)
{
  /*Stuff the function will provide
  double SJ[11]={23,4,2,1,4,3,5,3,2,5,12};
  double S[11]={2,41,2,1,2,4,5,9,6,8,12};
  double SS[11]={2,41,2,1,2,4,5,9,6,8,12};
  /*Actual stuff we need*/
  int i,j,sstar;
  double A,B,C,D,E,FF,seglen;
  /* coeffnew = (double *)calloc(5**nrows,sizeof(double)); /*take out this line to make it work */
  for(i = 0; i < *nrows; i++) {
    *(coeffnew+5*i)=*taustar;
    *(coeffnew+5*i+1)=*(coeffs+5*i);
    sstar=*(coeffnew+5*i+1);
    seglen=*taustar-sstar;
    A=(seglen+1)*(2*seglen+1)/(12*seglen* *sigsquared);
    B=(seglen*seglen-1)/(6*seglen* *sigsquared);
    C=(-1)/(seglen* *sigsquared)*(*(SJ+*taustar)-*(SJ+sstar)-sstar*(*(S+*taustar)-*(S+sstar)));
    D=seglen/2*log(2*M_PI* *sigsquared)+1/(2* *sigsquared)*(*(SS+*taustar)-*(SS+sstar));
    E=(-1)*C-1/ *sigsquared*(*(S+*taustar)-*(S+sstar));
    FF=(seglen-1)*(2*seglen-1)/(12*seglen* *sigsquared);
    if((FF==0) && (*(coeffs+5*i+2)==0)){
      if(B==0){
	*(coeffnew+5*i+4)=*(coeffs+5*i+4)+D+ *beta;
	*(coeffnew+5*i+3)=C;
	*(coeffnew+5*i+2)=A;
      }
      else{
	*(coeffnew+5*i+4)=(-E-*(coeffs+5*i+3))/B+ *beta;
	*(coeffnew+5*i+3)=0;
	*(coeffnew+5*i+2)=0;
      }
    }
    else{
      *(coeffnew+5*i+4)=*(coeffs+5*i+4)+D-(*(coeffs+5*i+3)+E)*(*(coeffs+5*i+3)+E)/(4*(*(coeffs+5*i+2)+FF))+ *beta;
      *(coeffnew+5*i+3)=C-(*(coeffs+5*i+3)+E)*B/(2*(*(coeffs+5*i+2)+FF));
      *(coeffnew+5*i+2)=A-B*B/(4*(*(coeffs+5*i+2)+FF));
    }
  }  

}
