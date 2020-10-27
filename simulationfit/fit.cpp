#include"func.hpp"
#include <time.h>
#include "nlopt.h"
using namespace std;

int main(int argc, char** argv){
  
  clock_t start,end;
  int n=grid*grid*bins*2;
  cout << n;
  double kinif[n];
  double t0 = 1e-5;
  double h0 = 1e-5;
  double prin = 2;
  for(int i=0;i<n;i++){
    kinif[i]=1e-7;
  }
  start = clock();
  nlopt_opt opt;
  opt = nlopt_create(NLOPT_LN_PRAXIS, n);
  nlopt_set_min_objective(opt, optimfuncnl, NULL);
  nlopt_set_xtol_rel(opt, 1e-5);
  //  nlopt_set_stopval(opt, sqrt(8./27.)+1e-3);
  // nlopt_set_stopval(opt, 1e-5);
  double minf;
  if (nlopt_optimize(opt, kinif, &minf) < 0) {
    printf("nlopt failed!\n");
}
else {
  r8vec_print ( n, kinif, "  Computed minimizer:" );
  //  printf("found minimum at f(%g,%g) = %0.10g\n", minf);
    }
  //  praxis(t0,h0,n,prin,kinif,optimfunc);
  end = clock();

  //  r8vec_print ( n, kinif, "  Computed minimizer:" );
  //cout << "  Function value = " << optimfunc( kinif, n ) << "\n";
  printf("all:%.4f[s]\n",(double)(end-start)/CLOCKS_PER_SEC);
  complex<double> a2mrfout[5][grid][grid][grid];
  complex<double> a2mrtout[5][grid][grid][grid];
  /*complex<double> *fp;
  fp = (complex<double> *)a2mrtout;
  use_func::read_4d_c("a2mt.dat",fp);
  for(int M=0;M<5;M++){
    cout << "true:" << endl;
    cout << a2mrtout[0][bins][bins][bins] << endl;
    }*/
  optimiresult(kinif,a2mrfout,a2mrtout);
}	 
