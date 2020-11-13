#include"func.hpp"
#include <time.h>
using namespace std;

int main(int argc, char** argv){
  
  make_ini m;
  use_func u;
  int length[3] = {grid,grid,bins};
  //  double knorm[grid][grid][bins]={};
  //  complex<double> k_y2m[5][grid][grid][grid];
  //  complex<double> r_y2m[5][grid][grid][grid];
  //  complex<double> r_box_harm[5][grid][grid][grid];
  //complex<double> k_box_harm[5][grid][grid][grid];
  complex<double> kinitial_true[grid][grid][bins];
  //  complex<double> QU_true[grid][grid][grid];
  char filename[] = "cambauto/test_transfer_z0.00_pig.csv";
  //  vector<double>  kt;
  //  vector<double>  transf;
  clock_t start,inter,end,fourieri,fouriere;
  cout << "start make initial condition"<< endl;
  start = clock();
  m.make_kspace(knorms,k_y2m);
  for(int M=0;M<5;M++){
    u.set_kturn(length,k_y2m[M]);
  }
  m.make_rspace(r_y2m);
  m.make_iniperturb(knorms,kinitial_true);
  m.read_transfer(filename,kt,transf);
  kt.erase(kt.begin());
  transf.erase(transf.begin());
  cout << "finish initial condition" <<endl;
  cout << "start the evolution of fluctuation" <<endl;
  inter = clock();
  make_trueQU mqu;
  //  complex<double> deltak[grid][grid][grid];  
  mqu.make_evoperturb(knorms,kinitial_true,kt,transf,k_y2m,k_box_harm);
  //  mqu.mluti_perturb_y2m(deltak,k_y2m,k_box_harm);

  fftw_complex *in ,*out;
  fftw_plan p;
  fourieri = clock();
  for(int M=0;M<5;M++){
    in = (fftw_complex *)k_box_harm[M];
    out =(fftw_complex *)r_box_harm[M];
    p = fftw_plan_dft_3d(grid,grid,grid,in,out, FFTW_BACKWARD, FFTW_ESTIMATE);
    fftw_execute(p);
    fftw_destroy_plan(p);
  }
  cout << "finish fftw" << endl;
  fouriere = clock();
  //  fftw_free(in); fftw_free(out);
  mqu.make_qu(QU_true,r_box_harm,r_y2m);
  //  mqu.make_QU(QU_true,a2m_true,r_y2m);
  end = clock();
  printf("mkkini:%.4f[s]\n",(double)(fourieri-inter)/CLOCKS_PER_SEC);
  printf("fourier:%.4f[s]\n",(double)(fouriere-fourieri)/CLOCKS_PER_SEC);
  printf("makeinitial:%.4f[s]\n",(double)(inter-start)/CLOCKS_PER_SEC);
  printf("makeperturb:%.4f[s]\n",(double)(end-inter)/CLOCKS_PER_SEC);
  printf("all:%.4f[s]\n",(double)(end-start)/CLOCKS_PER_SEC);

  /* fitting */
  make_fitQU fit(kt,transf);
  int n=grid*grid*bins*2;
  double kinif[n];
  double t0 = 0.00001;
  double h0 = 0.25;
  double prin = 0;
  for(int i=0;i<n;i++){
    kinif[i]=1e-7;
  }
  //  cout << knorms[5][3][3];
  cout << "  Function value = " << fit.optimized_fuction(kinif,n) << "\n";
  //  praxis(t0,h0,n,prin,kinif,fit.optimized_fuction);
}	 
