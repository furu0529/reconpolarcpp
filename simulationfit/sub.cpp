#include"func.hpp"
//#include<iostream>
//#include<string>
//#include<math.h>
#include"fftw3.h"
#include "ceres/ceres.h"
#include "tensorflow/c/c_api.h"
//#include "CppCubicSpline.hpp"
//#include <boost/math/interpolators/cubic_b_spline.hpp>
//#include <boost/math/interpolators/cardinal_cubic_b_spline.hpp>
//#include <boost/math/interpolators/barycentric_rational.hpp>
using namespace std;

int main(void){
  make_ini m;
  double ktheta[grid][grid][bins]={};
  double rtheta[grid][grid][grid]={};
  double kphi[grid][grid][bins]={};
  double rphi[grid][grid][grid]={};
  double knorm[grid][grid][bins]={};
  double rnorm[grid][grid][grid]={};
  //  complex<double> k_y2m_temp[grid][grid][grid];
  char filename[] = "cambauto/test_transfer_z0.00_pig.csv";
  vector<double>  kt;
  vector<double>  transf;

  //  double *ktheta;
  //  ktheta = malloc(sizeof(double)*grid*grid*bins);
  cout << "start";
  /* start make initial condition and array */
  m.make_kspace(ktheta,kphi,knorm);
  m.make_rspace(rtheta,rphi,rnorm);
  m.make_iniperturb(knorm,kinitial_true);
  //  m.make_harm(ktheta,kphi);
  //   m.make_kharm(ktheta,kphi,k_y2m[0],0);
  for(int M=0;M<5;M++){
    m.make_kharm(ktheta,kphi,k_y2m[M],M);
    m.set_harm_turn(k_y2m[M]);
    /*for(int i=0;i<grid;i++){
      for(int j=0;j<grid;j++){
    	for(int k=0;k<grid;k++){
	  k_y2m[M][i][j][k]= k_y2m_temp[i][j][k];
	}
      }
      }*/
  }
  m.make_rharm(rtheta,rphi,r_y2m);
  /*  cout
    << "Table extent is "
    << std::extent<decltype(k_y2m[0]),0>::value
    << "x"
    << std::extent<decltype(k_y2m),0>::value
    << std::endl;*/
  m.read_transfer(filename,kt,transf);
  /* end make initial condition and array */

  /* start true QU map */
  make_trueQU mqu;
  complex<double> deltak[grid][grid][grid];
  mqu.make_evoperturb(knorm,kinitial_true,kt,transf,deltak);
  mqu.set_perturb_turn(deltak);
  mqu.mluti_perturb_y2m(deltak,k_y2m);
  mqu.Fouriertransfer();
  mqu.make_a2m();
  mqu.make_QU();
  /* fftw_complex in[grid][grid][grid], out[grid][grid][grid];
  fftw_plan p;
  p = fftw_plan_dft_3d(grid,grid,grid,&in[grid][grid][grid],&out[grid][grid][grid], FFTW_FORWARD, FFTW_ESTIMATE);
  for(int M=0;M<5;M++){
    for(int i=0;i<grid;i++){
      for(int j=0;j<grid;j++){
    	for(int k=0;k<grid;k++){
	  in[i][j][k][0] = k_box_harm[M][i][j][k].real();
	  in[i][j][k][1] = k_box_harm[M][i][j][k].imag();
	}
      }
    }
  fftw_execute(p); 
  fftw_destroy_plan(p);
    for(int i=0;i<grid;i++){
      for(int j=0;j<grid;j++){
    	for(int k=0;k<grid;k++){
	  r_box_harm[M][i][j][k]= out[i][j][k][0]+out[i][j][k][1]*I;
	}
      }
    }
  }
  fftw_free(in); fftw_free(out);
  */
  complex<double> i(0.0,1.0);
  complex<double> j(3.0,0.0);
  double k=3.0;
  cout << i+j*k;
  /*  for(int i=0;i<30;i++){
  cout << kt[i] << " " << transf[i] <<endl;
  }*/

  //  cout << kt.size() << " " << transf.size();
  //  boost::math::barycentric_rational<double> inp(kt.data(), transf.data(),transf.size());
  
  // test the cppcubicspline
  /* cout<<"cpp spline sample"<<endl;
  vector<double> sx = {0,1,2,3};
  vector<double> sy = {2.7,6,5,6.5};

  CppCubicSpline cppCubicSpline(sy);
  vector<double> rx;
  vector<double> ry;
  for(double i=0.0;i<=3.2;i+=0.1){
    rx.push_back(i);
    ry.push_back(cppCubicSpline.Calc(i));
    }
  for(int i=0;i<30;i++){
    cout<< rx[i] << " " << ry[i] <<endl;
    }*/
}
