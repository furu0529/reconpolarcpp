#ifndef _FUNC_HPP_
#define _FUNC_HPP_
#include<complex>
#include<iostream>
#include<fstream>
#include<string>
#include<math.h>
#include <random>
#include<fstream>
#include "fftw3.h"
#include "global.hpp"
#include "linearinterp.hpp"
#include <boost/math/special_functions/spherical_harmonic.hpp>
#include "praxis.hpp"

using namespace std;

double optimfunc(double x[], int n);
double optimfuncnl(unsigned n, const double *x, double *grad, void *my_func_data);
void optimiresult(double x[],complex<double> (&a2mrf)[5][grid][grid][grid],complex<double> (&a2mrt)[5][grid][grid][grid]);

// 使用する関数及びnormalization factorなどを整理
class use_func{
private:
static constexpr double ns=0.9665;
static constexpr double As=2.089*1e-9;
static constexpr double kp=0.05;
  complex<double> spin;
public:
  static void scalar_power(double k,double &Pk) noexcept;
  static inline  void spinY2m(int m,double phi,double theta,complex<double>(&spin)) noexcept;
  static void set_kturn(int length[3],complex<double> (&k_box)[grid][grid][grid]) noexcept;
inline  static void set_kturn_conj(int length[3],complex<double> (&k_box)[grid][grid][grid]) noexcept;
  static double linearmy(double xi,double xi1,double yi,double yi1,double x);
  static double linearinterpmy(double x,vector<double> arr_x,vector<double>arr_y, int arr_length);
  static void writedata_1d(char name[], vector<double> yvec,int length,vector<double> vec);
  static void writedata_3d(char name[], double *box,int length[3],double *vec);
  static void writedata_3d_c(char name[], complex<double> *box,int length[3],double *vec);
  static void writedata_4d_c(char name[], complex<double> *box,int length[3],double *vec);
  static void read_4d_c(char name[], complex<double> *box);
};
//シミュレーションで用いる値のBOX(常に同じ値)を作成する。
class make_ini:public use_func{
 private:
 public:
  //double k_vec[grid][3],r_vec[grid][3];
  static  void make_kspace (double (&knorm)[grid][grid][bins],complex<double> (&k_y2m)[5][grid][grid][grid],double (&k_vec)[grid][3]);
  static  void make_rspace(complex<double> (&r_y2m)[5][grid][grid][grid],double (&r_vec)[grid][3]) noexcept;
  static  void make_iniperturb(double knorm[grid][grid][bins],complex<double> (&kinitial)[grid][grid][bins]);
  static void set_harm_turn_r(complex<double> (&r_box)[grid][grid][grid]);
  static void read_transfer(char name[],vector<double> &k,vector<double> &transf);
};
//真のQUmapの作製
class make_trueQU:public make_ini{
private:
public:
inline static void make_evoperturb(double knorm[grid][grid][bins],complex<double> kinitial[grid][grid][bins] ,vector<double> kt,vector<double>  transf,complex<double> k_y2m[5][grid][grid][grid],complex<double>(&a2mk)[5][grid][grid][grid]);
inline static void make_evoperturb_fast(double knorm[grid][grid][bins],const double *kinitial ,vector<double> kt,vector<double>  transf,complex<double> k_y2m[5][grid][grid][grid],complex<double>(&a2mk)[5][grid][grid][grid]);
inline  static void make_qu(complex<double> (&QU)[grid][grid][grid],complex<double> (&rbox)[5][grid][grid][grid],complex<double> r_y2m[5][grid][grid][grid]);
};
    
#endif // _FUNC_HPP_
