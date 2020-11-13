#include"func.hpp"
#include <time.h>
using namespace std;
int main(void){
  make_ini m;
  use_func u;
  int length[3] = {grid,grid,bins};
  int rlength[3] = {grid,grid,grid};
  double knorms[grid][grid][bins];
  double k_vec[grid][3],r_vec[grid][3];
  //double knorm[grid][grid][bins]={};
  //  complex<double> k_y2m[5][grid][grid][grid];
  complex<double> kinitial_true[grid][grid][bins];
  complex<double> QU_true[grid][grid][grid];
  double *fp, *fp2;
  complex<double> *fp3;
  char filename[] = "cambauto/test_transfer_z0.00_pig.csv";
  //char name[256];
  vector<double>  kt;
  vector<double>  transf;
  complex<double> r_y2m[5][grid][grid][grid];
  complex<double> k_y2m[5][grid][grid][grid];
  static  complex<double> r_box_harm[5][grid][grid][grid];
  static complex<double> k_box_harm[5][grid][grid][grid];
  double transftest[grid][grid][bins];
  clock_t start,end;
  start = clock();
  cout << "start make initial condition"<< endl;
  m.make_kspace(knorms,k_y2m,k_vec);
  for(int M=0;M<5;M++){
    u.set_kturn(length,k_y2m[M]);
  }
  fp = (double *)knorms;
  fp2 = (double *)k_vec;
  u.writedata_3d("knorms.dat",fp,length,fp2);
  fp3  = (complex<double> *)k_y2m;
  u.writedata_4d_c("kharm.dat",fp3,rlength,fp2);
  cout << "finish writing the sphere harmonic in kspace" << endl;

  m.make_rspace(r_y2m,r_vec);
  // check the distance and polar coordinate in real space
  fp2 = (double *)r_vec;
  fp3 = (complex<double> *)r_y2m;
  u.writedata_4d_c("rharm.dat",fp3,rlength,fp2);

  cout << "finish writing the sphere harmonic in rspace" << endl;
  
  m.make_iniperturb(knorms,kinitial_true);
  // check the initial perturbation
  fp2 = (double *)k_vec;
  fp3  = (complex<double> *)kinitial_true;
  u.writedata_3d_c("kinitial.dat",fp3,length,fp2);

  //check the initial curvature power spectrum
  double Pk[grid];
  ofstream outputfile0("inipk.dat");
  for(int i=0;i<grid;i++){
    u.scalar_power(k_vec[i][0],Pk[i]);
    outputfile0 << k_vec[i][0] << " " << Pk[i] << endl;
  }
  outputfile0.close();

  //check the normal random
  /*  random_device rnd;     // 非決定的な乱数生成器でシード生成機を生成
  mt19937 mt(rnd()); //  メルセンヌツイスターの32ビット版、引数は初期シード
  for(int i=0;i<grid;i++){
    normal_distribution<double> norm(0.0,sqrt(Pk[i]));
    cout << norm(mt) << endl;
    }*/
    
  // check the transfer function data
  m.read_transfer(filename,kt,transf);
  kt.erase(kt.begin());
  transf.erase(transf.begin());
  /*  ofstream outputfilet("transfercheck2.dat");
  for(int i=0;i<grid;i++){
    for(int j=0;j<grid;j++){
      for(int k=0;k<grid;k++){
	outputfilet << knorms[i][j][k] << " ";
	outputfilet << linearinterpmy(knorms[i][j][k],kt,transf,100) <<endl;
      }
    }
  }
  outputfilet.close(); */
  u.writedata_1d("transfer.dat",transf,kt.size(),kt);
  cout << "finish initial condition" <<endl;
  
  cout << "start the evolution of fluctuation" <<endl;
  make_trueQU mqu;
  
  // check the initial perturbation value in fourier space
  mqu.make_evoperturb(knorms,kinitial_true,kt,transf,k_y2m,k_box_harm);
  fp2 = (double *)k_vec;
  fp3 = (complex<double> *)k_box_harm;
  u.writedata_4d_c("a2mk.dat",fp3,rlength,fp2);
  
  fftw_complex *in ,*out;
  fftw_plan p;
  for(int M=0;M<5;M++){
    in = (fftw_complex *)k_box_harm[M];
    out =(fftw_complex *)r_box_harm[M];
    p = fftw_plan_dft_3d(grid,grid,grid,in,out, FFTW_BACKWARD, FFTW_ESTIMATE);
    fp3 =(complex<double> *)in;
    u.writedata_3d_c("a2min.dat",fp3,rlength,fp2);
    fftw_execute(p);
    fp3 =(complex<double> *)out;
    u.writedata_3d_c("a2mout.dat",fp3,rlength,fp2);
    fftw_destroy_plan(p);
  }

  // check the perturbation in real space after Fouriertransfer
  fp2 = (double *)r_vec;
  fp3 = (complex<double> *)r_box_harm;
  u.writedata_4d_c("a2mr.dat",fp3,rlength,fp2);

  for(int M=0;M<5;M++){
     m.set_harm_turn_r(r_box_harm[M]);
   }
   

  // check the QUvalue
  mqu.make_qu(QU_true,r_box_harm,r_y2m);
  fp3 = (complex<double> *)QU_true;
  u.writedata_3d_c("qu.dat",fp3,rlength,fp2);

  end = clock();
  printf("%.4f秒かかりました\n",(double)(end-start)/CLOCKS_PER_SEC);
  
}
