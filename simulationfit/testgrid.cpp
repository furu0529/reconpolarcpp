#include"func.hpp"
#include "ceres/ceres.h"
#include <time.h>
using namespace std;
int main(void){
  use_func u;
  make_ini m;
  int length[3] = {grid,grid,bins};
  double knorm[grid][grid][bins];
  double ***ktheta = new double**[grid];
  double ***kphi = new double**[grid];
  for (int i = 0; i < grid; i++) {
    ktheta[i] = new double*[grid];
    kphi[i] = new double*[grid];
    for (int j = 0; j < grid; j++) {
      ktheta[i][j] = new double[bins];
      kphi[i][j] = new double[bins];
    }
  }
  
  m.make_kspace(ktheta,kphi,knorm);/*
  m.make_iniperturb(knorm,kinitial_true);
  for(int M=0;M<5;M++){
    m.make_kharm(ktheta,kphi,k_y2m[M],M-2);
    u.set_kturn(length,k_y2m[M]);
  }
  for(int i=0;i<grid;i++){
    for(int j=0;j<grid;j++){
      delete[] ktheta[i][j];
      delete[] kphi[i][j];
    }
    delete[] ktheta[i];
    delete[] kphi[i];
  }
  delete[] ktheta;
  delete[] kphi;
  */
}
