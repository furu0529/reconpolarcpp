#include "func.hpp"
double optimfunc(double x[], int n){
  static int flag=0;
  static int length[3] = {grid,grid,bins};
  static double k_vec[grid][3],r_vec[grid][3];
  static double knorm[grid][grid][bins];
  static complex<double> k_y2m[5][grid][grid][grid];
  static complex<double> r_y2m[5][grid][grid][grid];
  static vector<double>  kt;
  static vector<double>  transf;
  static  complex<double> kinic[grid][grid][bins];
  static complex<double> kinit[grid][grid][bins];
  static complex<double> a2mkf[5][grid][grid][grid];
  static complex<double> a2mkt[5][grid][grid][grid];
  static complex<double> rboxf[5][grid][grid][grid];
  static complex<double> rboxt[5][grid][grid][grid];
  static complex<double> QUf[grid][grid][grid];
  static complex<double> QUt[grid][grid][grid];
  char filename[] = "cambauto/test_transfer_z0.00_pig.csv";
  fftw_complex *in ,*out;
  fftw_plan p;
  
   for(int i=0;i<grid;i++){
     for(int j=0;j<grid;j++){
       for(int k=0;k<bins;k++){
	 kinic[i][j][k]=x[i*grid*bins*2+j*bins*2+k*2]+x[i*grid*bins*2+j*bins*2+k*2+1]*I;
       }
     }
   }
   
   if(flag == 0){
     cout << "make_initial" <<endl;
     make_ini::make_kspace(knorm,k_y2m,k_vec);
     for(int M=0;M<5;M++){
       use_func::set_kturn(length,k_y2m[M]);
     }
     make_ini::make_rspace(r_y2m,r_vec);
     make_ini::read_transfer(filename,kt,transf);
     kt.erase(kt.begin());
     transf.erase(transf.begin());
     //cout << "make_true start"  <<endl;
     make_ini::make_iniperturb(knorm,kinit);
     kinit[0][0][0]=0.0;
     make_trueQU::make_evoperturb(knorm,kinit,kt,transf,k_y2m,a2mkt);
     for(int M=0;M<5;M++){
       in = (fftw_complex *)a2mkt[M];
       out =(fftw_complex *)rboxt[M];
       p = fftw_plan_dft_3d(grid,grid,grid,in,out, FFTW_BACKWARD, FFTW_ESTIMATE);
       fftw_execute(p);
       fftw_destroy_plan(p);
     }
     make_trueQU::make_qu(QUt,rboxt,r_y2m);
     ofstream outputfilet("opre.dat");
     for(int M=0;M<5;M++){
       outputfilet << "true" << endl;
       outputfilet << real(rboxt[M][bins][bins][bins])*(-4.0*M_PI)/sqrt(pow(2.0*M_PI*L,3.0)) << " ";
       outputfilet << imag(rboxt[M][bins][bins][bins])*(-4.0*M_PI)/sqrt(pow(2.0*M_PI*L,3.0)) << endl;
     }
     outputfilet.close();
     /*     ofstream outputfilet("a2mt.dat");
     for(int M=0;M<5;M++){
       for(int i=0;i<grid;i++){
	 for(int j=0;j<grid;j++){
	   for(int k=0;k<grid;k++){
	     outputfilet << real(rboxt[M][i][j][k])*(-4.0*M_PI)/sqrt(pow(2.0*M_PI*L,3.0)) << " ";
	     outputfilet << imag(rboxt[M][i][j][k])*(-4.0*M_PI)/sqrt(pow(2.0*M_PI*L,3.0)) << endl;
	   }
	 }
       }
     }
     outputfilet.close();*/
     //cout << "make_true end" <<endl;
     flag++;
     }
   /*else{
     cout << "already make" << endl;
     }*/
   make_trueQU::make_evoperturb(knorm,kinic,kt,transf,k_y2m,a2mkf);
  for(int M=0;M<5;M++){
    in = (fftw_complex *)a2mkf[M];
    out =(fftw_complex *)rboxf[M];
    p = fftw_plan_dft_3d(grid,grid,grid,in,out, FFTW_BACKWARD, FFTW_ESTIMATE);
    fftw_execute(p);
    fftw_destroy_plan(p);
  }
  /*  cout << QUf[3][3][3] << endl;
  for(int M=0;M<5;M++){
    cout << rboxf[M][3][3][3] << endl;
    }*/
  make_trueQU::make_qu(QUf,rboxf,r_y2m);
  /*for(int M=0;M<5;M++){
    cout << rboxf[M][3][3][3] << endl;
  }
  cout << QUf[3][3][3] << endl;*/
  double funcvalue=0.0;
  for(int i=0;i<grid;i++){
    for(int j=0;j<grid;j++){
      for(int k=0;k<grid;k++){
	//	if( (bins-i)*(bins-i) + (bins-j)*(bins-j) + (bins-k)*(bins-k) <= 3*3){
        funcvalue = funcvalue + pow(real(QUf[i][j][k])-real(QUt[i][j][k]),2);
        funcvalue = funcvalue + pow(imag(QUf[i][j][k])-imag(QUt[i][j][k]),2);
	//}
      }
    }
  }
  funcvalue = funcvalue-pow(real(QUf[bins][bins][bins])-real(QUt[bins][bins][bins]),2);
  funcvalue = funcvalue-pow(imag(QUf[bins][bins][bins])-imag(QUt[bins][bins][bins]),2);
  
  return funcvalue;
}

double optimfuncnl(unsigned n, const double *x, double *grad, void *my_func_data){
  static int flag=0;
  static int length[3] = {grid,grid,bins};
  static double k_vec[grid][3],r_vec[grid][3];
  static double knorm[grid][grid][bins];
  static complex<double> k_y2m[5][grid][grid][grid];
  static complex<double> r_y2m[5][grid][grid][grid];
  static vector<double>  kt;
  static vector<double>  transf;
  //static complex<double> kinic[grid][grid][bins];
  static complex<double> kinit[grid][grid][bins];
  static complex<double> a2mkf[5][grid][grid][grid];
  static complex<double> a2mkt[5][grid][grid][grid];
  static complex<double> rboxf[5][grid][grid][grid];
  static complex<double> rboxt[5][grid][grid][grid];
  //  static complex<double> QUf[grid][grid][grid];
  static complex<double> QUt[grid][grid][grid];
  static complex<double> a2m,qutemp;
  static  double con;
  char filename[] = "cambauto/test_transfer_z0_pig.csv";
  fftw_complex *in ,*out;
  fftw_plan p;
  
  if(flag == 0){
    cout << "make_initial" <<endl;
    make_ini::make_kspace(knorm,k_y2m,k_vec);
    for(int M=0;M<5;M++){
      use_func::set_kturn(length,k_y2m[M]);
    }
    make_ini::make_rspace(r_y2m,r_vec);
    make_ini::read_transfer(filename,kt,transf);
    kt.erase(kt.begin());
    transf.erase(transf.begin());
    cout << "make_true start"  <<endl;
    make_ini::make_iniperturb(knorm,kinit);
    kinit[0][0][0]=0.0;
    make_trueQU::make_evoperturb(knorm,kinit,kt,transf,k_y2m,a2mkt);
    for(int M=0;M<5;M++){
      in = (fftw_complex *)a2mkt[M];
      out =(fftw_complex *)rboxt[M];
      p = fftw_plan_dft_3d(grid,grid,grid,in,out, FFTW_BACKWARD, FFTW_ESTIMATE);
      fftw_execute(p);
      fftw_destroy_plan(p);
    }
    make_trueQU::make_qu(QUt,rboxt,r_y2m);
    cout << "make_true end" <<endl;
    con = (-4.0*M_PI)/sqrt(pow(2.0*M_PI*L,3.0));
    ofstream outputfilet("opre.dat");
    for(int M=0;M<5;M++){
      outputfilet << "true" << endl;
      outputfilet << real(rboxt[M][bins][bins][bins])*(-4.0*M_PI)/sqrt(pow(2.0*M_PI*L,3.0)) << " ";
      outputfilet << imag(rboxt[M][bins][bins][bins])*(-4.0*M_PI)/sqrt(pow(2.0*M_PI*L,3.0)) << endl;
    }
    outputfilet.close();
    flag++;
  }
  /*else{
    cout << "already make" << endl;
    }*/
  make_trueQU::make_evoperturb_fast(knorm,x,kt,transf,k_y2m,a2mkf);
  for(int M=0;M<5;M++){
    in = (fftw_complex *)a2mkf[M];
    out =(fftw_complex *)rboxf[M];
    p = fftw_plan_dft_3d(grid,grid,grid,in,out, FFTW_BACKWARD, FFTW_ESTIMATE);
    fftw_execute(p);
    fftw_destroy_plan(p);
  }
  double funcvalue=0.0;
  for(int i=0;i<grid;i++){
    for(int j=0;j<grid;j++){
      for(int k=0;k<grid;k++){
        qutemp = 0.0;
        for(int M=0;M<5;M++){
          a2m = rboxf[M][i][j][k]*con;
          qutemp += a2m*r_y2m[M][i][j][k];
        }
	funcvalue = funcvalue + pow(real(qutemp)-real(QUt[i][j][k]),2);
        funcvalue = funcvalue + pow(imag(qutemp)-imag(QUt[i][j][k]),2);
      }
    }
  }
  /* 中心を別計算 */
  qutemp=0.0;
  for(int M=0;M<5;M++){
    a2m = rboxf[M][bins][bins][bins]*con;
    qutemp += a2m*r_y2m[M][bins][bins][bins];
  }
  funcvalue = funcvalue-pow(real(qutemp)-real(QUt[bins][bins][bins]),2);
  funcvalue = funcvalue-pow(imag(qutemp)-imag(QUt[bins][bins][bins]),2); 
  return funcvalue;
}

  /*ofstream outputfilet("fittingcheck.dat");
  for(int i=0;i<grid;i++){
    //    outputfilet << kt[i] << " ";
    //    outputfilet << transf[i] << endl;
    for(int j=0;j<grid;j++){
      for(int k=0;k<grid;k++){
        //outputfilet << knorm[i][j][k] << " ";
	//outputfilet << k_y2m[0][i][j][k] << endl;
	//	outputfilet << r_y2m[1][i][j][k] << endl;
	//	outputfilet << kinic[i][j][k] << " "; 
	//outputfilet << a2mkf[1][i][j][k] << endl;
	outputfilet << QUt[i][j][k]-QUf[i][j][k] << endl;
	//	outputfilet << real(QUt[i][j][k]) << " ";
	//outputfilet << real(QUf[i][j][k]) <<endl;
	//        outputfilet << linearinterpmy(knorm[i][j][k],kt,transf,100) <<endl;
      }
    }
    }
    outputfilet.close();*/


void optimiresult(double x[],complex<double> (&a2mrf)[5][grid][grid][grid],complex<double> (&a2mrt)[5][grid][grid][grid]){
  static int length[3] = {grid,grid,bins};
  static double knorm[grid][grid][bins];
  static double k_vec[grid][3],r_vec[grid][3];
  static complex<double> k_y2m[5][grid][grid][grid];
  static complex<double> r_y2m[5][grid][grid][grid];
  static vector<double>  kt;
  static vector<double>  transf;
  static complex<double> kinic[grid][grid][bins];
  static complex<double> kinit[grid][grid][bins];
  static complex<double> a2mkf[5][grid][grid][grid];
  static complex<double> a2mkt[5][grid][grid][grid];
  char filename[] = "cambauto/test_transfer_z0.00_pig.csv";
  fftw_complex *in ,*out;
  fftw_plan p;
  
  for(int i=0;i<grid;i++){
    for(int j=0;j<grid;j++){
      for(int k=0;k<bins;k++){
	kinic[i][j][k]=x[i*grid*bins*2+j*bins*2+k*2]+x[i*grid*bins*2+j*bins*2+k*2+1]*I;
      }
    }
  }
  cout << "make_initial" <<endl;
  make_ini::make_kspace(knorm,k_y2m,k_vec);
  for(int M=0;M<5;M++){
    use_func::set_kturn(length,k_y2m[M]);
  }
  make_ini::make_rspace(r_y2m,r_vec);
  make_ini::read_transfer(filename,kt,transf);
  kt.erase(kt.begin());
  transf.erase(transf.begin());
  /* true a2mrを作製 */
  make_ini::make_iniperturb(knorm,kinit);
  make_trueQU::make_evoperturb(knorm,kinit,kt,transf,k_y2m,a2mkt);
  for(int M=0;M<5;M++){
    in = (fftw_complex *)a2mkt[M];
    out =(fftw_complex *)a2mrt[M];
    p = fftw_plan_dft_3d(grid,grid,grid,in,out, FFTW_BACKWARD, FFTW_ESTIMATE);
    fftw_execute(p);
    fftw_destroy_plan(p);
  }
  
  /* fitting a2mr作製 */
  make_trueQU::make_evoperturb(knorm,kinic,kt,transf,k_y2m,a2mkf);
  for(int M=0;M<5;M++){
    in = (fftw_complex *)a2mkf[M];
    out =(fftw_complex *)a2mrf[M];
    p = fftw_plan_dft_3d(grid,grid,grid,in,out, FFTW_BACKWARD, FFTW_ESTIMATE);
    fftw_execute(p);
    fftw_destroy_plan(p);
  }
  for(int i=0;i<grid;i++){
    for(int j=0;j<grid;j++){
      for(int k=0;k<grid;k++){
        for(int M=0;M<5;M++){
          a2mrt[M][i][j][k]=a2mrt[M][i][j][k]*(-4.0*M_PI)/sqrt(pow(2.0*M_PI*L,3.0));
	  a2mrf[M][i][j][k]=a2mrf[M][i][j][k]*(-4.0*M_PI)/sqrt(pow(2.0*M_PI*L,3.0));
        }
      }
    }
  }
  ofstream outputfilet("opre.dat",ios::app);
  for(int M=0;M<5;M++){
    outputfilet << "fit:" << endl;
    outputfilet << a2mrf[M][bins][bins][bins] << endl;
  }
}

void use_func::scalar_power(double k,double &Pk) noexcept{
  Pk = As*pow(k/kp,ns-1);
}
inline void use_func::spinY2m(int m,double phi,double theta,complex<double>(&spin)) noexcept{
  double spinre,spinim;
  switch(m){
  case -2:
    spinre = sqrt(5.0/64.0/M_PI)*(1.0+cos(theta))*(1.0+cos(theta))*cos(2.0*phi);
    spinim = -sqrt(5.0/64.0/M_PI)*(1.0+cos(theta))*(1.0+cos(theta))*sin(2.0*phi);
    break;
  case -1:
    spinre =  sqrt(5.0/16.0/M_PI)*sin(theta)*(1.0+cos(theta))*cos(phi);
    spinim =  -sqrt(5.0/16.0/M_PI)*sin(theta)*(1.0+cos(theta))*sin(phi);
    break;
  case 0:
    spinre = sqrt(15.0/32.0/M_PI)*sin(theta)*sin(theta);
    spinim = 0;
    break;
  case 1:
    spinre = sqrt(5.0/16.0/M_PI)*sin(theta)*(1.0-cos(theta))*cos(phi);
    spinim = sqrt(5.0/16.0/M_PI)*sin(theta)*(1.0-cos(theta))*sin(phi);
    break;
  case 2:
    spinre = sqrt(5.0/64.0/M_PI)*(1.0-cos(theta))*(1.0-cos(theta))*cos(2.0*phi);
    spinim = sqrt(5.0/64.0/M_PI)*(1.0-cos(theta))*(1.0-cos(theta))*sin(2.0*phi);
    break;
  }
  spin = spinre + spinim*I;
}

void use_func::set_kturn(int length[3],complex<double> (&k_box)[grid][grid][grid]) noexcept{
  /* check array of value */
  //  ofstream outputfile("klengar.dat");
  
  /* 空間上 */
  for(int i=1;i<length[0];i++){
    for(int j=1;j<length[1];j++){
      for(int k=1;k<length[2];k++){
	k_box[length[0]-i][length[1]-j][length[2]*2-k] = k_box[i][j][k];
      }
    }
  }    
  /* 平面上 */
  for(int i=1;i<length[0];i++){
    for(int j=1;j<length[1];j++){
      k_box[length[0]-i][length[1]-j][0]=k_box[i][j][0];
    }
    for(int j=1;j<length[2];j++){
      k_box[0][length[1]-i][length[2]*2-j]=k_box[0][i][j];
    }
  }
  for(int i=1;i<length[0];i++){
    for(int j=1;j<length[2];j++){
      k_box[length[1]-i][0][length[2]*2-j]=k_box[i][0][j];
    }
  }
  /* 軸上 */
  for(int i=1;i<length[0];i++){
    k_box[length[0]-i][0][0]=k_box[i][0][0];
  }
  for(int i=1;i<length[1];i++){
    k_box[0][length[1]-i][0]=k_box[0][i][0];
  }
  for(int i=1;i<length[2];i++){
    k_box[0][0][length[2]*2-i]=k_box[0][0][i];
  }
  //  outputfile.close();
}

void use_func::set_kturn_conj(int length[3],complex<double> (&k_box)[grid][grid][grid]) noexcept{
    /* 空間上 */
  for(int i=1;i<length[0];i++){
    for(int j=1;j<length[1];j++){
      for(int k=1;k<length[2];k++){
	k_box[length[0]-i][length[1]-j][length[2]*2-k] = conj(k_box[i][j][k]);
      }
    }
  }    
  /* 平面上 */
  for(int i=1;i<length[0];i++){
    for(int j=1;j<length[1];j++){
      k_box[length[0]-i][length[1]-j][0]=conj(k_box[i][j][0]);
    }
    for(int j=1;j<length[2];j++){
      k_box[0][length[1]-i][length[2]*2-j]=conj(k_box[0][i][j]);
    }
  }
  for(int i=1;i<length[0];i++){
    for(int j=1;j<length[2];j++){
      k_box[length[1]-i][0][length[2]*2-j]=conj(k_box[i][0][j]);
    }
  }
  /* 軸上 */
  for(int i=1;i<length[0];i++){
    k_box[length[0]-i][0][0]=conj(k_box[i][0][0]);
  }
  for(int i=1;i<length[1];i++){
    k_box[0][length[1]-i][0]=conj(k_box[0][i][0]);
  }
  for(int i=1;i<length[2];i++){
    k_box[0][0][length[2]*2-i]=conj(k_box[0][0][i]);
  }
}

void use_func::writedata_1d(char name[], vector<double> yvec,int length,vector<double> vec){
  ofstream outputfile(name);
  for(int i=0;i<length;i++){
    outputfile << vec[i] << " " << yvec[i] <<endl;
  }
  outputfile.close();
}

void use_func::writedata_3d(char name[], double *box,int length[3],double *vec){
  ofstream outputfile(name);
  for(int i=0;i<length[0];i++){
    for(int j=0;j<length[1];j++){
      for(int k=0;k<length[2];k++){
	outputfile << *(vec + i*3) << " ";
	outputfile << *(vec + j*3 +1) << " ";
	outputfile << *(vec + k*3 +2) << " ";
	outputfile << *(box+i*length[1]*length[2]+j*length[2]+k) <<endl;
      }
    }
  }
  outputfile.close();
}
 
void use_func::writedata_3d_c(char name[], complex<double> *box,int length[3],double *vec){
  ofstream outputfile(name);
  for(int i=0;i<length[0];i++){
    for(int j=0;j<length[1];j++){
      for(int k=0;k<length[2];k++){
	outputfile << *(vec + i*3) << " ";
	outputfile << *(vec + j*3 +1) << " ";
	outputfile << *(vec + k*3 +2) << " ";
	outputfile << real(*(box+i*length[1]*length[2]+j*length[2]+k)) << " ";
	outputfile << imag(*(box+i*length[1]*length[2]+j*length[2]+k)) << endl;
      }
    }
  }
  outputfile.close();
}

void use_func::writedata_4d_c(char name[], complex<double> *box,int length[3],double *vec){
  ofstream outputfile(name);
  for(int m=0;m<5;m++){
    //outputfile << "m=" << m-2 << endl;
    for(int i=0;i<length[0];i++){
      for(int j=0;j<length[1];j++){
	for(int k=0;k<length[2];k++){
	  outputfile << *(vec + i*3) << " ";
	  outputfile << *(vec + j*3 +1) << " ";
	  outputfile << *(vec + k*3 +2) << " ";
	  outputfile << real(*(box+m*length[0]*length[1]*length[2]
			  +i*length[1]*length[2]
			  +j*length[2]
			       +k)) << " ";
	  outputfile << imag(*(box+m*length[0]*length[1]*length[2]
			       +i*length[1]*length[2]
			       +j*length[2]
			       +k)) << endl;			       
	}
      }
    }
  }
  outputfile.close();
}

void use_func::read_4d_c(char name[],complex<double> *box){
  ifstream ifs(name);
  int i=0;
  double a[2];
  char str[256];
   if(ifs.fail()){
    cerr << "File do not exist.\n";
    exit(0);
  }
  while(ifs.getline(str,256)){
    sscanf(str,"%lf %lf",&a[0],&a[1]);
    *(box+i)=a[0]+a[1]*I;
    i++;
  }
  ifs.close();
}


double use_func::linearmy(double xi,double xi1,double yi,double yi1,double x){
  double y;
  y = yi + (yi1 - yi)* (x - xi)/(xi1 - xi);
  return y;
}

double use_func::linearinterpmy(double x,vector<double> arr_x,vector<double>arr_y, int arr_length){
  int i;
  for(i=1;i<arr_length;i++){
    if(arr_x[i]>=x) break;
  }
  return linearmy(arr_x[i-1],arr_x[i],arr_y[i-1],arr_y[i],x);
}


void make_ini::make_kspace (double (&knorm)[grid][grid][bins],complex<double> (&k_y2m)[5][grid][grid][grid],double (&k_vec)[grid][3]){
  double klist[grid];
  double kinner_x,kinner_y,kinner_z,k_theta,k_phi;
  // double k_vec[grid][3];
  for (int i=0;i<bins;i++){
    klist[i] = i*2.0*M_PI/L;
    klist[bins+i]= (-bins+i)*2.0*M_PI/L;
  }
  for (int i=0;i<grid;i++){
    for(int j=0;j<3;j++){
      k_vec[i][j] = klist[i];
    }
  }

  /* data check kphi ktheta */
  ofstream outputfile("kang.dat");
  /* 空間上 */
  for(int i=0;i<grid;i++){
    for(int j=0;j<grid;j++){
      for(int k=0;k<bins;k++){
	knorm[i][j][k]=sqrt(k_vec[i][0]*k_vec[i][0] +
			    k_vec[j][1]*k_vec[j][1] +
			    k_vec[k][2]*k_vec[k][2]);
	kinner_x = k_vec[i][0]/knorm[i][j][k];
	kinner_y = k_vec[j][1]/knorm[i][j][k];
	kinner_z = k_vec[k][2]/knorm[i][j][k];
	k_theta = acos(kinner_z);
	k_phi = acos(kinner_x/sin(k_theta));
	/* check field start */
	/*	outputfile << k_vec[i][0] << " ";
        outputfile << k_vec[j][1] << " ";
        outputfile << k_vec[k][2] << " ";
	outputfile << k_theta << " ";
	outputfile << k_phi << endl;*/
        /* check field end */
	for(int M=0;M<5;M++){
	  k_y2m[M][i][j][k]=conj(boost::math::spherical_harmonic(2,M-2,k_theta,k_phi));
	}
      }
    }
  }
  outputfile.close();

  /* y=0の場合、phiはxy座標系の偏角を表すため phi=0 or piに限定されるのでその補正 */
  ofstream outputfile2("kangcol.dat");
  for(int i=0;i<grid;i++){
    for(int j=0;j<1;j++){
      for(int k=0;k<bins;k++){
	knorm[i][j][k]=sqrt(k_vec[i][0]*k_vec[i][0] +
			    k_vec[j][1]*k_vec[j][1] +
			    k_vec[k][2]*k_vec[k][2]);
	kinner_x = k_vec[i][0]/knorm[i][j][k];
	kinner_y = k_vec[j][1]/knorm[i][j][k];
	kinner_z = k_vec[k][2]/knorm[i][j][k];
	k_theta = acos(kinner_z);
	if(kinner_x > 0.0){
	  k_phi = 0.0;
	}
	if(kinner_x < 0.0){
	  k_phi = 3.14159;
	}
	/* check field start */
	/*	outputfile2 << k_vec[i][0] << " ";
        outputfile2 << k_vec[j][1] << " ";
        outputfile2 << k_vec[k][2] << " ";
	outputfile2 << k_theta << " ";
	outputfile2 << k_phi << endl;*/
        /* check field end */
	for(int M=0;M<5;M++){
	  k_y2m[M][i][j][k]=conj(boost::math::spherical_harmonic(2,M-2,k_theta,k_phi));
	}
      }
    }
  }
  /* 原点の球面調和関数の値を0に設定 */
  for(int M=0;M<5;M++){
    k_y2m[M][0][0][0]=0.0;
    /* i,j,k=bins は共役を取る波数が存在しないため球面調和関数を0に設定*/
    for(int i=0;i<grid;i++){
      for(int j=0;j<grid;j++){
	k_y2m[M][i][j][bins]=0.0;
      }
      for(int j=0;j<bins;j++){
	k_y2m[M][i][bins][j]=0.0;
	k_y2m[M][bins][i][j]=0.0;
      }
    }	
  }
  outputfile2.close();
}

void make_ini::make_rspace(complex<double> (&r_y2m)[5][grid][grid][grid],double (&r_vec)[grid][3]) noexcept{
  double rlist[grid];
  double rnorm,rinner_x,rinner_y,rinner_z,r_theta,r_phi;
  complex<double> spin;
  //  double r_vec[grid][3];
  ofstream outputfile("rang.dat");
  for(int i=0;i<grid;i++){
    rlist[i]= L*i/(grid-1);
  }
  for(int i=0;i<grid;i++){
    for(int j=0;j<3;j++){
      r_vec[i][j] = rlist[i];
    }
  }
  for(int i=0;i<grid;i++){
    for(int j=0;j<grid;j++){
      for(int k=0;k<grid;k++){
	rnorm = sqrt(r_vec[i][0]*r_vec[i][0] +
		     r_vec[j][1]*r_vec[j][1] +
		     r_vec[k][2]*r_vec[k][2]);
	rinner_x = r_vec[i][0]/rnorm;
	rinner_y = r_vec[j][1]/rnorm;
	rinner_z = r_vec[k][2]/rnorm;
	r_theta = acos(rinner_z);
	r_phi = acos(rinner_x/sin(r_theta));
	/* check field start */
	outputfile << r_vec[i][0] << " ";
        outputfile << r_vec[j][1] << " ";
        outputfile << r_vec[k][2] << " ";
	outputfile << r_theta << " ";
	outputfile << r_phi << endl;
        /* check field end */
	for(int m=0;m<5;m++){
	  spinY2m(m-2,r_phi,r_theta,spin);
	  r_y2m[m][i][j][k] = spin;
	} 
      }
    }
  }
  for(int i=0;i<grid;i++){
    for(int j=0;j<1;j++){
      for(int k=0;k<grid;k++){
	rnorm = sqrt(r_vec[i][0]*r_vec[i][0] +
		     r_vec[j][1]*r_vec[j][1] +
		     r_vec[k][2]*r_vec[k][2]);
	rinner_x = r_vec[i][0]/rnorm;
	rinner_y = r_vec[j][1]/rnorm;
	rinner_z = r_vec[k][2]/rnorm;
	r_theta = acos(rinner_z);
	if(rinner_x >= 0.0){
	  r_phi = 0.0;
	}
	if(rinner_x < 0.0){
          r_phi = 3.14159;
        }
	for(int m=0;m<5;m++){
	  spinY2m(m-2,r_phi,r_theta,spin);
	  r_y2m[m][i][j][k] = spin;
	} 
      }
    }
  }

  for(int m=0;m<5;m++){
    r_y2m[m][0][0][0] = 0;
  }
	  
  outputfile.close();
}
void make_ini::make_iniperturb(double knorm[grid][grid][bins],complex<double> (&kinitial)[grid][grid][bins]){
  double Pk,random,ran_theta;
  random_device rnd;     // 非決定的な乱数生成器でシード生成機を生成
  mt19937 mt(rnd()); //  メルセンヌツイスターの32ビット版、引数は初期シード
  for(int i=0;i<grid;i++){
    for(int j=0;j<grid;j++){
      for(int k=0;k<bins;k++){
	random = (double)rand()/RAND_MAX;
	ran_theta= 2*M_PI*random;
	use_func::scalar_power(knorm[i][j][k],Pk);
	normal_distribution<double> norm(0,sqrt(Pk));
	kinitial[i][j][k]=complex<double>(norm(mt)*cos(ran_theta),norm(mt)*sin(ran_theta));
      }
    }
  }
}

void make_ini::set_harm_turn_r(complex<double> (&r_box)[grid][grid][grid]){
  for(int i=0;i<bins;i++){
    for(int j=0;j<bins;j++){
      for(int k=0;k<bins;k++){
	swap(r_box[i][j][k],r_box[bins+i][bins+j][bins+k]);
	swap(r_box[bins+i][j][k],r_box[i][bins+j][bins+k]);
	swap(r_box[i][bins+j][k],r_box[bins+i][j][bins+k]);
	swap(r_box[i][j][bins+k],r_box[bins+i][bins+j][k]);
      }
    }
  }
  /*  for(int i=0;i<grid;i++){
    for(int j=0;j<grid;j++){
      for(int k=0;k<bins;k++){
	swap(r_box[i][j][k],r_box[i][j][bins+k]);
      }
    }
    }*/
}

void make_ini::read_transfer(char name[],vector<double> &k,vector<double> &transf){
  ifstream ifs(name);
  int i=0;
  double a[14];
  char str[256];
  if(ifs.fail()){
    cerr << "File do not exist.\n";
    exit(0);
  }
  while(ifs.getline(str,256)){
    sscanf(str,"%lf   %lf   %lf   %lf   %lf   %lf    %lf   %lf   %lf   %lf   %lf   %lf   %lf   %lf"
	   ,&a[0],&a[1],&a[2],&a[3],&a[4],&a[5],&a[6],&a[7],&a[8],&a[9],&a[10],&a[11],&a[12],&a[13]);
    k.push_back(a[0]*ndh);
    transf.push_back(a[13]);
    //    cout << a[0] << " " << k[i] << endl;
    i++;
  }
}

void make_trueQU::make_evoperturb(double knorm[grid][grid][bins],complex<double> kinitial[grid][grid][bins],vector<double>  kt,vector<double>  transf,complex<double> k_y2m[5][grid][grid][grid],complex<double>(&a2mk)[5][grid][grid][grid]){
  complex<double> deltak;
  int length[3] = {grid,grid,bins};
  for(int i=0;i<grid;i++){
    for(int j=0;j<grid;j++){
      for(int k=0;k<bins;k++){
	deltak = sqrt(2.0*M_PI*M_PI*knorm[i][j][k])*linearinterpmy(knorm[i][j][k],kt,transf,200)*kinitial[i][j][k];
        for(int m=0;m<5;m++){
          a2mk[m][i][j][k]=deltak;
        }
      }
    }
  }
  for(int M=0;M<5;M++){
    use_func::set_kturn_conj(length,a2mk[M]);
  }
  for(int m=0;m<5;m++){
    for(int i=0;i<grid;i++){
      for(int j=0;j<grid;j++){
        for(int k=0;k<grid;k++){
          a2mk[m][i][j][k]=k_y2m[m][i][j][k]*a2mk[m][i][j][k];
        }
      }
    }
    a2mk[m][0][0][0] = 0.0;
  }
}


void make_trueQU::make_evoperturb_fast(double knorm[grid][grid][bins],const double *kinitial, vector<double>  kt,vector<double>  transf,complex<double> k_y2m[5][grid][grid][grid],complex<double>(&a2mk)[5][grid][grid][grid]){
  complex<double> deltak,kinic;
  double con;
  int length[3] = {grid,grid,bins};
  con = sqrt(2.0*M_PI*M_PI);
  for(int i=0;i<grid;i++){
    for(int j=0;j<grid;j++){
      for(int k=0;k<bins;k++){
	kinic = kinitial[i*grid*bins*2+j*bins*2+k*2]+kinitial[i*grid*bins*2+j*bins*2+k*2+1]*I;
	deltak = con*sqrt(knorm[i][j][k])*linearinterpmy(knorm[i][j][k],kt,transf,200)*kinic;
	for(int m=0;m<5;m++){
	  a2mk[m][i][j][k]=deltak;
	}
      }
    }
  }
  for(int M=0;M<5;M++){
    use_func::set_kturn_conj(length,a2mk[M]);
  }
  for(int m=0;m<5;m++){
    for(int i=0;i<grid;i++){
      for(int j=0;j<grid;j++){
        for(int k=0;k<grid;k++){
          a2mk[m][i][j][k]=k_y2m[m][i][j][k]*a2mk[m][i][j][k];
	}
      }
    }
    a2mk[m][0][0][0] = 0.0;
  }    
}

void make_trueQU::make_qu(complex<double> (&QU)[grid][grid][grid],complex<double> (&rbox)[5][grid][grid][grid],complex<double>r_y2m[5][grid][grid][grid]){
  complex<double> a2m;
  double con;
  con = (-4.0*M_PI)/sqrt(pow(2.0*M_PI*L,3.0));
  for(int i=0;i<grid;i++){
    for(int j=0;j<grid;j++){
      for(int k=0;k<grid;k++){
	QU[i][j][k] = 0.0;
	for(int M=0;M<5;M++){
	  a2m = rbox[M][i][j][k]*con;
	  QU[i][j][k] += a2m*r_y2m[M][i][j][k];
	}
      }
    }
  }
}
