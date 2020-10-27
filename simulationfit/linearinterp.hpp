//#ifndef _LINEARINTERP_HPP_
//#define _LINEARINTERP_HPP_
#include<iostream>

using namespace std;

static double linearmy(double xi,double xi1,double yi,double yi1,double x){
  double y;
  y = yi + (yi1 - yi)* (x - xi)/(xi1 - xi);
  return y;
}

static double linearinterpmy(double x,vector<double> arr_x,vector<double>arr_y, int arr_length){
  int i;
  for(i=1;i<arr_length;i++){
    if(arr_x[i]>=x) break;
  }
  return linearmy(arr_x[i-1],arr_x[i],arr_y[i-1],arr_y[i],x);
}

//#endif // _LINEARINTERP_HPP_
