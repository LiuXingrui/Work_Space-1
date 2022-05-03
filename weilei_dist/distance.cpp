#include <iostream>
#include "dist.h"
#include "lib.h"
#include "mm_read.h"
#include "mm_write.h"

#ifndef basic
#define basic
#include"BP_basic_func1.h"
#endif

#ifndef quantum
#define quantum
#include"BP_quantum_func1.h"
#endif

using namespace std;
using namespace itpp;

int main(int argc,char** argv){

  string filex=argv[1];
  string filez=argv[2];
  int nx,rx,nz,rz;
  bmat Hx0=read_matrix (nx,rx, filex);
  bmat Hz0=read_matrix (nz,rz, filez);
  bmat zero_mat1(rx,rz);
  zero_mat1.zeros();
  

  if((Hx0*transpose(Hz0))!=zero_mat1) {cout<<"Hx*Hz^T!=0, the two matrices donot match"<<endl;return 1;}


  GF2mat Hx(Hx0);
  GF2mat Hz(Hz0);
  int dx=common::quantum_dist_v2(Hx, Hz, 0);
  int dz=common::quantum_dist_v2(Hx, Hz, 1);

  cout<<"n="<<nx<<"  dx<="<<dx<<"  dz<="<<dz<<endl;
  


  return 0;
}
