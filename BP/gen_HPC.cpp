#include <iostream>
#include<cstring>
#include<string>
#include<vector>
#include <fstream>
#include <math.h> 
#include <random>
#include"BP_header.h"
#include <sstream>

using namespace std;
#include <itpp/itbase.h>
#include <itpp/itcomm.h>
using namespace itpp;


int main (int argc, char **argv){

  if (argc!=5){
    cout<<"need exactly 4 parameters, which are names for files stored H1, H2, Hx,Hz"<<endl;
    return 1;
  }

  
  string H1_file=argv[1];
  string H2_file=argv[2];
  string Hx_file=argv[3];
  string Hz_file=argv[4];  
  int n1,n2,r1,r2,k1,k2,n,k;
  
  bmat H1=read_matrix(n1,r1,H1_file);
  bmat H2=read_matrix(n2,r2,H2_file);
  
  bmat H1T=transpose(H1);
  bmat H2T=transpose(H2);
  
  k1=n1-r1;
  k2=n2-r2;
  //Identity matrices:
  bmat I_r1=eye_b (r1);
  bmat I_r2=eye_b (r2);
  bmat I_n1=eye_b (n1);
  bmat I_n2=eye_b (n2);

  //kronecker product:
  bmat Hx_left=kron(I_r2, H1);
  bmat Hx_right=kron(H2,I_r1);
  bmat Hz_left=kron(H2T,I_n1);
  bmat Hz_right=kron(I_n2,H1T);

  bmat Hx=merge_mat_hori(Hx_left,Hx_right);
  bmat Hz=merge_mat_hori(Hz_left,Hz_right);

  //test if they are orthogonal
  if (n1<10&&n2<10){
    cout<<Hx*transpose(Hz)<<endl;
  }

  write_matrix(Hx_file, Hx);
  write_matrix(Hz_file,Hz);


  return 0;

  



}
