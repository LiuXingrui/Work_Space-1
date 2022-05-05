#include <iostream>
#include<cstring>
#include<string>
#include<vector>
#include <fstream>
#include <math.h> 
#include <random>
#include"BP_basic_func1.h"
#include <sstream>

using namespace std;
#include <itpp/itbase.h>
#include <itpp/itcomm.h>
using namespace itpp;

//Hx=I prod H1,H2 prod I
//Hz=H2^T prod I, I prod H1^T
int main (int argc, char **argv){

  if (argc!=5&&argc!=6){
    cout<<"need 4 or 5 parameters, which are names for files stored H1, H2, Hx,Hz"<<endl;
    return 1;
  }

  
  string H1_file=argv[1];
  string H2_file=argv[2];
  string Hx_file=argv[3];
  string Hz_file=argv[4];  
  int n1,n2,r1,r2,k1,k2,n,k;
  
  bmat H1=read_matrix2(n1,r1,H1_file);
  bmat H2=read_matrix2(n2,r2,H2_file);
  
  bmat H1T=transpose(H1);
  bmat H2T=transpose(H2);
  
  // int thr_n=r1*n2+r2*n1;
  // int s1,s2;
  // s1=n1-r1;
  // s2=n2-r2;

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

  bmat Hx=merge_mat_hori2(Hx_left,Hx_right);
  bmat Hz=merge_mat_hori2(Hz_left,Hz_right);

  // n=Hx.rows();
  //if (n!=thr_n){cout<<"n!=thr_n, something is wrong"<<endl;}
  //test if they are orthogonal
  // if (n1<10&&n2<10){
  //    cout<<Hx*transpose(Hz)<<endl;
  //  }

  
  write_matrix2(Hx_file, Hx);
  write_matrix2(Hz_file,Hz);
  if (argc==6)
    {
      cout<<"H1= \n"<<H1<<"\n H1T=\n"<<H1T<<"\n H2=\n"<<H2<<"\nH2T=\n"<<H2T<<endl;
      cout<<"Hx=\n"<<Hx<<endl;
      cout<<"Hz=\n"<<Hz<<endl;

	}

  return 0;

  



}
