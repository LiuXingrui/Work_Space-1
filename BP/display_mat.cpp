#include <iostream>
#include<cstring>
#include<string>
#include<vector>
#include <fstream>
#include <math.h> 
#include <random>
#include <sstream>

#ifndef basic
#define basic
#include"BP_basic_func1.h"
#endif

#ifndef quantum
#define quantum
#include"BP_quantum_func1.h"
#endif

using namespace std;
#include <itpp/itbase.h>
#include <itpp/itcomm.h>
using namespace itpp;



int main(int argc, char **argv){
  if (argc!=2){cout<<"need exactly one parameter"<<endl;return 1;}

  string file_name=argv[1];
  int n,r;
  bmat H=read_matrix ( n,r, file_name);
  cout<<n<<"  "<<n-bmat_rank(H)<<endl;
  cout<<H<<endl;
  return 0;
}

