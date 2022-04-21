#include <iostream>
#include<cstring>
#include<string>
#include<vector>
#include <fstream>
#include <math.h> 
#include <random>

#ifndef basic
#define basic
#include"BP_basic_func1.h"
#endif

#include <sstream>

using namespace std;
#include <itpp/itbase.h>
#include <itpp/itcomm.h>
using namespace itpp;


int main(int argc, char **argv){
 
  string file_name;
  string file_name2;
  //get the parameters:
  
  if (argc != 3) { cout<<"need exactly 2 parameters: <H> <H^T>"<<endl;return 1;}
  file_name=argv[1];
  file_name2=argv[2];
  int n,r;
  
  bmat H= read_matrix ( n,r, file_name);
  //cout<<H<<endl;
  
  bmat HT=transpose(H);
  // cout<<HT<<endl;
     
  write_matrix(file_name2, HT);


  return 0;

}
