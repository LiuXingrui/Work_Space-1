#include <iostream>
#include<cstring>
#include<string>
#include<vector>
#include <fstream>
#include <math.h> 
#include <random>
#include"BP_header.h"
#include"BP_quantum_header.h"
#include <sstream>

using namespace std;
#include <itpp/itbase.h>
#include <itpp/itcomm.h>
using namespace itpp;


int main(){
  int n1,r1,n2,r2;
  string name1="testmat1";
  string name2="testmat3";
  string name3="merge_result1";
   
  bmat H1=read_matrix(n1,r1,name1 );
  bmat H2=read_matrix(n2,r2, name2);
  bmat m= merge_mat_hori(H1,H2);
  write_matrix(name3, m);

  return 0;
}
