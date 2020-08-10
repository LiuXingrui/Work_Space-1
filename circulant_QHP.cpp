#include <NTL/ZZ.h> 
#include <NTL/vector.h>
#include <iostream>
#include<cstring>
#include<string>
#include<vector>
#include <fstream>
#include <NTL/ZZXFactoring.h>
#include <NTL/ZZ_pXFactoring.h>
#include <NTL/ZZ_pEX.h>


using std::default_random_engine generator; 
using std::uniform_int_distribution<int>;
using std::cin;
using std::cout; 
using std::vector;
using namespace NTL;
vector<unsigned> rand_vec (const int n){
  static default_random_engine e;
  static uniform_int_distribution<int> u(0,1);
  vector<int> ret;
  for(size_t i=0;i<n;++i){
    ret.push_back(u(e));

  }
  return ret;
}


int main(int argc, char **argv){
vector<int> check_poly_a;
vector<int> check_poly_b;
 cout<<"test";

return 0;
}
