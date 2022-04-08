#include <iostream>
#include<cstring>
#include<string>
#include<vector>
#include <fstream>
#include <math.h> 
#include <random>
#include"BP_header.h"
#include"BP_quantum_header.h"

using namespace std;
#include <itpp/itbase.h>
#include <itpp/itcomm.h>
using namespace itpp;



int main(){
 GlobalRNG_randomize ();
  bvec cw(10);
  double p=0.1;
  vec pv(10);
  int num_of_1=0;
  // pro_dist(p,pv);
  // cout<<"pv is "<<pv<<endl;

  for (int i=0;i<10;i++)
    {
      pv(i)=p;
    }
  

  for (int i=0;i<100;i++)
    {
      cw.zeros();
      error_channel(cw,pv);
      for (int j=0;j<10;j++)
	{
	  if (cw(j)==1){num_of_1++;}
	}
    }

  cout<<"num_of_1= "<<num_of_1<<endl;

     
  



  return 0;
  
}
