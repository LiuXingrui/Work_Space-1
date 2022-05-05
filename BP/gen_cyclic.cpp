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
  int n;
  bvec h_poly(argc-3);
  string file_name;
  //get the parameters:
  
  if (argc >= 2)
    {
      file_name=argv[1];	
	if (argc>=3)
	  {
        istringstream argv2( argv[2] );	
	if ( argv2 >> n){}
	else
	  {
	    cout<<"n should be an int"<<endl;
	    return 1;
	  }
        for (int i=3;i<argc;i++)
	  {
	    istringstream argvi(argv[i]);
	    if(argvi>>h_poly(i-3)){}
	    else {cout<<"coefficients of check poly should be binary"<<endl;return 1;}
	  }
	  }
    }
  else
    {
      file_name="n6_cyc11";
      n=6;
      h_poly=" 1 1";
    }


  if(h_poly.size()>n){cout<<"error! check polynomial is longer than n!"<<endl;return 1;}

  int k=h_poly.size()-1;
  bmat H(n-k+1,n);
  H.zeros();

  for (int i=0;i<=n-k;i++)
    {
      int i2=0;
      for (int j=i;j<i+k+1;j++)
	{
	  H(i,j%n)=h_poly(i2);
	  i2++;
	}
    }


  write_matrix2(file_name, H);


  return 0;

}
