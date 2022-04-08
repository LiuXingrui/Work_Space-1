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



int main(int argc, char **argv){
  double p;
  int num_of_cws;
  string file_name;
  string file_name2;
  int lmax;
  //get the parameters:
  
  if (argc >= 3)
    {
        file_name=argv[1];
	file_name2=argv[2];
	if (argc>=4)
	  {
        istringstream argv3( argv[3] );	
	if ( argv3 >> p){}
	else
	  {
	    cout<<"p should be a double"<<endl;
	    return 1;
	  }
	if (argc>=5)
	  {
	    istringstream argv4( argv[4] );
	    if ( argv4 >> num_of_cws){}
	    else
	      {
		cout<<"num_of_cws should be an int"<<endl;
		return 1;
	      }
	    if (argc>=6)
	      {
		istringstream argv5( argv[5] );
		  if ( argv5 >> lmax){}
		  else
		    {
		      cout<<"lmax should be an int"<<endl;
		      return 1;
		    }
		if (argc>=7)
		  {
		    cout<<"more than 5 paremeters, example: my_prog Hx_file Hz_file p num_of_cws lmax"<<endl;
		    return 1;
		  }
	
	      }
	    //default parameters:
	    else
	      {
		lmax=20;
	      }
	  }
	else
	  {
	    num_of_cws=10000;
	    lmax=20;
	  }

	  }
	else
	  {	 
	    p=0.01;
	    num_of_cws=10000;
	    lmax=20;
	  }
    }
  else
    {
       file_name="";
       file_name2="";
       p=0.01;
       num_of_cws=10000;
       lmax=20;
    }
  GlobalRNG_randomize ();

  
  double num_iter=0.0; //for calculate average iterations
  int num_of_suc_dec=0;// number of successfully decoded results
  int n_valid_cws=0;  // number of  decoded results that are cordwords


  //read the parity check matrix:
  int n1,n2,n,k,r1,r2,r;
  bmat Hx=read_matrix ( n1,r1, file_name);
  bmat Hz=read_matrix ( n2,r2, file_name);

  if (n1!=n2){cout<<"nx!=nz, the two matrices donot match"<<endl;return 1}
  n=n1;
  r=r1+r2;
  k=n-r;
  bmat zero_mat1(r1,r2);
  zero_mat1.zeros();
  if(Hx*transpose(Hz)!=zero_mat1) {cout<<"Hx*Hz^T!=0, the two matrices donot match"<<endl; return 1}
  
