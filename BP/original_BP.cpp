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
//compile: g++   `pkg-config --cflags itpp` -o my_prog BP.cpp  BP1.cpp `pkg-config --libs itpp`
// format : ./my_prog file_name error_propability number_of_cordwords max_iteration
//the format for the parity check file: the first row is n n-k,  then the sparse matrix. Notice
// the first element stored in the file is 1 rathar than 0.

int main(int argc, char **argv){
  double p;
  int num_of_cws;
  string file_name;
  int lmax;
  //get the parameters:
  
  if (argc >= 2)
    {
        file_name=argv[1];
	if (argc>=3)
	  {
        istringstream argv2( argv[2] );	
	if ( argv2 >> p){}
	else
	  {
	    cout<<"p should be a double"<<endl;
	    return 1;
	  }
	if (argc>=4)
	  {
	    istringstream argv3( argv[3] );
	    if ( argv3 >> num_of_cws){}
	    else
	      {
		cout<<"num_of_cws should be an int"<<endl;
		return 1;
	      }
	    if (argc>=5)
	      {
		istringstream argv4( argv[4] );
		  if ( argv4 >> lmax){}
		  else
		    {
		      cout<<"lmax should be an int"<<endl;
		      return 1;
		    }
		if (argc>=6)
		  {
		    cout<<"more than 4 paremeters, example: my_prog file_name p num_of_cws lmax"<<endl;
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
       file_name="M2";
       p=0.01;
       num_of_cws=10000;
       lmax=20;
    }

  GlobalRNG_randomize ();
 
  double num_iter=0.0; //for calculate average iterations
  int num_of_suc_dec=0;// number of successfully decoded results
  int n_valid_cws=0;  // number of  decoded results that are cordwords


  //read the parity check matrix:
  int n;
  int k;;
  int &v=n;
  int r;
  int &c=r;
  int temp;
  int row_ind=0;
  int num_zero_cws=0;


  bmat H=read_matrix ( n,r,  file_name);
  
  nodes  checks[c];
  nodes  errors[v];
  BSC bsc(p);
  int E=0;
  //find the neighbourhoods of all nodes:
  initialize_checks (H, checks,  E);
  initialize_errors(H, errors);

  int er=0;  //er is the number of bits that are wrong after decoding
  
  for (int s=0;s<num_of_cws;s++)
    {     
      num_of_suc_dec= num_of_suc_dec+cla_decode( v,c,H, checks, errors,bsc, num_iter,  lmax, er,p);
    }
  
  cout<<"for p="<<p<<", there are total "<< num_of_suc_dec<<" successful decoding out of "<< num_of_cws<<" cws for a ["<<v<<", "<<v-c<<"] code"<<endl;
  // cout<<"and there are "<< n_valid_cws<<" decoding results are codewords"<<endl;
  cout<<"bit error rate after decoding is ";
  cout<<er/(1.0*n*num_of_cws)<<endl;
  // cout<< "number of zero errors:    "<<num_zero_cws<<endl;
  cout<<"average iterations:"<<endl;
  cout<<num_iter/num_of_suc_dec<<endl;
  return 0;
}



    
