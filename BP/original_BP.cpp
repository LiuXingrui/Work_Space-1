#include <iostream>
#include<cstring>
#include<string>
#include<vector>
#include <fstream>
#include <math.h> 
#include <random>
#include"BP1.h"

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
	    num_of_cws=1000;
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

  ifstream parity_check(file_name);
  string line;
  
  getline(parity_check, line);
  istringstream iss(line);
  
  if(  iss>>n>>r) {}
  else
    {
      cout<<"did not open the right parity check file"<<endl;
      return 1;
    }
  k=n-r;
  bmat H(n-k,n);

  while( getline(parity_check, line))
    {
      istringstream iss2(line);
      while (iss2>>temp){
	if (temp>=1&&temp<=n&&row_ind<r)
	  {
	    H(row_ind,temp-1)=1;
	  }
	else
	  {
	    cout<<"the format of the parity check is wrong, the first element is 1 rathar than 0"<<endl;
	    return 1;
	  }
	
      }
		  row_ind++;
    }
  parity_check.close();
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
      bvec real_eT(v);    //the transposed error vector, which is a row vector. 
      bmat zero_mat1(c,1); 
      bvec zero_vec(v);
      zero_vec.zeros();
      zero_mat1.zeros();
      
      real_eT=bsc(zero_vec);
      bmat real_e(v,1);  //the error vector which is a column vector,
      bmat zero_mat2(v,1);
      zero_mat2.zeros();
      
      for (int q=0;q<v;q++)
	{
	  real_e(q,0)=real_eT(q);
	}
      bmat syndrome=H*real_e;

      //is the syndrome a zero vector?
      if (syndrome==zero_mat1)
	{
	  n_valid_cws++;
	  if (real_eT==zero_vec)
	    {
	      num_of_suc_dec++;
	      num_zero_cws++;
	      // cout<<"success!  error=0"<<endl;     
	    }
	  else
	    {
	      er=er+ distance(zero_mat2, real_e, n);  
	      // cout<<"failure! error= some codeword"<<endl;
	    } 
	  continue;   
	}
   
      mat mcv(c,v);   // the messages from check nodes to variable nodes, which are p0/p1
      mat mvc(c,v);   // the messages from variable nodes to check nodes
      initialize_massages( mcv,mvc, H); //initialize to all-1 matrix
      bmat output_e(v,1);
      output_e.zeros();
      
      for (int l=1;l<=lmax;l++)
	{
	  //  cout<<l<<endl;
	  s_update(checks,errors, mcv,mvc,syndrome,p, c, v,output_e);
	  
	  if (H*output_e==syndrome)
	    {
	      num_iter=num_iter+l;
	      n_valid_cws++;
	      if(output_e==real_e)
		{
		  // cout<<"success! iteration number="<<l<<endl;
		  num_of_suc_dec++;
		  //  cout<<num_of_suc_dec<<endl;
		}
	      else
		{
		  er=er+ distance(output_e, real_e, n);	        
		}
	      break;	  
	    }
	  if(l==lmax){
	    num_iter=num_iter+l;
	    er=er+ distance(output_e, real_e, n);	 
	  }
	}
    }
  
  cout<<"for p="<<p<<", there are total "<< num_of_suc_dec<<" successful decoding out of "<< num_of_cws<<" cws for a ["<<v<<", "<<v-c<<"] code"<<endl;
  cout<<"and there are "<< n_valid_cws<<" decoding results are codewords"<<endl;
  cout<<"bit error rate after decoding is ";
  cout<<er/(1.0*n*num_of_cws)<<endl;
  cout<< "number of zero errors:    "<<num_zero_cws<<endl;
  cout<<"average iterations:"<<endl;
  cout<<num_iter/num_of_cws<<endl;
  return 0;
}



    
