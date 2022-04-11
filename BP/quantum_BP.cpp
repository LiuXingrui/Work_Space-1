#include <iostream>
#include<cstring>
#include<string>
#include<vector>
#include <fstream>
#include <math.h> 
#include <random>
#include"BP_basic_func1.h"
#include"BP_quantum_func1.h"
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

  
  double num_iter_x=0.0; //for calculate average iterations for e_x
  double num_iter_z=0.0;//for calculate average iterations for e_x
  int num_of_suc_dec=0;// number of successfully decoded results
  bool Hx_suc-false;
  int n_valid_cws=0;  // number of  decoded results that are cordwords


  //read the parity check matrix:
  int n1,n2,n,k,r1,r2,r;
  int&rx=r1;
  int&rz=r2;
  bmat Hx=read_matrix ( n1,r1, file_name);
  bmat Hz=read_matrix ( n2,r2, file_name);

  if (n1!=n2){cout<<"nx!=nz, the two matrices donot match"<<endl;return 1}
  n=n1;
  r=r1+r2;
  k=n-r;
  bmat zero_mat1(r1,r2);
  zero_mat1.zeros();
  if(Hx*transpose(Hz)!=zero_mat1) {cout<<"Hx*Hz^T!=0, the two matrices donot match"<<endl; return 1}

  
  nodes  xchecks[r1];//checks for Hx and Z errors
  nodes  zerrors[n];
  nodes  zchecks[r2];//checks for Hz and X errors
  nodes  xerrors[n];

  int E1=0;
  int E2=0;
  //find the neighbourhoods of all nodes:
  initialize_checks (Hx, xchecks,  E1);
  initialize_errors(Hx, zerrors);

  initialize_checks (Hz, zchecks,  E1);
  initialize_errors(Hz, xerrors);

  vec pxv(n);
  vec pzv(n);
  pro_dist( p, pxv);
  pro_dist( p, pzv);

  int er=0;  //er is the number of one-bit errors (x for 1, z for 1, y for 2) that are wrong after decoding 

   for (int s=0;s<num_of_cws;s++)
    {     
      bvec real_exT(n);    //the transposed error vector, which is a row vector.
      bvec real_ezT(n);
      real_exT.zeros();
      real_ezT.zeros();
      error_channel(real_exT, pxv);
      error_channel(real_exT, pzv);

      //if no error, break
      bvec zero_vec(n);
      zero_vec.zeros();
      if (real_exT==zero_vec&&real_ezT==zero_vec)
	{
	  num_of_suc_dec++;
	  break;
	}
      
      bmat real_ex(n,1);  //the error vector which is a column vector,
      bmat real_ez(n,1);
      bmat zero_mat1=(r1,1);
      bmat zero_mat2=(r2,1);
      zero_mat1.zeros();
      zero_mat2.zeros();
      
      for (int q=0;q<n;q++)
	{
	  real_ex(q,0)=real_exT(q);
	  real_ez(q,0)=real_ezT(q);
	}
      bmat x_syndrome=Hz*real_ex;
      bmat z_syndrome=Hx*real_ez;

      //is the syndrome a zero vector?
      if (x_syndrome==zero_mat2)
	{
	  n_valid_cws++;
	  if (real_eT==zero_vec)
	    {
	      num_of_suc_dec++;
	      
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
	  p_update(checks,errors, mcv,mvc,syndrome,p, c, v,output_e);
	  
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





  return 0;

}
  
