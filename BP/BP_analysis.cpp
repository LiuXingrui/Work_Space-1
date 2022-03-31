#include <iostream>
#include<cstring>
#include<string>
#include<vector>
#include <fstream>
#include <math.h> 
#include <random>
#include"BP_header.h"

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
		lmax=6;
	      }
	  }
	else
	  {
	    num_of_cws=10;
	    lmax=6;
	  }

	  }
	else
	  {	 
	    p=0.05;
	    num_of_cws=10;
	    lmax=6;
	  }
    }
  else
    {
       file_name="cycleMat1";
       p=0.05;
       num_of_cws=10;
       lmax=6;
    }
  GlobalRNG_randomize ();
 

  //read the parity check matrix:
  int n;
  int k;;
  int &v=n;
  int r;
  int &c=r;
  int temp;
  int row_ind=0;


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
  
  if (n>7||r>7) {cout<<"the size of the matrix is too large to analyze, the dimensions should less than 8"<<endl; parity_check.close();return 1;}
  
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
  int E=0;
  //find the neighbourhoods of all nodes:
  initialize_checks (H, checks,  E);
  initialize_errors(H, errors);

//only for 1 bit error -0331
    cout<<"lmax="<<lmax<<endl;
  for (int s=0;s<n;s++)
    {     
      bvec real_eT(v);    //the transposed error vector, which is a row vector. 
      bmat zero_mat1(c,1); 
      bvec zero_vec(v);
      zero_vec.zeros();
      zero_mat1.zeros();
      
      real_eT.zeros();
      real_eT(s)=1;
      cout<<"the error is \n"<<real_eT<<endl;
      bmat real_e(v,1);  //the error vector which is a column vector,
      bmat zero_mat2(v,1);
      zero_mat2.zeros();
      
      for (int q=0;q<v;q++)
	{
	  real_e(q,0)=real_eT(q);
	}
      bmat syndrome=H*real_e;
      cout<<"\n the syndrome is \n"<<syndrome<<endl;

        
      mat mcv(c,v);   // the messages from check nodes to variable nodes, which are p0/p1
      mat mvc(c,v);   // the messages from variable nodes to check nodes
      initialize_massages( mcv,mvc, H); //initialize to all-1 matrix
      bmat output_e(v,1);
      output_e.zeros();
      cout<<"\n for parallel:"<<endl;
    
      for (int l=1;l<=lmax;l++)
	{
	  cout<<"\n iteration "<<l<<": \n"<<endl;
	  p_update_a(checks,errors, mcv,mvc,syndrome,p, c, v,output_e);
	  cout<<"\n mcv: \n"<<mcv<<endl;
	  cout<<"\n  mvc: \n"<<mvc<<endl;
	  
	  
	  if (H*output_e==syndrome)
	    {
	      cout<<"get the right syndrome \n output_e: \n"<<output_e<<endl;
	      break;	  
	    }
	  if(l==lmax){ cout<<"failed \n output_e: \n"<<output_e<<endl;}

	}

      //for sequential
      initialize_massages( mcv,mvc, H); //initialize to all-1 matrix    
      output_e.zeros();
      cout<<"\n for sequential:"<<endl;
      for (int l=1;l<=lmax;l++)
	{
	  cout<<"\n iteration "<<l<<": \n"<<endl;
	  s_update_a(checks,errors, mcv,mvc,syndrome,p, c, v,output_e);
	  cout<<"\n mcv: \n"<<mcv<<endl;
	  cout<<"\n  mvc: \n"<<mvc<<endl;
	  
	  
	  if (H*output_e==syndrome)
	    {
	      cout<<"\n get the right syndrome output_e: \n"<<output_e<<endl;
	      break;	  
	    }
	  if(l==lmax){  cout<<"\n output_e: \n"<<output_e<<endl;cout<<"\n failed "<<endl;}

	}
    }
  

  return 0;
}



    
