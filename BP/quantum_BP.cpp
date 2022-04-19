#include <iostream>
#include<cstring>
#include<string>
#include<vector>
#include <fstream>
#include <math.h> 
#include <random>
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
	    num_of_cws=1000;
	    lmax=20;
	  }

	  }
	else
	  {	 
	    p=0.01;
	    num_of_cws=1000;
	    lmax=20;
	  }
    }
  else
    {
       file_name="testHx1";
       file_name2="testHz1";
       p=0.01;
       num_of_cws=1000;
       lmax=20;
    }
  GlobalRNG_randomize ();

  
  double num_iter=0.0; //for calculate average iterations for e_x
  int num_of_suc_dec=0;// number of successfully decoded results
  bool Hx_suc=false;
  bool Hz_suc=false;



  //read the parity check matrix:
  int n1,n2,n,k,r1,r2,r;
  int&rx=r1;
  int&rz=r2;

  bmat Hx=read_matrix ( n1,r1, file_name);
  bmat Hz=read_matrix ( n2,r2, file_name2);
    cout<<"rx="<<rx<<endl;
  cout<<"rz="<<rz<<endl;
  //are the parity check matrices right?
  if (n1!=n2){cout<<"nx!=nz, the two matrices donot match"<<endl;return 1;}  
  n=n1;
  r=r1+r2;
  int rankx=bmat_rank(Hx);
  int rankz=bmat_rank(Hz);
  // cout<<"rankx:"<<rankx<<endl;
  //cout<<"rankz:"<<rankz<<endl;
  k=n-rankx-rankz;
  bmat zero_mat1(r1,r2);
  zero_mat1.zeros();
  
  //  cout<<"Hx:"<<Hx<<endl;
  // cout<<"Hz:"<<Hz<<endl;
  // cout<<Hx*transpose(Hz)<<endl;
  if((Hx*transpose(Hz))!=zero_mat1) {cout<<"Hx*Hz^T!=0, the two matrices donot match"<<endl;return 1;}

  
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
      Hx_suc= quan_decode(Hx,Hz, xchecks,zerrors,pxv,num_iter,lmax,rankz);
     
      if (Hx_suc==true)
	{
	 
	  Hz_suc= quan_decode(Hz,Hx, zchecks,xerrors,pzv,num_iter,lmax,rankx);
	 
	  if (Hz_suc==true){num_of_suc_dec++;}
	}     
    }

   

   cout<<"for p around "<<p<<", there are total "<< num_of_suc_dec<<" successful decoding out of "<< num_of_cws<<" cws for a ["<<n<<", "<<k<<"] code"<<endl;
   cout<<"average iterations:"<<endl;
   cout<<num_iter/num_of_suc_dec/2<<endl;
   // cout<<"num of zero errors is about "<<pow(p,n)*num_of_cws<<endl;

  return 0;

}
  
