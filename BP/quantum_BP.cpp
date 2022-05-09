#include <iostream>
#include<cstring>
#include<string>
#include<vector>
#include <fstream>
#include <math.h> 
#include <random>
#include <sstream>


#ifndef basic
#define basic
#include"BP_basic_func1.h"
#endif

#ifndef quantum
#define quantum
#include"BP_quantum_func1.h"
#endif

using namespace std;
#include <itpp/itbase.h>
#include <itpp/itcomm.h>
using namespace itpp;

//GlobalRNG_reset (1);
int main(int argc, char **argv){
  GlobalRNG_randomize ();
  double pmax;
  double pmin;
  int num_of_cws;
  string file_name;
  string file_name2;
  string data_file;
  int lmax;
  
  if (argc!=8){cout<<" need 7 parameters: ./qBP  Hx_file Hz_file pmin pmax  num_of_cws lmax data_file"<<endl;return 1;}
  //get the parameters: 
   
  file_name=argv[1];
  file_name2=argv[2];
  data_file=argv[7];
 
  istringstream argv3( argv[3] );	
  if ( argv3 >> pmin){}
  else
    {
      cout<<"pmin should be a double"<<endl;
      return 1;
    }

  istringstream argv4( argv[4] );
  if ( argv4 >> pmax){}
  else
    {
      cout<<"pmax should be a double"<<endl;
      return 1;
    }
	   
	     
  istringstream argv5( argv[5] );
  if ( argv5 >> num_of_cws){}
  else
    {
      cout<<"num_of_cws should be an int"<<endl;
      return 1;
    }
	        
  istringstream argv6( argv[6] );
  if ( argv6 >> lmax){}
  else
    {
      cout<<"lmax should be an int"<<endl;
      return 1;
    }
		   		         	  
  double num_iter=0.0; //for calculate average iterations for e_x
  int num_of_suc_dec=0;// number of successfully decoded results
  int num_of_x_suc_dec=0;//number of Hx successfully decoded results
  bool Hx_suc=false;
  bool Hz_suc=false;



  //read the parity check matrix:
  int n1,n2,n,k,r1,r2,r;
  
 
  // bmat Hx2=read_matrix2( n1,r1, file_name);
  // bmat Hz2=read_matrix2( n2,r2, file_name2);
  //bmat zero_mat12(r1,r2);
  //cout<<Hx2*transpose(Hz2)<<endl;

  
  GF2mat Hx=read_matrix( n1,r1, file_name);
  GF2mat Hz=read_matrix( n2,r2, file_name2);
 
  //are the parity check matrices right?
  if (n1!=n2){cout<<"nx!=nz, the two matrices donot match"<<endl;return 1;}  
  n=n1;
  r=r1+r2;
  int rankx=GF2mat_rank(Hx);
  int rankz=GF2mat_rank(Hz);
  k=n-rankx-rankz;
  GF2mat zero_mat1(r1,r2);

 
  

  if((Hx*Hz.transpose())==zero_mat1) {}
  else{cout<<"Hx*Hz^T!=0, the two matrices donot match"<<endl;return 1;}
  //cout<<Hx*Hz.transpose()<<endl;
  
  GF2mat Gx=get_gen(Hx);
  //cout<<Gx*(Hx.transpose())<<endl;
  GF2mat Gz=get_gen(Hz);
  // cout<<Gz*(Hz.transpose())<<endl;
  
  nodes  xchecks[r1];//checks for Hx and Z errors
  nodes  zerrors[n];
  nodes  zchecks[r2];//checks for Hz and X errors
  nodes  xerrors[n];

  int E1=0;
  int E2=0;
  //find the neighbourhoods of all nodes:
  initialize_checks (Hx, xchecks,  E1);
  initialize_errors(Hx, zerrors);

  initialize_checks (Hz, zchecks,  E2);
  initialize_errors(Hz, xerrors);

  vec pxv(n);
  vec pzv(n);
     
    
  pro_dist( pmin,pmax, pxv);
  pro_dist( pmin,pmax, pzv);
    
  
  
  // int er=0;  //er is the number of one-bit errors (x for 1, z for 1, y for 2) that are wrong after decoding 
   for (int s=0;s<num_of_cws;s++)
    {
      Hx_suc= quan_decode(Hx,Gx, xchecks,zerrors,pxv,num_iter,lmax);
      // cout<<num_iter<<endl;
      if (Hx_suc==true)
	{
       
	  num_of_x_suc_dec++;
	  Hz_suc= quan_decode(Hz,Gz, zchecks,xerrors,pzv,num_iter,lmax);
	 
	  if (Hz_suc==true){num_of_suc_dec++;}
	}     
    }

   

   cout<<"for p in ( "<<pmin<<", "<<pmax<<"), there are total "<< num_of_suc_dec<<" successful decoding out of "<< num_of_cws<<" cws for a [["<<n<<", "<<k<<"]] code"<<endl;
   cout<<"average iterations:"<<endl;
   if((num_of_x_suc_dec+num_of_suc_dec)!=0)
     {
   cout<<num_iter/(num_of_x_suc_dec+num_of_suc_dec)<<endl;
     }
   else{cout<<0<<endl;}
   // cout<<"num of zero errors is about "<<pow(p,n)*num_of_cws<<endl;
 
   double midp=(pmax+pmin)/2;
      
   ofstream myfile;
   myfile.open (data_file,ios::app);
   myfile << n<<"  "<< 1.0*(num_of_cws-num_of_suc_dec)/num_of_cws<<"  "<<midp<<" "<<num_iter/(num_of_x_suc_dec+num_of_suc_dec)<<"  "<<num_of_suc_dec<<endl;
   myfile.close();

     
  return 0;

}
  
