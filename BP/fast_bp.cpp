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

int main(){
  GlobalRNG_randomize ();
  double p=0.01;
  int lmax=20;
  int num_of_cws=10000;
  int num_of_suc_dec=0;
  
  int n;
  int k;
  int r;
  int temp;
  int col_ind=0;
  ifstream parity_check("LDPC_Mat2.txt");
  string line;
  
  getline(parity_check, line);
  istringstream iss(line);
  iss>>n>>r;
  k=n-r;
   bmat H(n-k,n);
  //  bmat H="1 1 0;0 1 1";
  // bmat H= "1 0 0 1 1 1 0;  0 1 0 1 1 0 1;0 0 1 0 1 1 1";
  int  v=H.cols();
  int  c=H.rows();
     
  //   /*

    while( getline(parity_check, line)){
    istringstream iss2(line);
    while (iss2>>temp){
      H(col_ind,temp-1)=1;
   
  }

    col_ind++;
  }



    //	  */
     

  parity_check.close();

  nodes  checks[c];
  nodes  errors[v];
  BSC bsc(p);
  int E=0;
  initialize_checks (H, checks,  E);
  initialize_errors(H, errors);

  int er=0;
  for (int s=0;s<num_of_cws;s++){
 
  
    bvec real_eT(v);
    bmat zero_mat1(c,1);
    bvec zero_vec(v);
    zero_vec.zeros();
    zero_mat1.zeros();
    
   
    real_eT=bsc(zero_vec);
 
   
    bmat real_e(v,1);
    bmat zero_mat2(v,1);
    zero_mat2.zeros();
    
    for (int q=0;q<v;q++){

      real_e(q,0)=real_eT(q);
    }
    
    bmat syndrome=H*real_e;
     if (real_eT==zero_vec){
      num_of_suc_dec++;
      // cout<<"success!  error=0"<<endl;
      continue;
      
    }
     
  
     
    

    mat mcv(c,v);
    mat mvc(c,v);
    initialize_massages( mcv,mvc, H);
    bmat output_e(v,1);
    output_e.zeros();
      
    for (int l=1;l<=lmax;l++)
      {
      //  cout<<l<<endl;
      s_update(checks,errors, mcv,mvc,syndrome,p, c, v,output_e);
   
	if(output_e==real_e)
	  {
	  // cout<<"success! iteration number="<<l<<endl;
	  num_of_suc_dec++;
	  //  cout<<num_of_suc_dec<<endl;
	  break;
	  }
	  
      }
      
	    
      }
   


  
  
  cout<<"for p="<<p<<", there are total "<< num_of_suc_dec<<" successful decoding out of "<< num_of_cws<<" cws for a ["<<v<<", "<<v-c<<"] code"<<endl;
 





 
}



    
