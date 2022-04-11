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


void error_channel(bvec &cw, const vec &p){
  int temp=1e8;
  int temp2;
  bin one=1;
  if (cw.size()!=p.size())
    {cout<<"the size of p and cw do not match"<<endl;}

  else
    {
    for (int i=0;i<cw.size();i++)
      {    
	temp2=randi(1,temp);
	if(temp2<p[i]*temp)
	  {
	    cw[i]=cw[i]+one;
	  }	
      }
  }  
}


void pro_dist(double p, vec& pv){

  double minp=0.5*p;
  double maxp=1.5*p;
  int pvsize=pv.size();
  double temp;
  int temp2=1e8;

  for (int i=0;i<pvsize;i++)
    {
      temp=p*randi(1,temp2)/temp2;
      pv(i)=minp+temp;
    }
}

bool  quan_decode(int v,int c,const bmat &H,const bmat &H2,const nodes checks[],const nodes errors[],const vec &pv,double& num_iter, int lmax,int&num){

  bvec real_eT(v);    //the transposed error vector, which is a row vector.
  error_channel(real_eT, pv);
   

  //if no error, break
  bvec zero_vec(v);
  zero_vec.zeros();
  if (real_eT==zero_vec)
    {
	 
      return true;
    }
  
  bmat zero_mat1(c,1);
  bmat zero_mat2(v,1);
  zero_mat1.zeros();
  zero_mat2.zeros();
  bmat real_e(v,1);  //the error vector which is a column vector,
     
  for (int q=0;q<v;q++)
    {
      real_e(q,0)=real_eT(q);
    }
  bmat syndrome=H*real_e;

  //is the syndrome a zero vector?
  if (syndrome==zero_mat1)
	{
	  // is e a stablizer?
	  if (Q_inspan(real_eT,H2,c){return true;}        
	  else  
	    {
	      return false;
	      //er=er+ distance(zero_mat2, real_e, n);  
	      // cout<<"failure! error= some codeword"<<endl;
	    }  
	}
   
      mat mcv(c,v);   // the messages from check nodes to variable nodes, which are p0/p1
      mat mvc(c,v);   // the messages from variable nodes to check nodes
      initialize_massages( mcv,mvc, H); //initialize to all-1 matrix
      bmat output_e(v,1);
      output_e.zeros();
      
      for (int l=1;l<=lmax;l++)
	{
	  //  cout<<l<<endl;
	  quan_s_update(checks,errors, mcv,mvc,syndrome,pv, c, v,output_e);
	  
	  if (H*output_e==syndrome)
	    {
	      num_iter=num_iter+l;
	     
	      if(quan_check(output_e,real_e,H2,c,v))
		{
		  return true;
		  // cout<<"success! iteration number="<<l<<endl;
		  // num_of_suc_dec++;
		  //  cout<<num_of_suc_dec<<endl;
		}
	      else
		{
		  return false;
		  // er=er+ distance(output_e, real_e, n);	        
		}	    	  
	    }
	  if(l==lmax)
	    {
	    num_iter=num_iter+l;
	    return false;
	    // er=er+ distance(output_e, real_e, n);	 
	    }
	}
 }


  void quan_s_update(const nodes checks[],const nodes errors[],mat &mcv,mat& mvc,const bmat& syndrome,const vec &pv,int c, int v,  bmat& output_e){
  
    double ipr;

  for (int i=0;i<c;i++)
    {
      int ci_degree=checks[i].degree;
      
      	//update all v-to-ci massages first:
    for (int j=0;j<ci_degree;j++)
      {            
	int vnode=(checks[i].neighbors)(j);
        ipr=(1-pv[vnode])/pv[vnode];
	update_vj_to_ci(checks, errors,mcv, mvc,vnode,i, ipr);	
      }
    }
    //update c-to-vj massages:
  for (int j=0;j<v;j++)
    {
   int vj_degree=errors[j].degree;
   double  final_pr==(1-pv[j])/pv[j];

   for (int i=0;i<vj_degree;i++)
     {
       int cnode=(errors[j].neighbors)(i);
       update_ci_to_vj( checks, errors,mcv, mvc,cnode,j,syndrome(cnode,0));

       final_pr=final_pr*mcv(cnode,j);
      }   
       //  cout<<j<<"   "<<final_pr<<endl;
       output_e(j,0)=final_pr<1? 1:0;      
    }  
}

  //check if real_e is a stabilizer
bool Q_inspan(bvec &real_eT,bmat &H2,int r){
  int ori_rank=rank(H2);
  H2.ins_row (r, real_eT);
  int i=rank(H2);// the argument is const bmat&, so I guess itpp will creat a new matrix to do the Gaussian elimination, which will waste some time.
  H2.del_row ( r);
  
  if (i==ori_rank){return true;}
  else if (i==ori_ran+1){return false;}
  else {cout<<"some wrong happened with the Q_inspan func"<<endl;return false;}
  

}

//check if real_e-output_e is a stabilizer
 bool quan_check(bmat &output_e,bmat &real_e,bmat &H2,int r,int n){
   bvec difference(n);
   for (i=0;i<n;i++)
     {
       difference(i)=output_e(i,0)-real_e(i,0);
     }
   return Q_inspan(difference,H2,r);


  
 }
