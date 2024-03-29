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
  GF2mat  T;
  GF2mat  U;
  ivec    P; 

void error_channel(bvec &cw, const vec &p){
 
  double temp2;
  bin one=1;
  if (cw.size()!=p.size())
    {cout<<"the size of p and cw do not match"<<endl;}

  else
    {
    for (int i=0;i<cw.size();i++)
      {    
	temp2=randu();
	if(temp2<p[i])
	  {
	    cw[i]=cw[i]+one;
	  }	
      }
  }
  
}


void pro_dist(double pmin,double pmax, vec& pv){

 
  double pdiff=pmax-pmin;
  int pvsize=pv.size();
  double temp;
 

  for (int i=0;i<pvsize;i++)
    {
      temp=pdiff*randu();
      pv(i)=pmin+temp;
    }
}

bool  quan_decode(bmat &H, bmat &H2,const nodes checks[],const nodes errors[],const vec &pv,double& num_iter, int lmax,int rank2){
  int v=H.cols();
  int c=H.rows();
  // int r2=H2.rows();
  bvec real_eT(v);    //the transposed error vector, which is a row vector.
  real_eT.zeros();
  error_channel(real_eT, pv);
  // cout<<pv<<endl;

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
	  if (Q_inspan(real_eT,H2,rank2))
	    {
	      
	      return true;
	    }        
	  else  
	    {
	    
	      return false;
	      //er=er+ distance(zero_mat2, real_e, n);  
	      // cout<<"failure! real_error is in NS"<<endl;
	    }  
	}
  
      mat mcv(c,v);   // the messages from check nodes to variable nodes, which are p0/p1
      mat mvc(c,v);   // the messages from variable nodes to check nodes
      initialize_massages( mcv,mvc, H); //initialize to all-1 matrix
      bmat output_e(v,1);
      
      
      for (int l=1;l<=lmax;l++)
	{
	  //  cout<<l<<endl;
	  output_e.zeros();
	  quan_s_update(checks,errors, mcv,mvc,syndrome,pv, c, v,output_e);
	 
	  if (H*output_e==syndrome)
	    {
	      // num_iter=num_iter+l;
	     
	      if(quan_check(output_e,real_e,H2,c,rank2))
		{
		  cout<<"suc! Error wt= "<<distance(zero_mat2, real_e, v)<<endl;
		  num_iter= num_iter+l;
		  return true;
		  // cout<<"success! iteration number="<<l<<endl;
		  // num_of_suc_dec++;
		  //  cout<<num_of_suc_dec<<endl;
		}
	       else
	      	{
		   cout<<"failure! Error wt= "<<distance(zero_mat2, real_e, v)<<endl;
		  //cout<<"failure, e-e' is in NS"<<endl;
	      	  return false;
		  // er=er+ distance(output_e, real_e, n);	        
	      	}	    	  
	    }
	  
	}
     // num_iter=num_iter+lmax;
     // cout<<"failure, reach maximum iterations"<<endl;
       cout<<"failure and reach maximum iterations,  Error wt= "<<distance(zero_mat2, real_e, v)<<endl;
       return false;
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
   double  final_pr=(1-pv[j])/pv[j];

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
  bool Q_inspan(const bvec &real_eT,const bmat &H2, int ori_rank){
    bmat H2copy=H2;
    int r=H2.rows();
    H2copy.ins_row (r, real_eT);    
    int i=bmat_rank(H2copy);
    //  H2.del_row ( r);
  
    if (i==ori_rank){return true;}
    else if (i==ori_rank+1){return false;}
    else {cout<<"some wrong happened with the Q_inspan func"<<endl;return false;}
  

}

//check if real_e-output_e is a stabilizer
  bool quan_check(const bmat &output_e,const bmat &real_e,const bmat &H2, int n, int rank2){
    
   bvec difference(n);
   for (int i=0;i<n;i++)
     {
       difference(i)=output_e(i,0)+real_e(i,0);       
     }
    
   return Q_inspan(difference,H2,rank2);
 }



//get the rank of H and gaussian eliminate H:
int bmat_rank(const bmat& H){

  GF2mat Hp(H);
  return Hp.T_fact(T,U,P);	
}
