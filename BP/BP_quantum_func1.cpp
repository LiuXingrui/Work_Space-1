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

void error_channel(GF2mat &cw, const vec &p){
 
  double temp2;
  bin one=1;
  if (cw.cols()!=p.size())
    {cout<<"the size of p and cw do not match"<<endl;}

  else
    {
    for (int i=0;i<cw.cols();i++)
      {    
	temp2=randu();
	if(temp2<p[i])
	  {
	    cw.set(0,i,cw(0,i)+one);
	  }	
      }
  }
  
}


void error_channel2(GF2mat &error, int wt){
 
  double temp2;
  bin one=1;
  int n=error.cols();
  

    for (int i=0;i<wt;i++)
      {    
	temp2=randi(0,n-1);
	error.set(0,temp2,one);	      
  }
    if(weight(error)!=wt)
      
      { GF2mat error2(1,n);
	error=error2;
	//cout<<"!=wt, try again"<<endl;
	error_channel2(error,wt);
      }
}

void depolarizing(GF2mat &xerror,GF2mat &zerror, const vec &p){
  double temp1,temp2;
  bin one=1;
  int n=xerror.cols();

    for (int i=0;i<n;i++)
      {    
	temp2=randu();
	if(temp2<p[i])
	  {
	    temp1=randu();
	    if (temp1<0.25)
	      {
		xerror.set(0,i,xerror(0,i)+one);
	      }
	    else if (temp1>0.25&&temp1<0.5)
	      {
		 zerror.set(0,i,zerror(0,i)+one);
	      }
	    else if (temp1>0.5&&temp1<0.75)
	      {
		 zerror.set(0,i,zerror(0,i)+one);
		 xerror.set(0,i,xerror(0,i)+one);
	      }
	    
	  }	
      }
  
  

}


int weight(GF2mat &cw)
{
  int n=cw.cols();
  int wt=0;
  for (int i=0;i<n;i++)
    {
      if(cw(0,i)==1){wt++;}
    }
  return wt;
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


//here pavg and range are decode_p/decode_prange
bool  quan_decode(GF2mat &H, GF2mat &G,const nodes checks[],const nodes errors[],const vec &pv,const vec&pv_dec,double pavg,double range,double& num_iter, int lmax,int wt,int& max_fail, int&syn_fail,  int debug,double alpha){
  int v=H.cols();
  int c=H.rows();
  // int r2=H2.rows();
  GF2mat real_eT(1,v);    //the transposed error vector, which is a row vector.

  if (wt==0)

    {
  error_channel(real_eT, pv);
    }
  else
    {
      error_channel2(real_eT,wt);
    }

  // cout<<pv<<endl;

  //if no error, break
  GF2mat zero_rvec(1,v);
 
  if (real_eT==zero_rvec)
    {
	 
      return true;
    }
 
  GF2mat zero_mat1(c,1);
  GF2mat zero_mat2(v,1);
  GF2mat zero_rvec2(v-GF2mat_rank(H),1);
  GF2mat real_e(v,1);  //the error vector which is a column vector,

  
  for (int q=0;q<v;q++)
    {
      real_e.set(q,0,real_eT(0,q));
    }
  
  GF2mat syndrome=H*real_e;

  //is the syndrome a zero vector?
  if (syndrome==zero_mat1)
	{
	  // is e a stablizer?
	  if (G*real_e==zero_rvec2)
	    {
	      
	      return true;
	    }        
	  else  
	    {
	      syn_fail++;
	      return false;
	     
	    }  
	}
  
      mat mcv(c,v);   // the messages from check nodes to variable nodes, which are p0/p1
      mat mvc(c,v);   // the messages from variable nodes to check nodes
      mcv.zeros();
      mvc.zeros();
    
      initialize_massages( mcv,mvc, H); //initialize to all-1 matrix
      mat pre_mcv=mcv;
      mat pre_mvc=mvc;
      GF2mat output_e(v,1);
     
      
      for (int l=1;l<=lmax;l++)
	{
	  //  cout<<l<<endl;
	  if ((debug/2)%2==1)
	    {
	      quan_p_update(checks,errors, mcv,mvc,pre_mcv,pre_mvc,syndrome,pv_dec, c, v,output_e,alpha);
	    }
	  else
	    {
	      quan_s_update(checks,errors, mcv,mvc,syndrome,pv_dec, c, v,output_e,alpha);
	    }
	  if (H*output_e==syndrome)
	    {
	    
	      if(G*(output_e+real_e)==zero_rvec2)
		{
		  // cout<<"suc! Error wt= "<<distance(zero_mat2, real_e, v)<<endl;
		  num_iter= num_iter+l;
		  return true;
		
		}
	       else
	      	{
	          syn_fail++;
	      	  return false;
		  // er=er+ distance(output_e, real_e, n);	        
	      	}	    	  
	    }
	  
	}

      if (debug%2==1)
	{
	  
	  mcv.zeros();
	  mvc.zeros();
	  initialize_massages( mcv,mvc, H);
	  pre_mcv=mcv;
	  pre_mvc=mvc;
	  vec pv2(v);
	  pro_dist(pavg*(1-range),pavg*(1+range),pv2);
	  for (int l=1;l<=lmax;l++)
	    {
	   
	      if ((debug/2)%2==1)
		{
		  quan_p_update(checks,errors, mcv,mvc,pre_mcv,pre_mvc,syndrome,pv2, c, v,output_e,alpha);
		}
	      else
		{
		  quan_s_update(checks,errors, mcv,mvc,syndrome,pv2, c, v,output_e,alpha);
		}
	      
	      if ((debug/4)%2==1)
		{
		  cout<<"iter l: "<<l<<endl;
		  cout<<"\n real_e:\n"<<endl;
		  err_pos(errors,real_e);
		  cout<<"\n output_e:\n"<<endl;
		  err_pos(errors,output_e);
		  cout<<"\nmcv:\n"<<mcv<<"\n mvc \n"<<mvc<<endl;
		  //cout<<"\npremcv:\n"<<pre_mcv<<"\n premvc \n"<<pre_mvc<<endl;
		}
	      
	      if (H*output_e==syndrome)
		{
		  
		  if(G*(output_e+real_e)==zero_rvec2)
		    {
		      num_iter= num_iter+l;
		      return true;
		    }
		  else
		    {
		      // cout<<"\n syndrome is ok, but decoding fails:"<<endl;
		      syn_fail++;		
		      return false;      	        
		    }
 	    	  
		}	  
	     
	    }
	}
      max_fail++;
      return false;
 }




void quan_s_update(const nodes checks[],const nodes errors[],mat &mcv,mat& mvc,const GF2mat& syndrome,const vec &pv,int c, int v,  GF2mat& output_e,double alpha){
  
    double ipr;

  for (int i=0;i<c;i++)
    {
      int ci_degree=checks[i].degree;
      
      	//update all v-to-ci massages first:
    for (int j=0;j<ci_degree;j++)
      {            
	int vnode=(checks[i].neighbors)(j);
        ipr=(1-pv[vnode])/pv[vnode];
	update_vj_to_ci(checks, errors,mcv, mvc,vnode,i, ipr,alpha);	
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

       output_e.set(j,0,final_pr<1? 1:0);   
    }  
}

void quan_p_update(const nodes checks[],const nodes errors[],mat &mcv,mat& mvc,mat &pre_mcv,mat& pre_mvc,const GF2mat& syndrome,const vec &pv,int c, int v,  GF2mat& output_e,double alpha){
  
    double ipr;

  for (int i=0;i<c;i++)
    {
      int ci_degree=checks[i].degree;
      
      	//update all v-to-ci massages first:
    for (int j=0;j<ci_degree;j++)
      {            
	int vnode=(checks[i].neighbors)(j);
        ipr=(1-pv[vnode])/pv[vnode];
	update_vj_to_ci(checks, errors,pre_mcv, mvc,vnode,i, ipr,alpha);	
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
       update_ci_to_vj( checks, errors,mcv, pre_mvc,cnode,j,syndrome(cnode,0));

       final_pr=final_pr*mcv(cnode,j);
      }   
       //  cout<<j<<"   "<<final_pr<<endl;

       output_e.set(j,0,final_pr<1? 1:0);   
    }
  pre_mcv=mcv;
  pre_mvc=mvc;
}


//get the rank of H and gaussian eliminate H:
int GF2mat_rank(const GF2mat& H){

  GF2mat T,U;
  ivec P;
  return H.T_fact(T,U,P);	
}


GF2mat get_gen(const GF2mat &H){
  GF2mat HT=H.transpose();
  GF2mat T,U;
  ivec P;
  int Hrank= HT.T_fact(T,U,P);
  int r=H.rows();
  int n=H.cols();
  int k=n-r;


  GF2mat G=T.get_submatrix(Hrank,0,n-1,n-1);
  return G;


}


void err_pos(const nodes errors[],const GF2mat &error){
  int n=error.rows();
 
  
 
  // ivec pos; //an empty vector,ins() works well, 
  // int pos_size=0;


  for (int i=0;i<n;i++)
    {
      if (error(i,0)==1)
	{
	  // pos.ins(0,i);
	  // pos_size++;
	  cout<<i<<" neighbour checks: ";
	  cout<<errors[i].neighbors<<endl;
	}
    }
 
}
  
