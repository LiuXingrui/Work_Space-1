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

bool  quan_decode2(GF2mat &H, GF2mat &G,const nodes checks[],const nodes errors[],const vec &pv,double pmin,double pmax,double& num_iter, int lmax,int wt,int num_dec){
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
	    
	      return false;
	     
	    }  
	}
  
      mat mcv(c,v);   // the messages from check nodes to variable nodes, which are p0/p1
      mat mvc(c,v);   // the messages from variable nodes to check nodes
      mcv.zeros();
      mvc.zeros();
      initialize_massages( mcv,mvc, H); //initialize to all-1 matrix
      GF2mat output_e(v,1);
     
      
      for (int l=1;l<=lmax;l++)
	{
	  //  cout<<l<<endl;
	  quan_s_update(checks,errors, mcv,mvc,syndrome,pv, c, v,output_e);
	 
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
	   
	      	  return false;
		  // er=er+ distance(output_e, real_e, n);	        
	      	}	    	  
	    }
	  
	}
      if (num_dec==2){
       vec px(v);    
       pro_dist( pmin,pmax, px);
        for (int l=1;l<=lmax;l++)
	{
	  //  cout<<l<<endl;
	  quan_s_update(checks,errors, mcv,mvc,syndrome,pv, c, v,output_e);
	 
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
	   
	      	  return false;
		  // er=er+ distance(output_e, real_e, n);	        
	      	}	    	  
	    }
	  
	}
      }
       return false;
 }

bool  quan_decode(GF2mat &H, GF2mat &G,const nodes checks[],const nodes errors[],const vec &pv,double& num_iter, int lmax,int wt){
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
	    
	      return false;
	     
	    }  
	}
  
      mat mcv(c,v);   // the messages from check nodes to variable nodes, which are p0/p1
      mat mvc(c,v);   // the messages from variable nodes to check nodes
      mcv.zeros();
      mvc.zeros();
      initialize_massages( mcv,mvc, H); //initialize to all-1 matrix
      GF2mat output_e(v,1);
     
      
      for (int l=1;l<=lmax;l++)
	{
	  //  cout<<l<<endl;
	  quan_s_update(checks,errors, mcv,mvc,syndrome,pv, c, v,output_e);
	 
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
	   
	      	  return false;
		  // er=er+ distance(output_e, real_e, n);	        
	      	}	    	  
	    }
	  
	}
       mcv.zeros();
       mvc.zeros();
       initialize_massages( mcv,mvc, H);
        vec pv2(v);
       pro_dist(0.5*wt/v,1.5*wt/v,pv2);
        for (int l=1;l<=lmax;l++)
	{
	   
	  quan_s_update(checks,errors, mcv,mvc,syndrome,pv2, c, v,output_e);
	  
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
		  // syn_fail++;		
	      	  return false;      	        
	      	}
 	    	  
	    }	  
	     
	}
       

       return false;
 }


bool  quan_decode_ana(GF2mat &H, GF2mat &G,const nodes checks[],const nodes errors[],const  vec &pv,double& num_iter, int lmax,int wt,int &syn_fail,int &max_fail){
  int v=H.cols();
  int c=H.rows();
  GF2mat real_eT(1,v);    //the transposed error vector, which is a row vector.

  if (wt==0)

    {
  error_channel(real_eT, pv);
    }
  else
    {
      error_channel2(real_eT,wt);
    }

 

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
	      cout<<"\n syndrome is 0, but decode fails"<<endl;
	      syn_fail++;
	      return false;
	     
	    }  
	}
  
      mat mcv(c,v);   // the messages from check nodes to variable nodes, which are p0/p1
      mat mvc(c,v);   // the messages from variable nodes to check nodes
      mcv.zeros();
      mvc.zeros();
      initialize_massages( mcv,mvc, H); //initialize to all-1 matrix
      GF2mat output_e(v,1);
     
      
      for (int l=1;l<=lmax;l++)
	{
	  //  cout<<l<<endl;
	  quan_s_update(checks,errors, mcv,mvc,syndrome,pv, c, v,output_e);	 
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

       cout<<"\n reach maximum iterations:"<<endl;
		  
       cout<<"\n real_e:"<<endl;
       err_pos(errors,real_e);

       cout<<"\n output_e:"<<endl;
       err_pos(errors,output_e);
       
       
       cout<<"try to see what happens to the messages \n"<<endl;
       mcv.zeros();
       mvc.zeros();
       initialize_massages( mcv,mvc, H);
       // vec pv2(v);
       // pro_dist(0.5*wt/v,1.5*wt/v,pv2);
        for (int l=1;l<=lmax/3;l++)
	{
	   cout<<"iteration "<<l<<":\n"<<endl;
	  quan_s_update(checks,errors, mcv,mvc,syndrome,pv, c, v,output_e);
	   cout<<"\n output_e:"<<endl;
	   err_pos(errors,output_e);
	    cout<<"mcv:\n"<<mcv<<endl;
	   cout<<"\nmvc:\n"<<mvc<<endl;
	  if (H*output_e==syndrome)
	    {
	      /*
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
 */	    	  
	    }	  
	     
	}
       
	 max_fail++;
       return  false;
	// return false;
 }

  void quan_s_update(const nodes checks[],const nodes errors[],mat &mcv,mat& mvc,const GF2mat& syndrome,const vec &pv,int c, int v,  GF2mat& output_e){
  
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

       output_e.set(j,0,final_pr<1? 1:0);   
    }  
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
  
