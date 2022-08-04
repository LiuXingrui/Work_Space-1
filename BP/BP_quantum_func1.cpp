#include <iostream>
#include<cstring>
#include<string>
#include<vector>
#include <fstream>
#include <math.h> 
#include <random>
#include <sstream>
#include <stdlib.h> 

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

int error_channel(GF2mat &cw, const vec &p){

  double temp2;
  bin one=1;
  int temp=0;
  int r=cw.rows();
  if (cw.rows()!=p.size())
    {
      cout<<"the size of p and cw do not match"<<endl;
    }

  else
    {
    for (int i=0;i<r;i++)
      {

	temp2=randu();
	if(temp2<p[i])
	  {
	 
	    cw.set(i,0,cw(i,0)+one);
	    temp++;
	  }	
      }
    }
  return temp;
}


void error_channel2(GF2mat &error, int wt){
 
  double temp2;
  bin one=1;
  int n=error.rows();
  //cout<<1<<endl;

    for (int i=0;i<wt;i++)
      {    
	temp2=randi(0,n-1);
	error.set(temp2,0,one);	      
  }
    if (s_weight(error)!=wt)
      
      { GF2mat error2(n,1);
	error=error2;
	//cout<<"!=wt, try again"<<endl;
	error_channel2(error,wt);
      }
    //cout<<2<<endl;
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

int s_weight(GF2mat &s)
{
  int k=s.rows();
  int wt=0;
  for (int i=0;i<k;i++)
    {
      if(s(i,0)==1){wt++;}
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
bool  quan_decode(GF2mat &H, GF2mat &G,const nodes checks[],const nodes errors[],const vec &pv,const vec&pv_dec,double pavg,double range,double& num_iter, int lmax,int &wt,int& max_fail, int&syn_fail,  int debug, vec &LR,int rankH,int &OSD_suc,double alpha,double lambda){

 
  int wt_real_e=0;
  int v=H.cols();
  // int r2=H2.rows();
  GF2mat real_e(v,1); 

  if (wt==0)

    {


   wt_real_e=error_channel(real_e, pv);
     if (wt_real_e==0)
    {
	 
      return true;
    }

 
    }
  else
    {
      error_channel2(real_e,wt);
    }

  // cout<<pv<<endl;

  //if no error, break

 

  double pmin;
  double pmax;
  if (range>1)
    {
      pmin=0;
    }
  else
    {
    pmin=pavg*(1-range);
    }

  pmax=pavg*(1+range);
    
    GF2mat zero_rvec(1,v);
   
  int c=H.rows();
  LR.zeros();
  vec LR_avg=LR;
  GF2mat zero_mat1(c,1);//
  GF2mat zero_mat2(v,1);
  GF2mat zero_rvec2(v-rankH,1); //



  
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
      // mat pre_mcv=mcv;
      // mat pre_mvc=mvc;
      GF2mat output_e(v,1);
      // GF2mat output_e2=output_e;
     
      
      for (int l=1;l<=lmax;l++)
	{
	  //  cout<<l<<endl;
	  if ((debug/2)%2==1)
	    {
	      quan_p_update(checks,errors, mcv,mvc,syndrome,pv_dec, c, v,output_e,LR,alpha);
	    }
	  else
	    {
	      quan_s_update(checks,errors, mcv,mvc,syndrome,pv_dec, c, v,output_e,LR,alpha);
	    }
	  //cout<<"l="<<l<<endl;
	  // cout<<"LR="<<LR<<endl;
	  // LR_avg=LR*pow(LR_avg,lambda);
	  // LR_avg=0.9*LR_avg+LR;
	  //  for (int i=0;i<v;i++) {output_e2.set(i,0,LR_avg(i)<1? 1:0);}
	  if (H*output_e==syndrome)
	    {
	    
	      if(G*(output_e+real_e)==zero_rvec2)
		{
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
	  /*
	     if ((debug/4)%2==1)
	{
	  cout<<"BP fails for the first time, before OSD:"<<endl;
	  cout<<"output_e is:\n"<<endl;
	  err_pos(errors,output_e);
	  // cout<<"output_e2 is \n"<<endl;
	  // err_pos2(output_e2);
	  cout<<"real_e is \n"<<endl;
	  err_pos(errors,real_e);
	}
	  */
	  /*
	  else if (H*output_e2==syndrome)
	    {
	      if(G*(output_e2+real_e)==zero_rvec2)
		{
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
	  */
	  
	}
           if ((debug/8)%2==1)
	{
	  //cout<<"LR="<<LR<<endl;
	  //cout<<"mcv="<<mcv<<endl;
	  OSD(LR,H,syndrome,output_e,G,real_e);
	  //	if (OSD(LR,H,syndrome,output_e,G,real_e)==true);
	  if(G*(output_e+real_e)==zero_rvec2)
	    {
	      //   cout<<"OSD suc1.1"<<endl;
	      OSD_suc++;
	      return true;
	    }
	  /* else
	    {
	      OSD(LR_avg,H,syndrome,output_e2,G,real_e);
	      if(G*(real_e+output_e2)==zero_rvec2)
		{
		  //  cout<<"OSD suc1.2"<<endl;
		  OSD_suc++;
		  return true;
		}
	    }	 
	  */

	  
	
	  
	}
 
     
      if (debug%2==1)
	{
	  
	  mcv.zeros();
	  mvc.zeros();
	  initialize_massages( mcv,mvc, H);
	  //  pre_mcv=mcv;
	  //  pre_mvc=mvc;
	  vec pv2(v);
	  pro_dist(pmin,pmax,pv2);
	  for (int l=1;l<=lmax;l++)
	    {
	   
	      if ((debug/2)%2==1)
		{
		  quan_p_update(checks,errors, mcv,mvc,syndrome,pv2, c, v,output_e,LR,alpha);
		}
	      else
		{
		  quan_s_update(checks,errors, mcv,mvc,syndrome,pv2, c, v,output_e,LR,alpha);
		}
	      
	      //LR_avg=pow(LR,0.9)*pow(LR_avg,0.1);
	      // LR_avg=0.9*LR_avg+LR;
	      //  LR_avg=LR*pow(LR_avg,lambda);
	      //	for (int i=0;i<v;i++) {output_e2.set(i,0,LR_avg(i)<1? 1:0);}
	      /*
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
	      */
	      
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
	      /*
	       else if (H*output_e2==syndrome)
		 {
		   if(G*(output_e2+real_e)==zero_rvec2)
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
	      */
	    }
	}
      // GF2mat syndrome2=H*output_e+syndrome;
      //  cout<<" wt of s-s_output is"<<s_weight(syndrome2)<<endl;
      // cout<<"wt of s is"<<s_weight(syndrome)<<endl;
      //  cout<<"try OSD use s:"<<endl;
      /*
      if ((debug/4)%2==1)
	{
	  cout<<"BP fails the second time:"<<endl;
	  cout<<"output_e is:\n"<<endl;
	  err_pos(errors,output_e);
	  // cout<<"output_e2 is \n"<<endl;
	  // err_pos2(output_e2);
	  cout<<"real e wt="<<s_weight(real_e)<<endl;
	  cout<<"real_e is \n"<<endl;
	  err_pos(errors,real_e);
	}
      */
      
      if ((debug/8)%2==1)
	{
	  OSD(LR,H,syndrome,output_e,G,real_e);
	  //	if (OSD(LR,H,syndrome,output_e,G,real_e)==true);
	  if(G*(output_e+real_e)==zero_rvec2)
	    {
	      // cout<<"OSD suc2.1"<<endl;
	      OSD_suc++;
	      return true;
	    }
	 
	  /*
	  else
	    {
	      OSD(LR_avg,H,syndrome,output_e2,G,real_e);
	      if(G*(real_e+output_e2)==zero_rvec2)
		{
		  // cout<<"OSD suc2.2"<<endl;
		  OSD_suc++;
		  return true;
		}
	    }
	  */	 
	}
      // cout<<"OSD fail"<<endl;
      // cout<<"output_e is \n"<<endl;
      // err_pos2(output_e);

      //  cout<<"output_e2 is \n"<<endl;
      //  err_pos2(output_e2);
     
      // cout<<"OSD failed"<<endl;
       if ((debug/4)%2==1)
		  {
		    //cout<<"failed"<<endl;
		    //cout<<"output e wt="<<s_weight(output_e)<<endl;
		    // cout<<"show the position of errors: only for checkerboard codes"<<endl;
		    cout<<"output e is:\n"<<endl;
		    err_pos2(output_e);
		    // cout<<"output_e2 is \n"<<endl;
		    // err_pos2(output_e2);
		    //cout<<"real e wt="<<s_weight(real_e)<<endl;
		    cout<<"real e is \n"<<endl;
		 
		    err_pos2(real_e);
		    GF2mat sume=real_e+output_e;
		    //cout<<"output_e+real_e wt="<<s_weight(sume)<<endl;
		    cout<<"residual e is \n"<<endl;
		    err_pos2(sume);
		  }
      max_fail++;
      return false;
 }


void OSD(vec& LR,const GF2mat& H,const GF2mat& syndrome, GF2mat &output_e,const GF2mat& G,const GF2mat &real_e){ //r is the rank of H

  
  //get the first permutation that abs_LLR is in descending order
  int n=LR.length();
    // cout<<"LR is \n"<<LR<<endl;

  int r=H.rows();
  int k=n-r;
  vec LLR(n); //the sort function gives a ascending  order, LR=p0/p1
  for (int i=0;i<n;i++)
    {
      if (LR(i)>=0) 
	{
	  LLR(i)=log(LR(i));
	}
      else
	{

	  cout<<"negative LR!"<<endl;
	  cout<<LR(i)<<endl;
	}
    }
  //now we have H*perm1*perm1_inv*e=s, and def:H1=H*perm1
  ivec perm1=sort_index(LLR);
  //for (int i=0;i<n;i++){perm1(i)=i;}
  
  GF2mat H1=H;
  H1.permute_cols(perm1,0);
  GF2mat perm1_mat=col_permutation_matrix(perm1);
  GF2mat perm1_mat_inv=perm1_mat.inverse();
   
  GF2mat T;
  ivec perm2;
  GF2mat U;
  // now we have T*H1*perm2*perm2_inv*perm1_inv*e=T*s, def s1=T*s, U=T*H1*perm2, e1=perm2_inv*perm1_inv*e:

  int rankH=H1.T_fact(T,U,perm2);
  GF2mat perm2_mat=col_permutation_matrix(perm2);
  GF2mat perm2_mat_inv=perm2_mat.inverse();
      
  GF2mat syndrome1=T*syndrome;
    //U may be not a full rank matrix, need to delete some rows. So also need to delete some rows of syndrome1
  GF2mat H2=U.get_submatrix(0,0,rankH-1,n-1);
  GF2mat syndrome2=syndrome1.get_submatrix(0,0,rankH-1,0);
  
  
  // OSD-0: e1=(HS_inv*s2
  //                  0   )
  
  GF2mat H_S=H2.get_submatrix(0,0,rankH-1,rankH-1);
  GF2mat H_T=H2.get_submatrix(0,rankH,rankH-1,n-1);
  GF2mat HS_inv=H_S.inverse();

  

  GF2mat e_S=HS_inv*syndrome2;  
  GF2mat e_T(n-rankH,1);
  GF2mat new_e_S;
  int wt=s_weight(e_S)+1;
  //cout<<"original wt "<<wt<<endl;
  int temp1;
  int temp2=-1;
  int temp3=-1;
  int lambda=min(60,n-rankH);
  
  for (int i=0;i<n-rankH;i++)
    {
      //cout<<i<<": wt is "<<wt<<endl;
      GF2mat new_e_T=e_T;
      new_e_T.set(i,0,1);
      new_e_S=HS_inv*e_S+HS_inv*H_T*e_T;
      temp1=s_weight(new_e_S);
      // cout<<"another wt "<<temp1<<endl;
      if (temp1+1<wt)
	{
	  // cout<<"another e_S"<<endl;
	  // cout<<"old wt is "<<wt<<endl;
	  wt=temp1+1;
	 
	  // cout<<"new wt is "<<wt<<endl;
	  temp2=i;
	}
    }


  for (int i=0;i<lambda;i++)
    {
      for (int j=0;j<lambda;j++)
	{
	  GF2mat new_e_T=e_T;
	  new_e_T.set(i,0,1);
	  new_e_T.set(j,0,1);		  
	  new_e_S=HS_inv*e_S+HS_inv*H_T*e_T;
	  temp1=s_weight(new_e_S);
      // cout<<"another wt "<<temp1<<endl;
	  if (temp1+2<wt)
	    {
	      // cout<<"another e_S"<<endl;
	      wt=temp1+2;
	      // cout<<"new wt is "<<wt<<endl;

	      temp2=i;
	      temp3=j;
	    }
	}
    }

  
      if (temp2!=-1&&temp3==-1)
	{
	  e_T.set(0,temp2,1);
	  e_S=e_S+HS_inv*H_T*e_T;
	  // cout<<1<<endl;
	}
      else if (temp2!=-1&&temp3!=-1)
	{
	  e_T.set(0,temp2,1);
	  e_T.set(0,temp3,1);
	  e_S=e_S+HS_inv*H_T*e_T;
	  //cout<<2<<endl;
	}
      
      
      for (int i=0;i<rankH;i++){output_e.set(i,0,e_S(i,0));}
      for (int i=rankH;i<n;i++){output_e.set(i,0,e_T(i-rankH,0));}

	  
  if (H2*output_e==syndrome2){}
  else
    {
      cout<<"error H2*output_e!=syndrome2 output e is:"<<endl;
      err_pos2(output_e);
      cout<<"e_S is"<<endl;
      err_pos2(e_S);
      cout<<"e_T is"<<endl;
      err_pos2(e_T);
    }
 
  output_e=perm1_mat*perm2_mat*output_e;
   	  	  

   if (H*output_e==syndrome){}
   else
     {
       cout<<"error H*output_e!=syndrom"<<endl;
       //cout<<"syndrome is"<<syndrome<<endl;
       //cout<<"output_e is"<<output_e<<endl;
     }
  

  
}

GF2mat col_permutation_matrix(ivec & perm)
{
  int n=perm.length();
  GF2mat p(n,n);
  for ( int i=0;i<n;i++)
    {
      p.set(perm(i),i,1);
    }
  return p;


}



void quan_p_update(const nodes checks[],const nodes errors[],mat &mcv,mat& mvc,const GF2mat& syndrome,const vec &pv,int c, int v,  GF2mat& output_e, vec &LR,double alpha){
  
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

       final_pr=final_pr*pow(mcv(cnode,j),1.0/alpha);
      }   
       //  cout<<j<<"   "<<final_pr<<endl;
   LR(j)=final_pr;
   output_e.set(j,0,final_pr<1? 1:0);   
    }  
}

void quan_s_update(const nodes checks[],const nodes errors[],mat &mcv,mat& mvc,const GF2mat& syndrome,const vec &pv,int c, int v,  GF2mat& output_e, vec &LR,double alpha){
  
    double ipr;
    double final_pr;

   for (int j=0;j<v;j++)
      {
	 int vj_degree=errors[j].degree;
	 ipr=(1-pv[j])/pv[j];
	 final_pr=ipr;

	 //ci is the ith neibor of vj:
	 for (int i=0;i<vj_degree;i++)
	   {
	     //  update the  v_j to c_i massage:      
	     int ci=(errors[j].neighbors)(i);
	     update_vj_to_ci(checks, errors,mcv, mvc,j,ci, ipr);          
	   }
	 
	 for (int i=0;i<vj_degree;i++)
	   {     
	     int ci=(errors[j].neighbors)(i);
	     update_ci_to_vj( checks, errors,mcv, mvc,ci,j,syndrome(ci,0));	   // update all  c_i to v_j massages         
	     final_pr=final_pr*mcv(ci,j);      
	   } 
	 //  cout<<j<<"   "<<final_pr<<endl;
	 LR(j)=final_pr;
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


void err_pos2(const GF2mat &error){
  cout<<"\n";
  int n=error.rows();
  int d=sqrt(n);
  if (d*d==n)
    {
      for (int r=0;r<d;r++)
	{
	  cout<<r<<"th row: ";
	  for (int c=0;c<d;c++)	    
	    {
	      if (error(r*d+c,0)==1){cout<<1;}
	      else {cout<<".";}
	    }
	  cout<<endl;
	}
    }
    cout<<"\n";
 
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
  
/*
void OSD(vec& LR,const GF2mat& H,const GF2mat& syndrome,  const GF2mat &real_e1,GF2mat &output_e){ //r is the rank of H

  
  //get the first permutation that abs_LLR is in descending order
  GF2mat real_e=real_e1;
    int n=LR.length();
    // cout<<"LR is \n"<<LR<<endl;

    int r=H.rows();
  int k=n-r;
  vec nega_abs_LLR(n); //the sort function gives a ascending  order, we need descending order
  for (int i=0;i<n;i++)
    {
      if (LR(i)>=0) 
	{
	  nega_abs_LLR(i)=-abs(log(LR(i)));
	}
      else
	{

	  cout<<"something goes wrong in the OSD func"<<endl;
	}
    }
  //now we have H*perm1*perm1_inv*e=s, and def:H1=H*perm1
  ivec perm1=sort_index(nega_abs_LLR);
  
  //  cout<<"perm1 is\n"<<perm1<<endl;
  //  cout<<"H is \n"<<H<<endl;
  real_e.permute_rows(perm1,0);

  
  GF2mat H1=H;
  cout<<"H is"<<H<<endl;
  cout<<"perm1 is"<<perm1<<endl;
  H1.permute_cols(perm1,0);
  cout<<"after .per, H is"<<H1<<endl;
     GF2mat perm1_mat=col_permutation_matrix(perm1);
     GF2mat perm1_mat_inv=perm1_mat.inverse();

     GF2mat H11=H;
     cout<<"perm1_mat is"<<perm1_mat<<endl;
     H11=H11*perm1_mat;
     cout<<"after mul,H is"<<H11<<endl;
     
  //  cout<<"after perm1, H is"<<H1<<endl;
     if (H1*real_e==syndrome){cout<<" (H1*real_e==syndrome),ok"<<endl;}
   else {cout<<" (H1*real_e!=syndrome),error"<<endl;}
   
  GF2mat T;
  ivec perm2;
  GF2mat U;
  // now we have T*H1*perm2*perm2_inv*perm1_inv*e=T*s, def s1=T*s, U=T*H1*perm2, e1=perm2_inv*perm1_inv*e:
  //notice U!=TH1P, this function gives a wrong P, so we need to do it ourselves.
  int rankH=H1.T_fact(T,U,perm2);
     GF2mat perm2_mat=col_permutation_matrix(perm2);
     GF2mat perm2_mat_inv=perm2_mat.inverse();
    
  GF2mat H1copy=H1;
  H1copy.permute_cols(perm2,0);
  H1copy=T*H1copy;
    cout<<"H1copy="<<H1copy<<endl;
    cout<<"U is" <<U<<endl;
    
    GF2mat perm1_mat=mypermutation_matrix(perm1);
      GF2mat perm2_mat=mypermutation_matrix(perm2);
      GF2mat perm1_mat_inv=perm1_mat.inverse();
        GF2mat perm2_mat_inv=perm2_mat.inverse();
    
  // GF2mat U2=U;
  
  // U2.permute_cols(perm2,0);
  
  // cout<<"UP2^-1 should be TH1"<<U2<<endl;
  // cout<<"TH1 is"<<T*H1<<endl;

  // U2.permute_cols(perm1,0);
  
  // cout<<"UP2^-1P1^-1 should be TH"<<U2<<endl;
  // cout<<"TH is"<<T*H<<endl;
  
  //  cout<<"perm2 is"<<perm2<<endl; //perm2 should be close to Identity matrix
  //  cout<<"U is"<<U<<endl;
  
  GF2mat syndrome1=T*syndrome;
  real_e.permute_rows(perm2,0);
    //U may be not a full rank matrix, need to delete some rows. So also need to delete some rows of syndrome1
  GF2mat H2=U.get_submatrix(0,0,rankH-1,n-1);
    if (U*real_e==syndrome1){cout<<" (U*real_e==syndrome1),ok"<<endl;}
   else {cout<<" (U*real_e!=syndrome1),error"<<endl;}
     

 


  GF2mat syndrome2=syndrome1.get_submatrix(0,0,rankH-1,0);
    if (H2*real_e==syndrome2){cout<<" (H2*real_e==syndrome2),ok"<<endl;}
   else {cout<<" (H2*real_e!=syndrome2),error"<<endl;}
  
  // OSD-0: e1=(HS_inv*s2
  //                  0   )
  
  GF2mat H_S=H2.get_submatrix(0,0,rankH-1,rankH-1);
  GF2mat H_T=H2.get_submatrix(0,rankH,rankH-1,n-1);
  // cout<<"HS is"<<H_S<<endl;
  // cout<<"HT is"<<H_T<<endl;
  // cout<<"H2 is"<<H2<<endl;
 
  GF2mat HS_inv=H_S.inverse();

  // cout<<"Hs_inv is"<<HS_inv<<endl;
  //   cout<<"HS is"<<H_S<<endl;
  GF2mat e_S=HS_inv*syndrome2;
  
  GF2mat e_T(n-rankH,1);
  //   cout<<"eS is"<<e_S<<endl;
  //cout<<"eT is"<<e_T<<endl;

  for (int i=0;i<rankH;i++){output_e.set(i,0,e_S(i,0));}
   for (int i=rankH;i<n-1;i++){output_e.set(i,0,e_T(i-rankH,0));}
  

     if (H2*output_e==syndrome2){cout<<"H2*output_e==syndrome2 ok"<<endl;}
   else {cout<<"error H*2output_e!=syndrom2"<<endl;}
         if (U*output_e==syndrome1){cout<<"U*output_e==syndrome1 ok"<<endl;}
   else {cout<<"error U*output_e!=syndrom1"<<endl;}
  //   cout<<"output_e before perm is"<<output_e<<endl;
  // cout<<"real_e is"<<real_e<<endl;
   cout<<"H2*(real_e+output_e)"<<H2*(real_e+output_e)<<endl;
  //  cout<<"U*(real_e+output_e)"<<U*(real_e+output_e)<<endl;

  
  //    cout<<"perm1_mat_inv* perm1_mat"<<perm1_mat_inv* perm1_mat<<endl;
  //  cout<<"perm1_mat_inv*perm2_mat_inv*perm1_mat*perm2_mat"<<perm1_mat_inv*perm2_mat_inv*perm1_mat*perm2_mat<<endl;
  //e=perm1*perm2*e1:
 
  // output_e.permute_rows(perm1,0);
  GF2mat Ucopy=U;
     if (Ucopy*output_e==syndrome1){cout<<"before perm:Ucopy*output_e==syndrome1 ok"<<endl;}
   else {cout<<"error before perm:Ucopy*output_e!=syndrom1"<<endl;}
    
       cout<<"output_e before perm2 is"<<output_e<<endl;
         cout<< "Ucopy before perm2 is"<<Ucopy<<endl;
  
        cout<<"perm2 is"<<perm2<<endl;
       cout<<"perm2 mat is"<<perm2_mat<<endl;
     output_e=perm2_mat*output_e;
   
     Ucopy=Ucopy*perm2_mat_inv;
           cout<<"output_e after perm2 is"<<output_e<<endl;
         cout<< "Ucopy after perm2 is"<<Ucopy<<endl;
            if (Ucopy*output_e==syndrome1){cout<<"after perm2:Ucopy*output_e==syndrome1 ok"<<endl;}
   else {cout<<"error after perm2:Ucopy*output_e!=syndrom1"<<endl;}
	  
	  
  GF2mat T_inv=T.inverse();
 
     output_e=perm1_mat*output_e;
   
     Ucopy=Ucopy*perm1_mat_inv;
   if (Ucopy*output_e==syndrome1){cout<<"after perm1 TH*output_e==syndrome1 ok"<<endl;}
   else {cout<<"error after perm1 TH*output_e!=syndrom1"<<endl;}
   GF2mat syndrome1_pr=T_inv*syndrome1;
      if (syndrome1_pr==syndrome){cout<<"syndrome1_pr=syndrome) ok"<<endl;}
   else {cout<<"error syndrome1_pr!=syndrome)"<<endl;}
  H1copy=T_inv*Ucopy;

  cout<<"T_inv*Ucopy should be H"<<H1copy<<endl;
  cout<<"H is "<<H<<endl;

  // cout<<"permute reale and outpute"<<endl;
  // U=U*perm2_mat_inv*perm1_mat_inv;
  // cout<<"after permutation, u is"<<U<<endl;
  // cout<<"TH is"<<T*H<<endl;
  //   cout<<"U*(real_e+output_e)"<<U*(real_e+output_e)<<endl;
  //     cout<<"TH*(real_e+output_e)"<<T*H*(real_e+output_e)<<endl;
  //  cout<<"output_e after perm is"<<output_e<<endl;
  if (H*output_e==syndrome){cout<<"suc!"<<endl;}
  else {cout<<"error H*output_e!=syndrom"<<endl;}
  

  
}
*/
