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


//find the neighbors of check nodes  
void initialize_checks (const bmat &H, nodes checks[], int & E){
  int r=H.rows();
  int c=H.cols();
  
  for (int i=0; i<r;i++)
    {
    checks[i].degree=0;
    for (int j=0; j<c;j++)
      { 
      if (H(i,j)==1)
	{
	checks[i].degree++;
	(checks[i].neighbors).ins(0,j);
	E++;
      }
     
    }
    
  }
}


//find the neighbors of variable(error) nodes
void initialize_errors(const bmat &H, nodes errors[]){

    int r=H.rows();
    int c=H.cols();
    int index=0;
 
    for (int i=0; i<c;i++)
      {
      errors[i].degree=0;
      for (int j=0; j<r ;j++)
	{
	if (H(j,i)==1)
	  {
	  errors[i].degree++;
	  (errors[i].neighbors).ins(0,j);
	  index++;	
	}	
      }
    }
}

//initialize massages to all-one matrices
void initialize_massages(mat &mcv,mat& mvc, bmat &H){
    int r=H.rows();
    int c=H.cols();
    for (int i=0;i<r;i++)
      {
      for (int j=0;j<c;j++)
	{
	  if (H(i,j)==1){
	    mcv(i,j)=1;
	    mvc(i,j)=1;	 
	  }
	}	     
    }  
}


void p_update(const nodes checks[],const nodes errors[],mat &mcv,mat& mvc,const bmat& syndrome,double p,int c, int v,  bmat& output_e){
  
    double ipr=(1-p)/p;

  for (int i=0;i<c;i++)
    {
      int ci_degree=checks[i].degree;
      
      	//update all v-to-ci massages first:
    for (int j=0;j<ci_degree;j++)
      {            
	int vnode=(checks[i].neighbors)(j);
       
	update_vj_to_ci(checks, errors,mcv, mvc,vnode,i, ipr);	
      }
    }
    //update c-to-vj massages:
  for (int j=0;j<v;j++)
    {
   int vj_degree=errors[j].degree;
   double  final_pr=ipr;

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


void s_update(const nodes checks[],const nodes errors[],mat &mcv,mat& mvc,const bmat& syndrome,double p,int c, int v,  bmat& output_e){
  
    double ipr=(1-p)/p; //initial p0/p1
    double  final_pr;  //final p0/p1
    
    //fix j,  for every v_j, do the following:
    for (int j=0;j<v;j++)
      {
	 int vj_degree=errors[j].degree;
	 final_pr=ipr;

	 //ci is the ith neibor of vj:
	   for (int i=0;i<vj_degree;i++)
     {
       
       //  update the  v_j to c_i massage:      
        int ci=(errors[j].neighbors)(i);
	update_vj_to_ci(checks, errors,mcv, mvc,j,ci, ipr);

      // update all  c_i to v_k massages:          
	int ci_degree=checks[ci].degree;
	for (int k=0;k<ci_degree;k++)
	  {
	    int vk=(checks[ci].neighbors)(k);
	    update_ci_to_vj( checks, errors,mcv, mvc,ci,vk,syndrome(ci,0));	    
	  }
	final_pr=final_pr*mcv(ci,j);      
     } 
   	   //  cout<<j<<"   "<<final_pr<<endl;
	   output_e(j,0)=final_pr<1? 1:0;      
      }  
}


void update_ci_to_vj(const nodes checks[],const nodes errors[],mat &mcv,mat& mvc,int i,int j,bin s)
{
   int c_degree=checks[i].degree;
   double temp=1.0;

   for (int k=0;k<c_degree;k++)
     {
      int vk=(checks[i].neighbors)(k);
      if (vk!=j)
	{
	  temp=temp*(mvc(i,vk)-1)/(mvc(i,vk)+1);
	}
     }   
     mcv.set(i,j,s==0? (1+temp)/(1-temp):(1-temp)/(1+temp));
}

void update_vj_to_ci(const nodes checks[],const nodes errors[],mat &mcv,mat& mvc,int j,int i,double ipr){
  
   int v_degree=errors[j].degree;
   mvc.set(i,j,ipr);

   for (int k=0;k<v_degree;k++)
     {
      int ck=(errors[j].neighbors)(k);
      if (ck!=i)
	{
	  mvc.set(i,j,mvc(i,j)*mcv(ck,j));
	}
     }   
}

  
//the distance between 2 cws
int distance(const bmat &output_e, const bmat &real_e,int n){
  
  int er=0;
  for (int i=0;i<n;i++)
    {
    if( (output_e(i,0)!=real_e(i,0)))
      {
      er++;
      }
    }
  return er;
}
      

//print the positions of errors
void print_error_pos( const bmat &real_e,int n){
  
  int er=0;
  for (int i=0;i<n;i++)
    {
    if (real_e(i,0)!=0)
      {
	cout<<" "<<i<<" ";
      }
    }
}
  