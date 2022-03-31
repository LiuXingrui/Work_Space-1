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


void p_update_a(const nodes checks[],const nodes errors[],mat &mcv,mat& mvc,const bmat& syndrome,double p,int c, int v,  bmat& output_e){
  
    double ipr=(1-p)/p;

  for (int i=0;i<c;i++)
    {
      int ci_degree=checks[i].degree;
      
      	//update all v-to-ci massages first:
    for (int j=0;j<ci_degree;j++)
      {            
	int vnode=(checks[i].neighbors)(j);
       
	update_vj_to_ci(checks, errors,mcv, mvc,vnode,i, ipr);
	cout<<"v"<<vnode<<" to c"<<i<<":  "<<mvc(i,vnode)<<endl;
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
       cout<<"c"<<cnode<<" to v"<<j<<":  "<<mcv(cnode,j)<<endl;

       final_pr=final_pr*mcv(cnode,j);
      }   
       //  cout<<j<<"   "<<final_pr<<endl;
       output_e(j,0)=final_pr<1? 1:0;      
    }  
}


void s_update_a(const nodes checks[],const nodes errors[],mat &mcv,mat& mvc,const bmat& syndrome,double p,int c, int v,  bmat& output_e){
  
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
	     cout<<"v"<<j<<" to "<<ci<<":  "<<mvc(ci,j)<<endl;
	   }
	 
	 for (int i=0;i<vj_degree;i++)
	   {     
	     int ci=(errors[j].neighbors)(i);
	     update_ci_to_vj( checks, errors,mcv, mvc,ci,j,syndrome(ci,0));	   // update all  c_i to v_j massages
	     cout<<"c"<<ci<<" to v"<<j<<":  "<<mcv(ci,j)<<endl;
	     final_pr=final_pr*mcv(ci,j);      
	   } 
	 //  cout<<j<<"   "<<final_pr<<endl;
	 output_e(j,0)=final_pr<1? 1:0;      
      }  
}
