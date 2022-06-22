#include <iostream>
#include<cstring>
#include<string>
#include<vector>
#include <fstream>
#include <math.h> 
#include <random>

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


//find the neighbors of check nodes  
void initialize_checks (const GF2mat &H, nodes checks[], int & E){
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
void initialize_errors(const GF2mat &H, nodes errors[]){

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
void initialize_massages(mat &mcv,mat& mvc,const GF2mat &H){
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

//in fact, it is a sequential update schedule.
void p_update(const nodes checks[],const nodes errors[],mat &mcv,mat& mvc,const GF2mat& syndrome,vec& pv,int c, int v,  GF2mat& output_e){
  
   
  for (int i=0;i<c;i++)
    {
      int ci_degree=checks[i].degree;
      
      	//update all v-to-ci massages first:
    for (int j=0;j<ci_degree;j++)
      {            
	int vnode=(checks[i].neighbors)(j);
        double ipr=(1-pv[vnode])/pv[vnode];

	update_vj_to_ci(checks, errors,mcv, mvc,vnode,i, ipr);	
      }
    }
    //update c-to-vj massages:
  for (int j=0;j<v;j++)
    {
   int vj_degree=errors[j].degree;
   double  final_pr=(1-pv[j])/pv[j];;

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

//the parallel schedule
void real_p_update(const nodes checks[],const nodes errors[],mat &mcv,mat& mvc,mat &pre_mcv,mat& pre_mvc,const GF2mat& syndrome,double p,int c, int v,  GF2mat& output_e){
  
  double ipr=(1-p)/p;
  pre_mcv=mcv;
  pre_mvc=mvc;

  for (int i=0;i<c;i++)
    {
      int ci_degree=checks[i].degree;
      
      	//update all v-to-ci massages first:
    for (int j=0;j<ci_degree;j++)
      {            
	int vnode=(checks[i].neighbors)(j);
       
	update_vj_to_ci(checks, errors,pre_mcv, mvc,vnode,i, ipr);	
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
       update_ci_to_vj( checks, errors,mcv, pre_mvc,cnode,j,syndrome(cnode,0));

       final_pr=final_pr*mcv(cnode,j);
      }   
       //  cout<<j<<"   "<<final_pr<<endl;
       output_e.set(j,0,final_pr<1? 1:0);       
    }
 
}

//the sequential schedule from that paper
void s_update(const nodes checks[],const nodes errors[],mat &mcv,mat& mvc,const GF2mat& syndrome,double p,int c, int v,  GF2mat& output_e){
  
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
	   }
	 
	 for (int i=0;i<vj_degree;i++)
	   {     
	     int ci=(errors[j].neighbors)(i);
	     update_ci_to_vj( checks, errors,mcv, mvc,ci,j,syndrome(ci,0));	   // update all  c_i to v_j massages         
	     final_pr=final_pr*mcv(ci,j);      
	   } 
	 //  cout<<j<<"   "<<final_pr<<endl;
	 output_e.set(j,0,final_pr<1? 1:0);      
      }  
}


//update ci to vj massage:
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

//update vj to ci massage:
void update_vj_to_ci(const nodes checks[],const nodes errors[],mat &mcv,mat& mvc,int j,int i,double ipr,double alpha){
  
   int v_degree=errors[j].degree;
   mvc.set(i,j,ipr);

   for (int k=0;k<v_degree;k++)
     {
      int ck=(errors[j].neighbors)(k);
      if (ck!=i)
	{
	  mvc.set(i,j,mvc(i,j)*pow(mcv(ck,j),1.0/alpha));
	}
      else
	{
	   mvc.set(i,j,mvc(i,j)*pow(mcv(ck,j),1.0-1.0/alpha));
	}
     }   
}


int  cla_decode(int v,int c,const GF2mat &H,const nodes checks[],const nodes errors[],double& num_iter, int lmax,int & er,vec& pv){

  int n=v;
 
  //if no error, break
 
  
  GF2mat real_eT(1,v);    
  error_channel(real_eT, pv);
  GF2mat zero_vec(1,v);
  if (real_eT==zero_vec)
    {
      // cout<<"suc, no error"<<endl;
      return 1;
    }
  
  GF2mat zero_mat1(c,1);
  GF2mat zero_mat2(v,1);
  GF2mat real_e(v,1);  //the error vector which is a column vector,

     
  for (int q=0;q<v;q++)
    {
      real_e.set(q,0,real_eT(0,q));
    }
  GF2mat syndrome=H*real_e;

  //is the syndrome a zero vector?
  if (syndrome==zero_mat1)
	{
	  //cout<<"failure! error= some codeword"<<endl;	    
	  er=er+ distance(zero_mat2, real_e, n);
	  return 0;	 
		   
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
	  p_update(checks,errors, mcv,mvc,syndrome,pv, c, v,output_e);
	  
	  if (H*output_e==syndrome)
	    {
	      num_iter=num_iter+l;
	     
	      if(output_e==real_e)
		{
		  //cout<<"success! iteration number="<<l<<endl;
		  return 1;
		  
		  // num_of_suc_dec++;
		  //  cout<<num_of_suc_dec<<endl;
		}
	      else
		{
		  //cout<<"error!,get the wrong e"<<endl;
		  er=er+ distance(output_e, real_e, n);	  
		  return 0;		       
		}	    	  
	    }
	  if(l==lmax)
	    {
	      //cout<<"reach max iter"<<endl;
		er=er+ distance(output_e, real_e, n);
	      //num_iter=num_iter+l;
	      return 0;	   	 
	    }
	}
 }
  //the distance between 2 cws
int distance(const GF2mat &output_e, const GF2mat &real_e,int n){
  
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
void print_error_pos( const GF2mat &real_e,int n){
  
  int er=0;
  for (int i=0;i<n;i++)
    {
    if (real_e(i,0)!=0)
      {
	cout<<" "<<i<<" ";
      }
    }
}
  
//read a parity check matrix
GF2mat read_matrix (int& n,int &r, string & file_name){
  ifstream parity_check(file_name);
  string line;
  int row_ind=0;
  int temp;
  getline(parity_check, line);
  istringstream iss(line);
  
  if(  iss>>n>>r) {}
  else
    {
      cout<<"did not open the right parity check file"<<endl;
      
    }
 
  GF2mat H(r,n);
  while( getline(parity_check, line))
    {
      istringstream iss2(line);
      while (iss2>>temp){
	if (temp>=1&&temp<=n&&row_ind<r)
	  {
	    H.set(row_ind,temp-1,1);
	    // cout<< H(row_ind,temp-1)<<endl;
	  }
	else
	  {
	    cout<<"the format of the parity check is wrong, the first element is 1 rathar than 0"<<endl;
	    
	  }
	
      }
		  row_ind++;
    }
  parity_check.close();
  // cout<<H<<endl;
  return H;

  
}

bmat read_matrix2 (int& n,int &r, string & file_name){
  ifstream parity_check(file_name);
  string line;
  int row_ind=0;
  int temp;
  getline(parity_check, line);
  istringstream iss(line);
  
  if(  iss>>n>>r) {}
  else
    {
      cout<<"did not open the right parity check file"<<endl;
      
    }
 
  bmat H(r,n);
  H.zeros();
  while( getline(parity_check, line))
    {
      istringstream iss2(line);
      while (iss2>>temp){
	if (temp>=1&&temp<=n&&row_ind<r)
	  {
	    H(row_ind,temp-1)=1;
	    // cout<< H(row_ind,temp-1)<<endl;
	  }
	else
	  {
	    cout<<"the format of the parity check is wrong, the first element is 1 rathar than 0"<<endl;
	    
	  }
	
      }
		  row_ind++;
    }
  parity_check.close();
  // cout<<H<<endl;
  return H;

  
}

//write a matrix to a file
void write_matrix(string file_name, GF2mat &H){

  ofstream Hx;
  Hx.open (file_name,ios::trunc);
  int n=H.cols();
  int r=H.rows();


  Hx<<n<<" "<<r<<endl;

  for (int j=0;j<r;j++)
    {
      for (int i=0; i<n;i++)
	{
	  if (H(j,i)!=0)
	    {
	      Hx<<i+1<<" ";
	    }
	  
	}
      Hx<<endl;
    }
  

  Hx.close();

}

void write_matrix2(string file_name, bmat &H){

  ofstream Hx;
  Hx.open (file_name,ios::trunc);
  int n=H.cols();
  int r=H.rows();


  Hx<<n<<" "<<r<<endl;

  for (int j=0;j<r;j++)
    {
      for (int i=0; i<n;i++)
	{
	  if (H(j,i)!=0)
	    {
	      Hx<<i+1<<" ";
	    }
	  
	}
      Hx<<endl;
    }
  

  Hx.close();

}

//merge 2 matrices horizontally, [H1,H2]
bmat merge_mat_hori2(const bmat &left,const bmat &right){

  if (left.rows()!=right.rows())
    {
      cout<<"the 2 matrices cannot merge horizontally"<<endl;
      bmat error(1,1);
   
      return error;
    }

  else
    {
      int r=left.rows();
      int c=left.cols()+right.cols();
      int c1=left.cols();
      int c2=right.cols();
      bmat m(r,c);

      for (int i=0;i<r;i++)
	{
	  for (int j1=0;j1<c1;j1++){
	    m(i,j1)=left(i,j1);
	  }
	  for (int j2=0;j2<c2;j2++){
	    m(i,c1+j2)=right(i,j2);
	  }
	  
	}
       return m;
    }

 

}

GF2mat merge_mat_hori(const GF2mat &left,const GF2mat &right){

  if (left.rows()!=right.rows())
    {
      cout<<"the 2 matrices cannot merge horizontally"<<endl;
      GF2mat error(1,1);
   
      return error;
    }

  else
    {
      int r=left.rows();
      int c=left.cols()+right.cols();
      int c1=left.cols();
      int c2=right.cols();
      GF2mat m(r,c);

      for (int i=0;i<r;i++)
	{
	  for (int j1=0;j1<c1;j1++){
	    m.set(i,j1,left(i,j1));
	  }
	  for (int j2=0;j2<c2;j2++){
	    m.set(i,c1+j2,right(i,j2));
	  }
	  
	}
       return m;
    }

 

}
