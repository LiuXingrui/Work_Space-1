#include <iostream>
#include<cstring>
#include<string>
#include<vector>
#include <fstream>
#include <math.h> 
#include <random>
#include"BP_basic_func1.h"

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
void initialize_massages(mat &mcv,mat& mvc,const bmat &H){
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

//the parallel schedule
void real_p_update(const nodes checks[],const nodes errors[],mat &mcv,mat& mvc,mat &pre_mcv,mat& pre_mvc,const bmat& syndrome,double p,int c, int v,  bmat& output_e){
  
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
       output_e(j,0)=final_pr<1? 1:0;      
    }
 
}

//the sequential schedule from that paper
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
	   }
	 
	 for (int i=0;i<vj_degree;i++)
	   {     
	     int ci=(errors[j].neighbors)(i);
	     update_ci_to_vj( checks, errors,mcv, mvc,ci,j,syndrome(ci,0));	   // update all  c_i to v_j massages         
	     final_pr=final_pr*mcv(ci,j);      
	   } 
	 //  cout<<j<<"   "<<final_pr<<endl;
	 output_e(j,0)=final_pr<1? 1:0;      
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


int  cla_decode(int v,int c,const bmat &H,const nodes checks[],const nodes errors[],BSC &bsc,double& num_iter, int lmax,int & er,double p){

  int n=v;
  bvec real_eT(v);    //the transposed error vector, which is a row vector.
 
  
  //if no error, break
  bvec zero_vec(v);
  zero_vec.zeros();
  real_eT=bsc(zero_vec);
  if (real_eT==zero_vec)
    {
      // cout<<"zero"<<endl;
      return 1;
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
	  er=er+ distance(zero_mat2, real_e, n);
	  return 0;	 
	  // cout<<"failure! error= some codeword"<<endl;	   
	}
   
      mat mcv(c,v);   // the messages from check nodes to variable nodes, which are p0/p1
      mat mvc(c,v);   // the messages from variable nodes to check nodes
      initialize_massages( mcv,mvc, H); //initialize to all-1 matrix
      bmat output_e(v,1);
      output_e.zeros();
      
      for (int l=1;l<=lmax;l++)
	{
	  //  cout<<l<<endl;
	  p_update(checks,errors, mcv,mvc,syndrome,p, c, v,output_e);
	  
	  if (H*output_e==syndrome)
	    {
	      num_iter=num_iter+l;
	     
	      if(output_e==real_e)
		{
		  return 1;
		  // cout<<"success! iteration number="<<l<<endl;
		  // num_of_suc_dec++;
		  //  cout<<num_of_suc_dec<<endl;
		}
	      else
		{
		  er=er+ distance(output_e, real_e, n);	  
		  return 0;		       
		}	    	  
	    }
	  if(l==lmax)
	    {
	    er=er+ distance(output_e, real_e, n);
	    num_iter=num_iter+l;
	    return 0;	   	 
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
  
//read a parity check matrix
bmat read_matrix (int& n,int &r, string & file_name){
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

  while( getline(parity_check, line))
    {
      istringstream iss2(line);
      while (iss2>>temp){
	if (temp>=1&&temp<=n&&row_ind<r)
	  {
	    H(row_ind,temp-1)=1;
	  }
	else
	  {
	    cout<<"the format of the parity check is wrong, the first element is 1 rathar than 0"<<endl;
	    
	  }
	
      }
		  row_ind++;
    }
  parity_check.close();

  return H;

  
}

//write a matrix to a file
void write_matrix(string file_name, bmat &H){

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
bmat merge_mat_hori(const bmat &left,const bmat &right){

  if (left.rows()!=right.rows())
    {
      cout<<"the 2 matrices cannot merge horizontally"<<endl;
      bmat error(1,1);
      error.zeros();
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
