#include <iostream>
#include<cstring>
#include<string>
#include<vector>
#include <fstream>
#include <math.h> 
#include <random>
//#include<eigen3/Eigen/Dense>

using namespace std;
#include <itpp/itbase.h>
#include <itpp/itcomm.h>
using namespace itpp;



class  nodes{
public:
  ivec neighbors;
  int degree;


};


  
void initialize_checks (const bmat &H, nodes checks[], int & E){
  int r=H.rows();
  int c=H.cols();
  for (int i=0; i<r;i++){
    checks[i].degree=0;
    for (int j=0; j<c;j++){
     
      if (H(i,j)==1){
	checks[i].degree++;
	(checks[i].neighbors).ins(0,j);
	E++;
	
      }
     
    }
 
    //    cout<<"check  "<<i<<"  degree:"<<checks[i].degree<<endl;
    //   cout<<"check  "<<i<<"  neighbors:"<<checks[i].neighbors<<endl;
    
  }
}

void initialize_errors(const bmat &H, nodes errors[]){

    int r=H.rows();
    int c=H.cols();
    int index=0;
  


    for (int i=0; i<c;i++){
      errors[i].degree=0;
      for (int j=0; j<r ;j++){

	if (H(j,i)==1){
	  errors[i].degree++;
	  (errors[i].neighbors).ins(0,j);
	  //  massages[index].v=i;
	  // massages[index].c=j;
	  // massages[index].Mcv=1;
	  // massages[index].Mvc=1;
	  // massages[index].pre_Mcv=1;
	  // massages[index].pre_Mvc=1;
	  //	  cout<<massages[index].v<<"  to c:   "<<massages[index].c<<endl;
	  index++;
	
	}
	
      }
      //    cout<<"errors  "<<i<<"  degree:"<<errors[i].degree<<endl;
      //   cout<<"errors  "<<i<<"  neighbors:"<<errors[i].neighbors<<endl;

    }

    //  cout<<index<<endl;
}

void initialize_massages(mat &mcv,mat& mvc, bmat &H){
    int r=H.rows();
    int c=H.cols();
    for (int i=0;i<r;i++){
      for (int j=0;j<c;j++){
	  if (H(i,j)==1){
	    mcv(i,j)=1;
	    mvc(i,j)=1;
	 
	  }
	}
	     
    }
    
  
}


void update(const nodes checks[],const nodes errors[],mat &mcv,mat& mvc,const bmat& syndrome,double p,int c, int v,  bmat& output_e,int lll){
    double ipr=(1-p)/p;
    double l_p=log(ipr);
    //  cout<<"vc:"<<endl;
  for (int i=0;i<c;i++)
    {
      int tempdegree=checks[i].degree;
    for (int j=0;j<tempdegree;j++)
      {
      
      
      int vnode=(checks[i].neighbors)(j);
	//update all v-to-c massage first:
	mvc.set(i,vnode,ipr);
	int tempdegree2=errors[vnode].degree;
      for (int k1=0;k1<tempdegree2;k1++)
	{
	  
	int temp=(errors[vnode].neighbors)(k1);
	//	cout<<"temp="<<temp<<endl;
	//	cout<<"i-"<<i<<endl;
	if (temp!=i)

	  {//int temp33=mvc(i,vnode);
	  
	  mvc.set(i,vnode,mvc(i,vnode)*mcv(temp,vnode));
	  //  if (mvc(i,vnode)<-100){
	  //	    cout<<mvc(i,vnode)<<endl;
	  //  cout<<"got a inf: " <<i<<" , "<<vnode<<"\n"<<endl;
	  //    cout<<"before "<<temp33<<endl;
	  //    cout<<"l="<<lll<<endl;
	  //   }
	  // cout<<pre_mcv<<endl;
	  //  cout<<"k1="<<k1<<endl;
	  // cout<<pre_mcv(k1,j)<<endl;
	  
	
	}
      }
      //  cout<<i<<"  "<<vnode<<"  "<<mvc(i,vnode)<<endl;
      }
    }
  
  // cout<<"cv:"<<endl;
  for (int j=0;j<v;j++)
    {
   int tempdegree=errors[j].degree;
   double  final_pr=ipr;
   // double final_lp=lp;
   for (int i=0;i<tempdegree;i++)
     {
  
        int cnode=(errors[j].neighbors)(i);
   

      //update c-to-v massage:
      double temp2=1.0;
      
      int tempdegree2=checks[i].degree;
      for (int k2=0;k2<tempdegree2;k2++)
	{
	//	cout<<checks[i].degree<<endl;
	int temp=(checks[cnode].neighbors)(k2);
	
	if (temp!=j)

	  {
	    temp2=temp2*(mvc(cnode,temp)-1)/(mvc(cnode,temp)+1);
	  
	  //	  if (tanh(mvc(cnode,temp)/2)<=-1){
	    
	  //  cout<<"wrong:"<<tanh(mvc(cnode,temp)/2)<<endl;
	  //   cout<<mvc(cnode,temp)/2<<endl;
	  //    cout<<"l="<<lll<<endl;
	  //   }
	    //  cout<<pre_mvc<<endl;
	    // cout<<pre_mvc(i,k2);
	    // cout<<"temp2="<<temp2<<endl;
	  
	 }
      }
      //    if (temp2<=-1){
      //	cout<<temp2<<"wrong temp2"<<endl;
      //  }
      if (syndrome(cnode,0)==0)
	{

	  mcv.set(cnode,j,(1+temp2)/(1-temp2));
      }

      else
	{
          mcv.set(cnode,j,(1-temp2)/(1+temp2));
      }
      //  cout<<cnode<<"  "<<j<<"   "<<mcv(cnode,j)<<endl;
      final_pr=final_pr*mcv(cnode,j);
    
    
      }
   
   if (final_pr<1)
     {//cout<<final_l<<" final l"<<endl;
       // cout<<"?????????????"<<endl;
       cout<<j<<"   "<<final_pr<<endl;
       output_e(j,0)=1;
       //  cout<<final_pr<<endl;
     }
   else {
     output_e(j,0)=0;
   }
      
    }
 
   
   


  
  
  
}

void decode(mat& mcv, nodes errors[],bmat& e,double p,int v){
  double l_p=log((1-p)/p);
  for (int i=0;i<v;i++){
    double final_l=l_p;
    for (int j=0;j<errors[i].degree;j++){
      final_l=final_l+mcv((errors[i].neighbors)(j),i);

    }
    if (final_l<0){
      e(i,0)=1;
    }

    

  }
  




}

int n_bit_error(const bmat &output_e, const bmat &real_e,int n){
  int er=0;

  for (int i=0;i<n;i++){

    if( (output_e(i,0)!=real_e(i,0))){
      er++;
    }
  }

  return er;
}
      


void where_error( const bmat &real_e,int n){
  int er=0;

  for (int i=0;i<n;i++){

    if (real_e(i,0)!=0){
	cout<<" "<<i<<" ";
    }
  }


}
  




  

int main(){
  GlobalRNG_randomize ();
  int num_iter_end=0;
  double p=0.01;
  int lmax=20;
  int num_of_cws=10000;
  int num_of_suc_dec=0;
  int n_valid_cws=0;
  
  int n;
  int k;
  int r;
  int temp;
  int col_ind=0;
  int num_zero_cws=0;
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
      //  cout<<temp<<"  ";
  }
    // cout<<"\n"<<endl;
    col_ind++;
  }

    //  Sparse_Mat<bin> SH(H);
    //   cout<<SH<<endl;

    //	  */
     
  //  cout<<1<<endl;
  parity_check.close();
  //  cout<<2<<endl;
 
  // cout<<"H:"<<H<<endl;

  nodes  checks[c];
  nodes  errors[v];
  BSC bsc(p);
  int E=0;
  initialize_checks (H, checks,  E);
  initialize_errors(H, errors);

  int er=0;
  for (int s=0;s<num_of_cws;s++){
   
    //  bmat H=randb(3,5);
  
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
    
    
   
    if (syndrome==zero_mat1){
      n_valid_cws++;
     if (real_eT==zero_vec){
      num_of_suc_dec++;
      num_zero_cws++;
      // cout<<"success!  error=0"<<endl;
      
    }
     else {
       er=er+ n_bit_error(zero_mat2, real_e, n);
       //  cout<<"error wt: "<<n_bit_error(zero_mat2, real_e, n)<<endl;
       //  cout<<"er:  "<<er<<endl;
       // cout<<"failure! error= some codeword"<<endl;
     }
    
     continue;
     
    }
    //   cout<<syndrome<<"syndrome"<<endl;
  
  

    mat mcv(c,v);
    mat mvc(c,v);
      initialize_massages( mcv,mvc, H);
    

  

    bmat output_e(v,1);
    output_e.zeros();
      
    for (int l=1;l<=lmax;l++){
      //  cout<<l<<endl;
      update(checks,errors, mcv,mvc,syndrome,p, c, v,output_e,l);
      // decode(mcv,  errors,output_e, p, v);
    
      if (H*output_e==syndrome){
	 n_valid_cws++;
	if(output_e==real_e){
	  // cout<<"success! iteration number="<<l<<endl;
	  num_of_suc_dec++;
	  //  cout<<num_of_suc_dec<<endl;
	}
	else{
	  er=er+ n_bit_error(output_e, real_e, n);
	  //  cout<<n_bit_error(output_e, real_e, n)<<endl;
	    //     cout<<"error wt: "<<n_bit_error(zero_mat2, real_e, n)<<endl;
	    //    cout<<"er:  "<<er<<endl;
	  //  cout<<"Failure! get an output e such that H*e=s, but this is not the right error"<<endl;
	}
        
	//	cout<<l<<endl;
	  break;
	  
      }
      
    
	//update pre_massage:
	   if(l==lmax){
	     num_iter_end++;
	     er=er+ n_bit_error(output_e, real_e, n);
	     cout<<"\n real error:"<<endl;
	     where_error(real_e,n);
	     cout<<"\n output e"<<endl;
	     where_error(output_e,n);
	     //  cout<<n_bit_error(output_e, real_e, n)<<endl;
	     // cout<<"error wt: "<<n_bit_error(zero_mat2, real_e, n)<<endl;
	     //  cout<<"er:  "<<er<<endl;
	     // cout<<"failure, iterations finished and did not get an e "<<endl;
	     //    cout<<real_e<<endl;
      }
	    

      }
   
    
  
      //  cout<<"H:"<<H<<endl;
    // cout<<"syndrome: "<<syndrome<<endl;
    // cout<<"outpute:"<<output_e<<endl;
    //  cout<<"reale:"<<real_e<<endl;
  }
  
  cout<<"for p="<<p<<", there are total "<< num_of_suc_dec<<" successful decoding out of "<< num_of_cws<<" cws for a ["<<v<<", "<<v-c<<"] code"<<endl;
  cout<<"and there are "<< n_valid_cws<<" decoding results are codewords"<<endl;

  cout<<"bit error rate after decoding is ";
  cout<<er/(1.0*n*num_of_cws)<<endl;
  
  //  /*  
  //  cout<< num_iter_end<<endl;
  cout<< "number of real_error=0:"<<num_zero_cws<<endl;



  //  */
  
}



    
    

 
