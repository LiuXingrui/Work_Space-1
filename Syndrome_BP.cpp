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
/*
void initialize_massages( &mcv,mat& mvc, bmat &H){
    int r=H.rows();
    int c=H.cols();
    for (int i=0;i<r;i++){
      for (int j=0;j<c;j++){
	  if (H(i,j)==1){
	    mcv(i,j)=0;
	    mvc(i,j)=0;
	 
	  }
	}
	     
    }
    
  
}
*/

   void update(nodes checks[],nodes errors[],Sparse_Mat<double> &mcv,Sparse_Mat<double>& mvc,const bvec& syndrome,double p,int c, int v, const bmat &H){
  for (int i=0;i<c;i++){
    for (int j=0;j<v;j++){
      
      if (H(i,j)==1){
	
	//update all v-to-c massage first:
	mvc.set(i,j,log((1-p)/p));
      for (int k1=0;k1<errors[j].degree;k1++){
	int temp=(errors[j].neighbors)(k1);
	//	cout<<"temp="<<temp<<endl;
	//	cout<<"i-"<<i<<endl;
	if (temp!=i){
	  mvc.set(i,j,mvc(i,j)+mcv(temp,j));
	  // cout<<pre_mcv<<endl;
	  //  cout<<"k1="<<k1<<endl;
	  // cout<<pre_mcv(k1,j)<<endl;
	  
	
	}
      }
      }
    }
  }

  for (int j=0;j<v;j++){
    double belief_j=log((1-p)/p);
   for (int i=0;i<c;i++){
   
      
      if (H(i,j)==1){

      //update c-to-v massage:
      double temp2=1;
      for (int k2=0;k2<checks[i].degree;k2++){
	//	cout<<checks[i].degree<<endl;
	int temp=(checks[i].neighbors)(k2);
	
	if (temp!=j){
	  temp2=temp2*tanh(mvc(i,temp)/2);
	    //  cout<<pre_mvc<<endl;
	    // cout<<pre_mvc(i,k2);
	    // cout<<"temp2="<<temp2<<endl;
	  
	 }
      }
      
      if (syndrome(i)==0){

	mcv.set(i,j,2*atanh(temp2));
      }

      else{
	mcv.set(i,j,-2*atanh(temp2));
      }


      
    
    
      }
      
    }
 
   
   
  }


  
  
  
}

void decode(Sparse_Mat<double>& mcv, nodes errors[],bmat& e,double p,int v){
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

int main(){
  GlobalRNG_randomize ();

  int lmax=1000;
  int n;
  int k;
  int temp;
  int col_ind=0;

  ifstream parity_check("LDPC_Mat1.txt");

  string line;
  getline(parity_check, line);
  istringstream iss(line);
  iss>>n>>k;

  int v=n;
  int c=n-k;
   bmat H(n-k,n);
    // bmat H="1 1 0;0 1 1";
  // bmat H= "1 0 0 1 1 1 0;  0 1 0 1 1 0 1;0 0 1 0 1 1 1";
     v=H.cols();
     c=H.rows();
     
     //     /*

    while( getline(parity_check, line)){
    istringstream iss2(line);
    while (iss2>>temp){
      H(col_ind,temp)=1;
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

  for (int s=0;s<10;s++){
    double p=0.01;
    //  bmat H=randb(3,5);
  
    bvec syndrome;
    
    nodes  checks[c];
    nodes  errors[v];
    bmat syndromeT(c,1);
    int E=0;
    initialize_checks (H, checks,  E);
    Sparse_Mat<double> mcv(c,v);
    Sparse_Mat<double> mvc(c,v);
    mcv.zeros();
    mvc.zeros();
    initialize_errors(H, errors);
  

       syndrome=randb(c);
        for (int i=0;i<c;i++){

	 
	 syndromeT(i,0)=syndrome(i);

	 
       }

      // syndrome="0 1";
	bmat e(v,1);
      e.zeros();
      
      for (int l=1;l<=lmax;l++){
	//  cout<<l<<endl;
	update(checks,errors, mcv,mvc,syndrome,p, c, v,H);
	  decode(mcv,  errors, e, p, v);
	//	if (syndrome=="0 1"){
	
	//	  cout<<"c to v:"<<mcv<<"\n<<"<<endl;
	//	  cout<<"v to c:"<<mvc<<"\n<<"<<endl;
	//	}
	
	if(H*e==syndromeT){
	  cout<<"success! iteration number="<<l-1<<endl;
	  	  break;
	}
    
	//update pre_massage:
	   if(l==lmax){

	cout<<"failure"<<endl;
      }

      }
   
    
  
      //  cout<<"H:"<<H<<endl;
      //  cout<<"syndrome: "<<syndrome<<endl;
      // cout<<"e:"<<e<<endl;
  }
  
}



    
    

 

