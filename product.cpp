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
// using namespace Eigen;




void  burst_error_channel( bmat &m,double p1)
{
  int r=m.rows();
  int c=m.cols();
  int my_random_number;
  bin ONE=1;
  vec p(5);

  p.zeros();
  p[2]=p1*0.5;
  p[3]=p[2]*0.5;
  p[4]=p[3]*0.5;


  
  int max_sidelength=4;
  for (int b = 2; b <= max_sidelength; ++b){
    
    for (int r1=0;r1<=r-1-b;++r1){
      for (int c1=0;c1<=c-1-b;++c1){
	my_random_number=randi(1,1000);
	if(my_random_number<=1000*p[b])
	  {
	     for (int r2=0;r2<=b-1;++r2)
	       {
	       for (int c2=0;c2<=b-1;++c2)
	       {
		 m(r2,c2)=m(r2,c2)+randb();



	       }
	     }
	  }
	 }

 }
 }

  }

void product_encoder (const bmat &mas, int k1,bmat &cw, int L2,Reed_Solomon &RS, Hamming_Code & HC){



   for (int i1=0;i1<k1;i1++){

     cw.set_row(i1,RS.encode(mas.get_row(i1)));
  }
   //   cout<<cw<<"RSe\n"<<endl;

     bmat submat=cw(0,k1-1,0,L2-1);
     //  cout<<submat<<"sub\n"<<endl;
  for (int i2=0;i2<L2;i2++){

     cw.set_col(i2,HC.encode(submat.get_col(i2)));

     //  cout<<mas.get_col(i2)<<"\n"<<endl;
 }
  //  cout<<cw<<"HCe\n"<<endl;
}


void product_decoder (int L1,int L2,bmat &cor_cw,bmat &rec_mas,Reed_Solomon &RS, Hamming_Code & HC){

  int k1=rec_mas.rows();
  int ak2=rec_mas.cols();
  bvec temp_vec(ak2);
  bvec is_valid="0";
 
  //cout<<-1<<endl;
   for (int i1=0;i1<L2;i1++){
    
     cor_cw.set_col(i1,HC.decode(cor_cw.get_col(i1)));
     
   }
   //   cout<<cor_cw<<"HCd\n"<<endl;
 
  for (int i2=0;i2<k1;i2++){
   temp_vec= RS.decode(cor_cw.get_row(i2));
   //   cout<<temp_vec<<endl;
    //   cout<<3<<endl;
    rec_mas.set_row(i2,temp_vec);
   }
}

void random_channel ( bmat &cw,bmat &cor_cw1,const int L1,BSC &bsc){

  for (int i=0;i<L1;++i){
    cor_cw1.set_row(i,bsc(cw.get_row(i)));
  }

}

  


void product_simulation (const int n1,const int n2,const int k1,const int k2,const int a,const int m1,const double pmin,const double pmax,int num_of_p,const double num_of_cws){

  int n,k,d,t,q,num_lbits_RS,L1,L2,i1,i2,i3,d1,d2;
  double num_of_cor_cws=0;
  int temp=num_of_cws;
  vec   counter(num_of_cws);
  double x;
  double pintv=(pmax-pmin)/(num_of_p-1);
  // cout<<pintv<<endl;
  RNG_randomize();

   n=n1*n2*a;
   //   cout<<floor(n*9/10)<<"\n"<<floor(14/15*n)<<endl;


   // cout<<k<<endl;
   // cout<<"n1="<<n1<<"  n2="<<n2<<endl;
 
   
   k=k1*k2*a;
   // cout<<"  "<<k<<endl;
    
 
	  
	  
  d2=n2-k2+1;
  d1=3;
  t=(d2-1)/2;

   L1=n1;
  L2=a*n2;
  n=L1*L2;

 Array<bmat> cw(num_of_cws),cor_cw(num_of_cws),mas(num_of_cws),rec_mas(num_of_cws);
  bmat tempmat1(L1,L2),tempmat2(k1,a*k2);
  tempmat1.zeros();
  tempmat2.zeros();

    
  Reed_Solomon RS(a,t);
  Hamming_Code  HC(m1);



  //cout<<"k1="<<k1<<"  k2="<<k2<<"  a="<<a<<endl;
  cout<<"  n_1="<<n1<<", k_1="<<k1<<", d_1= "<<d1<<", n_2="<<n2<<", k_2="<<k2<<", d_2= "<<d2<<", n="<<n<<", k="<<k<<" Product code:ListPlot[{"<<endl;

 for (double p=pmin;p<=pmax+0.000000001;p=p+pintv){
     BSC bsc(p);
     num_of_cor_cws=0;
     counter.zeros();
            #pragma omp parallel for
  for (int i=0;i<temp;++i){
    mas(i)=randb(k1,a*k2);
    cw(i)=tempmat1;
    cor_cw(i)=tempmat1;
    rec_mas(i)=tempmat2;
   
  // cout<<mas<<"\n"<<endl;
  
    product_encoder (mas(i), k1,cw(i),  L2,RS, HC);
    random_channel (cw(i),cor_cw(i),L1,bsc);
    burst_error_channel( cor_cw(i), p);
	  
    product_decoder ( L1, L2,cor_cw(i),rec_mas(i),RS, HC);
  // cout<<rec_mas<<endl;
    if (rec_mas(i)!=mas(i)){

    counter(i)=1.0;
    //    cout<<counter(i);
	       
  }

  
  }
  for (int i=0;i<num_of_cws;i++){
    // cout<<counter(i)<<endl;
    num_of_cor_cws=num_of_cor_cws+counter(i);
    //  cout<<num_of_cor_cws<<endl;
    

  }
    cout<<"{"<<p<<", "<<num_of_cor_cws/num_of_cws<<"},";
 }
 double ndou=n*1.0;
 double kdou=k*1.0;
 double rate=kdou/ndou;
 cout<<"},Joined->True,AxesLabel->{RBER,UBER},PlotLabel->\"["<<n<<", "<<k<<", Rate="<<rate<<"] Product code\"]\n"<<endl;
  

}

int main()
{
  double pmin=pow(10,-6);
  //  cout<<pmin<<endl;
  
  double pmax=5*pow(10,-3);
  // cout<<pmax<<endl;
  double num_of_cws=100.0;
  int num_of_p=20;

  
  product_simulation (127,15,120,13,4,7,pmin,pmax,num_of_p,num_of_cws);
    product_simulation (63,15,57,13,4,6,pmin,pmax,num_of_p,num_of_cws);
      product_simulation (127,31,120,29,5,7,pmin,pmax,num_of_p,num_of_cws);
        product_simulation (31,15,26,13,4,5,pmin,pmax,num_of_p,num_of_cws);
    //  cout<<1<<endl;
  
  return 0;
}
