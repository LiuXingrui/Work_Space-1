#include <iostream>
#include<cstring>
#include<string>
#include<vector>
#include <fstream>
#include <math.h> 
#include <random>

using namespace std;
#include <itpp/itbase.h>
#include <itpp/itcomm.h>
using namespace itpp;


bmat bvec_to_bmatrix(const bvec &v,int r,int c)
//convert bvec to bmat, if v.size<r*c, add 0s, if v.size>r*c, disregrad those entries
{

  //cout<<"vsize="<<v.size()<<endl;

  if (v.size()>r*c){

    cout<<"error in bvec_to_bmat: size of v is larger than r*c"<<endl;
    return 0;
  }
  bmat m(r,c);
  m.zeros();
  int index=0;
  for (int i = 0; i <= r-1 && index<v.size(); ++i)
    {
      for (int j=0;j<=c-1&& index<v.size() ;++j)
	{
	  m(i,j)=v[index];
	  index++;
	  // cout<<index<<endl;

	}

      
    }

  
  return m;
  
}

bvec bmat_to_bvector(const bmat &m,int s)
{

  int r=m.rows();
  int c=m.cols();
  if (s>r*c) {
    cout<<"error in bmat_to_bvec: s is larger than r*c"<<endl;
    return 0;
  }
  int k=0;
  int i,j;
  bvec v(s);
  v.zeros();
  int index=0;
  while(index< s){
    i=k/c;
    j=k%c;
    v[index]= m(i,j);
	  index++;
	  k++;

   }

      
   

  return v;
}


bmat burst_error_channel(const bmat &m,double p2,double p3,double p4)
//bmat burst_error_channel  (const bmat &m,double p2)
{
  int r=m.rows();
  int c=m.cols();
  int my_random_number;
  bmat m1=m;
  bin ONE=1;
  vec p(5);

  p.zeros();
  p[2]=p2;
  p[3]=p3;
  p[4]=p4;


  
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
		 m1(r2,c2)=m1(r2,c2)+randb();



	       }
	     }


	  }
	
  
 }


    }


  }

  return m1;

}

int main()
{
  int n,d,t,q,num_bits,num_lbits,num_pbits,L1,L2,i1,i2,i3;
  double n_cor_cws=0;
  double num_of_cws=0;
  double p_i,p_f;
  bvec mas,cw, cor_cw1,cor_cw4,rec_cw,dec_m;
  bmat cor_cw2,cor_cw3;

 RNG_randomize();
 //num_of_cws=200.0;
 num_of_cws=1;
 for (int m=10;m<=13;m++){
// for (int m=7;m<=7;m++){
 
 
   //m at least be 7, or i3 will be zero and the prog will run forever
  n=pow(2,m)-1;
  q=n+1;
 
   L1=m*(pow(2,ceil_i(m/2-log(m)/log(2))));
   L2=pow(2,m-ceil_i(m/2-log(m)/log(2)));
   

  i1=n*14/15;
  i2=n*9/10;
  i3=n/90;
  
   for (int k=i2;k<=i1;k=i3+k){
     // for (int d=1;d<=2;d++){
     
    if (k%2==0){
      k=k+1;

  }
          // there is a" bug",  when decoding [n=2^m-1,k] code,  if k is even, then the length of the decoded massage is mk+m (it should be mk), I guess the reason is for k and k+1,  t=(d-1)/2 may be the same, and the prog choose the k+1 case.
    //so we choose k always be odd

 
  d=n-k+1;
  t=(d-1)/2;


  
  num_lbits=m*k;
  num_pbits=m*n;
  Reed_Solomon RS(m,t);
  cout<<"\n"<<endl;

  // cout<<"d="<<d<<endl;
  p_i=0.001;
  p_f=0.01;
  // p_f=0.0015;
    cout<<"for n="<<n<<", k="<<k<<", d= "<<d<<" Reed-Solomon code"<<endl;
     cout<<"{ p1, p2,p3,p4, error rare after decoding}=";
    
    // cout<<"raw bit error rate=";
    //  for (double p=p_i;p<=p_f; p=p+0.001){
    //  cout<<p<<",";
       
    //}
    
       #pragma omp parallel for
       for (int i=1;i<=10; i++){
	 double p=0.001*i;
    for (double p2=p_i;p2<=0.5*p_f;p2=p2+0.001){
      for (double p3=p_i;p3<=0.3*p_f;p3=p3+0.001){
	for (double p4=p_i;p4<=0.2*p_f;p4=p4+0.001){
	

  BSC bsc(p);
  

  n_cor_cws=0;
  for (int i=0;i<num_of_cws;i++)
  {

     mas.clear();
    cw.clear();
    cor_cw1.clear();
    cor_cw2.clear();
    cor_cw3.clear();
    rec_cw.clear();
    dec_m.clear();
    mas=randb(num_lbits);
    

    
    //    cout<<"massize:  "<<mas.size()<<endl;
    cw=RS.encode(mas);
    //  cout<<n*m<<"  l of cw:  "<<cw.size()<<endl;
  
     cor_cw1=bsc(cw);

    
    cor_cw2=bvec_to_bmatrix(cor_cw1,L1,L2);
    
    cor_cw3=burst_error_channel(cor_cw2,p2,p3,p4);
    rec_cw=bmat_to_bvector(cor_cw3,num_pbits);
    
    //  if(cw!=rec_cw){
    //    cout<<"error";
    //    }
    
    dec_m=RS.decode(rec_cw);

    //  if (dec_m.size()>mas.size()){

    // cout<<"size different"<<endl;
    // }
    
      if (dec_m!=mas)
    {
    n_cor_cws++;
    
  }
    
  }
  
  cout<<"{"<<p<<", "<<p2<<", "<<p3<<", "<<p4<<", "<<n_cor_cws/num_of_cws<<"},";



 }
		       }
	    }
		 }
  

 }

 }

  return 0;
}
