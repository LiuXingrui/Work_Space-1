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
  bmat m(r,c);
  m.zeros();
  int index=0;
  for (int i = 0; i <= r-1; ++i)
    {
      for (int j=0;j<=c-1 &&((i+1)*r+j+1)<=v.size();++j)
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


bmat burst_error_channel(const bmat &m,double initial_p)
{
  int r=m.rows();
  int c=m.cols();
  int my_random_number;
  bmat m1=m;
  bin ONE=1;
  double p=initial_p;


  
  int max_sidelength=4;
  for (int b = 2; b <= max_sidelength; ++b){
    p=p*0.1;
    for (int r1=0;r1<=r-1-b;++r1){
      for (int c1=0;c1<=c-1-b;++c1){
	my_random_number=randi(1,1000);
	if(my_random_number<=1000*p)
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
  int n,k,t,q,num_bits,num_lbits,num_pbits,L1,L2,i1,i2,i3;
  double n_cor_cws=0;
  double num_of_cws=0;
  bvec mas,cw, cor_cw1,cor_cw4,rec_cw,dec_m;
  bmat cor_cw2,cor_cw3;

 RNG_randomize();
num_of_cws=10.0;
 for (int m=10;m<=13;m++){
  n=pow(2,m)-1;
  q=n+1;
 
   L1=m*(pow(2,ceil_i(m/2+log(m)/log(2))));
   L2=pow(2,floor_i(m/2+log(m)/log(2)));
   

  i1=n/15;
  i2=n/10;
  i3=n/90;
  
   for (int d=i1;d<=i2;d=i3+d){
     // for (int d=1;d<=2;d++){
  
  t=(d-1)/2;
  k=n-d+1;
  num_lbits=m*k;
  num_pbits=m*n;
  Reed_Solomon RS(m,t);

  cout<<"d="<<d<<endl;
  for (double p=0.001;p<=0.001;p=p+0.001){
  

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
    cw=RS.encode(mas);

    
   cor_cw1=bsc(cw);
    cor_cw2=bvec_to_bmatrix(cor_cw1,L1,L2);
    
    cor_cw3=burst_error_channel(cor_cw2,p);
    rec_cw=bmat_to_bvector(cor_cw3,num_pbits);
    dec_m=RS.decode(rec_cw);
  if (dec_m!=mas)
    {
    n_cor_cws++;
  }
    
  }

  cout<<"for n="<<n<<", k="<<k<<", d= "<<d<<" Reed-Solomon code, raw bit error rate p="<<p<<", (b+1)*(b+1) burst error rate=(0.1)^b*p, max_b=4, the error rate is "<<n_cor_cws/num_of_cws<<"\n"<<endl;



 }

  

 }

 }

  return 0;
}
