#include <iostream>
#include<cstring>
#include<string>
#include<vector>
#include <fstream>
#include <math.h> 
#include <random>
#include <sstream>
#include <stdlib.h> 
using namespace std;
#include <itpp/itbase.h>
#include <itpp/itcomm.h>
using namespace itpp;

bool OSD(vec& LR,const GF2mat& H,int r,const GF2mat& syndrome){ //r is the rank of H
  cout<<"LR is"<<LR<<endl;
  cout<<"H is"<<H<<endl;
  
  cout<<"syndrome is"<<syndrome<<endl;
  
  //get the first permutation that abs_LLR is in descending order
  int n=LR.length();
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
	  cout<<"error:get a negative Likelyhood-ratio in OSD function"<<endl;
	  return false;
	}
    }

  cout<<"ok, no negative likelyhood-ratia, please comment those debug lines before runninng on the cluster"<<endl;

  cout<<"nega_abs_LLR is"<<nega_abs_LLR<<endl;
  
  ivec perm1=sort_index(nega_abs_LLR);
  cout<<"perm1 is"<<perm1<<endl;
  GF2mat H1=H;
  GF2mat HH=H;
  H1.permute_cols(perm1,1);
  cout<<"after perm1,1, H1 is"<<H1<<endl;
  HH.permute_cols(perm1,0);
    cout<<"after perm1,0, H1 is"<<HH<<endl;
  
  //get the second permutation that gives first rank(H) linear_independent cols, the result is H2, which is Hs in that paper
  GF2mat H2(r,r);
  bvec zero_vec(r);
  ivec perm2(n);
  ivec perm2_copy(n);
  zero_vec.zeros();
  H2.set_col(0,H1.get_col(0));
  int j1=1;// the current col of H2
  int j2=1; //the current col of H1

  for (int i=0;i<n;i++){perm2(i)=i;}

  while (j1<r)
    {
      H2.set_col(j1,H1.get_col(j2));
      	  cout<<j1<<H2<<endl;
	  cout<<"H2 rank="<<H2.row_rank()<<endl;
      if (H2.row_rank()==j1+1)
	{
	  j1++;
	  j2++;

	}
      else
	{
	  H2.set_col(j1,zero_vec);
	  perm2_copy=perm2;
	  for (int i=j2;i<n-1;i++)
	    {
	      perm2(i)=perm2_copy(i+1);	      
	    }
	  perm2(n-1)=perm2_copy(j2);
	  cout<<"perm2 is "<<perm2<<endl;
	  j2++;  
	}      
    }
  cout<<"perm2 is"<<perm2<<endl;
  cout<<"H2 is"<<H2<<endl;
  H1.permute_cols(perm2,1); // then H1= [H2,HT]
  cout<<"after perm2,H1 is"<<H1<<endl;
  
  GF2mat HT=H1.get_submatrix(0,r,r-1,n-1);
  cout<<"HT is"<<HT<<endl;
  GF2mat H2_inv=H2.inverse();
  cout<<"H2_inv="<<H2_inv<<endl;
  
  GF2mat e_S=H2_inv*syndrome;
  GF2mat e_T(n-r,1);
  return true;
  
}




int main(){

  vec LR=" 5 7 9 1 ";
  bmat H1="0 0 1 1; 0 1 1 0; 1 1 0 0;";
  int r=3;
  GF2mat H(H1);
  GF2mat G(H1);
  bmat s="1;1;1";
  GF2mat syndrome(s);
  OSD( LR, H, r, syndrome);


  return 0;
}
