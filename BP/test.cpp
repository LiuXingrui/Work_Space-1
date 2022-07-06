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

bool OSD(vec & LR,const GF2mat H){
  int n=LR.length();
  vec nega_abs_LLR(n); //the sort function gives a ascending  order, we need descending order
  for (int i=0;i<n;i++)
    {
      nega_abs_LLR(i)=-abs(log(LR(i)));
    }

  cout<<"nega_abs_LLR is"<<nega_abs_LLR<<endl;
  ivec perm_vec=sort_index(nega_abs_LLR);


  cout<<"perm_vec is"<<perm_vec<<endl;


  GF2mat H1=H;
  H1.permute_cols(perm_vec,1);
  cout<<H1<<endl;



}




int main(){
    cout<<"log1 is "<<log(1)<<endl;
  cout<<"abs log1 is "<<abs(log(1))<<endl;
  vec a=" 2 0 -1 ";
  bmat H1="0 0 0; 0 0 1; 0 1 1";
  GF2mat H(H1);
  cout<<"H= "<<H<<endl;
  cout<<"a is "<<a<<endl;
  OSD(a,H);
 

  return 0;
}
