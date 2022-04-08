#include <iostream>
#include<cstring>
#include<string>
#include<vector>
#include <fstream>
#include <math.h> 
#include <random>
#include"BP_quantum_header.h"
#include <sstream>

using namespace std;
#include <itpp/itbase.h>
#include <itpp/itcomm.h>
using namespace itpp;


void error_channel(bvec &cw, const vec &p){
  int temp=1e8;
  int temp2;
  bin one=1;
  if (cw.size()!=p.size())
    {cout<<"the size of p and cw do not match"<<endl;}

  else
    {
    for (int i=0;i<cw.size();i++)
      {    
	temp2=randi(1,temp);
	if(temp2<p[i]*temp)
	  {
	    cw[i]=cw[i]+one;
	  }	
      }
  }  
}


void pro_dist(double p, vec& pv){

  double minp=0.5*p;
  double maxp=1.5*p;
  int pvsize=pv.size();
  double temp;
  int temp2=1e8;

  for (int i=0;i<pvsize;i++)
    {
      temp=p*randi(1,temp2)/temp2;
      pv(i)=minp+temp;
    }
}
