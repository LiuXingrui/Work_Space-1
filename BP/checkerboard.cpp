#include <iostream>
#include<cstring>
#include<string>
#include<vector>
#include <fstream>
#include <math.h> 
#include <random>
#include <sstream>
#include <stdlib.h> 

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

int main(int argc, char **argv)
{

  if (argc!=4)
    {
      cout<<"wrong input parameters:  ./a d Hx_file Hz_file"<<endl;
      return 1;
    }

  int d;
  istringstream argv1( argv[1] );
  if ( argv1 >> d){}
  else
    {
      cout<<"d should be an int"<<endl;
      return 1;
    }
  if (d<3)
    {
      cout<<"d should larger or equal to 3"<<endl;
      return 1;
    }

  if (d%2==0)
    {
      cout<<"d should be odd"<<endl;
      return 1;
    }

  int n=d*d;   //9
  int k=1;
  int r=(n-k)/2;   //8
  GF2mat Hx(r,n);
  GF2mat Hz(r,n);
  int rx=0;
  int rz=0;
  
  int t1,t2,t3,t4;

  //in the middle:
  for (int i=1;i<d;i++)
    {
      t1=(i-1)*d;
      t2=t1+1;
      t3=t1+d;
      t4=t2+d;
      while (t1<i*d-1)
	{

	  if (i%2==1)
	    {
	      
	      Hz.set(rz,t1,1);
	      Hz.set(rz,t2,1);
	      Hz.set(rz,t3,1);
	      Hz.set(rz,t4,1);
	      cout<<"Hz row"<<rz<<"  "<<t1<<"  "<<t2<<"  "<<t3<<"  "<<t4<<endl;
	      rz++;
	      t1++;t2++;t3++;t4++;
	      Hx.set(rx,t1,1);
	      Hx.set(rx,t2,1);
	      Hx.set(rx,t3,1);
	      Hx.set(rx,t4,1);
	      cout<<"Hx row"<<rx<<"  "<<t1<<"  "<<t2<<"  "<<t3<<"  "<<t4<<endl;
	      t1++;t2++;t3++;t4++;
	      rx++;
	    }
	  else
	    {
	      Hx.set(rx,t1,1);
	      Hx.set(rx,t2,1);
	      Hx.set(rx,t3,1);
	      Hx.set(rx,t4,1);
	      cout<<"Hx row"<<rx<<"  "<<t1<<"  "<<t2<<"  "<<t3<<"  "<<t4<<endl;
	      t1++;t2++;t3++;t4++;
	      rx++;
	      Hz.set(rz,t1,1);
	      Hz.set(rz,t2,1);
	      Hz.set(rz,t3,1);
	      Hz.set(rz,t4,1);
	      cout<<"Hz row"<<rz<<"  "<<t1<<"  "<<t2<<"  "<<t3<<"  "<<t4<<endl;
	      rz++;
	      t1++;t2++;t3++;t4++;
	    }
	}     
    }

  
  //upper side:
  t1=0;
  t2=1;
  while(t1<d-1)
    {
      Hx.set(rx,t1,1);
      Hx.set(rx,t2,1);
      cout<<"Hx row"<<rx<<"  "<<t1<<"  "<<t2<<endl;
      rx++;
      t1=t1+2;
      t2=t2+2;
    }

  //bottom side:
  t1=(d-1)*d+1;
  t2=t1+1;
  while (t2<n)
    {
      Hx.set(rx,t1,1);
      Hx.set(rx,t2,1);
      cout<<"Hx row"<<rx<<"  "<<t1<<"  "<<t2<<endl;
      rx++;
      t1=t1+2;
      t2=t2+2;

    }

  if (rx!=r){cout<<"error! rx!=r!"<<endl;return 1;}

  //left side:
  t1=d;
  t2=2*d;
  while (t2<(d-1)*d+1)
    {
      Hz.set(rz,t1,1);
      Hz.set(rz,t2,1);
      cout<<"Hz row"<<rz<<"  "<<t1<<"  "<<t2<<endl;
      rz++;
      t1=t1+2*d;
      t2=t2+2*d;
    }
  
  //right side:
  t1=d-1;
  t2=2*d-1;
  while (t2<(d-1)*d)
    {
      Hz.set(rz,t1,1);
      Hz.set(rz,t2,1);
      cout<<"Hz row"<<rz<<"  "<<t1<<"  "<<t2<<endl;
      rz++;
      t1=t1+2*d;
      t2=t2+2*d;
    }
  if (rz!=r){cout<<"error! rz!=r!"<<endl;return 1;}


  GF2mat zero_mat1(rx,rz);
  if((Hx*Hz.transpose())==zero_mat1) {}
  else{cout<<"Hx*Hz^T!=0, the two matrices donot match"<<endl;cout<<"Hx:\n"<<Hx<<" \n Hz:\n"<<Hz<<endl;return 1;}
  cout<<"Hx:\n"<<Hx<<endl;
  cout<<"Hz:\n"<<Hz<<endl;
  
  string data_file1=argv[2]; 
  string data_file2=argv[3];
  ofstream myfile1;
  write_matrix(data_file1, Hx);
  write_matrix(data_file2, Hz);
  return 0;

}
