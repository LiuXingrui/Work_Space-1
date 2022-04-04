#include <iostream>
#include<cstring>
#include<string>
#include<vector>
#include <fstream>
#include <math.h> 
#include <random>
#include <sstream>

using namespace std;


int main(int argc, char **argv){
  int L;
  int k=2;

  if (argc!=3){

    cout<<"please input exactly 2  parameter which is the side length of the surface code"<<endl;
    return 1;
  }
  
  string file_name=argv[1]; 
  istringstream argv2(argv[2]);
  if (argv2>>L){}
  else {cout<<"the input is not a int"<<endl; return 1;}
 
  int n=2*L*L;
  int r=(n-k)/2;
  ofstream Hx;
  Hx.open (file_name,ios::trunc);



  Hx<<n<<" "<<r<<endl;

  for (int j=0;j<L;j++)
    {
      for (int i=1; i<=L;i++)
	{
	  Hx<<i+2*L*j<<" "<<i%L+1+2*L*j<<" "<<((i-L+n)%n+2*L*j-1+n)%n+1<<" "<<(i+L)%n+2*L*j<<endl;
	}
    }
  

  Hx.close();
  return 0;
}
