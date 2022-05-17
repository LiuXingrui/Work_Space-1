#include <iostream>
#include<cstring>
#include<string>
#include<vector>
#include <fstream>
#include <math.h>

using namespace std;
#include <itpp/itbase.h>
#include <itpp/itcomm.h>
using namespace itpp;

void write_matrix(string file_name, GF2mat &H){

  ofstream Hx;
  Hx.open (file_name,ios::trunc);
  int n=H.cols();
  int r=H.rows();

  Hx<<n<<" "<<r<<endl;

  for (int j=0;j<r;j++)
    {
      for (int i=0; i<n;i++)
	{
	  if (H(j,i)!=0)
	    {
	      Hx<<i+1<<" ";
	    }
	  
	}
      Hx<<endl;
    }
  
  Hx.close();
}



GF2mat read_mtx(string file_name){
  ifstream parity_check(file_name);
  string line;
  int row_ind=0;
  int cols,rows,entries;
  
  for (int i=0;i<4;i++)
    {
      getline(parity_check, line);
    }
  
  getline(parity_check, line);
  istringstream iss(line);
  if (iss>>rows>>cols>>entries){}

  GF2mat H(rows,cols);
  int temp;
  for (int i=0;i<entries;i++)
    {
       getline(parity_check, line);
       istringstream iss(line);
       int row_ind,col_ind;
       iss>>row_ind>>col_ind;
       iss>>temp;
       H.set(row_ind-1,col_ind-1,1);
    }
  return H;
}
  
int main(int argc, char** argv){
  string mtx_file=argv[1];
  string conv_file=argv[2];

  GF2mat H;
  H=read_mtx(mtx_file);
  write_matrix(conv_file,H);


  return 0;

 }

 
