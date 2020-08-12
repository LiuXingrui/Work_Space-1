 

#include <iostream>
#include<cstring>
#include<string>
#include<vector>
#include <fstream>
#include<random>
#include <chrono>



using std::cin;
using std::cout; 
using std::vector;
using std::endl;

vector<unsigned> random_binary_poly (const int n,const int w){
  static std::default_random_engine e;
  e.seed(std::chrono::system_clock::now().time_since_epoch().count());
  static std::uniform_int_distribution<unsigned> u(0,n-1);
  vector<unsigned> index;
  vector<unsigned> ret;
  unsigned temp;
  bool bool_temp=true;

  //randomly generate w positions where the value is 1

  for(size_t i=0;i<w;++i){
    temp=u(e);
    //check if there are two same positions generated
    for(auto it=index.begin();it!=index.end();++it){
      if(*it==temp){
	bool_temp=false;
      }
    }

    if(bool_temp==true){
      index.push_back(temp);
    }
    else {
      --i;
      bool_temp=true;
    }
  }
  bool_temp=true;
  for(int i=0;i<n;++i){
    for (auto it=index.begin();it!=index.end();++it)
      if (*it==i){
	bool_temp=false;

      }
    if(bool_temp==false){
      ret.push_back(1);
      bool_temp=true;
  
    }
    else{
      ret.push_back(0);
    }
  }
  return ret;
}


int main(int argc, char **argv){
  const int n1=6;
  const int n2=6;
  const int w1=4;
  const int w2=4;


  vector<unsigned> check_poly_1=random_binary_poly(n1,w1);
  vector<unsigned> check_poly_2=random_binary_poly(n2,w2);

  for (int i=0;i<10;++i){
    for (auto it=check_poly_1.begin();it!=check_poly_1.end();++it){
      cout<<*it;

    }  
    cout<<endl;
    check_poly_1=random_binary_poly(n1,w1);
}
 cout<<"test";

return 0;
}
