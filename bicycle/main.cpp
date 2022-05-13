#include <iostream>
#include <fstream>
#include<cmath>
#include <itpp/itbase.h>
#include <chrono>
#include "func.hpp" 

using namespace std;
//To run: make bike
//        ./bike.out

int main(int argc, char ** argv){

    std::chrono::high_resolution_clock::time_point start = std::chrono::high_resolution_clock::now();

    int tempcof, rank, d, Ad;
    long coff, g_deg;
    long verbose = 0;
    //To change: l = n/2
    //int l = 31;
    //
    std::string file_name1=argv[1];
    std::string file_name2=argv[2];
    std::istringstream argv3( argv[3] );
    int l;
    
    if (argv3>>l){}
    else{std::cout<<"l need to be an int"<<std::endl;return 1;}

  
    
     std::istringstream argv4( argv[4] );
    int a_wt;
    if (argv4>>a_wt){}
    else{std::cout<<"a_len need to be an int"<<std::endl;return 1;}

    std::istringstream argv5( argv[5] );
    int b_wt;
    if (argv5>>b_wt){}
    else{std::cout<<"b_len need to be an int"<<std::endl;return 1;}

    
    //  if (argc!=a_len+b_len+6){std::cout<<"the number of input parameters is not right, example: ./bike.out Hx_file Hz_file l a_len b_len a0 a1 ... b0 b1 ..., you input"<<argc-1<<"parameters,  but a_len+b_len+5="<<a_len+b_len+5<<std::endl;return 1;}
    itpp::bin tempcof2;
    itpp::GF2mat Aa(l,l), Bb(l, l), TAa(l,l), TBb(l,l), Hx(l,2*l), Hz(l,2*l);
    itpp::bvec vtemp(l), a(l), b(l);

    //construct polinomial a
    //To change: a(x) = ...
    a.zeros();
   for (int i=0;i<a_wt;i++)
     {  int temp;
	 std::istringstream argvi( argv[i+6] );
	 if (argvi>>temp){a(temp)=1;}
	 else {std::cout<<"1-position for a need to be an int"<<std::endl;return 1;}
      }
    //
    Aa = circulant_mat(a); 

    

    //construct polynomial b
    //To change: b(x) = ...
    b.zeros();
     for (int i=0;i<b_wt;i++)
      {
	int temp;
	std::istringstream argvi( argv[i+6+a_wt] );
	if (argvi>>temp){b(temp)=1;}
	 else {std::cout<<"b(i) need to be binary"<<std::endl;return 1;}
      }
    //
    Bb = circulant_mat(b);

    //construct Hx
    for(int i = 0; i < l; i++){
        vtemp = Aa.get_col(i);
        Hx.set_col(i, vtemp);
    }
    for(int j = l; j < 2*l; j++){
        vtemp = Bb.get_col(j-l);
        Hx.set_col(j, vtemp);
    }

    //construct Hz
    TAa = Aa.transpose();
    TBb = Bb.transpose();
    for(int p = 0; p < l; p++){
        vtemp = TBb.get_col(p);
        Hz.set_col(p, vtemp);
    }
    for(int q = l; q < 2*l; q++){
        vtemp = TAa.get_col(q-l);
        Hz.set_col(q, vtemp);
    }
    ofstream Hx_file;
    Hx_file.open (file_name1,ios::trunc);

    ofstream Hz_file;
    Hz_file.open (file_name2,ios::trunc);
  
    //std::cout << "for xingrui:" << "\n";
    // std::cout << "Hx" << "\n";
    Hx_file<< 2*l << " " << l << "\n";

    for(int onei = 0; onei < l; onei++){
        for(int onej = 0; onej < 2*l; onej++){
            if(Hx(onei,onej) == 1){
	      Hx_file<< onej+1 << " ";
            }
        }
        Hx_file<<"\n"; 
    }
   
    //std::cout << "Hz" << "\n";
    Hz_file << 2*l << " " << l << "\n";

    for(int onek = 0; onek < l; onek++){
        for(int onel = 0; onel < 2*l; onel++){
            if(Hz(onek,onel) == 1){
                Hz_file<< onel+1 << " ";
            }
        }
        Hz_file<< "\n"; 
    }

    std::chrono::high_resolution_clock::time_point end = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double, std::milli> time_span = (end - start)/1000;
    std::cout << "Run-time " << time_span.count() << " seconds.\n";
    // std::cout << "Run-time " << time_span.count()/60 << " minutes.\n";
    // std::cout << "Run-time " << time_span.count()/3600 << " hours.\n";
    Hx_file.close();
    Hz_file.close();
    return 0;
}
