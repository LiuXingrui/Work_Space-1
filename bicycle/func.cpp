#include<iostream>
#include<cmath>
#include<itpp/itbase.h>
#include<fstream>


// write a number into binary form
itpp::bvec binary_num(int r1, int x1){

    itpp::bvec co(r1);
    int j1 = r1 - 1;

    co.zeros();

    while (x1 > 0) { 
            co(j1) = x1 % 2; 
            x1 = x1 / 2; 
            j1--; 
    }
    
    return co;

}

int weight(itpp::bvec a2){
    int n2 = a2.size();
    int w2 = 0;
    for (int i2 = 0; i2 < n2; i2++) {
        if (a2(i2) == 1){
            w2++;
        }
    }

    return w2;
}

itpp::GF2mat circulant_mat(itpp::bvec a3){
    int l3 = a3.size();
    itpp::bvec a3_temp(l3),a3_temp2(l3);
    itpp::GF2mat Aa3(l3,l3);
    for(int j3 = 0; j3 < l3; j3++){
        a3_temp(j3) = a3(l3-j3-1);    
    }
    Aa3.set_row(l3-1,a3_temp);

    for(int k3 = l3-2; k3 >= 0; k3--){
        for(int m3 = 0; m3 < l3; m3++){
            a3_temp2(m3) = a3_temp((m3+1)%l3);
        }
        a3_temp = a3_temp2;
        Aa3.set_row(k3, a3_temp);
    }

    return Aa3;
}
