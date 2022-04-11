#include <iostream>
#include<cstring>
#include<string>
#include<vector>
#include <fstream>
#include <math.h> 
#include <random>
#include <sstream>

using namespace std;
#include <itpp/itbase.h>
#include <itpp/itcomm.h>
using namespace itpp;

// a "bsc" channel that different qubit has different error probablity
void error_channel(bvec &cw, const vec &p);

//get a error probability vector whose components are random from 0.5p t0 1.5p
void pro_dist(double p, vec& pv);

//a quantum sequential update schedule: 
void quan_s_update(const nodes checks[],const nodes errors[],mat &mcv,mat& mvc,const bmat& syndrome,const vec &pv,int c, int v,  bmat& output_e);

//decode an all-X or an all-Z error for a CSS code
bool  quan_decode(int v,int c,const bmat &H,const bmat &H2,const nodes checks[],const nodes errors[],const vec &pv,double& num_iter, int lmax);

  //check if real_e is a stabilizer
bool Q_inspan(bvec &real_eT,bmat &H2,int r);

//check if real_e-output_e is a stabilizer
bool quan_check(bmat &output_e,bmat &real_e,bmat &H2,int r,int n);
