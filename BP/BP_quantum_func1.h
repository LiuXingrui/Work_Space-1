#include <iostream>
#include<cstring>
#include<string>
#include<vector>
#include <fstream>
#include <math.h> 
#include <random>
#include <sstream>

#ifndef basic
#define basic
#include"BP_basic_func1.h"
#endif

using namespace std;
#include <itpp/itbase.h>
#include <itpp/itcomm.h>
using namespace itpp;

// a "bsc" channel that different qubit has different error probablity
void error_channel(GF2mat &cw, const vec &p);

//get a error probability vector whose components are random from 0.5p t0 1.5p
void pro_dist(double pmin,double pmax, vec& pv);

//a quantum sequential update schedule: 
void quan_s_update(const nodes checks[],const nodes errors[],mat &mcv,mat& mvc,const GF2mat& syndrome,const vec &pv,int c, int v,  GF2mat& output_e);

//decode an all-X or an all-Z error for a CSS code
bool  quan_decode(GF2mat &H, GF2mat &H2,const nodes checks[],const nodes errors[],const vec &pv,double& num_iter, int lmax);

GF2mat get_gen(const GF2mat &H);
int GF2mat_rank(const GF2mat& H);
