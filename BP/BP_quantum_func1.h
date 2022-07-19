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
void error_channel2(GF2mat &error, int wt);
void depolarizing(GF2mat &xerror,GF2mat &zerror, const vec &p);
int weight(GF2mat &cw);
int s_weight(GF2mat &s);
//get a error probability vector whose components are random from 0.5p t0 1.5p
void pro_dist(double pmin,double pmax, vec& pv);

//a quantum sequential update schedule: 
void quan_p_update(const nodes checks[],const nodes errors[],mat &mcv,mat& mvc,const GF2mat& syndrome,const vec &pv,int c, int v,  GF2mat& output_e, vec &LR,double alpha=1);

void quan_s_update(const nodes checks[],const nodes errors[],mat &mcv,mat& mvc,const GF2mat& syndrome,const vec &pv,int c, int v,  GF2mat& output_e, vec &LR,double alpha=1);
//decode an all-X or an all-Z error for a CSS code
bool  quan_decode(GF2mat &H, GF2mat &H2,const nodes checks[],const nodes errors[],const vec &pv,const vec&pv_dec,double pavg,double range,double& num_iter, int lmax,int wt,int &max_fail, int &syn_fail, int debug, vec &LR,int &OSD_sec,double alpha=1,double lambda=1);
void OSD(vec& LR,const GF2mat& H,const GF2mat& syndrome, GF2mat &output_e,const GF2mat& G,const GF2mat &real_e);
GF2mat get_gen(const GF2mat &H);
int GF2mat_rank(const GF2mat& H);
GF2mat col_permutation_matrix(ivec & perm);



void err_pos(const nodes errors[],const GF2mat &error);
void err_pos2(const GF2mat &error);
