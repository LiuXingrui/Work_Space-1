#include <iostream>
#include<cstring>
#include<string>
#include<vector>
#include <fstream>
#include <math.h> 
#include <random>

using namespace std;
#include <itpp/itbase.h>
#include <itpp/itcomm.h>
using namespace itpp;



class  nodes{
public:
  ivec neighbors;
  int degree;


};


//find the neighbors of check nodes  
void initialize_checks (const GF2mat &H, nodes checks[], int & E);

//find the neighbors of variable(error) nodes
void initialize_errors(const GF2mat &H, nodes errors[]);

//initialize massages to all-one matrices
void initialize_massages(mat &mcv,mat& mvc,const GF2mat &H);

//do a  parallel update 
void p_update(const nodes checks[],const nodes errors[],mat &mcv,mat& mvc,const GF2mat& syndrome,vec &pv,int c, int v,  GF2mat& output_e);

//do a  sequential update 
void s_update(const nodes checks[],const nodes errors[],mat &mcv,mat& mvc,const GF2mat& syndrome,double p,int c, int v,  GF2mat& output_e);

//for printing some information for analysis
void p_update_a(const nodes checks[],const nodes errors[],mat &mcv,mat& mvc,const GF2mat& syndrome,double p,int c, int v,  GF2mat& output_e);

void real_p_update(const nodes checks[],const nodes errors[],mat &mcv,mat& mvc,mat &pre_mcv,mat& pre_mvc,const GF2mat& syndrome,double p,int c, int v,  GF2mat& output_e);

void s_update_a(const nodes checks[],const nodes errors[],mat &mcv,mat& mvc,const GF2mat& syndrome,double p,int c, int v,  GF2mat& output_e);

//the distance between 2 cws
int distance(const GF2mat &output_e, const GF2mat &real_e,int n);

//print the positions of errors
void print_error_pos( const GF2mat &real_e,int n);

void update_ci_to_vj(const nodes checks[],const nodes errors[],mat &mcv,mat& mvc,int i,int j,bin s);

void update_vj_to_ci(const nodes checks[],const nodes errors[],mat &mcv,mat& mvc,int j,int i,double ipr,double alpha=1);

GF2mat read_matrix (int& n,int &r, string & file_name);
bmat read_matrix2 (int& n,int &r, string & file_name);

GF2mat merge_mat_hori(const GF2mat &left,const GF2mat &right);
bmat merge_mat_hori2(const bmat &left,const bmat &right);

void write_matrix(string file_name, GF2mat &H);
void write_matrix2(string file_name, bmat &H);

int  cla_decode(int v,int c,const GF2mat &H,const nodes checks[],const nodes errors[],double& num_iter, int lmax,int & er,vec & pv);
