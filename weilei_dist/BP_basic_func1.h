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
void initialize_checks (const bmat &H, nodes checks[], int & E);

//find the neighbors of variable(error) nodes
void initialize_errors(const bmat &H, nodes errors[]);

//initialize massages to all-one matrices
void initialize_massages(mat &mcv,mat& mvc,const bmat &H);

//do a  parallel update 
void p_update(const nodes checks[],const nodes errors[],mat &mcv,mat& mvc,const bmat& syndrome,double p,int c, int v,  bmat& output_e);

//do a  sequential update 
void s_update(const nodes checks[],const nodes errors[],mat &mcv,mat& mvc,const bmat& syndrome,double p,int c, int v,  bmat& output_e);

//for printing some information for analysis
void p_update_a(const nodes checks[],const nodes errors[],mat &mcv,mat& mvc,const bmat& syndrome,double p,int c, int v,  bmat& output_e);

void real_p_update(const nodes checks[],const nodes errors[],mat &mcv,mat& mvc,mat &pre_mcv,mat& pre_mvc,const bmat& syndrome,double p,int c, int v,  bmat& output_e);

void s_update_a(const nodes checks[],const nodes errors[],mat &mcv,mat& mvc,const bmat& syndrome,double p,int c, int v,  bmat& output_e);

//the distance between 2 cws
int distance(const bmat &output_e, const bmat &real_e,int n);

//print the positions of errors
void print_error_pos( const bmat &real_e,int n);

void update_ci_to_vj(const nodes checks[],const nodes errors[],mat &mcv,mat& mvc,int i,int j,bin s);

void update_vj_to_ci(const nodes checks[],const nodes errors[],mat &mcv,mat& mvc,int j,int i,double ipr);

bmat read_matrix (int& n,int &r, string & file_name);

bmat merge_mat_hori(const bmat &left,const bmat &right);

void write_matrix(string file_name, bmat &H);

int  cla_decode(int v,int c,const bmat &H,const nodes checks[],const nodes errors[],BSC &bsc,double& num_iter, int lmax,int & er,double p);
