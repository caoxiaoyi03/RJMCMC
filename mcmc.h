#ifndef MCMC_H
#define MCMC_H
#include"Matrix.h"
#include"ClusterFlags.h"
#include <iostream>
#include <vector>
#include "MCMCEnv.h"


class mcmc
{
private:

  typedef std::vector<std::pair<unsigned long, unsigned long long> > SplitMergeSet;

  //// the number of time points
  //long n_timepoints; 
  //// G
  //long n_genes;
  //// number of data
  //long num;
  //// Data number in each time point
  //std::vector <long> NT;
  // C
/*  std::vector<long> clusternum;
 */ // w: row relates to time , column relates to cluster.
  //matrix weight;
  //// tao
  //std::vector<bool> tao;
  //// z: row relates to time , column relates to cluster.
  //ClusterFlags z;
  //// Differentiation tree
  //matrixl tree;
  //// Report home of each cluster vector
  //matrixl clusterhome;
  //// X: Data
  //matrix Data; 
public:
	MCMCEnv env;			// may use pointer for further speeding up, need to check
	//mcmc(long n_timepointsdata, long n_genesdata, long numdata, const std::vector<long> & NTdata);
	mcmc(const MCMCEnv &env);
	~mcmc();
	//long gettime() const;
	//long getgene() const;
	//long getnum() const;
	//std::vector<long> getNT() const;
	//std::vector<long> getclusternum() const;
	//matrix getweight() const;
	//std::vector<bool> gettao() const;
	//long getfgenenum() const;
	//ClusterFlags getz() const;
	//matrixl gettree() const;
	//matrixl getclusterhome()const;
	//matrix getData() const;
	//void zw_reflag(ClusterFlags &x_z, matrix &x_weight);

	////matrix getData();
	//double Kc_weight(long flag)const;
	//matrix reData(long flag, bool tao_NULL) const;
	//std::vector<long> indextoflag(long x)const;
	//matrix datacov(long flag, bool tao_NULL, const matrix &mu1, const matrix &mu0)const;
	//
	//double density(const matrix &mu1, const matrix &mu0, const matrix &sigma0, const matrix &sigma1, const double beta0, 
	//	const double beta1, const double deltar) const;
	//double tao_density(const std::vector<bool> &newtao, const matrix &mu1, const matrix &mu0, const matrix &sigma0, 
	//	const matrix &sigma1, const double beta0, const double beta1,const double deltar) const;

	//matrix tao_reData(const std::vector<bool> &newtao, const long flag, const bool tao_NULL)const;
	//matrix tao_datacov(const std::vector<bool> &newtao, const long flag, const bool tao_NULL, const matrix &mu1, 
	//	const matrix &mu0)const;

	//double z_Kc_weight(const ClusterFlags &newz, const long flag)const;
	//matrix z_reData(const ClusterFlags &newz, const long flag, const bool tao_NULL)const;
	//matrix z_datacov(const ClusterFlags &newz, const long flag, const bool tao_NULL, const matrix &mu1, 
	//	const matrix &mu0)const;
	//double z_density(const ClusterFlags &newz, const matrix &mu1, const matrix &mu0, const matrix &sigma0, 
	//	const matrix &sigma1, const double beta0, const double beta1,const double deltar)const;

	//double zw_Kc_weight(const ClusterFlags &newz, const matrix &newweight, const long flag)const;
	//matrix zw_reData(const ClusterFlags &newz, const long flag, const bool tao_NULL)const;
	//matrix zw_datacov(const ClusterFlags &newz, const long flag, const bool tao_NULL, const matrix &mu1, 
	//	const matrix &mu0)const;
	//double zw_density(const ClusterFlags &newz, const matrix &newweight, const matrix &mu1, const matrix &mu0, 
	//const matrix &sigma0, const matrix &sigma1, const double beta0, const double beta1,const double deltar)const;

	//void tao_sample(const matrix &mu1, const matrix &mu0, const matrix &sigma0, const matrix &sigma1, const double beta0, 
	//	const double beta1, const double deltar);
	//void weight_sample(const double alpha, const double shape, const matrix &mu1, const matrix &mu0, const matrix &sigma0,
	//	const matrix &sigma1, const double beta0, const double beta1, const double deltar);
	//void z_sample(const matrix &mu1, const matrix &mu0, const matrix & sigma0, const matrix &sigma1, const double beta0, 
	//	const double beta1, const double deltar);

	////void flag_split(ClusterFlags &one, MCMCEnv &env, double &P_alloc, bool &NULLset, double &f_ui,  
	////	long time, long clusterindex)const;

	////void flag_merge(ClusterFlags &one, matrix &two, matrixl &newtree, matrixl &newchome,
	////    double &P_alloc, double &f_ui, long time, long clustermother, 
	////	long clusterson)const;

	////void flag_birth(ClusterFlags &one, matrix &two, matrixl &newtree, matrixl &newchome,  double &weightratio,
	////	double &jacobi, double &f_wstar, const long time, const long clusterindex, const double alpha)const;
	////
	////void flag_death(ClusterFlags &one, matrix &two, matrixl &newtree, matrixl &newchome,
	////double &jacobi, double &weightratio, double &f_wstar, long time, long clusterindex1, long clusterindex2, 
	////const double alpha)const;


	//void tail_birth(const long time, const long clusterindex,
	//ClusterFlags &newz, matrix &newweight, const double alpha)const;

	//void tail_death(const long time, const long clusterindex1, const long clusterindex2,
	//    ClusterFlags &newz, matrix &newweight, const double alpha)const;


    //long EmptyCluster(const ClusterFlags &one, const matrix &two, const long clusterindex)const;
	//long f_tdatanum(const long time)const;
	//long f_tidatanum(const long time, const long loc)const;
	//matrix f_clustermean(const long time, const long clusterindex, const ClusterFlags &one, const matrix &two)const;
	//double f_clusterdist(const long time1, const long clusterindex1, const long time2, const long clusterindex2,
	//const ClusterFlags &one, const matrix &two)const;
	//long f_similarity(const long time, const long clustermother, const long clusterson, const bool NULLset,
	//const matrixl splitsets,
	//const ClusterFlags &one, const matrix &two) const;
	//double f_pchos(const long time, const long clustermother, const long clusterson, const bool NULLset,
	//const matrixl splitsets,
	//const ClusterFlags &one, const matrix &two)const;


	//matrixl f_splitset(const ClusterFlags &newz, const matrix &newweight, const matrixl &newtree,
	//	const matrixl &newchome)const;
	//SplitMergeSet f_splitset(const MCMCEnv &env) const;
	//SplitMergeSet f_mergeset(const ClusterFlags &newz, const matrix &newweight, const matrixl &newtree)const;
	//SplitMergeSet f_deathset(const ClusterFlags &newz, const matrix &newweight, const matrixl &newtree)const;
	//SplitMergeSet f_taildeathset(const ClusterFlags &newz, const matrix &newweight, const matrixl &newtree)const;
	//SplitMergeSet f_tailbirthset(const ClusterFlags &newz, const matrix &newweight, const matrixl &newtree)const;

	//double AP_split(const long time, const long clusterindex, const bool NULLset, const matrixl &splitsets, 
	//	const long splitnum, const double P_alloc, const double f_ui,
	//	const ClusterFlags &newz, 
	//	const matrix &newweight, const double bk,
	//	const matrix &mu1, const matrix &mu0, const matrix &sigma0, const matrix &sigma1, const double beta0, 
	//	const double beta1, const double deltar, const double alpha)const;
	//double AP_merge(const long time, const long clustermother, const long clusterson, const matrixl &splitsetsaftermerge,
	//	 const long mergesetsrows,
	//     const long splitnumaftermerge, const double P_alloc, const double f_ui,
	//     const ClusterFlags &newz, const matrix &newweight, const double bk,
	//     const matrix &mu1, const matrix &mu0, const matrix &sigma0, const matrix &sigma1, 
	//     const double beta0, const double beta1, const double deltar, const double alpha)const;
	//double AP_birth(const long time, const long clusterindex, const long splitnumafterbirth, 
	//	 const double weightratio, const double jacobi, const double f_wstar, 
	//     const ClusterFlags &newz, const matrix &newweight, const double bk,
	//     const matrix &mu1, const matrix &mu0,
	//     const matrix &sigma0, const matrix &sigma1, const double beta0, const double beta1, const double deltar,
	//     const double alpha)const;
	//double AP_death(const long time, const long clustermother, const long clusterson,
	//     const long emptynumbeforedeath, const double f_wstar, const double weightratio, const double jacobi,
	//     const ClusterFlags &newz, const matrix &newweight, const double bk, 
	//     const matrix &mu1, const matrix &mu0, const matrix &sigma0, const matrix &sigma1,
	//     const double beta0, const double beta1, const double deltar, const double alpha)const;
	//double AP_tailbirth( 
	//     const ClusterFlags &newz, const matrix &newweight,
	//     const matrix &mu1, const matrix &mu0,
	//     const matrix &sigma0, const matrix &sigma1, const double beta0, const double beta1, const double deltar,
	//     const double alpha)const;

	//double AP_taildeath(
	//	 const ClusterFlags &newz, const matrix &newweight,
	//     const matrix &mu1, const matrix &mu0,
	//     const matrix &sigma0, const matrix &sigma1, const double beta0, const double beta1, const double deltar,
	//     const double alpha)const;


	//void SplitMerge_move(double bk, ClusterFlags &newz, matrix &newweight, matrixl &newtree, matrixl &newchome,
	//     const matrix &mu0, const matrix &mu1, const matrix &sigma0, 
	//     const matrix &sigma1, const double beta0, const double beta1,const double deltar, const double alpha,
	//	 bool &move, double &probability);
	bool SplitMerge_move(double &probability, bool isForest = false);
	bool SplitMerge_move_test(double &probability, bool isForest = false);
	//void BirthDeath_move(double bk, ClusterFlags &newz, matrix &newweight, matrixl &newtree, matrixl &newchome,
	//     const matrix &mu0, const matrix &mu1, const matrix &sigma0, 
	//     const matrix &sigma1, const double beta0, const double beta1,const double deltar, const double alpha,
	//	 bool &move, double &probability);
	bool BirthDeath_move(double &probability, bool isForest = false);
	bool TailBirthDeath_move(double &probability, bool isForest = false);



	//void C_sample(const matrix &mu0, const matrix &mu1, const matrix &sigma0, 
	//	const matrix &sigma1, const double beta0, const double beta1,const double deltar, const double alpha,
	//	const double shape, const double bk, const double ranc);

	//double treesummary()const;
	//double treesummary(const matrixl &newtree, const matrixl &newchome, const matrix &newweight,
	//	const ClusterFlags &newz)const;

    //void parameter_initialization(matrix &mu1, matrix &mu0, matrix &sigma1, matrix &sigma0, 
	   // double &beta1, double &beta0, double &deltar, double &alpha, double &shape, double &bk, double &ranc)const;




	//double AP_split(const long time, const long clusterindex, const bool NULLset, 
	//const double f_ui,
	//const ClusterFlags &newz, const matrix &newweight, const double bk,
	//const matrix &mu1, const matrix &mu0, 
	//const matrix &sigma0, const matrix &sigma1, const double beta0,	const double beta1, const double deltar,
	//const double alpha)const;



};

std::ostream & operator<<(std::ostream &os, const mcmc &mcm);



#endif