// Corrections:11/20/2011-11/30/2011
// 1. P_alloc is corrected in 'flag_split' function.
// 2. Kc_weight, z_Kc_weight, zw_Kc_weight are corrected
// 3. ran_num is corrected in order to get correct results from function z_sample.
// 4. K1 in parameter_initialization function should be 1%~10% of data variances to increase data effect.
// 5. beta1 should be set a large value in order to increase data effect.
// 6. mu1 is set to mean value of all data set.


#include "mcmc.h"
#include "ClusterFlags.h"
#include "Random.h"
#include <iostream>
#include <fstream>
#include <string>
#include <cmath>
#include <vector>
#include <stdexcept>

#define PI 3.1415926
#define LARGENUM 1e+300
#define SMALLNUM 1e-200

using namespace std;


//long n_timepoints; 
// // G
// long n_genes;
// // number of data
// long num;
// // Data number in each time point
// std::vector <long> NT;
// // C
// std::vector<long> clusternum(n_timepoints);
// // w
// matrix weight(n_timepoints, n_timepoints);
// // tao
// std::vector<bool> tao(n_genes);
// // z
// ClusterFlags z(n_timepoints, NT);
// // X: Data
// matrix Data(num, n_genes); 

//mcmc::mcmc(long n_timepointsdata, long n_genesdata, long numdata, const std::vector<long> & NTdata)
//	:n_timepoints(n_timepointsdata), n_genes(n_genesdata), num(numdata),NT(NTdata),clusternum(n_timepointsdata), 
//	weight(n_timepointsdata, n_timepointsdata), tao(n_genesdata), z(n_timepointsdata, NTdata), Data(numdata, n_genesdata){
//		/*for(int i=0; i<NTdata.size(); i++)
//		NT.at(i) = NTdata.at(i);*/
//}
//

mcmc::mcmc(const MCMCEnv &p_env): env(p_env) {
};


mcmc::~mcmc()
{

}

//
//long mcmc::gettime() const{
//	return n_timepoints;
//}
//long mcmc::getgene() const{
//	return n_genes;
//}
//long mcmc::getnum() const{
//	return num;
//}
//std::vector<long> mcmc::getNT() const{
//	return NT;
//}
//std::vector<long> mcmc::getclusternum() const{
//	return clusternum;
//}
//matrix mcmc::getweight() const{
//	return weight;
//}
//std::vector<bool> mcmc::gettao() const{
//	return tao;
//}
//
//long mcmc::getfgenenum() const{
//	int i,z = 0;
//	for(i=0; i<tao.size(); i++)
//		z += tao.at(i);
//	return z;
//}
//
//
//ClusterFlags mcmc::getz() const{
//	return z;
//}
//
//matrixl mcmc::gettree() const{
//	return tree;
//}
//
//
//matrixl mcmc::getclusterhome()const{
//	return clusterhome;
//}
//
//
//matrix mcmc::getData() const{
//	return Data;
//}

//// This function is used to reorder flag z and the relative weight
//// Note that if z is reorderred, weight must also be changed.
//void mcmc::zw_reflag(ClusterFlags &x_z, matrix &x_weight){
//
//
//
//}


//// Note that, if there is no flag appearing in the cluster list z, there will be a 'out of range' bug appearing.
//double mcmc::Kc_weight(long flag)const{
//	long kc = weight.getEarliestTime(flag);
//	long i;
//	double w=1.0;
//	// Note that i is not the real time point, its begginning point should be kc-1.
//	for(i=kc-1; i<n_timepoints; i++)
//		w *= pow(weight.getValue(i, flag-1),z.tflagnum(flag, i+1));
//
//
//	return w;
//}
//



//// In this function, cluster inicator z is assumbed to be reflaged using 'reflag' function.
//// tao_NULL = 1, delete tao=0, otherwise delete tao=1.
//matrix mcmc::reData(long flag, bool tao_NULL)const{
//	matrix one(Data);
//	std::vector<std::vector<long> > mclust = z.onezerocluster(flag);
//	long i,j,k;
//	long length = 0;
//	for(i=0; i<mclust.size(); i++){
//		for(j=0; j<mclust.at(i).size(); j++)
//			for(k=0; k<one.cols(); k++)
//				one.setValue(length+j, k, Data.getValue(length+j,k)*(tao_NULL? tao.at(k):(1-tao.at(k)))*(tao_NULL?mclust.at(i).at(j):1));
//
//		length += j;
//	}
//
//	return one;
//
//}
//
//matrix mcmc::datacov(long flag, bool tao_NULL, const matrix &mu1, const matrix &mu0)const{
//	matrix one = reData(flag, tao_NULL);
//	matrix two=(tao_NULL==true)?mu1:mu0;
//	
//	long x = (tao_NULL==true)?z.flagnum(flag):num;
//
//	matrix mid;
//	matrix one_mean = (one.datasum()*(1.0/(double)x));
//	std::vector<long> xx;
//
//	matrix cov_NULL(Data.cols(),Data.cols());
//	int i,j;
//
//	for(i=0; i<Data.rows(); i++)
//	{
//		xx = indextoflag(i);
//		if((z.getFlag(xx.at(0),xx.at(1))==flag)||(!tao_NULL)){
//			mid = (one.getRow(i)-one_mean);
//		} else{
//			mid = one.getRow(i);
//		}
//
//		//std::cout<<"Matrix transpose: "<<std::endl;
//		//std::cout<<mid<<std::endl;
//		cov_NULL += (mid.transpose()*mid);
//	}
//
//	mid = two - one_mean;
//	/*std::cout<<"mu0: "<<std::endl;
//	std::cout<<two<<std::endl;
//	std::cout<<"mean: "<<std::endl;
//	std::cout<<one_mean<<std::endl;*/
//	two = (mid.transpose()*mid);
//
//
//
//	return cov_NULL.rbind(two);
//}
//

//std::vector<long> mcmc::indextoflag(long x)const{
//	long i=0;
//	long length = NT.at(i);
//	std::vector<long> one(2);
//	if(x<0) {
//		throw(std::out_of_range("Index must be positive!!!"));
//	}
//	while(x>=length){
//		i++;
//		length += NT.at(i);
//	}
//	one.at(0) = i;
//	one.at(1) = x-length+NT.at(i);
//
//	return one;
//}


//
//
//double mcmc::density(const matrix &mu1, const matrix &mu0, const matrix &sigma0, const matrix &sigma1, const double beta0, 
//	const double beta1, const double deltar)const{
//
//		double den = 1.0;
//		long flag;
//		long i,j;
//		long clusternum = z.getMaxFlag();
//		long cclusternum;
//		double gammaprod = 1.0;
//		double K0=1.0, K1=1.0;
//		long dim = n_genes;
//		matrix S0, S1, mid;
//
//		long taonum=0;
//		for(i=0; i<dim; i++)
//			taonum += tao.at(i);
//		/*std::cout<<"taonum"<<std::endl;
//		std::cout<<taonum<<std::endl;*/
//		if((taonum>0)&&(taonum<n_genes)){
//			for(j=1; j<=(n_genes-taonum); j++)
//				gammaprod *= gammaf(0.5*(num+deltar+n_genes-taonum-j))/gammaf(0.5*(deltar+n_genes-taonum-j));
//
//			K0 = pow(beta0*num+1,-0.5*(n_genes-taonum))*gammaprod;
//			/*std::cout<<"K0"<<std::endl;
//			std::cout<<K0<<std::endl;*/
//			gammaprod = 1.0;
//
//			mid = datacov(1, 0, mu1, mu0);
//			/* std::cout<<"mid"<<std::endl;
//			std::cout<<mid<<std::endl;*/
//			S0 = mid.getBlock(1,dim,1,dim)
//				+ mid.getBlock(dim+1,2*dim,1,dim)*(num/(beta0*num+1));
//			/* std::cout<<"S0:"<<std::endl;
//			std::cout<<S0<<std::endl;*/
//			den = K0 * pow(fabs(sigma0.determinant()),0.5*(deltar+n_genes-taonum-1))
//				*pow(fabs((sigma0+S0).determinant()),-0.5*(num+deltar+n_genes-taonum-1));
//
//			for(flag=1; flag<=clusternum; flag++){
//				cclusternum = z.flagnum(flag);
//				if(cclusternum !=0 ){
//					for(j=1; j<=taonum; j++)
//						gammaprod *= gammaf(0.5*(cclusternum+deltar+taonum-j))/gammaf(0.5*(deltar+taonum-j));
//
//					K1 = Kc_weight(flag)*pow(beta1*z.flagnum(flag)+1,-(double)taonum/2.0)*gammaprod;
//					//std::cout<<"K1"<<std::endl;
//					//std::cout<<K1<<std::endl;
//					gammaprod = 1.0;
//
//
//
//					mid = datacov(flag, 1, mu1, mu0);
//					//std::cout<<"S1mid:"<<std::endl;
//					//std::cout<<mid<<std::endl;
//					S1 = mid.getBlock(1,dim,1,dim) 
//						+mid.getBlock(dim+1,2*dim,1,dim)*(cclusternum/(beta1*cclusternum+1));
//					/* std::cout<<"S1:"<<std::endl;
//					std::cout<<S1<<std::endl;*/
//
//					den *= K1*pow(fabs(sigma1.determinant()),0.5*(deltar+taonum-1))
//						*pow(fabs((sigma1+S1).determinant()),-0.5*(cclusternum+deltar+taonum-1));
//				}
//
//			}
//		} else if(taonum==0){
//			for(j=1; j<=(n_genes-taonum); j++)
//				gammaprod *= gammaf(0.5*(num+deltar+n_genes-taonum-j))/gammaf(0.5*(deltar+n_genes-taonum-j));
//
//			K0 = pow(beta0*num+1,-0.5*(n_genes-taonum))*gammaprod;
//			/*  std::cout<<"K0"<<std::endl;
//			std::cout<<K0<<std::endl;*/
//			gammaprod = 1.0;
//
//			mid = datacov(1, 0, mu1, mu0);
//			/*   std::cout<<"mid"<<std::endl;
//			std::cout<<mid<<std::endl;*/
//			S0 = mid.getBlock(1,dim,1,dim)
//				+ mid.getBlock(dim+1,2*dim,1,dim)*(num/(beta0*num+1));
//			/* std::cout<<"S0:"<<std::endl;
//			std::cout<<S0<<std::endl;*/
//			den = K0 * pow(fabs(sigma0.determinant()),0.5*(deltar+n_genes-taonum-1))
//				*pow(fabs((sigma0+S0).determinant()),-0.5*(num+deltar+n_genes-taonum-1));
//
//
//		} else{
//
//
//			for(flag=1; flag<=clusternum; flag++){
//				cclusternum = z.flagnum(flag);
//				if(cclusternum != 0){
//					for(j=1; j<=taonum; j++)
//						gammaprod *= gammaf(0.5*(cclusternum+deltar+taonum-j))/gammaf(0.5*(deltar+taonum-j));
//
//					/*std::cout<<"gammaprod"<<std::endl;
//					std::cout<<gammaprod<<std::endl;*/
//
//					K1 = Kc_weight(flag)*pow(beta1*z.flagnum(flag)+1,-(double)taonum/2.0)*gammaprod;
//					
//					
//					gammaprod = 1.0;
//
//
//
//					mid = datacov(flag, 1, mu1, mu0);
//					/*std::cout<<"S1mid:"<<std::endl;
//					std::cout<<mid<<std::endl;*/
//					
//					S1 = mid.getBlock(1,dim,1,dim) 
//						+mid.getBlock(dim+1,2*dim,1,dim)*(cclusternum/(beta1*cclusternum+1));
//					/*std::cout<<"density_K1"<<std::endl;
//					std::cout<<K1<<std::endl;
//					std::cout<<"density_Q0"<<std::endl;
//					std::cout<<pow(fabs(sigma1.determinant()),0.5*(deltar+taonum-1))<<std::endl;
//					std::cout<<"density_Q0+S1"<<std::endl;
//					std::cout<<pow(fabs((sigma1+S1).determinant()),-0.5*(cclusternum+deltar+taonum-1))<<std::endl;
//					*/
//					den *= K1*pow(fabs(sigma1.determinant()),0.5*(deltar+taonum-1))
//						*pow(fabs((sigma1+S1).determinant()),-0.5*(cclusternum+deltar+taonum-1));
//				}
//
//			}
//
//
//		}
//
//		return den;
//}
//
//
//
//
//matrix mcmc::tao_reData(const std::vector<bool> &newtao, const long flag, const bool tao_NULL)const{
//	matrix one(Data);
//	std::vector<std::vector<long> > mclust = z.onezerocluster(flag);
//	long i,j,k;
//	long length = 0;
//	for(i=0; i<mclust.size(); i++){
//		for(j=0; j<mclust.at(i).size(); j++)
//			for(k=0; k<one.cols(); k++)
//				one.setValue(length+j, k, Data.getValue(length+j,k)*(tao_NULL? newtao.at(k):(1-newtao.at(k)))*(tao_NULL?mclust.at(i).at(j):1));
//
//		length += j;
//	}
//
//	return one;
//}
//
//
//
//matrix mcmc::tao_datacov(const std::vector<bool> &newtao, const long flag, const bool tao_NULL, const matrix &mu1, 
//	const matrix &mu0)const{
//
//		matrix one = tao_reData(newtao, flag, tao_NULL);
//		matrix two=(tao_NULL==true)?mu1:mu0;
//		long x = (tao_NULL==true)?z.flagnum(flag):num;
//		matrix mid;
//		matrix one_mean = (one.datasum()*(1.0/(double)x));
//		std::vector<long> xx;
//
//		matrix cov_NULL(Data.cols(),Data.cols());
//		int i,j;
//		for(i=0; i<Data.rows(); i++)
//		{
//			xx = indextoflag(i);
//			if((z.getFlag(xx.at(0),xx.at(1))==flag)||(!tao_NULL)){
//				mid = (one.getRow(i)-one_mean);
//			} else{
//				mid = one.getRow(i);
//			}
//
//			//std::cout<<"Matrix transpose: "<<std::endl;
//			//std::cout<<mid<<std::endl;
//			cov_NULL += (mid.transpose()*mid);
//		}
//
//		mid = two - one_mean;
//		two = (mid.transpose()*mid);
//		return cov_NULL.rbind(two);
//
//
//}
//
//
//double mcmc::tao_density(const std::vector<bool> &newtao, const matrix &mu1, const matrix &mu0, const matrix &sigma0, 
//	const matrix &sigma1, const double beta0, const double beta1,const double deltar)const{
//		double den = 1.0;
//		long flag;
//		long i,j;
//		long clusternum = z.getMaxFlag();
//		long cclusternum;
//		double gammaprod = 1.0;
//		double K0, K1;
//		long dim = n_genes;
//		matrix S0, S1, mid;
//
//		long taonum=0;
//		for(i=0; i<dim; i++)
//			taonum += newtao.at(i);
//		if(newtao.size()==tao.size()){
//			if((taonum>0)&&(taonum<n_genes)){
//				for(j=1; j<=(n_genes-taonum); j++)
//					gammaprod *= gammaf(0.5*(num+deltar+n_genes-taonum-j))/gammaf(0.5*(deltar+n_genes-taonum-j));
//
//				K0 = pow(beta0*num+1,-0.5*(n_genes-taonum))*gammaprod;
//				/*   std::cout<<"K0"<<std::endl;
//				std::cout<<K0<<std::endl;*/
//				gammaprod = 1.0;
//
//				mid = tao_datacov(newtao, 1, 0, mu1, mu0);
//				//std::cout<<"mid"<<std::endl;
//				//std::cout<<mid<<std::endl;
//				S0 = mid.getBlock(1,dim,1,dim)
//					+ mid.getBlock(dim+1,2*dim,1,dim)*(num/(beta0*num+1));
//				//std::cout<<"S0:"<<std::endl;
//				//std::cout<<S0<<std::endl;
//				den = K0 * pow(fabs(sigma0.determinant()),0.5*(deltar+n_genes-taonum-1))
//					*pow(fabs((sigma0+S0).determinant()),-0.5*(num+deltar+n_genes-taonum-1));
//
//				for(flag=1; flag<=clusternum; flag++){
//					cclusternum = z.flagnum(flag);
//					if(cclusternum != 0){
//						for(j=1; j<=taonum; j++)
//							gammaprod *= gammaf(0.5*(cclusternum+deltar+taonum-j))/gammaf(0.5*(deltar+taonum-j));
//
//						K1 = Kc_weight(flag)*pow(beta1*z.flagnum(flag)+1,-(double)taonum/2.0)*gammaprod;
//						//std::cout<<"K1"<<std::endl;
//						//std::cout<<K1<<std::endl;
//						gammaprod = 1.0;
//
//
//
//						mid = tao_datacov(newtao, flag, 1, mu1, mu0);
//						//std::cout<<"S1mid:"<<std::endl;
//						//std::cout<<mid<<std::endl;
//						S1 = mid.getBlock(1,dim,1,dim) 
//							+mid.getBlock(dim+1,2*dim,1,dim)*(cclusternum/(beta1*cclusternum+1));
//						//std::cout<<"S1:"<<std::endl;
//						//std::cout<<S1<<std::endl;
//
//						den *= K1*pow(fabs(sigma1.determinant()),0.5*(deltar+taonum-1))
//							*pow(fabs((sigma1+S1).determinant()),-0.5*(cclusternum+deltar+taonum-1));
//					}
//
//				}
//			} else if(taonum==0){
//				for(j=1; j<=(n_genes-taonum); j++)
//					gammaprod *= gammaf(0.5*(num+deltar+n_genes-taonum-j))/gammaf(0.5*(deltar+n_genes-taonum-j));
//
//				K0 = pow(beta0*num+1,-0.5*(n_genes-taonum))*gammaprod;
//				//std::cout<<"K0"<<std::endl;
//				//std::cout<<K0<<std::endl;
//				gammaprod = 1.0;
//
//				mid = tao_datacov(newtao, 1, 0, mu1, mu0);
//				//std::cout<<"mid"<<std::endl;
//				//std::cout<<mid<<std::endl;
//				S0 = mid.getBlock(1,dim,1,dim)
//					+ mid.getBlock(dim+1,2*dim,1,dim)*(num/(beta0*num+1));
//				//std::cout<<"S0:"<<std::endl;
//				//std::cout<<S0<<std::endl;
//				den = K0 * pow(fabs(sigma0.determinant()),0.5*(deltar+n_genes-taonum-1))
//					*pow(fabs((sigma0+S0).determinant()),-0.5*(num+deltar+n_genes-taonum-1));
//
//
//			} else{
//
//
//				for(flag=1; flag<=clusternum; flag++){
//					cclusternum = z.flagnum(flag);
//					if(cclusternum != 0){
//						for(j=1; j<=taonum; j++)
//							gammaprod *= gammaf(0.5*(cclusternum+deltar+taonum-j))/gammaf(0.5*(deltar+taonum-j));
//
//						K1 = Kc_weight(flag)*pow(beta1*z.flagnum(flag)+1,-(double)taonum/2.0)*gammaprod;
//						//std::cout<<"K1"<<std::endl;
//						//std::cout<<K1<<std::endl;
//						gammaprod = 1.0;
//
//
//
//						mid = tao_datacov(newtao, flag, 1, mu1, mu0);
//						//std::cout<<"S1mid:"<<std::endl;
//						//std::cout<<mid<<std::endl;
//						S1 = mid.getBlock(1,dim,1,dim) 
//							+mid.getBlock(dim+1,2*dim,1,dim)*(cclusternum/(beta1*cclusternum+1));
//						//std::cout<<"S1:"<<std::endl;
//						//std::cout<<S1<<std::endl;
//
//						den *= K1*pow(fabs(sigma1.determinant()),0.5*(deltar+taonum-1))
//							*pow(fabs((sigma1+S1).determinant()),-0.5*(cclusternum+deltar+taonum-1));
//
//					}
//				}
//
//
//			}						   
//
//
//
//			return den;
//		} else{
//			throw(std::out_of_range("Tao length are differnt!!!"));
//			return NULL;
//		}
//
//
//
//}
//
//
//double mcmc::z_Kc_weight(const ClusterFlags &newz, const long flag)const{
//	long kc = weight.getEarliestTime(flag);
//	long i;
//	double w=1.0;
//	for(i=kc-1; i<n_timepoints; i++)
//		w *= pow(weight.getValue(i, flag-1),newz.tflagnum(flag, i+1));
//	return w;
//}
//
//
//double mcmc::zw_Kc_weight(const ClusterFlags &newz, const matrix &newweight, const long flag)const{
//	long kc = newweight.getEarliestTime(flag);
//	long i;
//	double w=1.0;
//	for(i=kc-1; i<n_timepoints; i++){
//		w *= pow(newweight.getValue(i, flag-1),newz.tflagnum(flag, i+1));
//	}
//	return w;
//}
//
//matrix mcmc::z_reData(const ClusterFlags &newz, const long flag, const bool tao_NULL)const{
//	matrix one(Data);
//	std::vector<std::vector<long> > mclust = newz.onezerocluster(flag);
//	long i,j,k;
//	long length = 0;
//	for(i=0; i<mclust.size(); i++){
//		for(j=0; j<mclust.at(i).size(); j++)
//			for(k=0; k<one.cols(); k++)
//				one.setValue(length+j, k, Data.getValue(length+j,k)*(tao_NULL? tao.at(k):(1-tao.at(k)))*(tao_NULL?mclust.at(i).at(j):1));
//
//		length += j;
//	}
//
//	return one;
//
//}
//
//
//matrix mcmc::zw_reData(const ClusterFlags &newz, const long flag, const bool tao_NULL)const{
//	matrix one(Data);
//	std::vector<std::vector<long> > mclust = newz.onezerocluster(flag);
//	long i,j,k;
//	long length = 0;
//	for(i=0; i<mclust.size(); i++){
//		for(j=0; j<mclust.at(i).size(); j++)
//			for(k=0; k<one.cols(); k++)
//				one.setValue(length+j, k, Data.getValue(length+j,k)*(tao_NULL? tao.at(k):(1-tao.at(k)))*(tao_NULL?mclust.at(i).at(j):1));
//
//		length += j;
//	}
//
//	return one;
//
//}
//
//
//matrix mcmc::z_datacov(const ClusterFlags &newz, const long flag, const bool tao_NULL, const matrix &mu1, 
//	const matrix &mu0)const{
//		matrix one = z_reData(newz,flag, tao_NULL);
//		matrix two=(tao_NULL==true)?mu1:mu0;
//		long x = (tao_NULL==true)?newz.flagnum(flag):num;
//		matrix mid;
//		matrix one_mean = (one.datasum()*(1.0/(double)x));
//		std::vector<long> xx;
//
//	
//
//		matrix cov_NULL(Data.cols(),Data.cols());
//		int i,j;
//		for(i=0; i<Data.rows(); i++)
//		{
//			xx = indextoflag(i);
//			if((newz.getFlag(xx.at(0),xx.at(1))==flag)||(!tao_NULL)){
//				mid = (one.getRow(i)-one_mean);
//			} else{
//				mid = one.getRow(i);
//			}
//
//			//std::cout<<"Matrix transpose: "<<std::endl;
//			//std::cout<<mid<<std::endl;
//			cov_NULL += (mid.transpose()*mid);
//		}
//
//		mid = two - one_mean;
//		/*std::cout<<"mu0: "<<std::endl;
//		std::cout<<two<<std::endl;
//		std::cout<<"mean: "<<std::endl;
//		std::cout<<one_mean<<std::endl;*/
//		two = (mid.transpose()*mid);
//		return cov_NULL.rbind(two);
//}
//
//
//matrix mcmc::zw_datacov(const ClusterFlags &newz, const long flag, const bool tao_NULL, const matrix &mu1, 
//	const matrix &mu0)const{
//		matrix one = zw_reData(newz,flag, tao_NULL);
//		matrix two=(tao_NULL==true)?mu1:mu0;
//		long x = (tao_NULL==true)?newz.flagnum(flag):num;
//		matrix mid;
//		matrix one_mean = (one.datasum()*(1.0/(double)x));
//		std::vector<long> xx;
//
//		matrix cov_NULL(Data.cols(),Data.cols());
//		int i,j;
//		for(i=0; i<Data.rows(); i++)
//		{
//			xx = indextoflag(i);
//			if((newz.getFlag(xx.at(0),xx.at(1))==flag)||(!tao_NULL)){
//				mid = (one.getRow(i)-one_mean);
//			} else{
//				mid = one.getRow(i);
//			}
//
//			//std::cout<<"Matrix transpose: "<<std::endl;
//			//std::cout<<mid<<std::endl;
//			cov_NULL += (mid.transpose()*mid);
//		}
//
//		mid = two - one_mean;
//		/*std::cout<<"mu0: "<<std::endl;
//		std::cout<<two<<std::endl;
//		std::cout<<"mean: "<<std::endl;
//		std::cout<<one_mean<<std::endl;*/
//		two = (mid.transpose()*mid);
//		return cov_NULL.rbind(two);
//}
//
//
//double mcmc::z_density(const ClusterFlags &newz, const matrix &mu1, const matrix &mu0, const matrix &sigma0, 
//	const matrix &sigma1, const double beta0, const double beta1,const double deltar)const{
//
//		double den = 1.0;
//		long flag;
//		long i,j;
//		long clusternum = newz.getMaxFlag();
//		long cclusternum;
//		double gammaprod = 1.0;
//		double K0, K1;
//		long dim = n_genes;
//		matrix S0, S1, mid;
//
//		long taonum=0;
//		for(i=0; i<dim; i++)
//			taonum += tao.at(i);
//
//
//		if((taonum>0)&&(taonum<n_genes)){
//			for(j=1; j<=(n_genes-taonum); j++)
//				gammaprod *= gammaf(0.5*(num+deltar+n_genes-taonum-j))/gammaf(0.5*(deltar+n_genes-taonum-j));
//
//			K0 = pow(beta0*num+1,-0.5*(n_genes-taonum))*gammaprod;
//
//			gammaprod = 1.0;
//
//			mid = z_datacov(newz, 1, 0, mu1, mu0);
//
//			S0 = mid.getBlock(1,dim,1,dim)
//				+ mid.getBlock(dim+1,2*dim,1,dim)*(num/(beta0*num+1));
//
//			den = K0 * pow(fabs(sigma0.determinant()),0.5*(deltar+n_genes-taonum-1))
//				*pow(fabs((sigma0+S0).determinant()),-0.5*(num+deltar+n_genes-taonum-1));
//
//			for(flag=1; flag<=clusternum; flag++){
//				cclusternum = newz.flagnum(flag);
//				if(cclusternum != 0){
//					for(j=1; j<=taonum; j++)
//						gammaprod *= gammaf(0.5*(cclusternum+deltar+taonum-j))/gammaf(0.5*(deltar+taonum-j));
//
//					K1 = z_Kc_weight(newz, flag)*pow(beta1*newz.flagnum(flag)+1,-(double)taonum/2.0)*gammaprod;
//
//					gammaprod = 1.0;
//
//
//
//					mid = z_datacov(newz, flag, 1, mu1, mu0);
//
//					S1 = mid.getBlock(1,dim,1,dim) 
//						+mid.getBlock(dim+1,2*dim,1,dim)*(cclusternum/(beta1*cclusternum+1));
//
//
//					den *= K1*pow(fabs(sigma1.determinant()),0.5*(deltar+taonum-1))
//						*pow(fabs((sigma1+S1).determinant()),-0.5*(cclusternum+deltar+taonum-1));
//
//				}
//			}
//		} else if(taonum==0){
//			for(j=1; j<=(n_genes-taonum); j++)
//				gammaprod *= gammaf(0.5*(num+deltar+n_genes-taonum-j))/gammaf(0.5*(deltar+n_genes-taonum-j));
//
//			K0 = pow(beta0*num+1,-0.5*(n_genes-taonum))*gammaprod;
//
//			gammaprod = 1.0;
//
//			mid = z_datacov(newz, 1, 0, mu1, mu0);
//
//			S0 = mid.getBlock(1,dim,1,dim)
//				+ mid.getBlock(dim+1,2*dim,1,dim)*(num/(beta0*num+1));
//
//			den = K0 * pow(fabs(sigma0.determinant()),0.5*(deltar+n_genes-taonum-1))
//				*pow(fabs((sigma0+S0).determinant()),-0.5*(num+deltar+n_genes-taonum-1));
//
//
//		} else{
//
//
//			for(flag=1; flag<=clusternum; flag++){
//				cclusternum = newz.flagnum(flag);
//				if(cclusternum != 0){
//					for(j=1; j<=taonum; j++)
//						gammaprod *= gammaf(0.5*(cclusternum+deltar+taonum-j))/gammaf(0.5*(deltar+taonum-j));
//
//					K1 = z_Kc_weight(newz, flag)*pow(beta1*newz.flagnum(flag)+1,-(double)taonum/2.0)*gammaprod;
//
//					gammaprod = 1.0;
//
//
//
//					mid = z_datacov(newz, flag, 1, mu1, mu0);
//
//					
//
//					S1 = mid.getBlock(1,dim,1,dim) 
//						+mid.getBlock(dim+1,2*dim,1,dim)*(cclusternum/(beta1*cclusternum+1));
//
//
//					den *= K1*pow(fabs(sigma1.determinant()),0.5*(deltar+taonum-1))
//						*pow(fabs((sigma1+S1).determinant()),-0.5*(cclusternum+deltar+taonum-1));
//
//				}
//			}
//
//
//		}
//
//		return den;
//}
//
//
//double mcmc::zw_density(const ClusterFlags &newz, const matrix &newweight, const matrix &mu1, const matrix &mu0, 
//	const matrix &sigma0, const matrix &sigma1, const double beta0, const double beta1,const double deltar)const{
//
//		double den = 1.0;
//		long flag;
//		long i,j;
//		long clusternum = newz.getMaxFlag();
//		long cclusternum;
//		double gammaprod = 1.0;
//		double K0, K1;
//		long dim = n_genes;
//		matrix S0, S1, mid;
//
//		long taonum=0;
//		for(i=0; i<dim; i++)
//			taonum += tao.at(i);
//
//
//		if((taonum>0)&&(taonum<n_genes)){
//			for(j=1; j<=(n_genes-taonum); j++)
//				gammaprod *= gammaf(0.5*(num+deltar+n_genes-taonum-j))/gammaf(0.5*(deltar+n_genes-taonum-j));
//
//			K0 = pow(beta0*num+1,-0.5*(n_genes-taonum))*gammaprod;
//
//			gammaprod = 1.0;
//
//			mid = zw_datacov(newz, 1, 0, mu1, mu0);
//
//			S0 = mid.getBlock(1,dim,1,dim)
//				+ mid.getBlock(dim+1,2*dim,1,dim)*(num/(beta0*num+1));
//
//			den = K0 * pow(fabs(sigma0.determinant()),0.5*(deltar+n_genes-taonum-1))
//				*pow(fabs((sigma0+S0).determinant()),-0.5*(num+deltar+n_genes-taonum-1));
//
//			for(flag=1; flag<=clusternum; flag++){
//				cclusternum = newz.flagnum(flag);
//				if(cclusternum != 0){
//					for(j=1; j<=taonum; j++)
//						gammaprod *= gammaf(0.5*(cclusternum+deltar+taonum-j))/gammaf(0.5*(deltar+taonum-j));
//
//					K1 = zw_Kc_weight(newz, newweight, flag)*pow(beta1*newz.flagnum(flag)+1,-(double)taonum/2.0)*gammaprod;
//
//					gammaprod = 1.0;
//
//
//
//					mid = zw_datacov(newz, flag, 1, mu1, mu0);
//
//					S1 = mid.getBlock(1,dim,1,dim) 
//						+mid.getBlock(dim+1,2*dim,1,dim)*(cclusternum/(beta1*cclusternum+1));
//
//
//					den *= K1*pow(fabs(sigma1.determinant()),0.5*(deltar+taonum-1))
//						*pow(fabs((sigma1+S1).determinant()),-0.5*(cclusternum+deltar+taonum-1));
//
//				}
//			}
//		} else if(taonum==0){
//			for(j=1; j<=(n_genes-taonum); j++)
//				gammaprod *= gammaf(0.5*(num+deltar+n_genes-taonum-j))/gammaf(0.5*(deltar+n_genes-taonum-j));
//
//			K0 = pow(beta0*num+1,-0.5*(n_genes-taonum))*gammaprod;
//
//			gammaprod = 1.0;
//
//			mid = zw_datacov(newz, 1, 0, mu1, mu0);
//
//			S0 = mid.getBlock(1,dim,1,dim)
//				+ mid.getBlock(dim+1,2*dim,1,dim)*(num/(beta0*num+1));
//
//			den = K0 * pow(fabs(sigma0.determinant()),0.5*(deltar+n_genes-taonum-1))
//				*pow(fabs((sigma0+S0).determinant()),-0.5*(num+deltar+n_genes-taonum-1));
//
//
//		} else{
//
//			for(flag=1; flag<=clusternum; flag++){
//				cclusternum = newz.flagnum(flag);
//				if(cclusternum != 0){
//					for(j=1; j<=taonum; j++)
//						gammaprod *= gammaf(0.5*(cclusternum+deltar+taonum-j))/gammaf(0.5*(deltar+taonum-j));
//					
//					/*std::cout<<gammaprod<<std::endl;
//					std::cout<<pow(beta1*((double)cclusternum)+1.0,-(double)taonum/2.0)<<std::endl;
//					std::cout<<zw_Kc_weight(newz, newweight, flag)<<std::endl;*/
//					K1 = zw_Kc_weight(newz, newweight, flag)*pow(beta1*((double)cclusternum)+1.0,-(double)taonum/2.0)*gammaprod;
//				
//					gammaprod = 1.0;
//
//
//
//					mid = zw_datacov(newz, flag, 1, mu1, mu0);
//					/*std::cout<<"mid"<<std::endl;
//					std::cout<<mid<<std::endl;*/
//					S1 = mid.getBlock(1,dim,1,dim) 
//						+mid.getBlock(dim+1,2*dim,1,dim)*(cclusternum/(beta1*cclusternum+1));
//					/*std::cout<<"zwdensity_K1"<<std::endl;
//					std::cout<<K1<<std::endl;
//					std::cout<<"zwdensity_Q0"<<std::endl;
//					std::cout<<pow(fabs(sigma1.determinant()),0.5*(deltar+taonum-1))<<std::endl;
//					std::cout<<"zwdensity_Q0+S1"<<std::endl;
//					std::cout<<pow(fabs((sigma1+S1).determinant()),-0.5*(cclusternum+deltar+taonum-1))<<std::endl;*/
//					den *= K1*pow(fabs(sigma1.determinant()),0.5*(deltar+taonum-1))
//						*pow(fabs((sigma1+S1).determinant()),-0.5*(cclusternum+deltar+taonum-1));
//
//				}
//			}
//
//
//		}
//
//		return den;
//}
//
//
//
//// Note that in each step, cluster index should be reassigned using function reflag.
//void mcmc::tao_sample(const matrix &mu1, const matrix &mu0, const matrix &sigma0, const matrix &sigma1, const double beta0, 
//	const double beta1, const double deltar)
//{
//	std::vector<bool> newtao(n_genes);
//	long i;
//	long accepted = 1;
//	long coin1, coin2, coin3;
//	bool tao_mid;
//	long tao_length = newtao.size();
//	double acceptedp = 0.0;
//	for(i=0; i<tao_length; i++)
//		newtao.at(i) = tao.at(i);
//	coin1 = ran_ber(0.5);
//	if(coin1 == 1){
//		coin2 = ran_iunif(0,tao_length-1);
//		newtao.at(coin2) = !newtao.at(coin2); //1. Change value.
//	} else{
//		coin2 = ran_iunif(0,tao_length-1);
//		coin3 = ran_iunif(0,tao_length-2);
//		if(coin3>=coin2){
//			coin3++;
//		}
//		tao_mid = newtao.at(coin2);
//		newtao.at(coin2) = newtao.at(coin3);
//		newtao.at(coin3) = tao_mid;          //or 2. Swap;
//	}
//
//	acceptedp = f_fmin(1.0, tao_density(newtao, mu1, mu0, sigma0, sigma1, beta0, beta1, 
//		deltar)/density(mu1, mu0, sigma0, sigma1, beta0, beta1, deltar));
//
//	accepted = ran_ber(acceptedp);
//	if(accepted == 1){
//		for(i=0; i<tao_length; i++)
//			tao.at(i) = newtao.at(i);
//	} 
//
//	/*std::cout<<"Tao"<<std::endl;
//	for(i=0; i<tao_length; i++)
//	std::cout<<tao.at(i)<<' ';
//	std::cout<<std::endl;*/
//
//
//}
//
//// 2011/11/10 General 
//// alpha is part of the shape parameter, shape is the scale parameter in gamma distribution.
//void mcmc::weight_sample(const double alpha, const double shape, const matrix &mu1, const matrix &mu0, const matrix &sigma0,
//	const matrix &sigma1, const double beta0, const double beta1, const double deltar)
//{
//
//	long t,i,flag,n;
//	long fnum;
//	double wnum=0.0;
//	double mid;
//	std::vector<long> midv;
//	for(t=1; t<n_timepoints; t++){
//		wnum = 0.0;
//		midv = weight.tclusterindex(t+1);
//		n = midv.size();
//		//std::cout<<"Clusternumber: "<<n<<std::endl;
//		for(i=0; i<n; i++){
//			flag = midv.at(i);
//			fnum = z.tflagnum(flag,t+1);
//			weight.setValue(t,flag-1,ran_gamma(alpha+fnum,shape));
//			wnum += weight.getValue(t,flag-1);
//		}
//
//		for(i=0; i<n; i++){
//			flag = midv.at(i);
//			mid = weight.getValue(t,flag-1);
//			weight.setValue(t,flag-1,mid/wnum);
//		}
//		//midv.clear();
//	}
//
//	/*std::cout<<"weight"<<std::endl;
//	std::cout<<weight<<std::endl;*/
//}
//
//
//
//
//
//
//
//void mcmc::z_sample(const matrix &mu1, const matrix &mu0, const matrix & sigma0, const matrix &sigma1, const double beta0, 
//	const double beta1, const double deltar){
//		long i,j,k,l,midnum;
//		std::vector<long> flagl;
//		std::vector<double> prob_z;
//		//ClusterFlags z_mid = z;
//		for(j=0; j<z.cols(0); j++)
//			z(0,j) = 1;
//		for(i=1; i<z.rows(); i++){
//
//			flagl = weight.tclusterindex(i+1);
//			for(j=0; j<z.cols(i); j++){
//				
//				for(k=0; k<flagl.size(); k++){
//				
//				    z(i,j) = flagl.at(k);
//					prob_z.push_back(z_density(z, mu1, mu0, sigma0, sigma1, beta0, beta1, deltar));
//					// optimize?
//					
//			    }
//				z(i,j) = ran_num(prob_z,flagl);
//				//z_mid(i,j) = z(i,j);
//				/*std::cout<<"Data"<<"("<<i<<","<<j<<")"<<std::endl;
//				 for(l=0; l<flagl.size(); l++)
//				std::cout<<prob_z.at(l)<<' ';
//				std::cout<<std::endl;*/
//				/*if(i==3&&j==0){
//				   for(l=0; l<flagl.size(); l++)
//				       std::cout<<prob_z.at(l)<<' ';
//				   std::cout<<std::endl;
//				}*/
//				
//				prob_z.clear();
//			}	 
//			//flagl.clear();
//		}
//		// Do we need to reflag the cluster indicators?
//		//z.reflag();
//		/*std::cout<<"Sampled z"<<std::endl;
//		std::cout<<z<<std::endl;*/
//}
//




// Return data number before a time point
//long mcmc::f_tdatanum(const long time)const{
//	long i,s = 0;
//
//	for(i=0; i<time; i++)
//		s += NT.at(i);
//	return s;
//}

//// Return data number before a time point and a location number;
//// Note that time and loc is the real reporting number.
//long mcmc::f_tidatanum(const long time, const long loc)const{
//	long i,s;
//	s = f_tdatanum(time-1);
//	return s+loc;
//}
//
//// Return average data value for a cluster
//// Here we use average to measure the similarity of two clusters.
//matrix mcmc::f_clustermean(const long time, const long clusterindex, const ClusterFlags &one, const matrix &two)const {
//	matrix gdata(1,n_genes);
//	long i,j,s=0,n=0;
//	for(i=time-1; i<n_timepoints; i++){
//		for(j=0; j<one.cols(i); j++){
//			if(one(i,j)==clusterindex){
//				gdata += Data.getRow(f_tidatanum(i+1,j+1)-1);
//				n++;
//			}
//		}
//	}
//
//	gdata *= (1/(double)n);
//
//	return gdata;
//
//
//
//}
//


// Return distance between two clusters
//double mcmc::f_clusterdist(const long time1, const long clusterindex1, const long time2, const long clusterindex2,
//	const ClusterFlags &one, const matrix &two)const{
//		double c=0.0;
//		matrix a,b;
//		long i;
//		a = this->f_clustermean(time1, clusterindex1, one, two);
//		b = this->f_clustermean(time2, clusterindex2, one, two);
//		/*std::cout<<a<<std::endl;
//		std::cout<<b<<std::endl;*/
//		for(i=0; i<n_genes; i++){
//			c +=pow(a(0,i)-b(0,i),2);
//		}
//		return c;
//
//}
//

// Return an integar indicating the relationship of the two splitted clusters
// 0: a cluster is empty.
// 1: The relationship is not closest.
// 2: The relationship is closest.
// NULLset tells clusterson is empty or not.
//long mcmc::f_similarity(const long time, const long clustermother, const long clusterson, const bool NULLset,
//	const matrixl splitsets,
//	const ClusterFlags &one, const matrix &two) const{
//		long n = splitsets.rows();
//		long i,m=2;
//		double dist_ms, dist_mi;
//		if(NULLset==true){
//			return 0;
//		}else{
//			dist_ms = this->f_clusterdist(time,clustermother,time,clusterson,one,two);
//			for(i=0; i<n; i++){
//				dist_mi = this->f_clusterdist(time,clustermother,splitsets(i,0),splitsets(i,1),one,two);
//				if(dist_ms>dist_mi){
//					m = 1;
//					break;
//				}
//			}
//		}
//		return m;
//}
//

//
//// P_chos is designed different from JASA2005. 
//double mcmc::f_pchos(const long time, const long clustermother, const long clusterson, const bool NULLset,
//	const matrixl splitsets,
//	const ClusterFlags &one, const matrix &two)const{
//
//		long G_star = two.clusternumber();
//		long G1_star = one.nonemptynum();
//		long G0_star = G_star - G1_star;
//		long n;
//		double s = 0.0;
//		/*std::cout<<"Three chosen probabilities."<<std::endl;
//		std::cout<<1.0/((double)(G_star*G1_star))<<std::endl;
//		std::cout<<1.0/(double)G_star<<std::endl;
//		std::cout<<2.0/(double)G_star<<std::endl;*/
//		if(NULLset){
//			return 1.0/((double)(G_star*G1_star));
//		}else{
//			n = this->f_similarity(time,clustermother,clusterson,NULLset,splitsets,one,two);
//
//			s = (n==2)?(2.0/(double)G_star):(1.0/(double)G_star);
//
//			return s;
//			//return s/(double)time;
//			//return pow(s,time-1);
//			//return pow(s,time);
//		}
//		
//
//}











// This function must use twice before and after 'flag_split' function.
// Return the clusters that can do split move.
// The result structure is matrix with the row element is (time, clusterindex).
// We consider all cluster vector except the one which already has next generation.
// Furthermore, the start cluster in the splitable cluster vector must be nonempty.
//matrixl mcmc::f_splitset(const ClusterFlags &newz, const matrix &newweight, const matrixl &newtree,
//	const matrixl &newchome)const{
//mcmc::SplitMergeSet mcmc::f_splitset(const MCMCEnv &env) const {
//
//		long i,j;
//		lvector splitone;
//		//Collect nonempty cluster vector.
//		std::vector<long> flags = newz.getTtoTindicator(n_timepoints);
//
//		long midf,midt,mids,sr = 0;
//		long fl = flags.size();
//		for(i=0; i<fl; i++){
//			midf = flags.at(i);
//			midt = newz.getEarliestTime(midf);
//			
//
//			for(j=1; j<n_timepoints; j++){ 
//
//				mids = (j==midt-1)?(newchome(0,midf-1)-1):(midf-1);
//
//				if((newtree(j-1,mids)==0)&&(newz.flaginT(midf,j+1))){
//				   splitone.push_back(j+1);
//				   splitone.push_back(midf);
//				   sr++;
//				}
//			}
//			
//		}
//
//
//		return matrixl(sr,2,splitone);
//}
//
//
//
//// Accept probability for split move, the input cluster index must be nonempty at least in the beginning.
//// newz is new cluster indicator after split move
//// newweight is new weight matrix after split move
//// This function follows a random sample from 'splitset' result.
//double mcmc::AP_split(const long time, const long clusterindex, const bool NULLset, const matrixl &splitsets,
//	const long splitnum, const double P_alloc, const double f_ui,
//	const ClusterFlags &newz, const matrix &newweight, const double bk,
//	const matrix &mu1, const matrix &mu0, 
//	const matrix &sigma0, const matrix &sigma1, const double beta0,	const double beta1, const double deltar,
//	const double alpha)const{
//		double densityratio, modelratio;
//		double weightratio = 1.0;
//		double P_chos = 1.0;
//		double AP = 1.0;
//		double dbratio = 1.0;
//		double jacobi = 1.0;
//		long i;
//		if(splitnum == 0){
//			return 0.0;
//		}else{
//			densityratio = this->zw_density(newz,newweight,mu1,mu0,sigma0,sigma1,beta0,beta1,deltar)/
//				this->density(mu1,mu0,sigma0,sigma1,beta0,beta1,deltar);
//			for(i=time-1; i<n_timepoints; i++){
//				// If the first cluster is an empty cluster, directly return 0.0.
//				if(weight(i,clusterindex-1)>SMALLNUM){
//					weightratio *= pow(newweight(i,clusterindex-1),alpha-1)*pow(newweight(i,clusterindex),alpha-1)
//						/pow(weight(i,clusterindex-1),alpha-1)/betaf(alpha,((double)weight.tclusternumber(i+1))*alpha);
//					jacobi *= weight(i,clusterindex-1);
//				}
//			}
//			// The ratio of f(Tree_new) and f(Tree_old).
//			// If we set the same unchange and change probability, then it is equal to 1.0; 
//			modelratio = 1.0;
//			// The ratio of d_new and b_old, here dbratio = 1.0, because d_new = b_old = 0.5.
//			dbratio = (1.0-bk)/bk;
//			// Parameter 'splitnum': The number of cluster vector that can be used to do split move.
//			//splitnum = 1.0;//
//			// A function must be written to summarize the number of individuals in new and old clusters
//			//P_alloc = 1.0;
//			// Probabilities must be designed to describe the merge clusters.
//			P_chos = this->f_pchos(time,clusterindex,clusterindex+1,NULLset,splitsets,newz,newweight);
//
//
//			AP = densityratio*weightratio*modelratio*dbratio*P_chos*((double)splitnum)/f_ui/P_alloc*jacobi;
//
//			/*std::cout<<"Split"<<std::endl;
//			std::cout<<"densityratio"<<std::endl;
//			std::cout<<densityratio<<std::endl;
//			std::cout<<"weightratio"<<std::endl;
//			std::cout<<weightratio<<std::endl;
//			std::cout<<"P_chos"<<std::endl;
//			std::cout<<P_chos<<std::endl;
//			std::cout<<"splitnum"<<std::endl;
//			std::cout<<splitnum<<std::endl;
//			std::cout<<"f_ui"<<std::endl;
//			std::cout<<f_ui<<std::endl;
//			std::cout<<"P_alloc"<<std::endl;
//			std::cout<<P_alloc<<std::endl;
//			std::cout<<"Jacobi"<<std::endl;
//			std::cout<<jacobi<<std::endl;*/
//			//std::cout<<"Split accepted probability: "<<std::endl;
//			return f_fmin(AP,1.0);
//		}
//}
//
//
//
//
//
//
//// Generate an empty cluster
//// Different from JASA 2005, to keep the tree structure, this empty cluster is generated from a branch.
//// This function is almost the same to flag_split except the clusterflags are not changed.
//void mcmc::flag_birth(ClusterFlags &one, matrix &two, matrixl &newtree, matrixl &newchome, double &weightratio,
//	double &jacobi,
//	double &f_wstar, const long time, const long clusterindex, const double alpha)const{
//		long i,j,k,n,timec;
//		double w_star;
//
//		jacobi = 1.0;
//		f_wstar = 1.0;
//		weightratio = 1.0;
//		if((time>=2)&&(time<=one.rows())&&(time>0)&&(one.flaginT(clusterindex, time))){
//			/*Justify whether clusterindex is an empty cluster.*/
//
//
//			for(i=1; i<=n_timepoints; i++){
//				if(i<time){
//					two(i-1,clusterindex) = 0.0;
//					for(j=clusterindex+1; j<two.cols(); j++)
//						two(i-1,j) = weight(i-1,j-1);
//					for(k=0; k<one.cols(i-1); k++){ 
//						one(i-1,k) = (z(i-1,k)<=clusterindex)?z(i-1,k):(z(i-1,k)+1);
//					}
//
//
//				}else{
//					/*weight*/
//					timec = weight.tclusternumber(i);
//					w_star = ran_beta(1.0,(double)timec);
//					jacobi *= pow(1-w_star,timec-1);
//					f_wstar *= betad(w_star, 1.0, (double)timec);
//					weightratio *= (pow(w_star, alpha-1)*pow(1-w_star, 
//						(double)timec*(alpha-1))/betaf(alpha, (double)timec*alpha));
//
//					/* This is a good choice for resetting the cluster indicators and weight parameter */
//					//two(i-1,clusterindex-1) = (1.0 - w_star) * weight(i-1,clusterindex-1);
//					two(i-1,clusterindex) = w_star;
//
//					for(j=0; j<two.cols(); j++){
//						if(j < clusterindex){
//							two(i-1,j) = (1.0 - w_star)*weight(i-1,j);
//						}else{
//							if(j>clusterindex){
//								two(i-1,j) = (1.0 - w_star)*weight(i-1,j-1);
//							}
//						}
//
//					}
//					/*weight relocated.*/
//					for(k=0; k<one.cols(i-1); k++){
//						one(i-1,k) = (z(i-1,k)<=clusterindex)?z(i-1,k):(z(i-1,k)+1);
//					}
//
//
//				}
//
//			}
//
//
//			// Tree reset
//			// This step is implemented according to the weight setting.
//			newtree(time-2,clusterindex-1) = 1;
//			for(i=0; i<n_timepoints; i++){
//				newtree(i,clusterindex) = 0;
//				for(j=clusterindex+1; j<two.cols(); j++){
//					newtree(i,j) = tree(i,j-1);
//				}
//			}
//
//			newchome(0,clusterindex) = clusterindex;
//			for(i=clusterindex+1; i<newchome.cols(); i++){
//				newchome(0,i) = (clusterhome(0,i-1)<=clusterindex)?(clusterhome(0,i-1)):(clusterhome(0,i-1)+1);
//			}
//
//
//			/*std::cout<<"Birth"<<std::endl;
//			std::cout<<"weightratio"<<std::endl;
//			std::cout<<weightratio<<std::endl;
//			std::cout<<"Jacobi"<<std::endl;
//			std::cout<<jacobi<<std::endl;
//			std::cout<<"f_wstar"<<std::endl;
//			std::cout<<f_wstar<<std::endl;*/
//
//
//		}else{
//			std::cout<<"This is an empty cluster!"<<std::endl;
//			std::cout<<"This cluster can not be splitted!"<<std::endl;
//		}
//
//
//}
//
//
//
//
//








// Return -1 if cluster does not exist.
//long mcmc::EmptyCluster(const ClusterFlags &one, const matrix &two, const long clusterindex)const{
//	long i,j;
//	long s = 0;
//	bool x = false;
//	for(i=0; i<n_timepoints; i++){
//		if(two(i,clusterindex-1)>SMALLNUM)
//			s++;
//	}
//
//	if(s == 0)
//	{
//		//std::cout<<"This cluster does not exist and can not be deleted!"<<std::endl;
//		return -1;
//	}else{
//		s = 0;
//		for(i=0; i<n_timepoints; i++){
//			if(one.flaginT(clusterindex,i+1)){
//
//				return 0;
//
//			}
//
//		}
//
//	}
//
//
//	return 1;
//}
//

//
//
//// Only nonempty cluster vectors are considered.
//// Return (time, clustermother, clusterson)
//matrixl mcmc::f_mergeset(const ClusterFlags &newz, const matrix &newweight, const matrixl &newtree)const{
//
//	//Collect nonempty cluster vector.
//	std::vector<long> flags = newz.getTtoTindicator(n_timepoints);
//	std::vector<long> freeflags;
//	std::vector<long> nonfreeflags;
//	lvector mergeone;
//	long i,j,k,s=0;
//	long time1, time2, time0;
//	long sr=0;
//	long n,m,l;
//
//	n = flags.size();
//	if(n != 1){
//		for(i=0; i<n; i++){
//			s=0;
//			for(j=0; j<n_timepoints; j++){
//				if(newtree(j,flags.at(i)-1)!=0){
//					s++;
//					nonfreeflags.push_back(flags.at(i));
//					break;
//				}
//			}
//			if(s==0){
//				freeflags.push_back(flags.at(i));
//			}
//
//		}
//		m = freeflags.size();
//		l = n-m;
//
//
//		for(i=0; i<m; i++){
//			time1 = newweight.getEarliestTime(freeflags.at(i));
//			for(k=i+1; k<m; k++){
//				time0 = newz.getEarliestTime(freeflags.at(k));
//				if(time1<=time0){
//					mergeone.push_back(time0);
//					mergeone.push_back(freeflags.at(i));
//					mergeone.push_back(freeflags.at(k));
//					sr++;
//				}else{
//					mergeone.push_back(time1);
//					mergeone.push_back(freeflags.at(k));
//					mergeone.push_back(freeflags.at(i));
//					sr++;
//				}
//			}
//
//
//			for(j=0; j<l; j++){
//				time2 = newweight.getEarliestTime(nonfreeflags.at(j));
//				if(time2<=time1){
//					mergeone.push_back(time1);
//					mergeone.push_back(nonfreeflags.at(j));
//					mergeone.push_back(freeflags.at(i));
//					sr++;
//				}
//
//			}
//		}
//	}
//	return matrixl(sr,3,mergeone);
//
//}
//
//
//
//
//
//matrixl mcmc::f_deathset(const ClusterFlags &newz, const matrix &newweight, const matrixl &newtree)const{
//	std::vector<long> emptyflags;
//	std::vector<long> nonemptyflags;
//	lvector deathsets;
//	long i,j,ec,n,m0,m1,time0, time1;
//	long sr = 0;
//	n = newweight.cols();
//	for(i=0; i<n; i++){
//		ec = this->EmptyCluster(newz, newweight, i+1);
//		if(ec == 1){
//			emptyflags.push_back(i+1);
//		}else{
//			if(ec==0){
//				nonemptyflags.push_back(i+1);
//			}
//		}
//	}
//	m0 = emptyflags.size();
//	m1 = nonemptyflags.size();
//	for(i=0; i<m1; i++){
//		time0 = newweight.getEarliestTime(nonemptyflags.at(i));
//		for(j=0; j<m0; j++){
//			time1 = newweight.getEarliestTime(emptyflags.at(j));
//			if((time0<=time1)&&(newtree(time1-1,emptyflags.at(j)-1)==0)){
//				deathsets.push_back(time1);
//				deathsets.push_back(nonemptyflags.at(i));
//				deathsets.push_back(emptyflags.at(j));
//				sr++;
//			}
//		}
//	}
//
//	return matrixl(sr,3,deathsets);
//
//}
//
//
//
//
//
//
//
//
//
//double mcmc::AP_merge(const long time, const long clustermother, const long clusterson, 
//	const matrixl &splitsetsaftermerge, const long mergesetsrows,
//	const long splitnumaftermerge, const double P_alloc, const double f_ui,
//	const ClusterFlags &newz, const matrix &newweight, const double bk,
//	const matrix &mu1, const matrix &mu0, const matrix &sigma0, const matrix &sigma1, 
//	const double beta0, const double beta1, const double deltar, const double alpha)const{
//
//		double densityratio,modelratio,bdratio; 
//		double weightratio = 1.0;
//		double jacobi = 1.0;
//		double P_chos;
//		long i;
//		double AP;
//		bool NULLset = false;
//
//		if((mergesetsrows==0)||(splitnumaftermerge==0)){
//			return 0.0;
//		}else{
//			densityratio = this->zw_density(newz,newweight,mu1,mu0,sigma0,sigma1,beta0,beta1,deltar)/
//				this->density(mu1,mu0,sigma0,sigma1,beta0,beta1,deltar);
//			/*std::cout<<"newz"<<std::endl;
//			std::cout<<newz<<std::endl;
//			std::cout<<"newweight"<<std::endl;
//			std::cout<<newweight<<std::endl;
//			std::cout<<"NEW DENSITY"<<std::endl;
//			std::cout<<zw_density(newz,newweight,mu1,mu0,sigma0,sigma1,beta0,beta1,deltar)<<std::endl;
//			std::cout<<"OLD DENSITY"<<std::endl;
//			std::cout<<density(mu1,mu0,sigma0,sigma1,beta0,beta1,deltar)<<std::endl;*/
//            
//			if(clustermother<clusterson){
//			   for(i=time-1; i<n_timepoints; i++){
//				   if((weight(i,clustermother-1)>SMALLNUM)&&(weight(i,clusterson-1)>SMALLNUM)){
//				      weightratio *= pow(newweight(i,clustermother-1),alpha-1)*
//					      betaf(alpha,((double)newweight.tclusternumber(i+1))*alpha)/pow(weight(i,clustermother-1),alpha-1)/
//					      pow(weight(i,clusterson-1),alpha-1);
//				      jacobi *= 1/newweight(i,clustermother-1);///  ?????!!!!!!!!  ///
//			
//				   }
//
//			   }
//			}else{// If clustermother>clusterson, it will be abstracted by 1 after merge moving.
//				for(i=time-1; i<n_timepoints; i++){
//				   if((weight(i,clustermother-1)>SMALLNUM)&&(weight(i,clusterson-1)>SMALLNUM)){
//				      weightratio *= pow(newweight(i,clustermother-2),alpha-1)*
//					      betaf(alpha,((double)newweight.tclusternumber(i+1))*alpha)/pow(weight(i,clustermother-1),alpha-1)/
//					      pow(weight(i,clusterson-1),alpha-1);
//				      jacobi *= 1/newweight(i,clustermother-2);///  ?????!!!!!!!!  ///
//			
//				   }
//
//			   }
//			  
//			}
//			modelratio = 1.0;
//			bdratio = bk/(1.0-bk);
//
//			P_chos = this->f_pchos(time,clustermother,clusterson,NULLset,splitsetsaftermerge,newz,newweight);
//
//		}
//
//
//		AP = densityratio*weightratio*modelratio*bdratio*P_alloc/((double)splitnumaftermerge)*f_ui/P_chos*jacobi;
//
//		    /*std::cout<<"Merge"<<std::endl;
//			std::cout<<"densityratio"<<std::endl;
//			std::cout<<densityratio<<std::endl;
//			std::cout<<"weightratio"<<std::endl;
//			std::cout<<weightratio<<std::endl;
//			std::cout<<"P_chos"<<std::endl;
//			std::cout<<P_chos<<std::endl;
//			std::cout<<"splitnumaftermerge"<<std::endl;
//			std::cout<<splitnumaftermerge<<std::endl;
//			std::cout<<"f_ui"<<std::endl;
//			std::cout<<f_ui<<std::endl;
//			std::cout<<"P_alloc"<<std::endl;
//			std::cout<<P_alloc<<std::endl;
//			std::cout<<"Jacobi"<<std::endl;
//			std::cout<<jacobi<<std::endl;
//			std::cout<<"Merge accepted probability: "<<std::endl;*/
//
//		return f_fmin(AP,1.0);
//}
//
//
//
//
//
//double mcmc::AP_death(const long time, const long clustermother, const long clusterson,
//	const long emptynumbeforedeath, const double f_wstar, const double weightratio, const double jacobi,
//	const ClusterFlags &newz, const matrix &newweight, const double bk, 
//	const matrix &mu1, const matrix &mu0, const matrix &sigma0, const matrix &sigma1,
//	const double beta0, const double beta1, const double deltar, const double alpha)const{
//
//		double densityratio,modelratio,bdratio;
//		double AP;
//
//		if(emptynumbeforedeath == 0){
//			return 0.0;
//		}else{
//			densityratio = this->zw_density(newz,newweight,mu1,mu0,sigma0,sigma1,beta0,beta1,deltar)/
//				this->density(mu1,mu0,sigma0,sigma1,beta0,beta1,deltar);
//			modelratio = 1.0;
//			bdratio = bk/(1.0-bk);
//		}
//
//
//		AP = densityratio*weightratio*modelratio*bdratio*f_wstar*emptynumbeforedeath*jacobi;
//
//
//		return f_fmin(AP,1.0);
//}
//
//
//
//
//
//
//


//
//
//
//
//matrixl mcmc::f_tailbirthset(const ClusterFlags &newz, const matrix &newweight, const matrixl &newtree)const{
//	lvector tailsets;
//	std::vector<long> clusterindex = newz.getTtoTindicator(n_timepoints);
//	long i,j,k,n,id,t;
//	long sr = 0;
//	long s = 0;
//	n = clusterindex.size();
//	for(i=0; i<n; i++){
//		id = clusterindex.at(i);
//		t = newz.getEarliestTime(id);
//		for(j = n_timepoints-1; j>t-2; j--){
//			if((newweight(j,id-1)>SMALLNUM)&&(j==n_timepoints-1)){
//				break;
//			}else{
//				if((newweight(j,id-1)>SMALLNUM)&&(j<n_timepoints-1)){
//					tailsets.push_back(j+2);
//					tailsets.push_back(id);
//					sr++;
//					break;
//				}
//			}
//		}
//
//	}
//
//	return matrixl(sr,2,tailsets);
//}
//



// The chosen cluster should begin with a nonempty cluster and end with series of empty cluster.
// All zero-weight cluster will be set nonzero weight
// ####OOOOOO
//void mcmc::tail_birth(const long time, const long clusterindex, 
//	ClusterFlags &newz, matrix &newweight, const double alpha)const{
//		long i,j;
//		double w_star;
//
//
//		for(i=time-1; i<n_timepoints; i++){
//
//			newweight(i,clusterindex-1) = w_star;
//			for(j=0; j<clusterindex-1; j++)
//				newweight(i,j) = weight(i,j)*(1.0-w_star);
//			for(j=clusterindex; j<newweight.cols(); j++)
//				newweight(i,j) = weight(i,j)*(1.0-w_star);
//		}
//
//
//
//
//}
//
//
//double mcmc::AP_birth(const long time, const long clusterindex, const long splitnumafterbirth, 
//	const double weightratio, const double jacobi, const double f_wstar, 
//	const ClusterFlags &newz, const matrix &newweight, const double bk,
//	const matrix &mu1, const matrix &mu0,
//	const matrix &sigma0, const matrix &sigma1, const double beta0, const double beta1, const double deltar,
//	const double alpha)const{
//		double densityratio, modelratio;
//		double AP = 1.0;
//		double dbratio = 1.0;
//		long i;
//		if(splitnumafterbirth == 0){
//			return 0.0;
//		}else{
//			densityratio = this->zw_density(newz,newweight,mu1,mu0,sigma0,sigma1,beta0,beta1,deltar)/
//				this->density(mu1,mu0,sigma0,sigma1,beta0,beta1,deltar);
//
//			// The ratio of f(Tree_new) and f(Tree_old).
//			// If we set the same unchange and change probability, then it is equal to 1.0; 
//			modelratio = 1.0;
//			// The ratio of d_new and b_old, here dbratio = 1.0, because d_new = b_old = 0.5.
//			dbratio = (1.0-bk)/bk;
//
//
//			AP = densityratio*weightratio*modelratio*dbratio*jacobi/(double)splitnumafterbirth/f_wstar;
//
//			/*std::cout<<"densityratio"<<std::endl;
//			std::cout<<densityratio<<std::endl;
//			std::cout<<"weightratio"<<std::endl;
//			std::cout<<weightratio<<std::endl;
//			std::cout<<"splitnumafterbirth"<<std::endl;
//			std::cout<<splitnumafterbirth<<std::endl;
//			std::cout<<"f_wstar"<<std::endl;
//			std::cout<<f_wstar<<std::endl;
//			std::cout<<"Jacobi"<<std::endl;
//			std::cout<<jacobi<<std::endl;*/
//			return f_fmin(AP,1.0);
//		}
//}
//
//
//double mcmc::AP_tailbirth(
//	const ClusterFlags &newz, const matrix &newweight,
//	const matrix &mu1, const matrix &mu0,
//	const matrix &sigma0, const matrix &sigma1, const double beta0, const double beta1, const double deltar,
//	const double alpha)const{
//		double AP;
//
//
//		AP = this->zw_density(newz,newweight,mu1,mu0,sigma0,sigma1,beta0,beta1,deltar)/
//			this->density(mu1,mu0,sigma0,sigma1,beta0,beta1,deltar);
//
//		return f_fmin(AP,1.0);
//
//}
//


//
//
//// Return (time, clustermother, deathtail)
//matrixl mcmc::f_taildeathset(const ClusterFlags &newz, const matrix &newweight, const matrixl &newtree)const{
//	lvector tailsets;
//	std::vector<long> t_cluster;
//	std::vector<long> t_tails, t_nontails;
//	long i,j,k,l,id,t;
//	long sr = 0;
//	long s = 0;
//	long n,m;
//
//	for(i=0; i<n_timepoints; i++){
//		t_cluster = newz.getTindicator(i+1);
//		for(j=0; j<t_cluster.size(); j++){
//			id = t_cluster.at(j);
//			t = newz.getEarliestTime(id);
//			if((newweight(i,id-1)>SMALLNUM)&&(i>=t)){
//				for(k=i; k<n_timepoints; k++){
//					if(newtree(k,id-1) != 0){
//						s++;
//						t_nontails.push_back(id);
//						break;
//					}
//				}
//				if(s==0){
//					t_tails.push_back(id);
//
//				}
//				s = 0;
//			}else{
//				t_nontails.push_back(id);
//			}
//		}
//		// Find merge sets;
//		n = t_tails.size();
//		m = t_nontails.size();
//
//
//
//		for(j=0; j<n-1; j++){
//			for(k=j+1; k<n; k++){
//				tailsets.push_back(i+1);
//				tailsets.push_back(t_tails.at(k));
//				tailsets.push_back(t_tails.at(j));
//				sr++;
//			}
//
//		}
//		for(j=0; j<n; j++){
//			for(l=0; l<m; l++){
//				tailsets.push_back(i+1);
//				tailsets.push_back(t_nontails.at(l));
//				tailsets.push_back(t_tails.at(j));
//				sr++;
//			}
//		}
//
//
//
//		t_tails.clear();
//		t_nontails.clear();
//
//	}
//
//
//
//
//	return matrixl(sr,3,tailsets);
//} 
//
//
//
//
//// 11/02/2011 new parameter 'newtree' added.
//// clusterindex 1 and 2 should not be empty clusters 
//// One of these two clusters must be merged.
//// clusterindex1 corresponds to clustermother, clusterindex2 corresponds to clusterson
//void mcmc::flag_merge(ClusterFlags &one, matrix &two, matrixl &newtree, matrixl &newchome,
//	double &P_alloc, double &f_ui, long time, long clusterindex1, 
//	long clusterindex2)const{
//		long mergeflag, leftflag;
//		long midnum;
//		long mergenum, leftnum;
//		long i,j,k;
//		double ut;
//
//		P_alloc = 1.0;
//		f_ui = 1.0;
//		if((time<=one.rows())&&(time>=0)&&(clusterindex2>=2)){
//			// Mother
//			leftflag = clusterindex1;
//			// Son
//			mergeflag = clusterindex2;
//			
//
//			for(i=1; i<=n_timepoints; i++){
//				
//
//				for(j=mergeflag-1; j<two.cols()-1; j++){
//					two(i-1,j) = weight(i-1,j+1);
//				}
//
//
//				if(leftflag<mergeflag){
//					two(i-1, leftflag-1) = weight(i-1, leftflag-1) + weight(i-1, mergeflag-1);
//				}else{
//					two(i-1, leftflag-2) = weight(i-1, leftflag-1) + weight(i-1, mergeflag-1);
//				}
//
//
//
//				if((i>=time)&&(weight(i-1, leftflag-1)>SMALLNUM)&&(weight(i-1, mergeflag-1)>SMALLNUM)){
//					midnum = (leftflag<mergeflag)?1:2;
//					ut = weight(i-1,leftflag-1)/two(i-1,leftflag-midnum);
//					leftnum = z.tflagnum(leftflag,i);
//					mergenum = z.tflagnum(mergeflag,i);
//					P_alloc *= pow(ut,leftnum)*pow(1-ut,mergenum);
//					f_ui *= betad(ut,2.0,2.0);
//				}
//
//				for(k=0; k<one.cols(i-1); k++){
//					
//					if(z(i-1,k) != mergeflag){
//						one(i-1,k) = (z(i-1,k) < mergeflag)?z(i-1,k):(z(i-1,k)-1);
//					} else{
//						one(i-1,k) = (leftflag<mergeflag)?leftflag:(leftflag-1); 
//					}
//					// cluster indicator reloaded.
//				}
//				
//			}
//
//			
//
//
//			newtree(time-2,clusterhome(0,mergeflag-1)-1) = 0;
//			
//			
//
//			for(i=0; i<n_timepoints; i++){
//				
//				for(j=mergeflag-1; j<two.cols()-1; j++){
//					newtree(i,j) = tree(i,j+1);
//				}
//			}
//
//			//if(leftflag<mergeflag){
//				for(i=mergeflag; i<newchome.cols()-1; i++){
//					newchome(0,i-1) = (clusterhome(0,i)<=mergeflag)?(clusterhome(0,i)):(clusterhome(0,i)-1);
//				}
//
//			//}else{
//			//	for(i=mergeflag; i<newchome.cols()-1; i++){
//			//		newchome(0,i-1) = (clusterhome(0,i)<=mergeflag)?(clusterhome(0,i)):(clusterhome(0,i)-1);
//			//	}
//			//}
//
//			/*std::cout<<"new"<<std::endl;
//			std::cout<<one<<std::endl;
//			std::cout<<two<<std::endl;
//			std::cout<<newtree<<std::endl;
//			std::cout<<newchome<<std::endl;
//			std::cout<<P_alloc<<std::endl;
//			std::cout<<f_ui<<std::endl;*/
//
//		} else{
//			std::cout<<"These two clusters can not be merged!"<<std::endl;
//			std::cout<<"Maybe empty cluster or maybe the time is not one of its starting time point!"<<std::endl;
//		}
//
//}
//
//
//
//
//
//// Death of an empty cluster
//// clusterindex2 corresponds to an empty cluster vector, clusterindex1 corresponds to a nonempty cluster vector.
//void mcmc::flag_death(ClusterFlags &one, matrix &two, matrixl &newtree, matrixl &newchome,
//	double &jacobi, double &weightratio, double &f_wstar, long time, long clusterindex1, long clusterindex2,
//	const double alpha)const{
//
//		long mergeflag, leftflag;
//		long mergenum, leftnum;
//		long i,j,k;
//		double midw,midn;
//
//		jacobi = 1.0;
//		weightratio = 1.0;
//		f_wstar = 1.0;
//
//		if((time<=one.rows())&&(time>=0)&&(clusterindex2>=2)){
//			// Mother
//			leftflag = clusterindex1;
//			// Son
//			mergeflag = clusterindex2;
//
//
//			for(i=0; i<n_timepoints; i++){
//				midw = weight(i,mergeflag-1);
//				if(midw > SMALLNUM){
//					if(i>=time-1){
//						midn = two.tclusternumber(i+1);
//						f_wstar *= betad(midw,1.0,(double)midn);
//						weightratio *= (betaf(alpha,(double)(midn-1)*alpha)/pow(midw,(alpha-1.0))
//							/pow(1.0-midw,(double)(midn-1)*(alpha-1.0)));
//						jacobi *= pow(1.0-midw,-1.0*(double)(midn-2));
//					}
//
//
//					for(j=0; j<mergeflag-1; j++){
//						two(i,j) = weight(i,j)/(1.0-midw);
//					}
//					for(j=mergeflag-1; j<(two.cols()-1); j++){
//						two(i,j) = weight(i,j+1)/(1.0-midw);
//					}
//
//				}else{
//					for(j=mergeflag-1; j<(two.cols()-1); j++)
//						two(i,j) = weight(i,j+1);
//				}
//
//
//
//				for(k=0; k<one.cols(i); k++){
//					one(i,k) = (z(i,k) <= mergeflag)?z(i,k):(z(i,k)-1);  
//				}
//				// cluster indicator reloaded.
//			}
//
//
//
//
//
//			newtree(time-2,clusterhome(0,mergeflag-1)-1) = 0;
//			for(i=0; i<n_timepoints; i++){
//				for(j=mergeflag-1; j<two.cols()-1; j++){
//					newtree(i,j) = tree(i,j+1);
//				}
//			}
//
//			if(leftflag<mergeflag){
//				for(i=mergeflag; i<newchome.cols()-1; i++){
//					newchome(0,i-1) = (clusterhome(0,i)<=mergeflag)?(clusterhome(0,i)):(clusterhome(0,i)-1);
//				}
//
//			}else{
//				for(i=mergeflag; i<newchome.cols()-1; i++){
//					newchome(0,i-1) = (clusterhome(0,i)<=mergeflag)?(clusterhome(0,i)):(clusterhome(0,i)-1);
//				}
//			}
//
//
//
//
//
//			/*std::cout<<"Death"<<std::endl;
//			std::cout<<"weightratio"<<std::endl;
//			std::cout<<weightratio<<std::endl;
//			std::cout<<"Jacobi"<<std::endl;
//			std::cout<<jacobi<<std::endl;
//			std::cout<<"f_wstar"<<std::endl;
//			std::cout<<f_wstar<<std::endl;*/
//
//		} else{
//			std::cout<<"These two clusters can not be merged!"<<std::endl;
//			std::cout<<"Maybe empty cluster or maybe the time is not one of its starting time point!"<<std::endl;
//		}
//
//}
//
//
//
//// 'clusterindex' is actually the parameter 'clustermother'.
//// After running flag_split, 'clusterson' is clusterindex+1.
//// 11/01/2011 Add new parameter 'newtree'.
//// Split!!! Return a resampled flag item, weight and parameter
//// time is the actual stage point,clusterindex is the flag which we want to summarize.
//// Note that one and two should be put into the parameter part of flag_split.
//// IMPORTANT! clusternumber should be comparied between time and time-1 before using this function.
//// one corresponds to a copy of cluster indicator z, two is a copy of weight matrix.
//// one can not be the cluster indiator z itself, the same is two.
//// After splitting the original flags, there may not be any order existing in the final flags.
//// Parameter 'time' means split move happens in stage 'time', while this move should be reported in 'time -1' row
//// in parameter 'newtree', so time must bigger than 1 (>=2).
//// 
//// 2011/12/26
//// Rewritten the data structure
//
//void mcmc::flag_split(ClusterFlags &one, matrix &two, matrixl &newtree, matrixl &newchome, 
//	double &P_alloc, bool &NULLset,
//	double &f_ui, const long time, 
//	const long clusterindex)const{
//		//ClusterFlags one = z;
//		//matrix two = weight;
//		long i,j,k,n,s=0,maxflag;
//		double ut;
//		double midpalloc = 1.0;
//		f_ui = 1.0;
//		if((time>=2)&&(time<=one.rows())&&(one.flaginT(clusterindex, time))){
//			/*Justify whether clusterindex is an empty cluster.*/
//
//			maxflag = z.getMaxFlag();
//			//std::cout<<"maxflag"<<std::endl;
//			//std::cout<<maxflag<<std::endl;
//			for(i=1; i<=n_timepoints; i++){
//				if(i<time){
//					two(i-1,clusterindex) = 0.0;
//					for(j=clusterindex+1; j<two.cols(); j++)
//						two(i-1,j) = weight(i-1,j-1);
//
//
//					for(k=0; k<one.cols(i-1); k++){
//						// cluster indicator reloaded. 
//						one(i-1,k) = (z(i-1,k)<=clusterindex)?z(i-1,k):(z(i-1,k)+1);
//
//					}
//				}else{
//
//					/*weight*/
//					ut = ran_beta(2.0,2.0);
//
//					f_ui *= (weight(i-1,clusterindex-1) > SMALLNUM)?betad(ut,2.0,2.0):1.0;
//					/* This is a good choice for resetting the cluster indicators and weight parameter */
//					two(i-1,clusterindex-1) = ut * weight(i-1,clusterindex-1);
//					two(i-1,clusterindex) = (1-ut) * weight(i-1,clusterindex-1);
//
//					for(j=clusterindex+1; j<two.cols(); j++){
//						two(i-1,j) = weight(i-1,j-1);
//					}
//					
//					/*weight relocated.*/
//					for(k=0; k<one.cols(i-1); k++){
//						if((weight(i-1,clusterindex-1) > SMALLNUM)&&(z(i-1,k)==clusterindex)){
//							one(i-1,k) = ran_nber(clusterindex, maxflag+1, ut);	// optimize? clusterindex+1 instead of maxflag+1
//							if(one(i-1,k) == clusterindex){
//							   //one(i-1,k) = one(i-1,k);
//							   midpalloc *= ut;
//							}else{
//							   one(i-1,k) = clusterindex + 1;
//							   midpalloc *= (1-ut);
//							   s++;
//							}			
//						} else{
//							one(i-1,k) = (z(i-1,k)<clusterindex)?z(i-1,k):(z(i-1,k)+1);
//						}
//
//						// cluster indicator reloaded.
//					}
//
//
//				}
//
//			}
//
//			P_alloc = midpalloc;
//			// Tree reset
//			// This step is implemented according to the weight setting.
//			newtree(time-2,clusterindex-1) = 1;
//			for(i=0; i<n_timepoints; i++){
//				newtree(i,clusterindex) = 0;
//				for(j=clusterindex+1; j<two.cols(); j++){
//					newtree(i,j) = tree(i,j-1);
//				}
//			}
//
//			newchome(0,clusterindex) = clusterindex;
//			for(i=clusterindex+1; i<newchome.cols(); i++){
//				newchome(0,i) = (clusterhome(0,i-1)<=clusterindex)?(clusterhome(0,i-1)):(clusterhome(0,i-1)+1);
//			}
//
//
//			NULLset = (s==0)?true:false;
//			/*std::cout<<"new"<<std::endl;
//			std::cout<<one<<std::endl;
//			std::cout<<two<<std::endl;
//			std::cout<<newtree<<std::endl;
//			std::cout<<newchome<<std::endl;
//			std::cout<<P_alloc<<std::endl;
//			std::cout<<NULLset<<std::endl;*/
//
//		}else{
//			std::cout<<"This is an empty cluster!"<<std::endl;
//			std::cout<<"This cluster can not be splitted!"<<std::endl;
//		}
//}
//
//


//
//// Nonempty cluster vector clusterindex will be set 0 weights from timebegin;
//// To make tail death, we combine it with the nearest cluster vector.
//void mcmc::tail_death(const long time, const long clusterindex1, const long clusterindex2,
//	ClusterFlags &newz, matrix &newweight, const double alpha)const{
//		long i,j,k;
//		long leftflag, mergeflag;
//
//
//		// Mother
//		leftflag = clusterindex1;
//		// Son
//		mergeflag = clusterindex2;
//
//
//		for(i=time; i<=n_timepoints; i++){
//
//			newweight(i-1, mergeflag-1) = 0.0;
//			newweight(i-1, leftflag-1) = weight(i-1, leftflag-1) + weight(i-1, mergeflag-1);
//
//
//			for(k=0; k<newz.cols(i-1); k++){
//				if(z(i-1,k) == mergeflag){
//					newz(i-1,k) = leftflag;  
//				}
//				// cluster indicator reloaded.
//			}
//
//		}
//
//}
//
//
//
//double mcmc::AP_taildeath(const ClusterFlags &newz, const matrix &newweight, 
//	const matrix &mu1, const matrix &mu0,
//	const matrix &sigma0, const matrix &sigma1, const double beta0, const double beta1, const double deltar,
//	const double alpha)const{
//		double AP;
//
//		AP = this->zw_density(newz,newweight,mu1,mu0,sigma0,sigma1,beta0,beta1,deltar)/
//			this->density(mu1,mu0,sigma0,sigma1,beta0,beta1,deltar);
//		return f_fmin(AP,1.0);
//
//}
//



//void mcmc::SplitMerge_move(double bk, ClusterFlags &newz, matrix &newweight, matrixl &newtree, matrixl &newchome,
//	const matrix &mu0, const matrix &mu1, const matrix &sigma0, 
//	const matrix &sigma1, const double beta0, const double beta1,const double deltar, const double alpha, 
//	bool &move, double &probability){
bool mcmc::SplitMerge_move(double &probability) {
	bool move_chosen;
	MCMCEnv::TreeSet n_splitset;
	MCMCEnv::TreeMergeSet n_mergeset;
	MCMCEnv::TreeSet::value_type chos_splitset;
	MCMCEnv::TreeMergeSet::value_type chos_mergeset;
	MCMCEnv newEnv(env);
	double P_alloc = 1.0;
	bool NULLset = false;
	double f_ui = 1.0;

	move_chosen = ran_ber(MCMCEnv::bk);

	if(move_chosen){
		n_splitset = newEnv.getSplitSet();
		unsigned long long childID;

		//std::cout<<"Split chosen " << n_splitset.size() <<std::endl;
		// Justify there is indeed cluster set that can be used to do split move.
		if(n_splitset.size()){
			chos_splitset = n_splitset[ran_iunif(0, n_splitset.size() - 1)];
			//std::cout<<"Split set chosen"<<std::endl;
			// Note that splitnum is the splitted cluster vector before splitting.
			childID = newEnv.flagSplit(P_alloc, NULLset, f_ui, chos_splitset);
			// Splitset after split move, it can be used to compute G*,G1*,G0* in P_chos.
			n_splitset = newEnv.getSplitSet();
			//P_chos = this->f_pchos(chos_splitset(0,0),chos_splitset(0,1),chos_splitset(0,1)+1,NULLset,n_splitset,
			//                       newz,newweight);
			probability = newEnv.apSplit(chos_splitset, NULLset, n_splitset, 
				n_splitset.size(), P_alloc, f_ui, env, childID);
			if(ran_unif(0.0,1.0) < probability) {
				env = newEnv;
			}
			//std::cout<<"split"<<std::endl;
			/*
			std::cout<<ranmite<<std::endl;
			std::cout<<z<<std::endl;
			std::cout<<weight<<std::endl;
			std::cout<<tree<<std::endl;
			std::cout<<clusterhome<<std::endl;*/
		}
	} else {// Merge move is chosen
		//std::cout<<"Merge chosen"<<std::endl;
		n_mergeset = newEnv.getMergeSet();
		if(n_mergeset.size()){
			chos_mergeset = n_mergeset[ran_iunif(0, n_splitset.size() - 1)];
			//cout<<"Merge set chosen"<<endl;
			//cout<<chos_mergeset<<endl;
			newEnv.flagMerge(P_alloc, f_ui, chos_mergeset);
			//cout<<"Merge flag complete"<<endl;
			// cluster vectors can be used to do split move after merge.
			n_splitset = newEnv.getSplitSet();
			probability = newEnv.apMerge(chos_mergeset, n_splitset, 
				n_mergeset.size(), n_splitset.size(), P_alloc, f_ui, env);
			//cout<<"APmerge complete, prob = " << probability <<endl;
			//cout<<APmerge<<endl;
			if(ran_unif(0.0,1.0) < probability){
				env = newEnv;
			}
			/*std::cout<<"merge"<<std::endl;
			std::cout<<ranmite<<std::endl;
			std::cout<<z<<std::endl;
			std::cout<<weight<<std::endl;
			std::cout<<tree<<std::endl;
			std::cout<<clusterhome<<std::endl;*/
		}
	}
	return move_chosen;
}

bool mcmc::SplitMerge_move_test(double &probability) {
	bool move_chosen;
	MCMCEnv::TreeSet n_splitset;
	MCMCEnv::TreeMergeSet n_mergeset;
	MCMCEnv::TreeSet::value_type chos_splitset;
	MCMCEnv::TreeMergeSet::value_type chos_mergeset;
	env.testNumberSet();
	MCMCEnv newEnv(env);
	double P_alloc = 1.0;
	bool NULLset = false;
	double f_ui = 1.0;

	
	chos_splitset = MCMCEnv::TreeSet::value_type(1, 1);
	cout << "Mergemove_test" << env.calcLogDensity() << " " << newEnv.calcLogDensity() << endl;
			//std::cout<<"Split set chosen"<<std::endl;
			// Note that splitnum is the splitted cluster vector before splitting.
			//newEnv.flagSplitTest(P_alloc, NULLset, f_ui, chos_splitset);
			// Splitset after split move, it can be used to compute G*,G1*,G0* in P_chos.
		n_mergeset = newEnv.getMergeSet();
		if(n_mergeset.size()){
			chos_mergeset = n_mergeset[ran_iunif(0, n_splitset.size() - 1)];
			//cout<<"Merge set chosen"<<endl;
			//cout<<chos_mergeset<<endl;
			newEnv.flagMerge(P_alloc, f_ui, chos_mergeset);
			//cout<<"Merge flag complete"<<endl;
			// cluster vectors can be used to do split move after merge.
			n_splitset = newEnv.getSplitSet();
			probability = newEnv.apMerge(chos_mergeset, n_splitset, 
				n_mergeset.size(), n_splitset.size(), P_alloc, f_ui, env);
			//cout<<"APmerge complete, prob = " << probability <<endl;
			//cout<<APmerge<<endl;
			cerr << "oldEnv" << endl << env;
			env.writeStatus(cerr);

			cerr << "newEnv" << endl << newEnv;
			newEnv.writeStatus(cerr);

			if(ran_unif(0.0,1.0) < probability){
				env = newEnv;
			}

			
			cerr << probability << endl;
			throw(0);
			/*std::cout<<"merge"<<std::endl;
			std::cout<<ranmite<<std::endl;
			std::cout<<z<<std::endl;
			std::cout<<weight<<std::endl;
			std::cout<<tree<<std::endl;
			std::cout<<clusterhome<<std::endl;*/
		}
			//n_splitset = newEnv.getSplitSet();
			////P_chos = this->f_pchos(chos_splitset(0,0),chos_splitset(0,1),chos_splitset(0,1)+1,NULLset,n_splitset,
			////                       newz,newweight);
			//probability = newEnv.apMerge(chos_splitset, NULLset, n_splitset, 
			//	n_splitset.size(), P_alloc, f_ui, env);
			//if(ran_unif(0.0,1.0) < probability) {
			//	env = newEnv;
			//}
			//cerr << probability << endl;
			//throw(0);
			//std::cout<<"split"<<std::endl;
			/*
			std::cout<<ranmite<<std::endl;
			std::cout<<z<<std::endl;
			std::cout<<weight<<std::endl;
			std::cout<<tree<<std::endl;
			std::cout<<clusterhome<<std::endl;*/
	return move_chosen;
}

bool mcmc::BirthDeath_move(double &probability) {

	bool move_chosen;
	MCMCEnv::TreeSet n_splitset;
	MCMCEnv::TreeMergeSet n_deathset;
	MCMCEnv::TreeSet::value_type chos_splitset;
	MCMCEnv::TreeMergeSet::value_type chos_deathset;
	MCMCEnv newEnv(env);
	long splitnumafterbirth, emptynumbeforedeath;
	double weightratio = 1.0;
	double jacobi = 1.0;
	double f_wstar = 1.0;

	move_chosen = ran_ber(MCMCEnv::bk);
	if(move_chosen) { 
		// Empty cluster birth.
		//std::cout<<"BirthBirthBirth"<<std::endl;
		n_splitset = newEnv.getSplitSet();
		// Justify there is indeed cluster set that can be used to do birth move.
		if(n_splitset.size()){
			chos_splitset = n_splitset[ran_iunif(0, n_splitset.size() - 1)];
			//cout<<"Birth set chosen"<<endl;
			newEnv.flagBirth(chos_splitset, weightratio, jacobi, f_wstar);
			// Number of empty sets after birth				
			splitnumafterbirth = newEnv.getDeathSet().size();
			probability = newEnv.apBirth(chos_splitset, splitnumafterbirth,
				weightratio, jacobi, f_wstar, env);
			//std::cout<<probability<<std::endl;
			if(ran_unif(0.0, 1.0) < probability){
				env = newEnv;
			}
			/*std::cout<<ranmite<<std::endl;
			std::cout<<z<<std::endl;
			std::cout<<weight<<std::endl;
			std::cout<<tree<<std::endl;
			std::cout<<clusterhome<<std::endl;*/
		}
	} else {
		// Empty cluster death.
		//std::cout<<"DeathDeathDeath"<<std::endl;
		n_deathset = newEnv.getDeathSet();
		emptynumbeforedeath = n_deathset.size();
		if(emptynumbeforedeath){
			chos_deathset = n_deathset[ran_iunif(0, n_splitset.size() - 1)];
			//cout<<"Death set chosen"<<endl;
			newEnv.flagDeath(chos_deathset,jacobi, weightratio, f_wstar);
			//cout<<"Flagdeath complete"<<endl;
			probability = newEnv.apDeath(chos_deathset, emptynumbeforedeath, f_wstar,
				weightratio, jacobi, env);
			//cout<<"ApDeath complete, probability = " << probability <<endl;
			//std::cout<<probability<<std::endl;

			if(ran_unif(0.0, 1.0) < probability) {
				env = newEnv;
			}
		}

	}
	return move_chosen;
}

//
//
//
//// Notet that 'TailBirthDeath_move' is not a step for model selection or tree jumping, because
//// the number of parameters is not changed.
//void mcmc::TailBirthDeath_move(const double ranc, ClusterFlags &newz, matrix &newweight, matrixl &newtree,
//	matrixl &newchome,
//	const matrix &mu0, const matrix &mu1, const matrix &sigma0, 
//	const matrix &sigma1, const double beta0, const double beta1,const double deltar, const double alpha){
//
//		double move_chosen; 
//		matrixl n_tailbirthset, chos_tailbirthset;
//		matrixl n_taildeathset, chos_taildeathset;
//		double APtailbirth,APtaildeath;
//
//		double ranmite;
//
//
//		move_chosen = ran_ber(ranc);
//		if(move_chosen){ 
//			// Tail birth.
//			n_tailbirthset = this->f_tailbirthset(newz,newweight,newtree);
//
//			if(n_tailbirthset.rows()!=0){
//				chos_tailbirthset = n_tailbirthset.ran_getRow();
//				this->tail_birth(chos_tailbirthset(0,0),chos_tailbirthset(0,1),
//					newz,newweight,alpha);
//				APtailbirth = this->AP_tailbirth(newz,newweight,
//					mu1,mu0,sigma0,sigma1,beta0,beta1,deltar,alpha);
//				ranmite = ran_unif(0.0,1.0);
//				if(ranmite < APtailbirth){
//					z = newz;
//					weight = newweight;
//					tree = newtree;
//					clusterhome = newchome;
//				}else{
//				    newz = z;
//					newweight = weight;
//					newtree = tree;
//					newchome = clusterhome;
//				}
//			}
//
//		}else{
//			// Tail death.
//			n_taildeathset = this->f_taildeathset(newz,newweight,newtree);
//			if(n_taildeathset.rows()!=0){
//				chos_taildeathset = n_taildeathset.ran_getRow();
//				this->tail_death(chos_taildeathset(0,0),chos_taildeathset(0,1),chos_taildeathset(0,2),
//					newz,newweight,alpha);
//				APtaildeath = this->AP_taildeath(newz,newweight,mu1,mu0,sigma0,sigma1,beta0,beta1,deltar,alpha);
//				ranmite = ran_unif(0.0,1.0);
//				if(ranmite < APtaildeath){
//					z = newz;
//					weight = newweight;
//					tree = newtree;
//					clusterhome = newchome;
//				}else{
//				    newz = z;
//					newweight = weight;
//					newtree = tree;
//					newchome = clusterhome;
//				}
//			}
//
//
//		}
//}


//
//
//void mcmc::C_sample(const matrix &mu0, const matrix &mu1, const matrix &sigma0, 
//	const matrix &sigma1, const double beta0, const double beta1,const double deltar, const double alpha,
//	const double shape, const double bk, const double ranc)
//{
//	bool move;
//	MCMCEnv newenv(env);
//	ClusterFlags newz;
//	//matrix newweight;
//	//matrixl newtree;
//	//matrixl newchome;
//	double probability;
//
//
//	tao_sample(mu1,mu0,sigma0,sigma1,beta0,beta1,deltar);
//
//	weight_sample(alpha,shape,mu1,mu0,sigma0,sigma1,beta0,beta1,deltar);
//
//	z_sample(mu1,mu0,sigma0,sigma1,beta0,beta1,deltar);
//
//
//
//	newz = z;
//	//newweight = weight;
//	//newtree = tree;
//	//newchome = clusterhome;
//
//
//	SplitMerge_move(bk,newz,newweight,newtree,newchome,mu0,mu1,sigma0,sigma1,beta0,beta1,deltar,alpha,move,
//		probability);
//	
//
//
//
//	BirthDeath_move(bk,newz,newweight,newtree,newchome,mu0,mu1,sigma0,sigma1,beta0,beta1,deltar,alpha,move,
//		probability);
//
//
//
//	//TailBirthDeath_move(ranc,newz,newweight,newtree,newchome,mu0,mu1,sigma0,sigma1,beta0,beta1,deltar,alpha);
//	//newz = z;
//	//newweight = weight;
//	//newtree = tree;
//	//newchome = clusterhome;
//
//}


// Because parameters of tree and clusterhome must be considerred at one time in order to 
// generate a value summarizing a tree result, it is inconvient for us to set a longitudal value
// for a tree.
//double mcmc::treesummary()const{
//   long n = tree.rows();
//   long m = tree.cols();
//   long i, j, ss=0, k;
//   long cc, ct, zt, chome;
//   double s = 0.0;
//   for(j = 0; j < m; j ++){
//	   for(i = 0; i < n; i ++){
//		   if(tree(i, j) != 0){
//			   for(k = 1; k < clusterhome.cols(); k++){
//				   if(clusterhome(0,k)==(j+1)){
//				      cc = k+1;
//                      // Does cluster cc come from cluster j+1?
//				      ct = weight.getEarliestTime(cc);
//				      // Is cc an empty cluster?
//				      zt = z.getEarliestTime(cc);
//
//				      if((ct == i+2)&&(zt != -1)){
//				         s += pow(10.0,i)*(1.0+0.5*(double)j);
//				      }
//				   }
//				   
//				   
//			   }
//	         
//		   }
//	   }
//
//
//	   
//   }
//   return s;
//}
//
//
//
//double mcmc::treesummary(const matrixl &newtree, const matrixl &newchome, const matrix &newweight, 
//	const ClusterFlags &newz)const{
//   long n = newtree.rows();
//   long m = newtree.cols();
//   long i, j, ss=0, k;
//   long cc, ct, zt, chome;
//   double s = 0.0;
//   for(j = 0; j < m; j ++){
//	   for(i = 0; i < n; i ++){
//		   if(newtree(i, j) != 0){
//			   //std::cout<<i<<","<<j<<std::endl;
//			   for(k = 1; k < newchome.cols(); k++){
//				   if(newchome(0,k)==(j+1)){
//				      cc = k+1;
//				      ct = newweight.getEarliestTime(cc);
//				      zt = newz.getEarliestTime(cc);
//				      if((ct == i+2)&(zt != -1)){
//				         s += pow(10.0,i)*(1.0+0.5*(double)i);
//				      }
//				   }
//				   
//			   }
//	         
//		   }
//	   }
//
//
//	   
//   }
//   return s;
//}
//
//
//
//
////std::ostream & operator<<(std::ostream &os, const mcmc &one) {
////	long i_t = one.gettime();
////	long i_g = one.getgene();
////	int i;
////	os<<"Time length:"<<std::endl;
////	os<<i_t;
////	os<<std::endl;
////	os<<"Gene number:"<<std::endl;
////	os<<i_g;
////	os<<std::endl;
////	os<<"Data number:"<<std::endl;
////	os<<one.getnum();
////	os<<std::endl;
////	os<<"Data number in each cluster:"<<std::endl;
////	for(i=0; i<i_t; i++)
////		os<<one.getNT().at(i)<<' ';
////	os<<std::endl;
////	os<<"Cluster number in each cluster:"<<std::endl;
////	for(i=0; i<i_t; i++)
////		os<<one.getclusternum().at(i)<<' ';
////	os<<std::endl;
////	os<<"Cluster proportions in each cluster: "<<std::endl;
////	os<<one.getweight();
////	os<<std::endl;
////	os<<"Featureed genes: "<<std::endl;
////	for(i=0; i<i_g; i++)
////		os<<one.gettao().at(i)<<' ';
////	os<<std::endl;
////	os<<"Cluster index: "<<std::endl;
////	os<<one.getz();
////	os<<std::endl;
////	os<<"Tree: "<<std::endl;
////	os<<one.gettree();
////	os<<std::endl;
////	os<<"The original data: "<<std::endl;
////	os<<one.getData();
////	os<<std::endl;
////
////
////	return os;
////}
////
//
//
//
//
//
//void mcmc::parameter_initialization(matrix &mu1, matrix &mu0, matrix &sigma1, matrix &sigma0, 
//	double &beta1, double &beta0, double &deltar, double &alpha, double &shape, double &bk, double &ranc)const{
//
//	double k1,k0;
//	// K1 and deltar are both used to control data effect, it should be 
//	// Note that K1 must be <= 1;
//	// With a large or small K1, the algorithm both seems converge very slowly.
//	k1 = 10.0;
//	k0 = 10.0;
//	mu1 = Data.matrix_mean();
//	mu0 = mu1;
//	sigma1 *= 1.0/k1;
//	sigma0 *= 1.0/k0;
//
//	// Beta controls the prior effect to the sampling probability
//	// A large Beta value should be set, if you want to decrease the prior effect.
//	beta1 = 1000.0;
//	beta0 = 1000.0;
//	alpha = 1.0;
//	// The data effect is controlled by deltar, if you want to increase data effect, a large deltar should be set.
//	// However the minimum of 
//	deltar = 3.0;
//	shape = 1.0;
//
//
//
//
//
//}

//
//
//// AP_split function for testing 
//double mcmc::AP_split(const long time, const long clusterindex, const bool NULLset, 
//	const double f_ui,
//	const ClusterFlags &newz, const matrix &newweight, const double bk,
//	const matrix &mu1, const matrix &mu0, 
//	const matrix &sigma0, const matrix &sigma1, const double beta0,	const double beta1, const double deltar,
//	const double alpha)const{
//		double densityratio, modelratio;
//		double weightratio = 1.0;
//		double P_chos = 1.0;
//		double AP = 1.0;
//		double dbratio = 1.0;
//		double jacobi = 1.0;
//		long splitnum=1;
//		double P_alloc;
//		long i;
//		if(splitnum == 0){
//			return 0.0;
//		}else{
//			densityratio = this->zw_density(newz,newweight,mu1,mu0,sigma0,sigma1,beta0,beta1,deltar)/
//				this->density(mu1,mu0,sigma0,sigma1,beta0,beta1,deltar);
//			for(i=time-1; i<n_timepoints; i++){
//				// If the first cluster is an empty cluster, directly return 0.0.
//				if(weight(i,clusterindex-1)>SMALLNUM){
//					weightratio *= pow(newweight(i,clusterindex-1),alpha-1)*pow(newweight(i,clusterindex),alpha-1)
//						/pow(weight(i,clusterindex-1),alpha-1)/betaf(alpha,((double)weight.tclusternumber(i+1))*alpha);
//					jacobi *= weight(i,clusterindex-1);
//				}
//			}
//			// The ratio of f(Tree_new) and f(Tree_old).
//			// If we set the same unchange and change probability, then it is equal to 1.0; 
//			modelratio = 1.0;
//			// The ratio of d_new and b_old, here dbratio = 1.0, because d_new = b_old = 0.5.
//			dbratio = (1.0-bk)/bk;
//			// Parameter 'splitnum': The number of cluster vector that can be used to do split move.
//			//splitnum = 1.0;//
//			// A function must be written to summarize the number of individuals in new and old clusters
//			//P_alloc = 1.0;
//			// Probabilities must be designed to describe the merge clusters.
//			P_chos = 2.0/3.0;
//			P_alloc = pow(0.5,10);
//			splitnum = 2;
//
//
//			AP = densityratio*weightratio*modelratio*dbratio*P_chos*((double)splitnum)/f_ui/P_alloc*jacobi;
//
//			/*std::cout<<"Split"<<std::endl;
//			std::cout<<"densityratio"<<std::endl;
//			std::cout<<densityratio<<std::endl;
//			std::cout<<"weightratio"<<std::endl;
//			std::cout<<weightratio<<std::endl;
//			std::cout<<"P_chos"<<std::endl;
//			std::cout<<P_chos<<std::endl;
//			std::cout<<"splitnum"<<std::endl;
//			std::cout<<splitnum<<std::endl;
//			std::cout<<"f_ui"<<std::endl;
//			std::cout<<f_ui<<std::endl;
//			std::cout<<"P_alloc"<<std::endl;
//			std::cout<<P_alloc<<std::endl;
//			std::cout<<"Jacobi"<<std::endl;
//			std::cout<<jacobi<<std::endl;
//			std::cout<<"Split accepted probability: "<<std::endl;*/
//			return f_fmin(AP,1.0);
//		}
//}
//


#undef PI
#undef LARGENUM
#undef SMALLNUM