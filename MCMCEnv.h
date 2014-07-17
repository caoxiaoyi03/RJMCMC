#ifndef MCMCENV_H
#define MCMCENV_H
#pragma once

#include "ClusterTree.h"
#include "ClusterFlags.h"
#include "Matrix.h"

#include <sstream>

#ifdef LINUX
#include <tr1/unordered_map>
#else
#include <unordered_map>
#endif


class MCMCEnv
{
public:
	typedef std::vector<bool> TaoSet;
	typedef std::tr1::unordered_map<unsigned long long, ClusterTree> TreeContainer;
	typedef std::vector<std::pair<unsigned long, unsigned long long> > TreeSet;
	// time and ID

	typedef std::vector<std::pair<unsigned long long, unsigned long long> > TreeMergeSet;
	// parent ID and child ID

private:
	typedef std::pair<unsigned long, unsigned long> TimeCoordinate;

	typedef std::map<unsigned long, unsigned long> SortedULContainer;
	typedef std::map<unsigned long long, char> RootIDContainer;		// char is a dummy holder

	TreeContainer MapClusterTrees;
	unsigned long long ID_max;
	RootIDContainer rootIDs;

	// Notice that the row of flags is time, column is sample
	//ClusterFlags Flags;
	// tao
	TaoSet Tao;

	static unsigned long NumOfTP;
	static unsigned long NumOfGenes;
	static unsigned long NumOfData;
	static ClusterTree::NumOfSampleContainer DataNumEachTP;	// number of data in each time point
	static SortedULContainer TPStartIndex, ReverseTPStartIndex;	// the starting index for each time point

	MCMCEnv(const TaoSet &tao, const ClusterFlags &flags, const matrix &data);

	friend class mcmc;

	// The mutable variants are the values intended to speed up the calculation process
	mutable bool UpToDate;		// if true, no calculation is needed

	mutable double LogDensityValue;
	mutable double LogDensityValueNonFlag;
	mutable double LogDensityValueFlag;

	static matrix mu1;
	static matrix mu0;
	static matrix sigma0;
	static matrix sigma1;
	static double beta0;
	static double beta1;
	static double deltar;
	static double alpha;
	static double shape;
	static double bk; // For split_merge move and birth_death move, we give equal probability.
	static double ranc; // Tail birth is given much more porbability.
	static double k1;
	static double k0;
	static double constSigma0;
	static double constSigma1;
	
	static bool buildForest;		// whether to build a forest or simply a tree

	const matrix &Data;

	TimeCoordinate IndexToCoor(unsigned long index) const;
	unsigned long CoorToIndex(const TimeCoordinate &corr) const;

public:
	ClusterFlags Flags;
	static MCMCEnv initMCMCEnv(unsigned long numOfTimePoints, unsigned long numOfGenes,
		unsigned long numOfData, const std::vector<unsigned long> &numOfDataInEachTimePoint,
		const TaoSet &tao, const ClusterFlags &flags, const matrix &p_data);
	MCMCEnv(const MCMCEnv &oldenv);		// duplicate the environment
	~MCMCEnv(void);

	MCMCEnv &operator=(const MCMCEnv &oldenv);

	static unsigned long getDataNum();

	matrix reData(unsigned long long flag, bool tao_NULL, const TaoSet &newTao, 
		std::vector<unsigned long> &taoIndices) const;
	matrix reData(unsigned long long flag, bool tao_NULL, const TaoSet &newTao, std::ostream &os) const;
	matrix datacov(unsigned long long flag, bool tao_NULL, const TaoSet &newTao) const;
	matrix datacov(unsigned long long flag, bool tao_NULL, const TaoSet &newTao, std::ostream &os) const;

	static void initializeParameters(const matrix &p_mu1, const matrix &p_mu0, const matrix &p_sigma1,
		const matrix &p_sigma0, double p_beta1, double p_beta0, double p_deltar, double p_alpha, 
		double p_shape, double p_bk, double p_ranc, bool bForest);

	inline ClusterTree &getTreeFromID(unsigned long long ID) {
		if(MapClusterTrees.find(ID) == MapClusterTrees.end()) {
			std::ostringstream ostr;
			ostr << "ID " << ID << " not in the list!";
			throw(std::logic_error(ostr.str()));
		}
		return MapClusterTrees.find(ID)->second;
	}

	inline const ClusterTree &getTreeFromID(unsigned long long ID) const {
		if(MapClusterTrees.find(ID) == MapClusterTrees.end()) {
			std::ostringstream ostr;
			ostr << "ID " << ID << " not in the list!";
			throw(std::logic_error(ostr.str()));
		}
		return MapClusterTrees.find(ID)->second;
	}

	//inline ClusterTree &getTreeHead();
	//inline const ClusterTree &getTreeHead() const;

	ClusterTree &createTree(unsigned long long parentID, unsigned long time);
	ClusterTree &createTree();

	void removeTree(unsigned long long TreeID);

	// get the set that can be split, merged, born or killed
	TreeSet getSplitSet(bool hasForest, bool forestOnly) const;
	//TreeSet getTreeSplitSet() const;		// this is to get the new tree (separate process)
	TreeMergeSet getMergeSet(bool hasForest, bool forestOnly) const;
	//TreeMergeSet getTreeMergeSet() const;
	TreeMergeSet getDeathSet(bool hasForest, bool forestOnly) const;

	TreeSet getTailBirthSet(bool hasForest, bool forestOnly) const;		// time actually doesn't matter here
	TreeSet getTailDeathSet(bool hasForest, bool forestOnly) const;		// time and ID

	// the return value is the new tree / branch ID
	unsigned long long flagSplit(double &P_alloc, bool &NULLset, double &f_ui,  
		const TreeSet::value_type &split_pair);

	void flagSplitTest(double &P_alloc, bool &NULLset, double &f_ui,  
		const TreeSet::value_type &split_pair);
	//void flagRootSplit(double &P_alloc, bool &NULLset, double &f_ui, unsigned long long FromID);
	
	void flagMerge(double &P_alloc, double &f_ui, const TreeMergeSet::value_type &merge_pair);
	//void flagRootMerge(double &P_alloc, double &f_ui, const TreeMergeSet::value_type &merge_pair);
	
	void flagBirth(const TreeSet::value_type &splitPair,
		double &weightratio, double &jacobi, double &f_wstar);
	void flagBirth(const TreeSet::value_type &splitPair);

	void flagTailBirth(const TreeSet::value_type &splitPair,
		double &weightratio, double &jacobi, double &f_wstar);
	
	void flagDeath(const TreeMergeSet::value_type &merge_pair, 
		double &jacobi, double &weightratio, double &f_wstar);

	void flagTailDeath(const TreeSet::value_type &merge_pair,
		double &jacobi, double &weightratio, double &f_wstar);

	void testNumberSet();

	long getNoZeroWeights(unsigned long time) const;

	double apSplit(const TreeSet::value_type &splitPair, const bool NULLset, const TreeSet &splitsets, 
		const long splitnum, const double P_alloc, const double f_ui,
		const MCMCEnv &oldenv, unsigned long long newTreeID) const;
	double apMerge(const TreeMergeSet::value_type &mergePair, const TreeSet &splitSetsAfterMerge,
		const long mergesetsrows, const long splitnumaftermerge, const double P_alloc, const double f_ui,
		const MCMCEnv &oldenv) const;
	double apBirth(const TreeSet::value_type &splitPair, const long splitnumafterbirth,
		const double weightratio, const double jacobi, const double f_wstar,
		const MCMCEnv &oldenv) const;
	double apDeath(const TreeMergeSet::value_type &mergePair, const long emptynumbeforedeath,
		const double weightratio, const double jacobi, const double f_wstar, 
		const MCMCEnv &oldenv) const;

	double apTailBirth(const TreeSet::value_type &splitPair, const long splitnumafterbirth,
		const double weightratio, const double jacobi, const double f_wstar,
		const MCMCEnv &oldenv) const;
	double apTailDeath(const TreeSet::value_type &mergePair, const long emptynumbeforedeath,
		const double weightratio, const double jacobi, const double f_wstar, 
		const MCMCEnv &oldenv) const;

	double apRootSplit(unsigned long long FromID, const bool NULLset, const TreeSet &splitsets, 
		const long splitnum, const double P_alloc, const double f_ui,
		const MCMCEnv &oldenv) const;
	double apRootMerge(const TreeMergeSet::value_type &mergePair, const TreeSet &splitSetsAfterMerge,
		const long mergesetsrows, const long splitnumaftermerge, const double P_alloc, const double f_ui,
		const MCMCEnv &oldenv) const;
	double apRootBirth(const long splitnumafterbirth,
		const double weightratio, const double jacobi, const double f_wstar,
		const MCMCEnv &oldenv) const;
	double apRootDeath(const TreeMergeSet::value_type &mergePair, const long emptynumbeforedeath,
		const double weightratio, const double jacobi, const double f_wstar, 
		const MCMCEnv &oldenv) const;

	double calcLogDensity() const;
	double calcLogDensity(const TaoSet &newTao) const;

	long taoCount(const TaoSet &newTao) const;
	long taoCount() const;
	//void changeSigmaToNewTao(const TaoSet &newTao) const;
	double calcLogDensityNonFlagPart(long taonum) const;
	double calcLogDensityFlagPart(long taonum) const;
	double calcLogDensityNonFlagPart(const TaoSet &newTao, long taonum) const;
	double calcLogDensityFlagPart(const TaoSet &newTao, long taonum) const;
	//void returnSigmaFromTao() const;

	double calcPChos(const TreeMergeSet::value_type &parentChildPair, const bool NULLset,
		const TreeSet &SplitSet) const;
	double calcSimilarity(unsigned long time, const ClusterTree &parent, const ClusterTree &child,
		const bool NULLset, const TreeSet &SplitSet) const;
	double calcClusterDist(const unsigned long time1, const ClusterTree &cluster1, 
		const long time2, const ClusterTree &cluster2) const;
	matrix calcClusterMean(const unsigned long time, const ClusterTree &cluster) const;

	unsigned long getClusterNumber() const;		// get the number of trees
	unsigned long getNonEmptyClusterNumber() const;	
	// get the number of trees that are not empty throughout the time
	unsigned long getClusterNumber(unsigned long time) const;	// get the number of trees
	unsigned long getNonEmptyClusterNumber(unsigned long time) const;	
	// get the number of trees that are not empty throughout the time

	double getWeightFromSample(unsigned long time, unsigned long sampleNum) const;

	void sampleWeight();
	void sampleTao();
	void sampleZ();

	double treeSummary() const;
	unsigned long getTreeNumbers() const;

	std::ostream &writeTree(std::ostream &os) const;
	std::ostream &writeTao(std::ostream &os) const;
	std::ostream &writeFlags(std::ostream &os) const;
	std::ostream &writeWeights(std::ostream &os) const;

	std::ostream &writeStatus(std::ostream &os) const;

	double & weight(unsigned long long ID, unsigned long time);

	std::ostream &writeSet(std::ostream &os, const TreeSet &set) const;
	std::ostream &writeSet(std::ostream &os, const TreeMergeSet &set) const;

	friend std::ostream &operator<<(std::ostream &os, const MCMCEnv &env);

};

std::ostream &operator<<(std::ostream &os, const MCMCEnv &env);

#endif