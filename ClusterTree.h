#ifndef CLUSTERTREE_H
#define CLUSTERTREE_H
#pragma once

#include <map>
#include <vector>
#include <iostream>

class ClusterTree;

// Notice that ClusterTree is a component of MCMCEnv (hence its construction is private)

class ClusterTree
{
	typedef std::vector<unsigned long long> ChildContainer;
	// note that the index (unsigned long) is the time point for the child to be born
	// can be NULL if no child is present

	typedef std::vector<double> WeightContainer;
	// index is the time point, value is the weights;

	typedef std::vector<unsigned long> NumOfSampleContainer;
	// index is the time point again

	//static unsigned long long ID_max;
	static unsigned long numTimePoints;

	ClusterTree(unsigned long long id);
	ClusterTree(ClusterTree &parent, unsigned long long id, unsigned long btime);

	friend class MCMCEnv;

	unsigned long numOfChildren;

public:
	//static TreeContainer MapClusterTree;
	ClusterTree(const ClusterTree &old);

	unsigned long long ID;
	unsigned long legacy_ID;		// used to ensure the tree term is the same as legacy

	unsigned long born_time;
	// the first time point this tree is born (appears)
	// therefore, only the root node has a born time of 0, others should be at least 1

	unsigned long long parentID;
	ChildContainer children;
	WeightContainer weights;
	NumOfSampleContainer samples;

	static bool setNumOfTimePoints(unsigned long num);
	static unsigned long getNumOfTimePoints();

	//const ClusterTree &getChild(unsigned long time) const;
	//ClusterTree &getChild(unsigned long time);

	unsigned long long getChildID(unsigned long time) const;

	ClusterTree &setChild(ClusterTree &child, unsigned long time);
	unsigned long long setChild(unsigned long long childID, unsigned long time);

	unsigned long long removeChild(unsigned long time);
	unsigned long long removeChild(unsigned long long childID);
	unsigned long long removeChild(ClusterTree &child);

	void fillWeights(double value);
	void fillWeights(WeightContainer values);
	// this is a warpped weight container (if values is shorter than weights)
	std::ostream &writeWeights(std::ostream &os) const;

	//double getWeight(unsigned long time);
	//void setWeight(unsigned long value);

	long getSampleNum(unsigned long time) const;
	long getSampleNum() const;

	double calcKcWeight() const;

	unsigned long getNumOfChildren() const;

	~ClusterTree(void);

};

#endif