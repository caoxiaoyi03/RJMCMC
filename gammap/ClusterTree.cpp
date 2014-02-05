#include "ClusterTree.h"

#include <iostream>
#include <stdexcept>
#include <cmath>
using namespace std;

//unsigned long long ClusterTree::ID_max = 0;
unsigned long ClusterTree::numTimePoints = 0;

//TreeContainer ClusterTree::MapClusterTree;

bool ClusterTree::setNumOfTimePoints(unsigned long num) {
	if(numTimePoints > 0) {
		cerr << "Time point number already defined, cannot change!" << endl;
		return false;
	}
	numTimePoints = num;
	return true;
}

unsigned long ClusterTree::getNumOfTimePoints() {
	return numTimePoints;
}


ClusterTree::ClusterTree(unsigned long long id)
	: numOfChildren(0), ID(id), legacy_ID(0),
	born_time(0), parentID(0), children(numTimePoints, NULL),
	weights(numTimePoints, 0.0), samples(numTimePoints, 0) {
}

ClusterTree::ClusterTree(ClusterTree &par, unsigned long long id, unsigned long btime)
	: numOfChildren(0), ID(id),
	legacy_ID(0), born_time(btime), parentID(par.ID), children(numTimePoints, NULL),
	weights(numTimePoints, 0.0), samples(numTimePoints, 0) {
}

ClusterTree::ClusterTree(const ClusterTree &old)
	: numOfChildren(old.numOfChildren), ID(old.ID),
	legacy_ID(old.legacy_ID), born_time(old.born_time), parentID(old.parentID), 
	children(old.children), weights(old.weights), samples(old.samples) {
}

ClusterTree::~ClusterTree(void)
{
}

//const ClusterTree &ClusterTree::getChild(unsigned long time) const;
//ClusterTree &ClusterTree::getChild(unsigned long time);

unsigned long long ClusterTree::getChildID(unsigned long time) const {
	return children.at(time);
}

ClusterTree &ClusterTree::setChild(ClusterTree &child, unsigned long time) {
	if(children.at(time)) {
		throw(std::logic_error("Child already born at this time!"));
	}
	children.at(time) = child.ID;
	numOfChildren++;
	return child;
}

unsigned long long ClusterTree::setChild(unsigned long long childID, unsigned long time) {
	if(children.at(time)) {
		throw(std::logic_error("Child already born at this time!"));
	}
	children.at(time) = childID;
	numOfChildren++;
	return childID;
}

unsigned long long ClusterTree::removeChild(unsigned long time) {
	unsigned long long oldChildID = children.at(time);
	if(!oldChildID) {
		throw(std::logic_error("There is no child at this time!"));
	}
	children.at(time) = 0;
	numOfChildren--;
	return oldChildID;
}

unsigned long long ClusterTree::removeChild(unsigned long long childID) {
	for(ChildContainer::iterator itor = children.begin();
		itor != children.end(); itor++) {
			if(*itor == childID) {
				*itor = 0;
				numOfChildren--;
				return childID;
			}
	}
	throw(std::logic_error("There is no child with this ID!"));
	return 0;
}

unsigned long long ClusterTree::removeChild(ClusterTree &child) {
	return removeChild(child.ID);
}

void ClusterTree::fillWeights(double value) {
	cout << weights.size();
	weights.assign(weights.size(), value);
}

void ClusterTree::fillWeights(WeightContainer values) {
	WeightContainer::iterator weight_itor = weights.begin();
	WeightContainer::const_iterator values_itor = values.begin();
	for(; weight_itor != weights.end(); weight_itor++) {
		if(values_itor == values.end()) {
			values_itor = values.begin();
		}
		*weight_itor = *values_itor;
	}
}

std::ostream &ClusterTree::writeWeights(std::ostream &os) const {
	os << "[" << ID << "]";
	for(WeightContainer::const_iterator itor = weights.begin();
		itor != weights.end(); itor++) {
			os << " " << *itor;
	}
	os << endl;
	return os;
}

long ClusterTree::getSampleNum(unsigned long time) const {
	return samples.at(time);
}

long ClusterTree::getSampleNum() const {
	long result = 0;
	for(NumOfSampleContainer::const_iterator itor = samples.begin();
		itor != samples.end(); itor++) {
			result += *itor;
	}
	return result;
}

// Note that, if there is no flag appearing in the cluster list z, there will be a 'out of range' bug appearing.
double ClusterTree::calcKcWeight() const {
	//long kc = weight.getEarliestTime(flag);
	double w=1.0;
	// Note that i is not the real time point, its begginning point should be kc-1.
	for(long i = born_time; i < numTimePoints; i++) {
		w *= pow(weights[i], (int) samples[i]);
	}	// notice that pow(0, 0) = 1 so it should be OK

	return w;
}

unsigned long ClusterTree::getNumOfChildren() const {
	return numOfChildren;
}