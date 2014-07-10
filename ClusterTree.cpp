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
	: numOfChildren(0), born_time(0), death_time(numTimePoints), 
	ID(id), legacy_ID(0),
	parentID(0), children(numTimePoints, 0),
	weights(numTimePoints, 0.0), samples(numTimePoints, 0) {
}

ClusterTree::ClusterTree(ClusterTree &par, unsigned long long id, unsigned long btime)
	: numOfChildren(0), born_time(btime), death_time(par.death_time), ID(id),
	legacy_ID(0), parentID(par.ID), children(numTimePoints, 0),
	weights(numTimePoints, 0.0), samples(numTimePoints, 0) {
}

ClusterTree::ClusterTree(const ClusterTree &old)
	: numOfChildren(old.numOfChildren), born_time(old.born_time), death_time(old.death_time), 
	ID(old.ID),	legacy_ID(old.legacy_ID), parentID(old.parentID), 
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

bool ClusterTree::setDeathTime(unsigned long time) {
	if(time < death_time) {
		// performing an early death
		if(time <= born_time || time < getEarliestDeathableTime()) {
			throw(std::logic_error("Death time is invalid or children and/or samples are present!"));
		} else {
			// valid death time
			death_time = time;
			return true;
		}
	} else {
		// is actually tail-birth
		if(time > numTimePoints || time <= born_time) {
			throw(std::logic_error("Death time is invalid!"));
		}
		death_time = time;
		return true;
	}
}

bool ClusterTree::tailBirth() {
	return setDeathTime(numTimePoints);
}

unsigned long ClusterTree::getEarliestDeathableTime() const {
	unsigned long current_deathable;
	for(unsigned long current_deathable = death_time; current_deathable > born_time; index--) {
		if(children.at(current_deathable - 1) || samples.at(current_deathable - 1)) {
			// there is child or sample, cannot call death
			break;
		} 
	}
	// current deathable index will be index + 1
	// if this value == death_time, then this tree cannot be tail-deathed
	return current_deathable;
}

ClusterTree &ClusterTree::setChild(ClusterTree &child, unsigned long time) {
	if(children.at(time)) {
		throw(std::logic_error("Child already born at this time!"));
	}
	if(death_time <= time) {
		throw(std::logic_error("Child cannot be born to a dead tree!");
	}
	children.at(time) = child.ID;
	numOfChildren++;
	return child;
}

unsigned long long ClusterTree::setChild(unsigned long long childID, unsigned long time) {
	if(children.at(time)) {
		throw(std::logic_error("Child already born at this time!"));
	}
	if(death_time <= time) {
		throw(std::logic_error("Child cannot be born to a dead tree!");
	}
	children.at(time) = childID;
	numOfChildren++;
	return childID;
}

unsigned long long ClusterTree::removeChild(unsigned long time) {
	if(death_time <= time) {
		throw(std::logic_error("Tree is already dead at that time!");
	}
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
	// cout << weights.size();
	weights.assign(weights.size(), 0.0);
	for(unsigned int index = born_time; index < death_time; index++) {
		weights.at(index) = value;
	}
}

void ClusterTree::fillWeights(WeightContainer values) {
	WeightContainer::const_iterator values_itor = values.begin();
	weights.assign(weights.size(), 0.0);
	for(unsigned int index = born_time; index < death_time; index++, values_itor++) {
		if(values_itor == values.end()) {
			values_itor = values.begin();
		}
		weights.at(index) = *values_itor;
	}
}

std::ostream &ClusterTree::writeWeights(std::ostream &os) const {
	os << "[" << ID << "]";
	for(unsigned int index = 0; index < numTimePoints; index++) {
		if(index < born_time || index >= death_time) {
			os << " -";
		} else {
			os << " " << weights.at(index);
		}
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
	for(long i = born_time; i < death_time; i++) {
		w *= pow(weights[i], (int) samples[i]);
	}	// notice that pow(0, 0) = 1 so it should be OK

	return w;
}

unsigned long ClusterTree::getNumOfChildren() const {
	return numOfChildren;
}
