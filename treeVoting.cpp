#include "string_manip_stl.h"
#include <vector>
#include <string>
#include <ifstream>
#include <algorithm>

using namespace std;

class Tree {
	unsigned long depth;
	unsigned long ID;
	vector<Tree> child;

public:
	Tree(unsigned long, unsigned long);
	Tree(const Tree& tree);

	int compare(const Tree& tree) const;
	bool operator<(const Tree& tree) const;
	bool operator==(const Tree& tree) const;
	bool operator>(const Tree& tree) const;
	bool operator!=(const Tree& tree) const;
	bool operator<=(const Tree& tree) const;
	bool operator>=(const Tree& tree) const;

	bool branch(unsigned long time, unsigned long ID);
	bool branchChild(unsigned long parentID, unsigned long time, unsigned long ID);
	bool removeChild(unsigned long);

}

Tree::Tree(unsigned long id, unsigned long newdepth): ID(id), depth(newdepth) {
	if(newdepth > 0) {
		child.push_back(Tree(id, newdepth - 1));
	}
}

Tree::Tree(const Tree& tree): ID(tree.ID), depth(tree.depth) {
	for(vector<Tree>::const_iterator itor = tree.child.begin();
		itor != tree.child.end(); itor++) {
			child.push_back(*itor);
	}
}

int Tree::compare(const Tree& tree) const {
	// return -1 if *this < tree
	// return 0 if *this == tree
	// return 1 if *this > tree
	if(child.size() != tree.child.size()) {
		return (child.size() < tree.child.size())? -1: 1;
	} else {
		if(depth <= 1) {
			return 0;	// should be equal
		} else {
			// check from the smallest child
			for(vector<Tree>::const_iterator itor = child.begin(), targetItor = tree.child.begin();
				itor != child.end() && targetItor != child.end(); itor++, targetItor++) {
					int result = itor->compare(*targetItor);
					if(result) {
						return result;
					} 
			}
			return 0;
		}
	}
}

bool Tree::operator==(const Tree& tree) const {
	return !(this->compare(tree));
}

bool Tree::operator<(const Tree& tree) const {
	return (this->compare(tree) == -1);
}

bool Tree::operator>(const Tree& tree) const {
	return (this->compare(tree) == 1);
}

bool Tree::operator!=(const Tree& tree) const {
	return !(*this == tree);
}

bool Tree::operator<=(const Tree& tree) const {
	return !(*this > tree);
}

bool Tree::operator>=(const Tree& tree) const {
	return !(*this < tree);
}

bool Tree::branch(unsigned long childdepth, unsigned long id) {
	// default parentID = ID
	return branch(ID, childdepth, id);
}

bool Tree::branchChild(unsigned long parentID, unsigned long childdepth, unsigned long id) {
	if(childdepth == depth - 1) {
		if(child.size() < 2 && ID == parentID) {
			// should branch here and can branch here
			child.push_back(Tree(id, childdepth));
			sort(child.begin(), child.end());
			return true;
		} else {
			return false;
		}
	} else {
		for(vector<Tree>::iterator itor = child.begin(); itor != child.end(); itor++) {
			if(itor->branchChild(parentID, childdepth, id)) {
				sort(child.begin(), child.end());
				return true;
			}
		}
		return false;
	}	
}

bool Tree::removeChild(unsigned long id) {
	for(vector<Tree>::iterator itor = child.begin(); itor != child.end(); itor++) {
		if(itor->ID == id) {
			child.erase(itor);
			sort(child.begin(), child.end());
			return true;
		}
		if(itor->removeChild(id)) {
			sort(child.begin(), child.end());
			return true;
		}
	}
	return false;
}
