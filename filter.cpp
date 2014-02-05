#include <iostream>
#include <sstream>
#include <fstream>
#include <set>
#include <string>
#include <cctype>
#include <vector>

using namespace std;

void executeFilter(unsigned long begin, unsigned long interval, const set<char> &cmd, const vector<string> &fileNames) {
	for(vector<string>::const_iterator itor = fileNames.begin(); itor != fileNames.end(); itor++) {
		ifstream fin(itor->c_str());
		if(!fin) {
			cerr << "ERROR: Cannot open file \"" << *itor << "\"." << endl;
			continue;
		} else {
			cout << *itor << endl << endl;
		}
	
		string line;
		while(getline(fin, line)) {
			if(line.find("STEP:") != string::npos) {
				// this is step info;
				string step;
				unsigned long currStep;
				istringstream istr(line);
				istr >> step >> currStep;
				if(currStep > begin && !((currStep - begin) % interval)) {
					// this is the step we wanted
					cout << line << endl << endl;
				
					while(getline(fin, line) && line.find("BirthDeath") == string::npos);
					if(cmd.find('w') != cmd.end()) {
						// output weight
						while(getline(fin, line) && line.find("Cluster proportions") == string::npos);
						if(fin) 
							cout << line << endl;
						while(getline(fin, line) && !line.empty()) {
							if(fin) 
								cout << line << endl;
						}
						cout << endl;
					}
					if(cmd.find('f') != cmd.end()) {
						// output feature
						while(getline(fin, line) && line.find("Feature") == string::npos);
						if(fin) 
							cout << line << endl;
						while(getline(fin, line) && !line.empty()) {
							if(fin) 
								cout << line << endl;
						}
						cout << endl;
					}
					if(cmd.find('c') != cmd.end()) {
						// output cluster index
						while(getline(fin, line) && line.find("Cluster index") == string::npos);
							if(fin) 
								cout << line << endl;
						while(getline(fin, line) && !line.empty()) {
							if(fin) 
								cout << line << endl;
						}
						cout << endl;
					}
					if(cmd.find('t') != cmd.end()) {
						// output tree structure
						while(getline(fin, line) && line.find("Tree:") == string::npos);
							if(fin) 
								cout << line << endl;
						while(getline(fin, line) && !line.empty()) {
							if(fin) 
								cout << line << endl;
						}
						cout << endl;
					}
				}
			}
		}
	}
}


int main(int argc, char * argv[]) {
	// usage ./filter [-wfct] -s <start> -i <interval> <output_name>
	// w: weight; f: feature; c: cluster index; t: tree structure
	if(argc < 2) {
		cout << "usage: ./filter [-wfct] <-s start> <-i interval> <output_name>" << endl
			<< "w: weight; f: feature; c: cluster index; t: tree structure" << endl;
	} else {
		set<char> commands;
		vector<string> fileNames;
		unsigned long start, interval;
		for(int i = 1; i < argc; i++) {
			string arg(argv[i]);
			if(arg[0] == '-') {
				// This is argument
				if(arg.length() < 2) {
					continue;
				}
				if(arg[1] == 's') {
					i++;
					if(i < argc) {
						string nextArg(argv[i]);
						istringstream istr(nextArg);
						istr >> start;
					}
				} else if(arg[1] == 'i') {
					i++;
					if(i < argc) {
						string nextArg(argv[i]);
						istringstream istr(nextArg);
						istr >> interval;
					}
				} else {
					for(int j = 1; j < arg.length(); j++) {
						commands.insert(tolower(arg[j]));
					}
				}
			} else {
				// it's filename
				fileNames.push_back(arg);
			}
		}
		cout << commands.size();
		executeFilter(start, interval, commands, fileNames);
	}

}