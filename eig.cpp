/*
 * eig.cpp
 *
 *  Created on: Jul 18, 2016
 *      Author: j4yan
 */

#include <iostream>
#include <sstream>
#include <fstream>
#include <string>
#include <cmath>

using namespace std;

template<typename T>
void BubbleSort(T* a, unsigned int n);

int main(int argc, char **argv) {
	ifstream fin;
  string fname = "eigs_real1.dat";
	// string fname = argv[1];
	stringstream ss;
	string line;
	int nLines = 0;
	fin.open(fname.c_str());
  if (!fin.is_open()) {
    return 0;
  }
	while (std::getline(fin, line)) {
		++nLines;
	}

	cout << "file-" << fname << ", nLines = " << nLines <<endl;
	fin.clear();
	fin.seekg(0, ios::beg);
	double* eig = new double[nLines];
	for (int i = 0; i < nLines; ++i) {
		getline(fin, line);
		ss.clear();
		ss.str(line);
		ss >> eig[i];
	}
	BubbleSort(eig, nLines);
	ofstream fout;
	ofstream fcn;
	fout.open("eig.curve");
	fcn.open("cn.dat");
	for (int i = 0; i < nLines; ++i) {
		fout << i << ", " << eig[i] << endl;
	}
	
	double eig_min = 1e10;
	double eig_max = 1e-10;
	for(int i=0; i<nLines; i++){
		eig_min = min(eig_min, fabs(eig[i]));
		eig_max = max(eig_max, fabs(eig[i]));
	}
	cout << "eig_min = " << eig[0] << ", eig_max = " << eig[nLines-1] <<endl;
	double cn = eig_max/eig_min;

	if (eig[0]*eig[nLines-1] < 0.0){
		cn *= -1.0;
	}

	fcn << eig[nLines-1] << " " << eig[0] << " " << cn << endl;
	cout << "cn = " << cn << endl;
	delete[] eig;

	fout.close();
	fcn.close();
	return 0;
}

template<typename T>
void BubbleSort(T* a, unsigned int n) {
	unsigned int flag = 1;	// set flag to 1 to start first pass
	T tmp;
	for (unsigned int i = 0; (i < n) && flag; i++) {
		flag = 0;
		for (unsigned int j = 0; j < (n - 1); j++) {
			if (a[j + 1] > a[j]) {
				tmp = a[j];
				a[j] = a[j + 1];
				a[j + 1] = tmp;
				flag = 1;
			}
		}
	}
	return;
}

