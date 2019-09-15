/*
 * functions.h
 *
 *  Created on: Mar 19, 2014
 *      Author: esakhaee
 */

#ifndef FUNCTIONS_H_
#define FUNCTIONS_H_
#include <vector>
#include <math.h>
#include <algorithm>
#include <iostream>
#include <fstream>
#include <iomanip>
#include <string>

#define NUM_VERTEX 4   //this is for face ambiguity set it as 8 for cube ambiguity

using namespace std;

class convolution {
public:

	convolution(double* means, double* bw, int n, double sampRate);
	virtual ~convolution();
	void findPDF();
	double decider (double isovalue);

private:
	vector <double> distrib1,distrib2,distrib3,distrib4;
	vector <double> binWidth;
	double sampRate;
	vector <double> kern;
	vector <double> pdf;
	double kernelWidth, startPT,endPT;

	int Heaviside (double value);
	void Cubickernel();
	void LinearKernel();
	void degree7Kernel();
	void printVector(vector <double> vect);
	void printVector2File(vector <double> vect);


};

bool MidPointWrapper(double* means, double* binWidths , int sampNum,  double isoValue , double samplingRate);
int readData (string filename, double* means, double* binWidths);

#endif /* FUNCTIONS_H_ */


