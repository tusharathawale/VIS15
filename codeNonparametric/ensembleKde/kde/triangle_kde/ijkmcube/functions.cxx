/*
 * functions.cpp
 *
 *  Created on: Mar 19, 2014
 *      Author: esakhaee
 */

#include "functions.h"

convolution::convolution(double* means, double* bw, int n, double sr){
        int i=0; int q=0;
        for (q=0 ; q<n ;q++){
                distrib1.push_back(means[q]/NUM_VERTEX);
        }
        binWidth.push_back(bw[i++]/NUM_VERTEX);
        for (q=0 ; q<n ;q++){
                distrib2.push_back(means[i*n+q]/NUM_VERTEX);
        }
        binWidth.push_back(bw[i++]/NUM_VERTEX);
        for (q=0 ; q<n ;q++){
                distrib3.push_back(means[i*n+q]/NUM_VERTEX);
        }
        binWidth.push_back(bw[i++]/NUM_VERTEX);
        for (q=0 ; q<n ;q++){
                distrib4.push_back(means[i*n+q]/NUM_VERTEX);
        }
        binWidth.push_back(bw[i++]/NUM_VERTEX);

        sampRate = sr;

}


convolution::~convolution() {
	//vectors are automatically cleaned as they go out of scope :)
}

bool MidPointWrapper(double* means, double* binwidth,int sampNum,double isoValue,double sampRate){

	convolution C (means, binwidth, sampNum, sampRate);
    C.findPDF();
    //cout << "distrib1[1] :" << C.distrib1[1]<<endl;
	double prob = C.decider(isoValue);
	cout << "Pr(X>" << isoValue << ")=" <<  prob << endl;
	return (prob>=0.5);
}

//===================================================================================
//				DIFFERENT KERNELS USING FINITE DIFFERENCING IN MATLAB
//=====================================================================================
void convolution::LinearKernel(){
	double h1 = binWidth[0];
	vector <double>x;

	for (double i=-h1; i<=h1; i+=sampRate){
		x.push_back(i);
	}

	double p;
	for (int j=0; j<(int)x.size(); j++){
		p = (Heaviside(h1 + x[j])*(h1 + x[j]) - x[j]*Heaviside(x[j]))/(h1*h1) - (Heaviside(x[j] - h1)*(h1 - x[j]) + x[j]*Heaviside(x[j]))/(h1*h1);

		if (p<1e-6) p=0;
		this->kern.push_back(p);
	}
}


void convolution::Cubickernel(){
	double h1 = binWidth[0];
	double h2 = binWidth[1];

	vector <double>x;

	for (double i=-(h1+h2); i<=(h1+h2); i+=sampRate){
		x.push_back(i);
	}

	double p;
	for (int j=0; j<(int)x.size(); j++){
		p =((double)((((Heaviside(h1 - h2 + x[j]) * pow(h1 - h2 + x[j],3.0)
				+ Heaviside(x[j] - h2) * pow(h2 - x[j],3.0)) / h1 +
				(Heaviside(x[j] - h2) * pow(h2 - x[j],3.0) - Heaviside(x[j] - h2 - h1)
						* pow(h1 + h2 - x[j],3.0)) / h1) / h1 /6 + ((pow(x[j],3.0)
								* Heaviside(x[j]) + Heaviside(x[j] - h1) * pow(h1 - x[j],3.0)) / h1 +
								(pow(x[j],3.0) * Heaviside(x[j]) -
										Heaviside(h1 + x[j]) * pow(h1 + x[j],3.0)) / h1) / h1 /6)*(1/(h2*h2))
										+ (((Heaviside(h2 - h1 + x[j]) *
												pow(h2 - h1 + x[j],3.0) - Heaviside(h2 + x[j]) *
												pow(h2 + x[j],3.0)) / h1 + (Heaviside(h1 + h2 + x[j]) *
														pow(h1 + h2 + x[j],3.0) - Heaviside(h2 + x[j]) *
														pow(h2 + x[j],3.0)) / h1) / h1 /6 +
														((pow(x[j],3.0) * Heaviside(x[j]) + Heaviside(x[j] - h1) *
																pow(h1 - x[j],3.0)) / h1 + (pow(x[j],3.0) * Heaviside(x[j]) -
																		Heaviside(h1 + x[j]) * pow(h1 + x[j],3.0)) / h1) / h1 /6) *
																		(1/(h2*h2)) ));


		if (p<1e-6) p=0;
		//p = p/(h1*h2);    // no need kernel is normalized by its own
		this->kern.push_back(p);
	}

}



void convolution::degree7Kernel(){
	double h1 = binWidth[0];
	double h2 = binWidth[1];
	double h3 = binWidth[2];
	double h4 = binWidth[3];

	vector <double> x;

	for (double i=-(h1+h2+h3+h4); i<=(h1+h2+h3+h4); i+=sampRate){
		x.push_back(i);
	}

	double p;
	for (int j=0; j<(int)x.size(); j++){

		p= ((double)((((((((Heaviside(h1 - h2 + x[j]) * pow(h1 - h2 + x[j], 7) + Heaviside(x[j] - h2) * pow(h2 - x[j], 7)) / h1 + (Heaviside(x[j] - h2) * pow(h2 - x[j], 7) - Heaviside(x[j] - h2 - h1) * pow(h1 + h2 - x[j], 7)) / h1) / h1 + ((pow(x[j], 7) * Heaviside(x[j]) + Heaviside(x[j] - h1) * pow(h1 - x[j], 7)) / h1 + (pow(x[j], 7) * Heaviside(x[j]) - Heaviside(h1 + x[j]) * pow(h1 + x[j], 7)) / h1) / h1) / h2 + (((Heaviside(h2 - h1 + x[j]) * pow(h2 - h1 + x[j], 7) - Heaviside(h2 + x[j]) * pow(h2 + x[j], 7)) / h1 + (Heaviside(h1 + h2 + x[j]) * pow(h1 + h2 + x[j], 7) - Heaviside(h2 + x[j]) * pow(h2 + x[j], 7)) / h1) / h1 + ((pow(x[j], 7) * Heaviside(x[j]) + Heaviside(x[j] - h1) * pow(h1 - x[j], 7)) / h1 + (pow(x[j], 7) * Heaviside(x[j]) - Heaviside(h1 + x[j]) * pow(h1 + x[j], 7)) / h1) / h1) / h2) / h2 + ((((Heaviside(h1 - h3 + x[j]) * pow(h1 - h3 + x[j], 7) + Heaviside(x[j] - h3) * pow(h3 - x[j], 7)) / h1 + (Heaviside(x[j] - h3) * pow(h3 - x[j], 7) - Heaviside(x[j] - h3 - h1) * pow(h1 + h3 - x[j], 7)) / h1) / h1 + ((Heaviside(h2 - h3 + x[j]) * pow(h2 - h3 + x[j], 7) + Heaviside(h2 - h1 - h3 + x[j]) * pow(h1 - h2 + h3 - x[j], 7)) / h1 - (Heaviside(h1 + h2 - h3 + x[j]) * pow(h1 + h2 - h3 + x[j], 7) - Heaviside(h2 - h3 + x[j]) * pow(h2 - h3 + x[j], 7)) / h1) / h1) / h2 + (((Heaviside(h1 - h3 + x[j]) * pow(h1 - h3 + x[j], 7) + Heaviside(x[j] - h3) * pow(h3 - x[j], 7)) / h1 + (Heaviside(x[j] - h3) * pow(h3 - x[j], 7) - Heaviside(x[j] - h3 - h1) * pow(h1 + h3 - x[j], 7)) / h1) / h1 - ((Heaviside(x[j] - h3 - h2) * pow(h2 + h3 - x[j], 7) + Heaviside(h1 - h2 - h3 + x[j]) * pow(h1 - h2 - h3 + x[j], 7)) / h1 + (Heaviside(x[j] - h3 - h2) * pow(h2 + h3 - x[j], 7) - Heaviside(x[j] - h2 - h3 - h1) * pow(h1 + h2 + h3 - x[j], 7)) / h1) / h1) / h2) / h2) / h3 + (((((Heaviside(h3 - h1 + x[j]) * pow(h3 - h1 + x[j], 7) - Heaviside(h3 + x[j]) * pow(h3 + x[j], 7)) / h1 + (Heaviside(h1 + h3 + x[j]) * pow(h1 + h3 + x[j], 7) - Heaviside(h3 + x[j]) * pow(h3 + x[j], 7)) / h1) / h1 + ((Heaviside(h3 - h2 + x[j]) * pow(h3 - h2 + x[j], 7) + Heaviside(h3 - h2 - h1 + x[j]) * pow(h1 + h2 - h3 - x[j], 7)) / h1 - (Heaviside(h1 - h2 + h3 + x[j]) * pow(h1 - h2 + h3 + x[j], 7) - Heaviside(h3 - h2 + x[j]) * pow(h3 - h2 + x[j], 7)) / h1) / h1) / h2 + (((Heaviside(h3 - h1 + x[j]) * pow(h3 - h1 + x[j], 7) - Heaviside(h3 + x[j]) * pow(h3 + x[j], 7)) / h1 + (Heaviside(h1 + h3 + x[j]) * pow(h1 + h3 + x[j], 7) - Heaviside(h3 + x[j]) * pow(h3 + x[j], 7)) / h1) / h1 - ((Heaviside(h2 - h1 + h3 + x[j]) * pow(h2 - h1 + h3 + x[j], 7) - Heaviside(h2 + h3 + x[j]) * pow(h2 + h3 + x[j], 7)) / h1 + (Heaviside(h1 + h2 + h3 + x[j]) * pow(h1 + h2 + h3 + x[j], 7) - Heaviside(h2 + h3 + x[j]) * pow(h2 + h3 + x[j], 7)) / h1) / h1) / h2) / h2 + ((((Heaviside(h1 - h2 + x[j]) * pow(h1 - h2 + x[j], 7) + Heaviside(x[j] - h2) * pow(h2 - x[j], 7)) / h1 + (Heaviside(x[j] - h2) * pow(h2 - x[j], 7) - Heaviside(x[j] - h2 - h1) * pow(h1 + h2 - x[j], 7)) / h1) / h1 + ((pow(x[j], 7) * Heaviside(x[j]) + Heaviside(x[j] - h1) * pow(h1 - x[j], 7)) / h1 + (pow(x[j], 7) * Heaviside(x[j]) - Heaviside(h1 + x[j]) * pow(h1 + x[j], 7)) / h1) / h1) / h2 + (((Heaviside(h2 - h1 + x[j]) * pow(h2 - h1 + x[j], 7) - Heaviside(h2 + x[j]) * pow(h2 + x[j], 7)) / h1 + (Heaviside(h1 + h2 + x[j]) * pow(h1 + h2 + x[j], 7) - Heaviside(h2 + x[j]) * pow(h2 + x[j], 7)) / h1) / h1 + ((pow(x[j], 7) * Heaviside(x[j]) + Heaviside(x[j] - h1) * pow(h1 - x[j], 7)) / h1 + (pow(x[j], 7) * Heaviside(x[j]) - Heaviside(h1 + x[j]) * pow(h1 + x[j], 7)) / h1) / h1) / h2) / h2) / h3) / h3 / 5040 - ((((((Heaviside(x[j] - h4 - h3) * pow(h3 + h4 - x[j], 7) + Heaviside(h1 - h3 - h4 + x[j]) * pow(h1 - h3 - h4 + x[j], 7)) / h1 + (Heaviside(x[j] - h4 - h3) * pow(h3 + h4 - x[j], 7) - Heaviside(x[j] - h3 - h4 - h1) * pow(h1 + h3 + h4 - x[j], 7)) / h1) / h1 + ((Heaviside(h2 - h3 - h4 + x[j]) * pow(h2 - h3 - h4 + x[j], 7) + Heaviside(h2 - h1 - h3 - h4 + x[j]) * pow(h1 - h2 + h3 + h4 - x[j], 7)) / h1 - (Heaviside(h1 + h2 - h3 - h4 + x[j]) * pow(h1 + h2 - h3 - h4 + x[j], 7) - Heaviside(h2 - h3 - h4 + x[j]) * pow(h2 - h3 - h4 + x[j], 7)) / h1) / h1) / h2 - (((Heaviside(x[j] - h3 - h4 - h2) * pow(h2 + h3 + h4 - x[j], 7) - Heaviside(x[j] - h2 - h3 - h4 - h1) * pow(h1 + h2 + h3 + h4 - x[j], 7)) / h1 + (Heaviside(x[j] - h3 - h4 - h2) * pow(h2 + h3 + h4 - x[j], 7) - Heaviside(h1 - h2 - h3 - h4 + x[j]) * pow(h2 - h1 + h3 + h4 - x[j], 7)) / h1) / h1 - ((Heaviside(x[j] - h4 - h3) * pow(h3 + h4 - x[j], 7) + Heaviside(h1 - h3 - h4 + x[j]) * pow(h1 - h3 - h4 + x[j], 7)) / h1 + (Heaviside(x[j] - h4 - h3) * pow(h3 + h4 - x[j], 7) - Heaviside(x[j] - h3 - h4 - h1) * pow(h1 + h3 + h4 - x[j], 7)) / h1) / h1) / h2) / h2 - ((((Heaviside(h1 - h4 + x[j]) * pow(h1 - h4 + x[j], 7) + Heaviside(x[j] - h4) * pow(h4 - x[j], 7)) / h1 + (Heaviside(x[j] - h4) * pow(h4 - x[j], 7) - Heaviside(x[j] - h4 - h1) * pow(h1 + h4 - x[j], 7)) / h1) / h1 + ((Heaviside(h2 - h4 + x[j]) * pow(h2 - h4 + x[j], 7) + Heaviside(h2 - h1 - h4 + x[j]) * pow(h1 - h2 + h4 - x[j], 7)) / h1 - (Heaviside(h1 + h2 - h4 + x[j]) * pow(h1 + h2 - h4 + x[j], 7) - Heaviside(h2 - h4 + x[j]) * pow(h2 - h4 + x[j], 7)) / h1) / h1) / h2 + (((Heaviside(h1 - h4 + x[j]) * pow(h1 - h4 + x[j], 7) + Heaviside(x[j] - h4) * pow(h4 - x[j], 7)) / h1 + (Heaviside(x[j] - h4) * pow(h4 - x[j], 7) - Heaviside(x[j] - h4 - h1) * pow(h1 + h4 - x[j], 7)) / h1) / h1 - ((Heaviside(x[j] - h4 - h2) * pow(h2 + h4 - x[j], 7) + Heaviside(h1 - h2 - h4 + x[j]) * pow(h1 - h2 - h4 + x[j], 7)) / h1 + (Heaviside(x[j] - h4 - h2) * pow(h2 + h4 - x[j], 7) - Heaviside(x[j] - h2 - h4 - h1) * pow(h1 + h2 + h4 - x[j], 7)) / h1) / h1) / h2) / h2) / h3 + (((((Heaviside(h2 + h3 - h4 + x[j]) * pow(h2 + h3 - h4 + x[j], 7) - Heaviside(h2 - h1 + h3 - h4 + x[j]) * pow(h2 - h1 + h3 - h4 + x[j], 7)) / h1 + (Heaviside(h2 + h3 - h4 + x[j]) * pow(h2 + h3 - h4 + x[j], 7) - Heaviside(h1 + h2 + h3 - h4 + x[j]) * pow(h1 + h2 + h3 - h4 + x[j], 7)) / h1) / h1 - ((Heaviside(h3 - h4 + x[j]) * pow(h3 - h4 + x[j], 7) + Heaviside(h3 - h1 - h4 + x[j]) * pow(h1 - h3 + h4 - x[j], 7)) / h1 - (Heaviside(h1 + h3 - h4 + x[j]) * pow(h1 + h3 - h4 + x[j], 7) - Heaviside(h3 - h4 + x[j]) * pow(h3 - h4 + x[j], 7)) / h1) / h1) / h2 - (((Heaviside(h3 - h2 - h4 + x[j]) * pow(h2 - h3 + h4 - x[j], 7) - Heaviside(h3 - h2 - h1 - h4 + x[j]) * pow(h1 + h2 - h3 + h4 - x[j], 7)) / h1 + (Heaviside(h1 - h2 + h3 - h4 + x[j]) * pow(h1 - h2 + h3 - h4 + x[j], 7) + Heaviside(h3 - h2 - h4 + x[j]) * pow(h2 - h3 + h4 - x[j], 7)) / h1) / h1 + ((Heaviside(h3 - h4 + x[j]) * pow(h3 - h4 + x[j], 7) + Heaviside(h3 - h1 - h4 + x[j]) * pow(h1 - h3 + h4 - x[j], 7)) / h1 - (Heaviside(h1 + h3 - h4 + x[j]) * pow(h1 + h3 - h4 + x[j], 7) - Heaviside(h3 - h4 + x[j]) * pow(h3 - h4 + x[j], 7)) / h1) / h1) / h2) / h2 - ((((Heaviside(h1 - h4 + x[j]) * pow(h1 - h4 + x[j], 7) + Heaviside(x[j] - h4) * pow(h4 - x[j], 7)) / h1 + (Heaviside(x[j] - h4) * pow(h4 - x[j], 7) - Heaviside(x[j] - h4 - h1) * pow(h1 + h4 - x[j], 7)) / h1) / h1 + ((Heaviside(h2 - h4 + x[j]) * pow(h2 - h4 + x[j], 7) + Heaviside(h2 - h1 - h4 + x[j]) * pow(h1 - h2 + h4 - x[j], 7)) / h1 - (Heaviside(h1 + h2 - h4 + x[j]) * pow(h1 + h2 - h4 + x[j], 7) - Heaviside(h2 - h4 + x[j]) * pow(h2 - h4 + x[j], 7)) / h1) / h1) / h2 + (((Heaviside(h1 - h4 + x[j]) * pow(h1 - h4 + x[j], 7) + Heaviside(x[j] - h4) * pow(h4 - x[j], 7)) / h1 + (Heaviside(x[j] - h4) * pow(h4 - x[j], 7) - Heaviside(x[j] - h4 - h1) * pow(h1 + h4 - x[j], 7)) / h1) / h1 - ((Heaviside(x[j] - h4 - h2) * pow(h2 + h4 - x[j], 7) + Heaviside(h1 - h2 - h4 + x[j]) * pow(h1 - h2 - h4 + x[j], 7)) / h1 + (Heaviside(x[j] - h4 - h2) * pow(h2 + h4 - x[j], 7) - Heaviside(x[j] - h2 - h4 - h1) * pow(h1 + h2 + h4 - x[j], 7)) / h1) / h1) / h2) / h2) / h3) / h3 / 5040) *(1/(h4*h4)) + (((((((Heaviside(h1 - h2 + x[j]) * pow(h1 - h2 + x[j], 7) + Heaviside(x[j] - h2) * pow(h2 - x[j], 7)) / h1 + (Heaviside(x[j] - h2) * pow(h2 - x[j], 7) - Heaviside(x[j] - h2 - h1) * pow(h1 + h2 - x[j], 7)) / h1) / h1 + ((pow(x[j], 7) * Heaviside(x[j]) + Heaviside(x[j] - h1) * pow(h1 - x[j], 7)) / h1 + (pow(x[j], 7) * Heaviside(x[j]) - Heaviside(h1 + x[j]) * pow(h1 + x[j], 7)) / h1) / h1) / h2 + (((Heaviside(h2 - h1 + x[j]) * pow(h2 - h1 + x[j], 7) - Heaviside(h2 + x[j]) * pow(h2 + x[j], 7)) / h1 + (Heaviside(h1 + h2 + x[j]) * pow(h1 + h2 + x[j], 7) - Heaviside(h2 + x[j]) * pow(h2 + x[j], 7)) / h1) / h1 + ((pow(x[j], 7) * Heaviside(x[j]) + Heaviside(x[j] - h1) * pow(h1 - x[j], 7)) / h1 + (pow(x[j], 7) * Heaviside(x[j]) - Heaviside(h1 + x[j]) * pow(h1 + x[j], 7)) / h1) / h1) / h2) / h2 + ((((Heaviside(h1 - h3 + x[j]) * pow(h1 - h3 + x[j], 7) + Heaviside(x[j] - h3) * pow(h3 - x[j], 7)) / h1 + (Heaviside(x[j] - h3) * pow(h3 - x[j], 7) - Heaviside(x[j] - h3 - h1) * pow(h1 + h3 - x[j], 7)) / h1) / h1 + ((Heaviside(h2 - h3 + x[j]) * pow(h2 - h3 + x[j], 7) + Heaviside(h2 - h1 - h3 + x[j]) * pow(h1 - h2 + h3 - x[j], 7)) / h1 - (Heaviside(h1 + h2 - h3 + x[j]) * pow(h1 + h2 - h3 + x[j], 7) - Heaviside(h2 - h3 + x[j]) * pow(h2 - h3 + x[j], 7)) / h1) / h1) / h2 + (((Heaviside(h1 - h3 + x[j]) * pow(h1 - h3 + x[j], 7) + Heaviside(x[j] - h3) * pow(h3 - x[j], 7)) / h1 + (Heaviside(x[j] - h3) * pow(h3 - x[j], 7) - Heaviside(x[j] - h3 - h1) * pow(h1 + h3 - x[j], 7)) / h1) / h1 - ((Heaviside(x[j] - h3 - h2) * pow(h2 + h3 - x[j], 7) + Heaviside(h1 - h2 - h3 + x[j]) * pow(h1 - h2 - h3 + x[j], 7)) / h1 + (Heaviside(x[j] - h3 - h2) * pow(h2 + h3 - x[j], 7) - Heaviside(x[j] - h2 - h3 - h1) * pow(h1 + h2 + h3 - x[j], 7)) / h1) / h1) / h2) / h2) / h3 + (((((Heaviside(h3 - h1 + x[j]) * pow(h3 - h1 + x[j], 7) - Heaviside(h3 + x[j]) * pow(h3 + x[j], 7)) / h1 + (Heaviside(h1 + h3 + x[j]) * pow(h1 + h3 + x[j], 7) - Heaviside(h3 + x[j]) * pow(h3 + x[j], 7)) / h1) / h1 + ((Heaviside(h3 - h2 + x[j]) * pow(h3 - h2 + x[j], 7) + Heaviside(h3 - h2 - h1 + x[j]) * pow(h1 + h2 - h3 - x[j], 7)) / h1 - (Heaviside(h1 - h2 + h3 + x[j]) * pow(h1 - h2 + h3 + x[j], 7) - Heaviside(h3 - h2 + x[j]) * pow(h3 - h2 + x[j], 7)) / h1) / h1) / h2 + (((Heaviside(h3 - h1 + x[j]) * pow(h3 - h1 + x[j], 7) - Heaviside(h3 + x[j]) * pow(h3 + x[j], 7)) / h1 + (Heaviside(h1 + h3 + x[j]) * pow(h1 + h3 + x[j], 7) - Heaviside(h3 + x[j]) * pow(h3 + x[j], 7)) / h1) / h1 - ((Heaviside(h2 - h1 + h3 + x[j]) * pow(h2 - h1 + h3 + x[j], 7) - Heaviside(h2 + h3 + x[j]) * pow(h2 + h3 + x[j], 7)) / h1 + (Heaviside(h1 + h2 + h3 + x[j]) * pow(h1 + h2 + h3 + x[j], 7) - Heaviside(h2 + h3 + x[j]) * pow(h2 + h3 + x[j], 7)) / h1) / h1) / h2) / h2 + ((((Heaviside(h1 - h2 + x[j]) * pow(h1 - h2 + x[j], 7) + Heaviside(x[j] - h2) * pow(h2 - x[j], 7)) / h1 + (Heaviside(x[j] - h2) * pow(h2 - x[j], 7) - Heaviside(x[j] - h2 - h1) * pow(h1 + h2 - x[j], 7)) / h1) / h1 + ((pow(x[j], 7) * Heaviside(x[j]) + Heaviside(x[j] - h1) * pow(h1 - x[j], 7)) / h1 + (pow(x[j], 7) * Heaviside(x[j]) - Heaviside(h1 + x[j]) * pow(h1 + x[j], 7)) / h1) / h1) / h2 + (((Heaviside(h2 - h1 + x[j]) * pow(h2 - h1 + x[j], 7) - Heaviside(h2 + x[j]) * pow(h2 + x[j], 7)) / h1 + (Heaviside(h1 + h2 + x[j]) * pow(h1 + h2 + x[j], 7) - Heaviside(h2 + x[j]) * pow(h2 + x[j], 7)) / h1) / h1 + ((pow(x[j], 7) * Heaviside(x[j]) + Heaviside(x[j] - h1) * pow(h1 - x[j], 7)) / h1 + (pow(x[j], 7) * Heaviside(x[j]) - Heaviside(h1 + x[j]) * pow(h1 + x[j], 7)) / h1) / h1) / h2) / h2) / h3) / h3 / 5040 + ((((((Heaviside(h4 - h1 + x[j]) * pow(h4 - h1 + x[j], 7) - Heaviside(h4 + x[j]) * pow(h4 + x[j], 7)) / h1 + (Heaviside(h1 + h4 + x[j]) * pow(h1 + h4 + x[j], 7) - Heaviside(h4 + x[j]) * pow(h4 + x[j], 7)) / h1) / h1 + ((Heaviside(h4 - h2 + x[j]) * pow(h4 - h2 + x[j], 7) + Heaviside(h4 - h2 - h1 + x[j]) * pow(h1 + h2 - h4 - x[j], 7)) / h1 - (Heaviside(h1 - h2 + h4 + x[j]) * pow(h1 - h2 + h4 + x[j], 7) - Heaviside(h4 - h2 + x[j]) * pow(h4 - h2 + x[j], 7)) / h1) / h1) / h2 + (((Heaviside(h4 - h1 + x[j]) * pow(h4 - h1 + x[j], 7) - Heaviside(h4 + x[j]) * pow(h4 + x[j], 7)) / h1 + (Heaviside(h1 + h4 + x[j]) * pow(h1 + h4 + x[j], 7) - Heaviside(h4 + x[j]) * pow(h4 + x[j], 7)) / h1) / h1 - ((Heaviside(h2 - h1 + h4 + x[j]) * pow(h2 - h1 + h4 + x[j], 7) - Heaviside(h2 + h4 + x[j]) * pow(h2 + h4 + x[j], 7)) / h1 + (Heaviside(h1 + h2 + h4 + x[j]) * pow(h1 + h2 + h4 + x[j], 7) - Heaviside(h2 + h4 + x[j]) * pow(h2 + h4 + x[j], 7)) / h1) / h1) / h2) / h2 - ((((Heaviside(h2 - h3 + h4 + x[j]) * pow(h2 - h3 + h4 + x[j], 7) - Heaviside(h2 - h1 - h3 + h4 + x[j]) * pow(h2 - h1 - h3 + h4 + x[j], 7)) / h1 + (Heaviside(h2 - h3 + h4 + x[j]) * pow(h2 - h3 + h4 + x[j], 7) - Heaviside(h1 + h2 - h3 + h4 + x[j]) * pow(h1 + h2 - h3 + h4 + x[j], 7)) / h1) / h1 - ((Heaviside(h4 - h3 + x[j]) * pow(h4 - h3 + x[j], 7) + Heaviside(h4 - h3 - h1 + x[j]) * pow(h1 + h3 - h4 - x[j], 7)) / h1 - (Heaviside(h1 - h3 + h4 + x[j]) * pow(h1 - h3 + h4 + x[j], 7) - Heaviside(h4 - h3 + x[j]) * pow(h4 - h3 + x[j], 7)) / h1) / h1) / h2 - (((Heaviside(h4 - h3 - h2 + x[j]) * pow(h2 + h3 - h4 - x[j], 7) - Heaviside(h4 - h2 - h3 - h1 + x[j]) * pow(h1 + h2 + h3 - h4 - x[j], 7)) / h1 + (Heaviside(h1 - h2 - h3 + h4 + x[j]) * pow(h1 - h2 - h3 + h4 + x[j], 7) + Heaviside(h4 - h3 - h2 + x[j]) * pow(h2 + h3 - h4 - x[j], 7)) / h1) / h1 + ((Heaviside(h4 - h3 + x[j]) * pow(h4 - h3 + x[j], 7) + Heaviside(h4 - h3 - h1 + x[j]) * pow(h1 + h3 - h4 - x[j], 7)) / h1 - (Heaviside(h1 - h3 + h4 + x[j]) * pow(h1 - h3 + h4 + x[j], 7) - Heaviside(h4 - h3 + x[j]) * pow(h4 - h3 + x[j], 7)) / h1) / h1) / h2) / h2) / h3 + (((((Heaviside(h4 - h1 + x[j]) * pow(h4 - h1 + x[j], 7) - Heaviside(h4 + x[j]) * pow(h4 + x[j], 7)) / h1 + (Heaviside(h1 + h4 + x[j]) * pow(h1 + h4 + x[j], 7) - Heaviside(h4 + x[j]) * pow(h4 + x[j], 7)) / h1) / h1 + ((Heaviside(h4 - h2 + x[j]) * pow(h4 - h2 + x[j], 7) + Heaviside(h4 - h2 - h1 + x[j]) * pow(h1 + h2 - h4 - x[j], 7)) / h1 - (Heaviside(h1 - h2 + h4 + x[j]) * pow(h1 - h2 + h4 + x[j], 7) - Heaviside(h4 - h2 + x[j]) * pow(h4 - h2 + x[j], 7)) / h1) / h1) / h2 + (((Heaviside(h4 - h1 + x[j]) * pow(h4 - h1 + x[j], 7) - Heaviside(h4 + x[j]) * pow(h4 + x[j], 7)) / h1 + (Heaviside(h1 + h4 + x[j]) * pow(h1 + h4 + x[j], 7) - Heaviside(h4 + x[j]) * pow(h4 + x[j], 7)) / h1) / h1 - ((Heaviside(h2 - h1 + h4 + x[j]) * pow(h2 - h1 + h4 + x[j], 7) - Heaviside(h2 + h4 + x[j]) * pow(h2 + h4 + x[j], 7)) / h1 + (Heaviside(h1 + h2 + h4 + x[j]) * pow(h1 + h2 + h4 + x[j], 7) - Heaviside(h2 + h4 + x[j]) * pow(h2 + h4 + x[j], 7)) / h1) / h1) / h2) / h2 - ((((Heaviside(h3 - h2 + h4 + x[j]) * pow(h3 - h2 + h4 + x[j], 7) - Heaviside(h3 - h2 - h1 + h4 + x[j]) * pow(h3 - h2 - h1 + h4 + x[j], 7)) / h1 + (Heaviside(h3 - h2 + h4 + x[j]) * pow(h3 - h2 + h4 + x[j], 7) - Heaviside(h1 - h2 + h3 + h4 + x[j]) * pow(h1 - h2 + h3 + h4 + x[j], 7)) / h1) / h1 + ((Heaviside(h3 - h1 + h4 + x[j]) * pow(h3 - h1 + h4 + x[j], 7) - Heaviside(h3 + h4 + x[j]) * pow(h3 + h4 + x[j], 7)) / h1 + (Heaviside(h1 + h3 + h4 + x[j]) * pow(h1 + h3 + h4 + x[j], 7) - Heaviside(h3 + h4 + x[j]) * pow(h3 + h4 + x[j], 7)) / h1) / h1) / h2 + (((Heaviside(h2 + h3 + h4 + x[j]) * pow(h2 + h3 + h4 + x[j], 7) - Heaviside(h2 - h1 + h3 + h4 + x[j]) * pow(h2 - h1 + h3 + h4 + x[j], 7)) / h1 - (Heaviside(h1 + h2 + h3 + h4 + x[j]) * pow(h1 + h2 + h3 + h4 + x[j], 7) - Heaviside(h2 + h3 + h4 + x[j]) * pow(h2 + h3 + h4 + x[j], 7)) / h1) / h1 + ((Heaviside(h3 - h1 + h4 + x[j]) * pow(h3 - h1 + h4 + x[j], 7) - Heaviside(h3 + h4 + x[j]) * pow(h3 + h4 + x[j], 7)) / h1 + (Heaviside(h1 + h3 + h4 + x[j]) * pow(h1 + h3 + h4 + x[j], 7) - Heaviside(h3 + h4 + x[j]) * pow(h3 + h4 + x[j], 7)) / h1) / h1) / h2) / h2) / h3) / h3 / 5040) *(1/(h4*h4))));

		if (p<1e-6) p=0;
		this->kern.push_back(p);
	}

}

//====================================================================================
//			FIND PDF OF THE MIDPOINT  - OR EACH VERTEX
//====================================================================================

void convolution::findPDF(){

	vector <double> shifts;
	int i,j,k,w,p;

	if(NUM_VERTEX==4){
		for (i=0 ; i<(int)(distrib1.size()) ; i++){
			for ( j=0 ; j<(int)(distrib2.size()) ; j++){
				for (k=0 ; k<(int)distrib3.size() ; k++){
					for (w=0 ; w<(int)distrib3.size() ; w++){
						shifts.push_back(distrib1[i]+distrib2[j]+distrib3[k]+distrib4[w]);
					}
				}
			}
		}
		this->degree7Kernel();
	}
	else if (NUM_VERTEX==1){
		for (i=0 ; i<(int)(distrib1.size()) ; i++){
			shifts.push_back(distrib1[i]);
		}
		this->LinearKernel();
	}

	kernelWidth = (kern.size()-1)*(sampRate/2) ; //this is half width of kernel in units
	startPT = *(std::min_element(shifts.begin(),shifts.end()))-kernelWidth;   //it is the starting point of the PDF
	endPT = *(std::max_element(shifts.begin(),shifts.end()))+kernelWidth;


	//cout << "start of the pdf = " << startPT << endl;
	//cout << "end of the pdf = " << endPT << endl;
	for (double q=startPT; q<=endPT ; q+=sampRate){
		pdf.push_back(0.0);
	}
	for (i=0; i<(int)shifts.size(); i++){
		for (p=0; p<(int)kern.size();p++){

			pdf[(int)((shifts[i]-kernelWidth-startPT)/sampRate)+p] =
					pdf[(int)((shifts[i]-kernelWidth-startPT)/sampRate)+p] + kern[p];

		}
	}
	this->printVector2File(pdf);
}

//====================================================================================
//			DICIDE   Pr(X>isovalue)
//====================================================================================

double convolution::decider (double isovalue){
	if(isovalue <= this->startPT)
		return 1.0;
	if (isovalue >= this->endPT)
		return 0.0;
	double prob = 0.0;
	for (int i=0 ; (int)((isovalue-startPT)/sampRate)+i < (int)pdf.size() ; i++){
		prob += pdf[(int)((isovalue-startPT)/sampRate)+i];
	}
	int norm_factor=1;
	for (int i = 0 ; i< NUM_VERTEX ; i++){
		norm_factor *= distrib1.size();   //assuming distrib1.size() = distrib2.size() = ...
	}
	//otherwise norm_factor = (distrib1.size()*distrib2.size()*distrib3.size()*distrib4.size());

	//cout << "divide by " << norm_factor << endl;
	return (prob*sampRate)/norm_factor;

}

// ==============================================================================
// helper functions
//===============================================================================
int convolution::Heaviside (double value){
	return (value>0);
}

void convolution::printVector(vector <double> vect){
	std::cout << "[";
	for (int i=0; i<(int)vect.size() ; i++){
		std::cout << setprecision(18)  << vect[i] << " ," ;
	}
	std::cout << "];" << std::endl;
}

void convolution::printVector2File(vector <double> vect){
	ofstream myfile ("mypdf.txt");
	if (myfile.is_open())
	{
		myfile << "[";
		for (int i=0; i<(int)vect.size() ; i++){
			myfile << setprecision(6)  << vect[i] << " ,\n" ;
		}
		myfile << "];\n";
		myfile.close();

	}
	else cout << "Unable to open file";
}
int readData (string filename, double means[], double binWidths []){
	ifstream ifile (filename.c_str());
	double mean , binWidth;
	if ( ifile.fail() ) {
		cout << "File not found: " << filename << endl;;
		return -1;
	}
	int q=0; int sampNum = 0;
	char ch = ifile.peek();
	for (int i=0 ; i< NUM_VERTEX ; i++){   //4 is number of vertices of the face we are looking at
		ifile.ignore(80, '\n');
		ifile.ignore(80, '\n');
		ifile.ignore(80, '\n');
		ifile.ignore(80, '\n');

		ifile.ignore(80, ':');
		while( ch != '\n'){
			ifile >> mean;
			means[q] = (mean/NUM_VERTEX);
			ch = ifile.peek();
			q++;
		}
		if(i==0) {sampNum = q;}
		ch = ' ';
		ifile.ignore(80, '\n');
		ifile.ignore(80, '\n');
		ifile.ignore(80, ':');
		ifile >> binWidth;
		binWidths[i]=(binWidth/NUM_VERTEX);

	}
	ifile.close();
	return sampNum;
}
