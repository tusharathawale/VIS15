/*
 * pdf.cpp
 *
 *  Created on: Mar 19, 2014
 *      Author: esakhaee
 */

#include <iostream>
#include "functions.h"

using namespace std;




int main(){

	double isoValue = 56;
	double samplingRate = 0.01;

	double means [4*27];
	double bw [4];
    // readData function is ONLY used to read the arrays from the provided file
	int sampNum = readData("face2.txt", means, bw);  //here sampNum = 27;

	//the prototype you need to call in your program
	//means is a double array of size say 27*4
	//bw is double array of size 5
	//sampNum denotes how many samples at each pdf (say 27)
	bool decision =  MidPointWrapper(means, bw , sampNum, isoValue , samplingRate);
}


