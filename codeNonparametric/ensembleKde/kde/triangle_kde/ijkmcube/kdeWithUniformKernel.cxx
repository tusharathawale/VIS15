//---------------------------------------------------------------------------------------
/*
	Title : Uncertainty Quantification in Linear Interpolation for Isosurface Extraction	
       Authors: Tushar Athawale, Alireza Entezari
        Date  : Jun 27, 2013.
*/
//---------------------------------------------------------------------------------------

#include "stdio.h"
#include <iostream>
#include "stdlib.h"
#include "math.h"
#include <algorithm>
#include "kdeWithUniformKernel.h"

# define infinity 10000
# define minus_infinity -10000
#define _USE_MATH_DEFINES

using namespace std;

// p is pdf over -infinity to infinity range. Function returns part of pdf just over range 0 and 1.
piecewise z_density_uniform :: getPdfOver0To1(piecewise p)
{

	//cout<<"p.pc[0].k3:"<<p.pc[0].k3;
	double epsilon = 1e-10;
	// if pieces are already lying in 0 to 1 range return the input pieces as they are
	 if((fabs(p.limits[0])<epsilon) && (fabs(p.limits[p.numPieces]-1)<epsilon))
		return p;

	piecewise p01, answer;

	int  k=0;       
	int lowestLimitGreaterThanZeroFlag = 0;
	int highestLimitLessThanOneFlag = 0;

	if(p.limits[0] > 0)
		lowestLimitGreaterThanZeroFlag = 1;
	if(p.limits[p.numPieces] < 1)
		highestLimitLessThanOneFlag = 1;

	// obtain number of pieces in the piecewise function 
	int piece_count = p.numPieces;

	// a1<p1<b1<p2<a2  .. so if 2 pieces, then 3 limits
	int limit_count = piece_count + 1;       

	// Determine limits that fall in [0,1]

	for (int i=0; i<(limit_count-1); i++) 
	{
		//if both limits are below 0 or above 1 don't do anything. skip to next limit
     
		// if any limits of piece lie in 0-1 range   
		if (((p.limits[i] >= 0) && (p.limits[i] <= 1)) || ((p.limits[i+1] >= 0) && (p.limits[i+1] <= 1)))
		{
			p01.limits[k] = p.limits[i];
			p01.limits[k+1] = p.limits[i+1];
			p01.pc[k] = p.pc[i];
			k++;
		}

		// if a piece limits completely overlap [0,1] such that 1 limit is less than 0 and other is greater than 1  
		if ((p.limits[i] < 0) && (p.limits[i+1] > 1))
		{
			p01.limits[k] = 0;
			p01.limits[k+1] = 1;
			p01.pc[k] = p.pc[i];
			k++; 
		}	  

	} 

	// if limit is less than 0 then set it to 0. That can happen for only 1st entry of the lm_in_0_1[10][2] array
	// if limit is greater than 1 then set it to 1. That can happen for only last entry of the lm_in_0_1[10][2] array

	if (p01.limits[0] < 0)
		p01.limits[0] = 0;
	if (p01.limits[k] > 1)
		p01.limits[k] = 1;

	p01.numPieces = k;

	// Add k1=0,k2=0 pieces at the ends if necessary
	if((lowestLimitGreaterThanZeroFlag == 1) && (highestLimitLessThanOneFlag!=1))
	{
		answer.limits[0] = 0;
		answer.pc[0].k1 = 0;
		answer.pc[0].k2 = 0;
		//answer.pc[0].k3 = denominator;
		answer.pc[0].k3 = p.pc[0].k3;
		for(int i=0; i<=k; i++)
		{
			answer.limits[i+1] = p01.limits[i];
		}
		for(int i=0; i<k; i++)
		{
			answer.pc[i+1] = p01.pc[i];
		}
		answer.numPieces = k + 1;
	}
	else if((lowestLimitGreaterThanZeroFlag != 1) && (highestLimitLessThanOneFlag==1))
	{
		for(int i=0; i<=k; i++)
		{
			answer.limits[i] = p01.limits[i];
		}
		for(int i=0; i<k; i++)
		{
			answer.pc[i] = p01.pc[i];
		}
		answer.limits[k+1] = 1;
		answer.pc[k].k1 = 0;
		answer.pc[k].k2 = 0;
		//answer.pc[k].k3 = denominator;
		answer.pc[k].k3 = p.pc[0].k3;
		answer.numPieces = k + 1;
	}
	else if((lowestLimitGreaterThanZeroFlag == 1) && (highestLimitLessThanOneFlag == 1))
	{
		answer.limits[0] = 0;
		answer.pc[0].k1 = 0;
		answer.pc[0].k2 = 0;
		//bug removed : following line added
		answer.pc[0].k3 = p.pc[0].k3;
		for(int i=0; i<=k; i++)
		{
			answer.limits[i+1] = p01.limits[i];
		}
		for(int i=0; i<k; i++)
		{
			answer.pc[i+1] = p01.pc[i];
		}
		answer.limits[k+2] = 1;
		answer.pc[k+1].k1 = 0;
		answer.pc[k+1].k2 = 0;	
		//answer.pc[k+1].k3 = denominator;
		answer.pc[k+1].k3 = p.pc[0].k3;
		answer.numPieces = k + 2;
	}
	// bug fix
	else
	{
		for(int i=0; i<=k; i++)
		{
			answer.limits[i] = p01.limits[i];
		}
		for(int i=0; i<k; i++)
		{
			answer.pc[i] = p01.pc[i];
		}
		answer.numPieces = k;	
	}
	
	/*for(int i=0;i<answer.numPieces;i++)
		cout<<" k1:"<<answer.pc[i].k1<<" k2:"<<answer.pc[i].k2<<" k3:"<<answer.pc[i].k3<<"\n";
	cout<<"\n";
		for(int i=0;i<=answer.numPieces;i++)
		cout<<answer.limits[i]<<" ";
	cout<<"\n";*/
	
	return answer;
	//return p01; 
}

double z_density_uniform :: getExpectedValOverSinglePiece(double k1, double k2, double k3, double a, double b)
{
	double expectation = 0;
	double epsilon = 1e-10;

	double c1 = (double)(k2/k3)*log((double)((b)/(a))); 
	double c2 = (double)(k1/k3)*log((double)((b-1)/(a-1))); 
	
	if(fabs(k1)<epsilon && fabs(k2)<epsilon)
		expectation = 0;

	else if(fabs(k1)<epsilon && fabs(k2)>epsilon)
		expectation = c1;

	else if(fabs(k1)>epsilon && fabs(k2)<epsilon)
		expectation = c2 - (double)(k1/k3)*((double)(1 / (b - 1))) + (double)(k1/k3)*((double)(1 / (a - 1)));
	
	else 
		expectation = c1 + c2 - (double)(k1/k3)*((double)(1 / (b - 1))) + (double)(k1/k3)*((double)(1 / (a - 1)));

	return expectation;
}

double z_density_uniform :: getEdgeCrossingProbability(double k1, double k2, double k3, double a, double b)
{
	double crossingProb = 0;
	double epsilon = 1e-10; 
	if((fabs(k1)<epsilon) && (fabs(k2)<epsilon))
		crossingProb = 0;
	else if((fabs(k1)<epsilon) && (fabs(k2)>epsilon))
	{
		crossingProb = (double)(k2/k3)*((double)(1/a)-(double)(1/b));
	}
	else if((fabs(k1)>epsilon) && (fabs(k2)<epsilon))
	{
		crossingProb = (double)(k1/k3)*((double)(1/(1-b))-(double)(1/(1-a)));	
	}
	else
	{ 			
		crossingProb = (double)(1/k3)*(k2*((double)(1/a) - (double)(1/b)) + k1*((double)(1/(a-1))- (double)(1/(b-1))));
	}
	return crossingProb;
}

double z_density_uniform :: getSecondMoment(double k1, double k2, double k3, double a, double b)
{
	double secondMoment = 0;
	double epsilon = 1e-10; 
	
	if(fabs(k1)<epsilon && fabs(k2)<epsilon)
		secondMoment = 0;

	else if(fabs(k1)<epsilon && fabs(k2)>epsilon)
		secondMoment = (double)(k2/k3)*(b-a);

	else if(fabs(k1)>epsilon && fabs(k2)<epsilon)
		secondMoment = (double)(k1/k3)*(b-a) + (double)(k1/k3)*((double)(1/(a-1))-(double)(1/(b-1))) + (double)((2*k1)/k3)*log((double)((b-1)/(a-1)));
	else 
		secondMoment = ((double)(k1/k3) + (double)(k2/k3))*(b-a) + (double)(k1/k3)*((double)(1/(a-1))-(double)(1/(b-1))) + (double)((2*k1)/k3)*log((double)((b-1)/(a-1)));

	return secondMoment;

}

void z_density_uniform :: Compute0To1(piecewise p, double* expt, double* cross_prob, double* var, double *second_moment, double *first_moment)
{
/*	for(int i=0;i<p.numPieces;i++)
		cout<<" k1:"<<p.pc[i].k1<<" k2:"<<p.pc[i].k2<<" k3:"<<p.pc[i].k3<<"\n";
	cout<<"\n";
		for(int i=0;i<=p.numPieces;i++)
		cout<<p.limits[i]<<" ";
	cout<<"\n";*/

	piecewise pdf01 = getPdfOver0To1(p);

	/*for(int i=0;i<pdf01.numPieces;i++)
		cout<<" k1:"<<pdf01.pc[i].k1<<" k2:"<<pdf01.pc[i].k2<<" k3:"<<pdf01.pc[i].k3<<"\n";
	cout<<"\n";
		for(int i=0;i<=pdf01.numPieces;i++)
		cout<<pdf01.limits[i]<<" ";
	cout<<"\n";*/

	// compute expectation (pdf which is unnormalized over range 0 to 1)
	double expectation_over_0_1 = 0;	
	for (int i=0; i<pdf01.numPieces; i++) 
	{				
		expectation_over_0_1 += getExpectedValOverSinglePiece(pdf01.pc[i].k1 , pdf01.pc[i].k2, pdf01.pc[i].k3, pdf01.limits[i], pdf01.limits[i+1]);
	} 

	// compute crossing probability 
	double crossProb_0_1 = 0;	
	for (int i=0; i<pdf01.numPieces; i++) 
	{				
		crossProb_0_1 += getEdgeCrossingProbability(pdf01.pc[i].k1 , pdf01.pc[i].k2, pdf01.pc[i].k3, pdf01.limits[i], pdf01.limits[i+1]);
	} 

	// compute variance 
	double sec_0_1 = 0;
	for (int i=0; i<pdf01.numPieces; i++) 
	{				
		sec_0_1 += getSecondMoment(pdf01.pc[i].k1 , pdf01.pc[i].k2, pdf01.pc[i].k3, pdf01.limits[i], pdf01.limits[i+1]);
	} 

	*expt = expectation_over_0_1;
	*cross_prob = crossProb_0_1;
	*var = (double)(((sec_0_1*crossProb_0_1) - (double)(expectation_over_0_1*expectation_over_0_1))/(crossProb_0_1*crossProb_0_1));
	*second_moment = sec_0_1;
	*first_moment = (double)(expectation_over_0_1/crossProb_0_1);
}

// Create piece by setting piece coefficients
// Set k1, k2, k3, k4 for piece type 0, Set k1 if piece is of type 1 or 2
piecewise z_density_uniform :: setPiece(piecewise P, int pieceIndex, double k1, double k2)
{                
	P.pc[pieceIndex].k1 = k1; 
	P.pc[pieceIndex].k2 = k2;
	P.pc[pieceIndex].k3 = denominator;
	return P; 
}

// 1) nonOverlapping Intervals
piecewise z_density_uniform :: nonOverlapping(double mu1, double delta1, double mu2, double delta2, double c) 
{	
	piecewise P; 
	int pieceID;

	// mu1-delta1 <= c <= mu1+delta1
	if ((c >= (mu1-delta1)) && (c <= (mu1+delta1))) 
	{
		// Vertex Order : SPQR			
		// specify number of pieces
		P.numPieces = 3;
		pieceID = 0;

		// set pieces

		// Set piece 0 coefficients                    
		P = setPiece(P, pieceID, -f4, f1);        
		pieceID++;           	

		// Set piece 1 coefficients
		P = setPiece(P, pieceID, (double)(4*delta2*(mu2-c)),0);   
		pieceID++;          	

		// Set piece 2 coefficients
		P = setPiece(P, pieceID, -f4, f2);    
		pieceID++;                 	

		// set limits
		P.limits[0] = slope_OS;  
		P.limits[1] = slope_OP;
		P.limits[2] = slope_OQ; 
		P.limits[3] = slope_OR;

	}    

	// mu1+delta1 < c <= m2-delta2           
	else if ((c > (mu1+delta1)) && (c <= (mu2-delta2))) 
	{
		
		// Vertex Order : PQSR		

		if(slope_OS > slope_OQ)
		{
			// specify number of pieces
			P.numPieces = 3;
			pieceID = 0; 


			// set pieces

			// Set piece 0 coefficients                    
			P = setPiece(P, pieceID, f3, -f1);        
			pieceID++;           	

			// Set piece 1 coefficients
			P = setPiece(P, pieceID, 0, (double)(4*delta1*(c-mu1)));   
			pieceID++;          	

			// Set piece 2 coefficients
			P = setPiece(P, pieceID, -f4, f2);    
			pieceID++;                 	

			// set limits
			P.limits[0] = slope_OP;  
			P.limits[1] = slope_OQ;
			P.limits[2] = slope_OS; 
			P.limits[3] = slope_OR;

		}

		else if (slope_OQ >= slope_OS)
		{			
			// Vertex Order : PSQR		

			// specify number of pieces
			P.numPieces = 3;
			pieceID = 0; 


			// set pieces

			// Set piece 0 coefficients                    
			P = setPiece(P, pieceID, f3, -f1);        
			pieceID++;           	

			// Set piece 1 coefficients
			P = setPiece(P, pieceID, (double)(4*delta2*(mu2-c)),0);   
			pieceID++;          	

			// Set piece 2 coefficients
			P = setPiece(P, pieceID, -f4, f2);    
			pieceID++;               

			// set limits
			P.limits[0] = slope_OP;  
			P.limits[1] = slope_OS;
			P.limits[2] = slope_OQ; 
			P.limits[3] = slope_OR;
		}                
	}

	// mu2-delta2 < c <= mu2+delta2
	else if ((c > (mu2-delta2)) && (c <= (mu2+delta2)))
	{				
				// Vertex Order : PQRS		
	
				// specify number of pieces
				P.numPieces = 3;
				pieceID = 0; 

				// set pieces

				// Set piece 0 coefficients                    
				P = setPiece(P, pieceID, f3, -f1);        
				pieceID++;           	

				// Set piece 1 coefficients
				P = setPiece(P, pieceID, 0, (double)(4*delta1*(c-mu1)));   
				pieceID++;          	

				// Set piece 2 coefficients
				P = setPiece(P, pieceID, f4, -f1);    
				pieceID++;       

				// set limits
				P.limits[0] = slope_OP;  
				P.limits[1] = slope_OQ;
				P.limits[2] = slope_OR; 
				P.limits[3] = slope_OS;  
	}

	return P; 
}


// 2) Overlapping Intervals
piecewise z_density_uniform :: overlapping(double mu1, double delta1, double mu2, double delta2, double c) 
{	
	piecewise P;
	int pieceID;

	// mu1-delta1 <= c <= mu2-delta2
	if ((c >= mu1-delta1) && (c <= mu2-delta2))
	{

		// Vertex Order : PQRS		

		// specify number of pieces
		P.numPieces = 5;
		pieceID = 0; 

		// set pieces

		// Set piece 0 coefficients                        
		P = setPiece(P, pieceID, -f4, f1);        
		pieceID++;          

		// Set piece 1 coefficients
		P = setPiece(P, pieceID, (double)(4*delta2*(mu2-c)),0);   
		pieceID++;       		  	

		// Set piece 2 coefficients
		P = setPiece(P, pieceID, -f4, f2);        
		pieceID++;    					      	

		// Set piece 3 coefficients                        
		P = setPiece(P, pieceID, 0, 0);   
		pieceID++;   

		// Set piece 4 coefficients
		P = setPiece(P, pieceID, -f4, f1);        
		pieceID++;       

		// set limits
		// slope_OS is -infinity
		P.limits[0] = -infinity;  
		P.limits[1] = slope_OP;
		P.limits[2] = slope_OQ; 
		P.limits[3] = slope_OR;
		P.limits[4] = slope_OS;
		P.limits[5] = infinity;
	}     

	// mu2-delta2 < c <= mu1+delta1 
	else if ((c > mu2-delta2) && (c <= mu1+delta1))
	{
		if (slope_OS < slope_OQ)
		{
			// Vertex Order : PSQR	

			// specify number of pieces
			P.numPieces = 5;
			pieceID = 0; 

			// set pieces

			// Set piece 0 coefficients                        
			P = setPiece(P, pieceID, f4, f1);        
			pieceID++;          

			// Set piece 1 coefficients
			P = setPiece(P, pieceID, (double)(2*((c-mu2)*(c-mu2)+delta2*delta2)),0);   
			pieceID++;       		  	

			// Set piece 2 coefficients
			P = setPiece(P, pieceID, f3, f1);        
			pieceID++;    					      	

			// Set piece 3 coefficients                        
			P = setPiece(P, pieceID, 0,(double)(2*((c-mu1)*(c-mu1)+delta1*delta1)));
			pieceID++;   

			// Set piece 4 coefficients
			P = setPiece(P, pieceID, f4, f1);        
			pieceID++;       

			// set limits

			// slope_OS is -infinity
			P.limits[0] = -infinity;  
			P.limits[1] = slope_OP;
			P.limits[2] = slope_OS; 
			P.limits[3] = slope_OQ;
			P.limits[4] = slope_OR;
			P.limits[5] = infinity;
		}

		else if (slope_OQ <= slope_OS) 
		{

			// Vertex Order : PQSR	

			// specify number of pieces
			P.numPieces = 5;
			pieceID = 0; 

			// set pieces

			// Set piece 0 coefficients                        
			P = setPiece(P, pieceID, f4, f1);        
			pieceID++;          

			// Set piece 1 coefficients
			P = setPiece(P, pieceID, (double)(2*((c-mu2)*(c-mu2)+delta2*delta2)),0);   
			pieceID++;       		  	

			// Set piece 2 coefficients
			P = setPiece(P, pieceID, f4, f2);        
			pieceID++;    					      	

			// Set piece 3 coefficients                        
			P = setPiece(P, pieceID, 0,(double)(2*((c-mu1)*(c-mu1)+delta1*delta1)));
			pieceID++;   

			// Set piece 4 coefficients
			P = setPiece(P, pieceID, f4, f1);        
			pieceID++;       

			// set limits
			// slope_OS is -infinity
			P.limits[0] = -infinity;  
			P.limits[1] = slope_OP;
			P.limits[2] = slope_OQ; 
			P.limits[3] = slope_OS;
			P.limits[4] = slope_OR;
			P.limits[5] = infinity;
		}		
	}

	// case mu1+delta1 < c <= mu2+delta2
	else if ((c > mu1+delta1) && (c <= mu2+delta2))
	{
			// Vertex Order : SPQR	

			// specify number of pieces
			P.numPieces = 5;
			pieceID = 0; 

			// set pieces

			// Set piece 0 coefficients                        
			P = setPiece(P, pieceID, f4, -f1);        
			pieceID++;          

			// Set piece 1 coefficients
			P = setPiece(P, pieceID, 0, 0);  
			pieceID++;       		  	

			// Set piece 2 coefficients
			P = setPiece(P, pieceID, f3, -f1);        
			pieceID++;    					      	

			// Set piece 3 coefficients          
			P = setPiece(P, pieceID, 0, (double)(4*delta1*(c-mu1)));                  
			pieceID++;   

			// Set piece 4 coefficients
			P = setPiece(P, pieceID, f4, -f1);        
			pieceID++;       

			// set limits
			// slope_OS is -infinity
			P.limits[0] = -infinity;  
			P.limits[1] = slope_OS;
			P.limits[2] = slope_OP; 
			P.limits[3] = slope_OQ;
			P.limits[4] = slope_OR;
			P.limits[5] = infinity;     
	}	

	return P;	
}

// 3a) Contained Intervals
piecewise z_density_uniform :: containedA(double mu1, double delta1, double mu2, double delta2, double c) 
{	
	piecewise P;
	int pieceID;

	// mu1-delta1 <= c <= mu2-delta2 
	if ((c >= (mu1-delta1)) && (c <= (mu2-delta2))) 
	{

		// Vertex Order : QRSP		

		// specify number of pieces
		P.numPieces = 5;
		pieceID = 0; 


		// set pieces

		// Set piece 0 coefficients       
		P = setPiece(P, pieceID, (double)(4*delta2*(mu2-c)), 0);                  
		pieceID++;          

		// Set piece 1 coefficients
		P = setPiece(P, pieceID, -f4, f2);  
		pieceID++;       		  	

		// Set piece 2 coefficients
		P = setPiece(P, pieceID, 0, 0);        
		pieceID++;    					      	

		// Set piece 3 coefficients          
		P = setPiece(P, pieceID, -f4, f1);                         
		pieceID++;   

		// Set piece 4 coefficients
		P = setPiece(P, pieceID, (double)(4*delta2*(mu2-c)), 0);   
		pieceID++;       

		// set limits
		// slope_OS is -infinity
		P.limits[0] = -infinity;  
		P.limits[1] = slope_OQ;
		P.limits[2] = slope_OR; 
		P.limits[3] = slope_OS;
		P.limits[4] = slope_OP;
		P.limits[5] = infinity;
	}	   

	// mu2-delta2 < c <= mu2+delta2 
	else if ((c > (mu2-delta2)) && (c <= (mu2+delta2)))
	{
		if ((slope_OS <= slope_OQ) && (slope_OR <= slope_OP))
		{

			// Vertex Order : SQRP	

			// specify number of pieces
			P.numPieces = 5;
			pieceID = 0; 

			// set pieces

			// Set piece 0 coefficients             
			P = setPiece(P, pieceID, (double)(2*((c-mu2)*(c-mu2)+delta2*delta2)),0);             
			pieceID++;          

			// Set piece 1 coefficients
			P = setPiece(P, pieceID, f3, f1);         
			pieceID++;       		  	

			// Set piece 2 coefficients
			P = setPiece(P, pieceID, 0,(double)(2*((c-mu1)*(c-mu1)+delta1*delta1)));    
			pieceID++;    					      	

			// Set piece 3 coefficients                        
			P = setPiece(P, pieceID, f4, f1);    
			pieceID++;   

			// Set piece 4 coefficients
			P = setPiece(P, pieceID, (double)(2*((c-mu2)*(c-mu2)+delta2*delta2)),0);             
			pieceID++;     

			// set limits
			// slope_OS is -infinity
			P.limits[0] = -infinity;  
			P.limits[1] = slope_OS;
			P.limits[2] = slope_OQ; 
			P.limits[3] = slope_OR;
			P.limits[4] = slope_OP;
			P.limits[5] = infinity;
		}

		else if ((slope_OS >= slope_OQ) && (slope_OR <= slope_OP))
		{
	
			// Vertex Order : QSRP
	
			// specify number of pieces
			P.numPieces = 5;
			pieceID = 0; 

			// set pieces

			// Set piece 0 coefficients             
			P = setPiece(P, pieceID, (double)(2*((c-mu2)*(c-mu2)+delta2*delta2)),0);             
			pieceID++;          

			// Set piece 1 coefficients
			P = setPiece(P, pieceID, f4, f2);         
			pieceID++;       		  	

			// Set piece 2 coefficients
			P = setPiece(P, pieceID, 0,(double)(2*((c-mu1)*(c-mu1)+delta1*delta1)));    
			pieceID++;    					      	

			// Set piece 3 coefficients                        
			P = setPiece(P, pieceID, f4, f1);    
			pieceID++;   

			// Set piece 4 coefficients
			P = setPiece(P, pieceID, (double)(2*((c-mu2)*(c-mu2)+delta2*delta2)),0);             
			pieceID++;     


			// set limits
			// slope_OS is -infinity
			P.limits[0] = -infinity;  
			P.limits[1] = slope_OQ;
			P.limits[2] = slope_OS; 
			P.limits[3] = slope_OR;
			P.limits[4] = slope_OP;
			P.limits[5] = infinity;
		}

		else if ((slope_OS <= slope_OQ) && (slope_OR >= slope_OP))
		{

			// Vertex Order : SQPR	
	
			// specify number of pieces
			P.numPieces = 5;
			pieceID = 0; 

			// set pieces

			// Set piece 0 coefficients             
			P = setPiece(P, pieceID, (double)(2*((c-mu2)*(c-mu2)+delta2*delta2)),0);             
			pieceID++;          

			// Set piece 1 coefficients
			P = setPiece(P, pieceID, f3, f1);         
			pieceID++;       		  	

			// Set piece 2 coefficients
			P = setPiece(P, pieceID, 0,(double)(2*((c-mu1)*(c-mu1)+delta1*delta1)));    
			pieceID++;    					      	

			// Set piece 3 coefficients                        
			P = setPiece(P, pieceID, f3, f2);    
			pieceID++;   

			// Set piece 4 coefficients
			P = setPiece(P, pieceID, (double)(2*((c-mu2)*(c-mu2)+delta2*delta2)),0);             
			pieceID++;     


			// set limits
			// slope_OS is -infinity
			P.limits[0] = -infinity;  
			P.limits[1] = slope_OS;
			P.limits[2] = slope_OQ; 
			P.limits[3] = slope_OP;
			P.limits[4] = slope_OR;
			P.limits[5] = infinity;
		}

		else if ((slope_OS >= slope_OQ) && (slope_OR >= slope_OP)) 
		{
		
			// Vertex Order : QSPR	

			// specify number of pieces
			P.numPieces = 5;
			pieceID = 0; 

			// set pieces

			// Set piece 0 coefficients             
			P = setPiece(P, pieceID, (double)(2*((c-mu2)*(c-mu2)+delta2*delta2)),0);             
			pieceID++;          

			// Set piece 1 coefficients
			P = setPiece(P, pieceID, f4, f2);         
			pieceID++;       		  	

			// Set piece 2 coefficients
			P = setPiece(P, pieceID, 0,(double)(2*((c-mu1)*(c-mu1)+delta1*delta1)));    
			pieceID++;    					      	

			// Set piece 3 coefficients                        
			P = setPiece(P, pieceID, f3, f2);    
			pieceID++;   

			// Set piece 4 coefficients
			P = setPiece(P, pieceID, (double)(2*((c-mu2)*(c-mu2)+delta2*delta2)),0);             
			pieceID++;     	


			// set limits

			// slope_OS is -infinity
			P.limits[0] = -infinity;  
			P.limits[1] = slope_OQ;
			P.limits[2] = slope_OS; 
			P.limits[3] = slope_OP;
			P.limits[4] = slope_OR;
			P.limits[5] = infinity;
		}	
	}	

	// mu2+delta2 < c <= mu1+delta1 
	else if ((c > (mu2+delta2)) && (c <= (mu1+delta1)))
	{

		// Vertex Order : SPQR

		// specify number of pieces
		P.numPieces = 5;
		pieceID = 0; 

		// set pieces

		// Set piece 0 coefficients       
		P = setPiece(P, pieceID, (double)(4*delta2*(c-mu2)), 0);                  
		pieceID++;          

		// Set piece 1 coefficients
		P = setPiece(P, pieceID, -f3, f1);  
		pieceID++;       		  	

		// Set piece 2 coefficients
		P = setPiece(P, pieceID, 0, 0);        
		pieceID++;    					      	

		// Set piece 3 coefficients          
		P = setPiece(P, pieceID, -f3, f2);                         
		pieceID++;   

		// Set piece 4 coefficients
		P = setPiece(P, pieceID, (double)(4*delta2*(c-mu2)), 0);   
		pieceID++;       

		// set limits
		// slope_OS is -infinity
		P.limits[0] = -infinity;  
		P.limits[1] = slope_OS;
		P.limits[2] = slope_OP; 
		P.limits[3] = slope_OQ;
		P.limits[4] = slope_OR;
		P.limits[5] = infinity;

	}   
	return P;	
}

// 3b) Contained Intervals
piecewise z_density_uniform :: containedB(double mu1, double delta1, double mu2, double delta2, double c) 
{	
	piecewise P;
	int pieceID;

	// mu2-delta2 <= c <= mu1-delta1
	if ((c >= (mu2-delta2)) && (c <= (mu1-delta1)))
	{		
			// Vertex Order : PQRS	

			// specify number of pieces
			P.numPieces = 5;
			pieceID = 0; 

			// set pieces

			// Set piece 0 coefficients       
			P = setPiece(P, pieceID, 0, (double)(4*delta1*(mu1-c)));                  
			pieceID++;          

			// Set piece 1 coefficients
			P = setPiece(P, pieceID, f3, -f2);  
			pieceID++;       		  	

			// Set piece 2 coefficients
			P = setPiece(P, pieceID, 0, 0);        
			pieceID++;    					      	

			// Set piece 3 coefficients          
			P = setPiece(P, pieceID, f4, -f2);                         
			pieceID++;   

			// Set piece 4 coefficients
			P = setPiece(P, pieceID, 0, (double)(4*delta1*(mu1-c)));
			pieceID++;       

			// set limits
			// slope_OS is -infinity
			P.limits[0] = -infinity;  
			P.limits[1] = slope_OP;
			P.limits[2] = slope_OQ; 
			P.limits[3] = slope_OR;
			P.limits[4] = slope_OS;
			P.limits[5] = infinity;	
	}      

	// mu1-delta1 < c <= mu1+delta1
	else if ((c > (mu1-delta1)) && (c <= (mu1+delta1)))
	{
		if ((slope_OP >= slope_OR) && (slope_OQ >= slope_OS))
		{

			// Vertex Order : RPSQ	

			// specify number of pieces
			P.numPieces = 5;
			pieceID = 0; 

			// set pieces

			// Set piece 0 coefficients             
			P = setPiece(P, pieceID, 0, (double)(2*((c-mu1)*(c-mu1)+delta1*delta1)));             
			pieceID++;          

			// Set piece 1 coefficients
			P = setPiece(P, pieceID, f4, f1);         
			pieceID++;       		  	

			// Set piece 2 coefficients
			P = setPiece(P, pieceID, (double)(2*((c-mu2)*(c-mu2)+delta2*delta2)),0);    
			pieceID++;    					      	

			// Set piece 3 coefficients                        
			P = setPiece(P, pieceID, f3, f1);    
			pieceID++;   

			// Set piece 4 coefficients
			P = setPiece(P, pieceID, 0, (double)(2*((c-mu1)*(c-mu1)+delta1*delta1)));               
			pieceID++;     


			// set limits
			// slope_OS is -infinity
			P.limits[0] = -infinity;  
			P.limits[1] = slope_OR;
			P.limits[2] = slope_OP; 
			P.limits[3] = slope_OS;
			P.limits[4] = slope_OQ;
			P.limits[5] = infinity;
		}

		else if ((slope_OP <= slope_OR) && (slope_OQ >= slope_OS))
		{

			// Vertex Order : PRSQ	

			// specify number of pieces
			P.numPieces = 5;
			pieceID = 0; 

			// set pieces

			// Set piece 0 coefficients             
			P = setPiece(P, pieceID, 0, (double)(2*((c-mu1)*(c-mu1)+delta1*delta1)));             
			pieceID++;          

			// Set piece 1 coefficients
			P = setPiece(P, pieceID, f3, f2);         
			pieceID++;       		  	

			// Set piece 2 coefficients
			P = setPiece(P, pieceID, (double)(2*((c-mu2)*(c-mu2)+delta2*delta2)),0);    
			pieceID++;    					      	

			// Set piece 3 coefficients                        
			P = setPiece(P, pieceID, f3, f1);    
			pieceID++;   

			// Set piece 4 coefficients
			P = setPiece(P, pieceID, 0, (double)(2*((c-mu1)*(c-mu1)+delta1*delta1)));               
			pieceID++;     
	
			// set limits
			// slope_OS is -infinity
			P.limits[0] = -infinity;  
			P.limits[1] = slope_OP;
			P.limits[2] = slope_OR; 
			P.limits[3] = slope_OS;
			P.limits[4] = slope_OQ;
			P.limits[5] = infinity;
		}

		else if ((slope_OP >= slope_OR) && (slope_OQ <= slope_OS)) 
		{

			// Vertex Order : RPQS	

			// specify number of pieces
			P.numPieces = 5;
			pieceID = 0; 

			// set pieces

			// Set piece 0 coefficients             
			P = setPiece(P, pieceID, 0, (double)(2*((c-mu1)*(c-mu1)+delta1*delta1)));             
			pieceID++;          

			// Set piece 1 coefficients
			P = setPiece(P, pieceID, f4, f1);         
			pieceID++;       		  	

			// Set piece 2 coefficients
			P = setPiece(P, pieceID, (double)(2*((c-mu2)*(c-mu2)+delta2*delta2)),0);    
			pieceID++;    					      	

			// Set piece 3 coefficients                        
			P = setPiece(P, pieceID, f4, f2);    
			pieceID++;   

			// Set piece 4 coefficients
			P = setPiece(P, pieceID, 0, (double)(2*((c-mu1)*(c-mu1)+delta1*delta1)));               
			pieceID++;     

			// set limits
			// slope_OS is -infinity
			P.limits[0] = -infinity;  
			P.limits[1] = slope_OR;
			P.limits[2] = slope_OP; 
			P.limits[3] = slope_OQ;
			P.limits[4] = slope_OS;
			P.limits[5] = infinity;
		}

		else if ((slope_OP <= slope_OR) && (slope_OQ <= slope_OS)) 
		{

			// Vertex Order : PRQS	

			// specify number of pieces
			P.numPieces = 5;
			pieceID = 0; 

			// set pieces

			// Set piece 0 coefficients             
			P = setPiece(P, pieceID, 0, (double)(2*((c-mu1)*(c-mu1)+delta1*delta1)));             
			pieceID++;          

			// Set piece 1 coefficients
			P = setPiece(P, pieceID, f3, f2);         
			pieceID++;       		  	

			// Set piece 2 coefficients
			P = setPiece(P, pieceID, (double)(2*((c-mu2)*(c-mu2)+delta2*delta2)),0);    
			pieceID++;    					      	

			// Set piece 3 coefficients                        
			P = setPiece(P, pieceID, f4, f2);    
			pieceID++;   

			// Set piece 4 coefficients
			P = setPiece(P, pieceID, 0, (double)(2*((c-mu1)*(c-mu1)+delta1*delta1)));               
			pieceID++;     
	
			// set limits
			// slope_OS is -infinity
			P.limits[0] = -infinity;  
			P.limits[1] = slope_OP;
			P.limits[2] = slope_OR; 
			P.limits[3] = slope_OQ;
			P.limits[4] = slope_OS;
			P.limits[5] = infinity;     						                              	
		}
	}       

	// mu1+delta1 < c <= mu2+delta2
	else if ((c > (mu1+delta1)) && (c <= mu2+delta2))
	{	
			// Vertex Order : RSPQ	

			// specify number of pieces
			P.numPieces = 5;
			pieceID = 0; 

			// set pieces

			// Set piece 0 coefficients       
			P = setPiece(P, pieceID, 0, (double)(4*delta1*(c-mu1)));                  
			pieceID++;          

			// Set piece 1 coefficients
			P = setPiece(P, pieceID, f4, -f1);  
			pieceID++;       		  	

			// Set piece 2 coefficients
			P = setPiece(P, pieceID, 0, 0);        
			pieceID++;    					      	

			// Set piece 3 coefficients          
			P = setPiece(P, pieceID, f3, -f1);                         
			pieceID++;   

			// Set piece 4 coefficients
			P = setPiece(P, pieceID, 0, (double)(4*delta1*(c-mu1)));
			pieceID++;       

			// set limits
			// slope_OS is -infinity
			P.limits[0] = -infinity;  
			P.limits[1] = slope_OR;
			P.limits[2] = slope_OS; 
			P.limits[3] = slope_OP;
			P.limits[4] = slope_OQ;
			P.limits[5] = infinity;
	}	 	   
	return P;	
}

// Get Distribution of Z = (c-X1)/(X2-X1) assuming endpoints have independent uniform noise
piecewise z_density_uniform::alpha_pdf(double mu1, double delta1, double mu2, double delta2, double c)
{
	piecewise P;         
	
	piecewise defaultP;
	
	// set default piece
	defaultP.numPieces = 1;
	defaultP.pc[0].k1 = 0;
	defaultP.pc[0].k2 = 0;
  defaultP.pc[0].k3 = 10;
  defaultP.limits[0] = 0;
  defaultP.limits[1] = 1;
    
	// Assume mu2 is always greater than mu1
	if (mu1 > mu2)
	{
		mean_temp = mu1;
		mu1 = mu2;
		mu2 = mean_temp;

		delta_temp = delta1;
		delta1 = delta2;
		delta2 = delta_temp;
	}   
	
	// If c is out of range [mu1-delta1,mu2+delta2] return default piece
	
	//if((c < mu1-delta1) || (c > mu2+delta2))
		//	return defaultP;

	// Major fix
	double rangemin, rangemax;
	if((mu1-delta1) < (mu2-delta2))
		rangemin = (mu1-delta1);
	else
		rangemin = (mu2-delta2);
	
	if((mu1+delta1) > (mu2+delta2))
		rangemax = (mu1+delta1);
	else
		rangemax = (mu2+delta2);
	
	if((c < rangemin) || (c > rangemax))
			return defaultP;

	// Precomputation of slopes and variables which are part of final PDF    
	if ((mu2 - mu1 + delta2 - delta1) != 0)
		slope_OP = (double)((c - mu1 - delta1) / (mu2 - mu1 + delta2 - delta1));

	if ((mu2 - mu1 + delta2 + delta1) != 0) 
		slope_OQ = (double)((c - mu1 + delta1) / (mu2 - mu1 + delta2 + delta1)); 

	if ((mu2 - mu1 - delta2 - delta1) != 0)
		slope_OS = (double)((c - mu1 - delta1) / (mu2 - mu1 - delta2 - delta1));		

	if ((mu2 - mu1 + delta1 - delta2) != 0)
		slope_OR = (double)((c - mu1 + delta1) / (mu2 - mu1 + delta1 - delta2));     	


	// Pre-define constants
	
	f1 = (mu1 + delta1 -c)*(mu1 + delta1 -c);
	f2 = (mu1 - delta1 -c)*(mu1 - delta1 -c);
	f3 = (mu2 + delta2 -c)*(mu2 + delta2 -c);
	f4 = (mu2 - delta2 -c)*(mu2 - delta2 -c);
	denominator = 8*delta1*delta2;


	// 1) Non-overlapping Intervals
	if ((mu2 - delta2) >= (mu1 + delta1))
	{            
		// latest bug fix (Sept 29, 2013)
		if((mu2-delta2) == (mu1+delta1))
		{
				slope_OS = infinity;					
		}		
		P = nonOverlapping(mu1, delta1, mu2, delta2, c); 
	} 

	// 2) Overlapping Intervals
	else if (((mu2-delta2) <= (mu1+delta1)) && ((mu2-delta2) > (mu1-delta1)) && ((mu2+delta2) > (mu1+delta1)))	
	{
		P = overlapping(mu1, delta1, mu2, delta2, c);          
	}	

	// 3a) Contained Intervals
	else if (((mu2-delta2) < (mu1+delta1)) && ((mu2-delta2) >= (mu1-delta1)) && ((mu2+delta2) <= (mu1+delta1)) && ((mu2+delta2) > (mu1-delta1))) 	     	
	{
		if((mu2+delta2) == (mu1+delta1))
		{
			if(c < mu1+delta1)
				slope_OP = infinity;					
		}		 
		P = containedA(mu1, delta1, mu2, delta2, c);          
	}	

	// 3b) Contained Intervals
	else if ((mu2-delta2) <= (mu1-delta1))
	{
		if((mu2-delta2) == (mu1-delta1))
		{
			slope_OR = -infinity;					
		}		

		P = containedB(mu1, delta1, mu2, delta2, c);            
	}		 

	return P; 	
} 

// Takes a piecewise function defined over 0 to 1 and adjusts its limits if mu1 > mu2. 
// After adjusting the piece limits, it is possible that quantities expected value, crossing probability
// and variance are not possible to compute in closed form.
piecewise z_density_uniform :: adjustPieceLimits(piecewise p)
{
	piecewise temp;	
	
	int pieceCount = p.numPieces;

	int limitCount = pieceCount + 1;

	double* limitDiff = new double[pieceCount];	

	// store difference between consecutive limits i.e., if limits are 0, 0.45, 0.8, 1 store difference as 0.45, 0.35, 0.2
	for(int i=0; i<limitCount; i++)

		limitDiff[i] = p.limits[i+1] - p.limits[i];

	// Now store pieces in reverse order and limits accordingly in a temporary variable e.g., if  p1, p2, p3 fall in 0, 0.45, 0.8, 1, now store
	// p3, p2, p1 in the limits 0, 0.2, 0.55, 1

	for(int i=0; i<pieceCount; i++)
	{		
		temp.pc[i] = p.pc[pieceCount-i-1];
	}

	for(int i=0; i<limitCount; i++)
	{
		if(i==0)
			temp.limits[i] = 0;
		else
			temp.limits[i] = temp.limits[i-1] + limitDiff[limitCount-i-1];		

	}
	// modify piece order and limits of original picewise function	
	for(int i=0; i<pieceCount; i++)
	{
		p.pc[i] = temp.pc[i];
	}

	for(int i=0; i<limitCount; i++)
	{
		p.limits[i] = temp.limits[i];
	}

	return p;
}

// set number of uniform distributions in kde1 
void z_density_uniform :: setNumUniformsInKde1(int a)
{
	numUniformsInKde1 = a;
}

// get number of uniform distributions in kde2 
int z_density_uniform :: getNumUniformsInKde1()
{
	return numUniformsInKde1;
}

// set number of uniform distributions in kde1 
void z_density_uniform :: setNumUniformsInKde2(int a)
{
	numUniformsInKde2 = a;
}

// get number of uniform distributions in kde2 
int z_density_uniform :: getNumUniformsInKde2()
{
	return numUniformsInKde2;
}

// add pieces in the range [low high]	
piece z_density_uniform :: add(piecewise P1, piecewise P2, double low, double high)
{
	piece pc1, pc2; 

	// store addition of pc1 and pc2
	piece pc3;

	for(int i=0; i< (P1.numPieces + 1); i++)
	{
		if((P1.limits[i]<=low) && (high<=P1.limits[i+1]))
		{
			pc1 = P1.pc[i];
			break;
		}
	}

	for(int i=0; i< (P2.numPieces + 1); i++)
	{
		if((P2.limits[i]<=low) && (high<=P2.limits[i+1]))
		{
			pc2 = P2.pc[i];
			break;
		}
	}

	// (pc1.k1 * z^2 + pc1.k2 * (1-z)^2)/(pc1.k3 * z^2* (1-z)^2) + (pc2.k1 * z^2 + pc2.k2 * (1-z)^2)/(pc2.k3 * z^2* (1-z)^2)

	// Normalize the density functuion
	/*double expected, c_prob, sig, second, first;

	Compute0To1(P1,&expected,&c_prob,&sig,&second,&first);

	pc1.k1 = (double)(pc1.k1/c_prob);
	pc1.k2 = (double)(pc1.k2/c_prob);

	Compute0To1(P2,&expected,&c_prob,&sig,&second,&first);

	pc2.k1 = (double)(pc2.k1/c_prob);
	pc2.k2 = (double)(pc2.k2/c_prob);*/
	

	// Calculate product of the denominators
	pc3.k3 = pc1.k3*pc2.k3;

	pc3.k1 = pc1.k1*pc2.k3 + pc2.k1*pc1.k3;

	pc3.k2 = pc1.k2*pc2.k3 + pc2.k2*pc1.k3;

	return pc3;

}

// Adds 2 piecewise functions   
piecewise z_density_uniform :: addTwoPiecewiseFunctions(piecewise P1, piecewise P2)
{
	// To store addition of two piecewise functions
	piecewise P3;

	// Number of different limits after addition
	int p1LimitCount, p2LimitCount;
	p1LimitCount = P1.numPieces + 1; 
	p2LimitCount = P2.numPieces + 1;	
	// move over array limits as if merging two arrays and while going add piecewise functions.

	int i=0,j=0;

	// To store limits of new piecewise function
	int k=0;

	double temp1,temp2;
	double epsilon = 1e-10;
	//cout<< "Limit order:";
	while((i< p1LimitCount) || (j< p2LimitCount))
	{
		// temp1 and temp2 to store immediate smallest limits among two arrays.
		// if its not very start of the loop
		if((i!=0) || (j!=0))
			temp1 = temp2;

		//temp2 keeps on storing the limits in increasing order
		if ((i< p1LimitCount) && (j< p2LimitCount))
		{
			if (P1.limits[i] < P2.limits[j])
			{
				temp2 = P1.limits[i];
				i++;
			}
			else
			{
				temp2 = P2.limits[j];	
				j++;
			}
			//n++;
		}
		else if (i < p1LimitCount)
		{
			temp2 = P1.limits[i];
			i++; 
			//n++;
		}

		else if (j < p2LimitCount)
		{
			temp2 = P2.limits[j];
			j++; 
			//n++;
		}
		//cout<<temp2<<" ";	

		// At very start
		if(k==0)
		{
			P3.limits[k] = temp2;
			k++;		
		}

		//if((i!=0) || (j!=0))
		else
		{
			if(temp1 < temp2)
			{
				//Add pieces in the range [temp1 temp2]
				P3.pc[k-1] = add(P1,P2,temp1,temp2);
				P3.limits[k] = temp2;
				k++;						
				// Add pieces in the range temp1 and temp2 and also store the piece limits	
			}
		}		
	}
	P3.numPieces = k-1;
	cout<<"\n";
	/*cout<<"Limits of piecewise addition:";
	for(int m=0; m<k;m++)
		cout<<P3.limits[m]<<" ";
	cout<<"\n";*/

	/*for(int m=0; m<=P1.numPieces;m++)
		cout<<P1.limits[m]<<" ";
	cout<<"\n";

	for(int m=0; m<P1.numPieces;m++)
		cout<<"k1:"<<P1.pc[m].k1<<" "<<"k2:"<<P1.pc[m].k2<<" "<<"k3:"<<P1.pc[m].k3<<"\n";
	cout<<"\n\n";

	for(int m=0; m<=P2.numPieces;m++)
		cout<<P2.limits[m]<<" ";
	cout<<"\n";

	for(int m=0; m<P2.numPieces;m++)
		cout<<"k1:"<<P2.pc[m].k1<<" "<<"k2:"<<P2.pc[m].k2<<" "<<"k3:"<<P2.pc[m].k3<<"\n";
	cout<<"\n\n";

	for(int m=0; m<=P3.numPieces;m++)
		cout<<P3.limits[m]<<" ";
	cout<<"\n";

	for(int m=0; m<P3.numPieces;m++)
		cout<<"k1:"<<P3.pc[m].k1<<" "<<"k2:"<<P3.pc[m].k2<<" "<<"k3:"<<P3.pc[m].k3<<"\n";
	cout<<"\n\n";*/

	return P3;

}


// Takes in distributions for individual uniform pairs and returns piecewise function for kde  
piecewise z_density_uniform :: getKdePdf(piecewise* P)
{
	int numKde1 = getNumUniformsInKde1();
	int numKde2 = getNumUniformsInKde2();
	int numPiecewiseFun = numKde1*numKde2;

	int totalLimits = 0;

	// count total number of limits from all piecewise functions
	for(int i=0; i<numPiecewiseFun; i++)
	{
		totalLimits += P[i].numPieces + 1;
	}
	
	// Create array to store limits from all piecewise functions
	double* lim = new double[totalLimits];

	int cntr = 0;

	for(int i=0; i<numPiecewiseFun; i++)
	{
		for(int j=0; j<=P[i].numPieces; j++)
		{
			lim[cntr] = P[i].limits[j];
			cntr++;
		}
	}

	cout<<"Limits array:"<<"\n";
	for(int i=0; i<totalLimits; i++)
	{
		cout<<lim[i]<<" ";
	}
	cout<<"\n";

	// sort the limits array
	std::sort(lim, lim + totalLimits);

	cout<<"Limits array:"<<"\n";
	for(int i=0; i<totalLimits; i++)
	{
		cout<<lim[i]<<" ";
	}
	cout<<"\n";

	// Remove duplicates from the sorted array
	int duplicateCount = 0;
	double epsilon = 1e-10;
	for(int i=0; i<totalLimits-1; i++)
	{
		if((lim[i+1] - lim[i])<epsilon)
		{
			duplicateCount++;
		}
	}

	cout<<"Duplicate Count:"<<duplicateCount<<"\n";

	// Create array to store sorted limits without duplicates
	double* limSorted = new double[totalLimits-duplicateCount];
	int ctr = 0;
	int ptr = 1;
	for(int i=0; i<totalLimits-1; i++)
	{
		if((lim[ptr] - lim[i])<epsilon)
		{	
			//If pointer reaches end copy the last element
			if(ptr == totalLimits-1)
			{
				limSorted[ctr] = lim[ptr];
				ctr++;
			}	
			ptr++;
		}
		else
		{
			limSorted[ctr] = lim[i];
			ctr++;
			//If pointer reaches end copy the both last element
			if(ptr == totalLimits-1)
			{
				limSorted[ctr] = lim[ptr];
				ctr++;
			}	
			ptr++;	
		}
	}

	cout<<"Sorted limits array without duplicates:"<<"\n";
	for(int i=0; i<ctr; i++)
	{
		cout<<limSorted[i]<<" ";
	}
	cout<<"\n";


	// Create a new piecewise function representing uncertainty in the linear interpolation assuming 
	// data at the end points is sampled from the kernel-based density.
	piecewise kdeout;
	kdeout.numPieces = ctr;
	for(int j=0; j<=kdeout.numPieces; j++)
	{
		kdeout.limits[j] = limSorted[j];
		cntr++;
	} 
}

// Piecewise function returned assuming data is sampled from a kernel density estimation
// Arrays of mu1, delta1, mu2, delta2 are passed. cardinality of mu1 and delta1 is same (similar for mu2 and delta2).
piecewise z_density_uniform :: kde_alpha_pdf(float* mu1, float* delta1, float* mu2, float* delta2, double c)
{
	int numKde1 = getNumUniformsInKde1();
	int numKde2 = getNumUniformsInKde2();

	cout<<"\nNumkde1:"<<numKde1;
	cout<<"\nNumkde2:"<<numKde2;
	// Store Piecewise density functions for each pair of uniform distributions in an array
	piecewise* P = new piecewise[numKde1*numKde2];

	// Store information whether mu1 > mu2
	int* flag = new int[numKde1*numKde2];


	piecewise temp;

	piecewise pieceAdd;

	for(int i=0; i<numKde1; i++)
	{
		for(int j=0; j<numKde2; j++)
		{
			temp = alpha_pdf(mu1[i], delta1[i], mu2[j], delta2[j], c);
	
			temp = getPdfOver0To1(temp);
		
			if(mu1[i] > mu2[j])
			{
				//temp = adjustPieceLimits(temp);
				flag[i*numKde1 + j] = 1;
			}
			else 
			{
				flag[i*numKde1 + j] = 0;
			}

			P[i*numKde1 + j] = temp;

			/*cout<<"\n"<<P[i*numKde1 + j].numPieces<<"\n";
			for(int k=0;k<=P[i*numKde1 + j].numPieces;k++)
				cout<<P[i*numKde1 + j].limits[k]<<" ";
			cout<<"\n";*/
			
		}
	}

	


/*	for(int j=0; j<=pieceAdd.numPieces; j++)
	{
		pieceAdd.pc[j].k1 = (double)(pieceAdd.pc[j].k1);
		pieceAdd.pc[j].k2 = (double)(pieceAdd.pc[j].k2);	
				
	}*/
/*	double expected, c_prob, sig, second, first;

	Compute0To1(P[0],&expected,&c_prob,&sig,&second,&first);

	for(int j=0; j<=P[0].numPieces; j++)
	{
		P[0].pc[j].k1 = (double)(P[0].pc[j].k1/c_prob);
		P[0].pc[j].k2 = (double)(P[0].pc[j].k2/c_prob);					
	}

	if(flag[0])
		P[0] = adjustPieceLimits(P[0]);

	cout<<"Expected Value is:"<<first<<"\n";
	cout<<"Crossing Probability is:"<<c_prob<<"\n";
	cout<<"Variance is:"<<sig<<"\n";

	Compute0To1(P[1],&expected,&c_prob,&sig,&second,&first);

	for(int j=0; j<=P[1].numPieces; j++)
	{
		P[1].pc[j].k1 = (double)(P[1].pc[j].k1/c_prob);
		P[1].pc[j].k2 = (double)(P[1].pc[j].k2/c_prob);					
	}

	if(flag[1])
		P[1] = adjustPieceLimits(P[1]);

	cout<<"Expected Value is:"<<first<<"\n";
	cout<<"Crossing Probability is:"<<c_prob<<"\n";
	cout<<"Variance is:"<<sig<<"\n";

	Compute0To1(P[2],&expected,&c_prob,&sig,&second,&first);

	for(int j=0; j<=P[2].numPieces; j++)
	{
		P[2].pc[j].k1 = (double)(P[2].pc[j].k1/c_prob);
		P[2].pc[j].k2 = (double)(P[2].pc[j].k2/c_prob);					
	}

	if(flag[2])
		P[2] = adjustPieceLimits(P[2]);

	cout<<"Expected Value is:"<<first<<"\n";
	cout<<"Crossing Probability is:"<<c_prob<<"\n";
	cout<<"Variance is:"<<sig<<"\n";

	Compute0To1(P[3],&expected,&c_prob,&sig,&second,&first);
	for(int j=0; j<=P[3].numPieces; j++)
	{
		P[3].pc[j].k1 = (double)(P[3].pc[j].k1/c_prob);
		P[3].pc[j].k2 = (double)(P[3].pc[j].k2/c_prob);					
	}

	if(flag[3])
		P[3] = adjustPieceLimits(P[3]);

	cout<<"Expected Value is:"<<first<<"\n";
	cout<<"Crossing Probability is:"<<c_prob<<"\n";
	cout<<"Variance is:"<<sig<<"\n";*/

	
	pieceAdd = addTwoPiecewiseFunctions(P[0],P[1]);
	//piecewise kdeAnswer = getKdePdf(P);

	pieceAdd = addTwoPiecewiseFunctions(pieceAdd,P[2]);
	pieceAdd = addTwoPiecewiseFunctions(pieceAdd,P[3]);
  double expected, c_prob, sig, second, first;
	Compute0To1(pieceAdd,&expected,&c_prob,&sig,&second,&first);

	cout<<"Expected Value is:"<<(double)(expected/c_prob)<<"\n";
	cout<<"Crossing Probability is:"<<c_prob<<"\n";
	cout<<"Variance is:"<<sig<<"\n\n\n";

	//return P1;
}

// Piecewise function returned assuming data is sampled from a kernel density estimation
// h1 and h2 are bandwidths or variances for each of samples in mu1 and mu2. 
/*double z_density_uniform :: kde_alpha_pdf_expected(float* mu1, double h1, float* mu2, double h2, double c)
{
	int numKde1 = getNumUniformsInKde1();
	int numKde2 = getNumUniformsInKde2();

	piecewise temp;

	double expected, c_prob, sig, second, first;

	double kde_expected = 0;

	double EPSILON = 0.000000001;

	for(int i=0; i<numKde1; i++)
	{
		for(int j=0; j<numKde2; j++)
		{
			temp = alpha_pdf(mu1[i], h1, mu2[j], h2, c);
	
			Compute0To1(temp,&expected,&c_prob,&sig,&second,&first);

			//cout<<"c_prob:"<<c_prob<<"\n";
		
			// What happens when crossing probability is 0
			if(mu1[i] > mu2[j])
			{
				if(fabs(c_prob)<EPSILON)
				{
					kde_expected = kde_expected + 0.5;
				}
				else
				{
					kde_expected = kde_expected + (1 - (double)(expected/c_prob));					
				}
			}	
			else
			{
				if(fabs(c_prob)<EPSILON)
				{
					kde_expected = kde_expected + 0.5;
				}
				else
				{
					kde_expected =  kde_expected + (double)(expected/c_prob);
				}
			}
			
		}
	}

	kde_expected = (double)(kde_expected/(numKde1*numKde2));

	cout<<"kde_expected:"<<kde_expected<<"\n";
	
	return kde_expected;
	
}*/


double z_density_uniform :: kde_alpha_pdf_expected(float* mu1, double h1, float* mu2, double h2, double c)
{
	int numKde1 = getNumUniformsInKde1();
	int numKde2 = getNumUniformsInKde2();

	piecewise temp;

	double expected, c_prob, sig, second, first;

	double kde_expected = 0;

	double total_crossing_prob = 0;
	
	double EPSILON = 0.000000001;

	for(int i=0; i<numKde1; i++)
	{
		for(int j=0; j<numKde2; j++)
		{
			temp = alpha_pdf(mu1[i], h1, mu2[j], h2, c);
	
			Compute0To1(temp,&expected,&c_prob,&sig,&second,&first);

			//cout<<"c_prob:"<<c_prob<<"\n";
		
			// What happens when crossing probability is 0 skip
			if(mu1[i] > mu2[j])
			{
				if(fabs(c_prob)<EPSILON)
				{
					// Means given pair of kernels doesn't contribute to the expected value
				}
				else
				{
					//kde_expected = kde_expected + (1 - expected); If pdf is not normalized 1 should be replace with crossing probability
					// Major update
					kde_expected = kde_expected + (c_prob - expected);
					total_crossing_prob = total_crossing_prob + c_prob;
				}
			}	
			else
			{
				if(fabs(c_prob)<EPSILON)
				{
					// Means given pair of kernels doesn't contribute to the expected value
				}
				else
				{
					kde_expected =  kde_expected + expected;
					total_crossing_prob = total_crossing_prob + c_prob;
				}
			}
			
		}
	}

  if(total_crossing_prob > 0)
  
		kde_expected = (double)(kde_expected/total_crossing_prob);

	else
		
		kde_expected = 0.5;
		
	//cout<<"kde_expected:"<<kde_expected<<"\n";
	
	return kde_expected;
}

// Piecewise function returned assuming data is sampled from a kernel density estimation
// h1 and h2 are bandwidths or variances for each of samples in mu1 and mu2. 
/*double z_density_uniform :: kde_alpha_pdf_variance(float* mu1, double h1, float* mu2, double h2, double c)
{

	int numKde1 = getNumUniformsInKde1();
	int numKde2 = getNumUniformsInKde2();
	
	

	piecewise temp;

	double expected, c_prob, sig, second, first;

	double kde_second_moment = 0;
		

	double kde_first_moment = kde_alpha_pdf_expected(mu1, h1, mu2, h2, c);

	
	
	double kde_variance = 0;

	// kde second moment can be obtained in following two ways when mu1 > mu2.

	// reverse the pdf using function adjustPieceLimits and compute second moment in the forward direction (may be closed form compution is not possible because of things like log of negative number etc.)

        // or


        // use the following formula with current pdf instead of reversing the pdf
	// second moment in reverse direction = second moment in forward direction + 1 - 2*first moment in forward direction


	double EPSILON = 0.000000001;

	for(int i=0; i<numKde1; i++)
	{
		for(int j=0; j<numKde2; j++)
		{
			temp = alpha_pdf(mu1[i], h1, mu2[j], h2, c);

			Compute0To1(temp,&expected,&c_prob,&sig,&second,&first);
		
			if(mu1[i] > mu2[j])
			{
					if(fabs(c_prob)<EPSILON)
					{
						kde_second_moment = kde_second_moment + 0.1;
					}
					else
					{
						kde_second_moment = kde_second_moment + (double)(second/c_prob) + 1 - 2*(double)(expected/c_prob);
					}
			}	
			else
			{
					if(fabs(c_prob)<EPSILON)
					{
						kde_second_moment = kde_second_moment + 0.1;
					}
					else
					{
						kde_second_moment =  kde_second_moment + (double)(second/c_prob);
					}
			}
			
		}
	}

	kde_variance  = (double)(kde_second_moment/(numKde1*numKde2)) - kde_first_moment*kde_first_moment;
	
	cout<<"kde_variance:"<<kde_variance<<"\n";

	return kde_variance;
}*/


double z_density_uniform :: kde_alpha_pdf_variance(float* mu1, double h1, float* mu2, double h2, double c)
{

	int numKde1 = getNumUniformsInKde1();
	int numKde2 = getNumUniformsInKde2();
	
	

	piecewise temp;

	double expected, c_prob, sig, second, first;

	double kde_second_moment = 0;
	double total_crossing_prob = 0;	

	double kde_first_moment = kde_alpha_pdf_expected(mu1, h1, mu2, h2, c);

	
	
	double kde_variance = 0;

	// kde second moment can be obtained in following two ways when mu1 > mu2.

	// reverse the pdf using function adjustPieceLimits and compute second moment in the forward direction (may be closed form compution is not possible because of things like log of negative number etc.)

        // or


        // use the following formula with current pdf instead of reversing the pdf
	// second moment in reverse direction = second moment in forward direction + 1 - 2*first moment in forward direction


	double EPSILON = 0.000000001;

	for(int i=0; i<numKde1; i++)
	{
		for(int j=0; j<numKde2; j++)
		{
			temp = alpha_pdf(mu1[i], h1, mu2[j], h2, c);

			Compute0To1(temp,&expected,&c_prob,&sig,&second,&first);			
		
			if(mu1[i] > mu2[j])
			{
					if(fabs(c_prob)<EPSILON)
					{
						// kde_second_moment = kde_second_moment + 0.1;
						// Means given pair of kernels doesn't determine variance
					}
					else
					{
						//If pdf is not normalized 1 should be replace with crossing probability	
						kde_second_moment = kde_second_moment + second + c_prob - 2*expected;
						total_crossing_prob = total_crossing_prob + c_prob;
					}
			}	
			else
			{
					if(fabs(c_prob)<EPSILON)
					{
						// kde_second_moment = kde_second_moment + 0.1;
						// Means given pair of kernels doesn't determine variance
					}
					else
					{
						kde_second_moment =  kde_second_moment + second;
						total_crossing_prob = total_crossing_prob + c_prob;
					}
			}
			
		}
	}

	if(total_crossing_prob > 0)
	
		kde_variance  = (double)(kde_second_moment/total_crossing_prob) - kde_first_moment*kde_first_moment;

	else
		
		kde_variance = 0;	
	//cout<<"kde_variance:"<<kde_variance<<"\n";

	return kde_variance;
}

/*int main()
{
	z_density_uniform* d = new z_density_uniform();
	piecewise p1;
	double expected, c_prob, sig, second, first;
	p1 = d->alpha_pdf(10,2,30,4,30);
	for(int i=0;i<p1.numPieces;i++)
	{
		cout<<"limits:"<<p1.limits[i]<<"-"<<p1.limits[i+1];
		cout<<" k1:"<<p1.pc[i].k1<<" k2:"<<p1.pc[i].k2<<" k3:"<<p1.pc[i].k3<<"\n";
	}
	cout<<"\n\n\n";
	
	d->Compute0To1(p1,&expected,&c_prob,&sig,&second,&first);
	cout<<"Expected Value is:"<<expected<<"\n";
	cout<<"Crossing Probability is:"<<c_prob<<"\n";
	cout<<"Variance is:"<<sig<<"\n";

	//mu1: 2.85714 delta1: 57.1429 mu2: 11.4286 delta2: 114.286
	
	// ratio distribution assuming data is uniformly distributed
	z_density_uniform* d = new z_density_uniform();

	// ratio distribution assuming data is sampled from kde
	z_density_uniform* kde_d = new z_density_uniform();
	
	// Allocate memory for mu and delta arrays
	float* mu1 = new float[10];
	float* mu2 = new float[10];
	float* delta1 = new float[10];
	float* delta2 = new float[10];

	// kernel 1
	mu2[0] = 10;
	delta2[0] = 2;
	mu2[1] = 9;
	delta2[1] = 2;*/

	/*mu1[0] = 51.2044;
	delta1[0] = 3.16207;
	mu1[1] = 56.6435;
	delta1[1] = 3.16207;
	mu1[2] = 57.0106;
	delta1[2] = 3.16207;
	mu1[3] = 62.8328;
	delta1[3] = 3.16207;
	mu1[4] = 52.6929;
	delta1[4] = 3.16207;
	mu1[5] = 58.9998;
	delta1[5] = 3.16207;
	mu1[6] = 62.0882;
	delta1[6] = 3.16207;
	mu1[7] = 51.0841;
	delta1[7] = 3.16207;
	mu1[8] = 51.5606;
	delta1[8] = 3.16207;
	mu1[9] = 63.7621;
	delta1[9] = 3.16207;*/

	// kernel 2
/*	mu1[0] = 20;
	delta1[0] = 3;
	mu1[1] = 27;
	delta1[1] = 3;*/

	/*mu2[0] = 47.5542;
	delta2[0] = 2.05051;
	mu2[1] = 51.9143;
	delta2[1] = 2.05051;
	mu2[2] = 51.7099;
	delta2[2] = 2.05051;
	mu2[3] = 46.1688;
	delta2[3] = 2.05051;
	mu2[4] = 42.431;
	delta2[4] = 2.05051;
	mu2[5] = 45.839;
	delta2[5] = 2.05051;
	mu2[6] = 47.3078;
	delta2[6] = 2.05051;
	mu2[7] = 46.1065;
	delta2[7] = 2.05051;
	mu2[8] = 47.79 ;
	delta2[8] = 2.05051;
	mu2[9] = 52.7583;
	delta2[9] = 2.05051;

	piecewise p1, p2;

	// Assuming data is uniformly distributed
	p1 = d->alpha_pdf(9,2,27,3,20);
	p2 = d->getPdfOver0To1(p1);

	//for(int i=0;i<p2.numPieces;i++)
	//	cout<<" k1:"<<p2.pc[i].k1<<" k2:"<<p2.pc[i].k2<<" k3:"<<p2.pc[i].k3<<"\n";
	//cout<<"\n";
	

	double expected, crossProb, var, second, first;
	d->Compute0To1(p2, &expected, &crossProb, &var, &second, &first);
	cout<<"Expected Value is:"<<expected<<"\n";
	cout<<"Crossing Probability is:"<<crossProb<<"\n";
	cout<<"Variance is:"<<var<<"\n";*/

	//p1 = d->alpha_pdf(mu1[0], delta1[0], mu2[0], delta2[0], 10);
	/*cout<<p1.numPieces<<"\n";
	for(int i=0;i<=p1.numPieces;i++)
		cout<<p1.limits[i]<<" ";
	cout<<"\n";

	cout<<p2.numPieces<<"\n";
	for(int i=0;i<=p2.numPieces;i++)
		cout<<p2.limits[i]<<" ";
	cout<<"\n";*/

	/*piecewise p3, p4;
	// Number of uniform distribution pairs
	kde_d->setNumUniformsInKde1(2);
	kde_d->setNumUniformsInKde2(2);
	kde_d->kde_alpha_pdf(mu1, delta1, mu2, delta2, 15);
	// Assuming data is kde sampled
	double expected_value = kde_d->kde_alpha_pdf_expected(mu1, 3, mu2, 2, 15);
	cout<<"Expected value:"<<expected_value<<"\n";
	double variance_value = kde_d->kde_alpha_pdf_variance(mu1, 3, mu2, 2, 15);
	cout<<"Variance value:"<<variance_value<<"\n";
	/*cout<<"\n"<<p3.numPieces<<"\n";
	for(int i=0;i<=p3.numPieces;i++)
		cout<<p3.limits[i]<<" ";
	cout<<"\n";*/

	// Important:::probabilities or expected values can not be computed in closed form if mu1 > mu2.
	// However for such case derive these probabilities by swapping mu1 and mu2 and corresponding deltas
	// and then finding distribution as done in the vis 2013 paper.
	/*d->Compute0To1(p3, &expected, &crossProb, &var);
	cout<<"Expected Value is:"<<expected<<"\n";
	cout<<"Crossing Probability is:"<<crossProb<<"\n";
	cout<<"Variance is:"<<var<<"\n";*/

//	return 0;
//}
