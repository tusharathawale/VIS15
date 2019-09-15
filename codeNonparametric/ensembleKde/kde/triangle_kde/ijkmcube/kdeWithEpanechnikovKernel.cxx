#include "stdio.h"
#include <iostream>
#include "stdlib.h"
#include "math.h"
//#include "cmath.h"
#include <algorithm>
#include "kdeWithEpanechnikovKernel.h"
#include <cfloat>
#include <limits> 
#include "epanechnikov_kernel_polynomial.h"

//# define infinity 10000
//# define minus_infinity -10000
#define _USE_MATH_DEFINES

using namespace std;

double e_mu1, e_delta1, e_mu2, e_delta2, e_c;


/*

                      b
             r -  -  -  -  -  - q
            /        /        /  
           /   P3   /   P2   /
          /        /  	    /
        c  -  -  - e -  -  - a
        /         /         /
       /    P4   /   P1    /
      /         /         /
      s-  -  -  d- - - - - p

% P1, P2, P3, P4 are different polynomials.*/

//triangular_kernel_polynomial* tri_poly;

epanechnikov_kernel_polynomial* epan_poly;

double e_eps = 1.0e-4;

// How to determine constant 1e-4? It is highly sensitive and can make some results unstable
// Test a>=b
inline int e_ge(double a, double b)
{
	if(a>b)
		return 1;
	if(fabs(a-b) < e_eps)
		return 1;
	return 0;
}

// Test a<=b
inline int e_le(double a, double b)
{
	if(a<b)
		return 1;
	if(fabs(a-b) < e_eps)
		return 1;	
	return 0;
}

// Test a==b
inline int e_eq(double a, double b)
{
	if(fabs(a-b) < e_eps)
		return 1;
	return 0;
}

double e_integrate_pdf(double low, double high, int piecenum)
{

	// assuming low is always greater than 0 and  low<=high

	// bug1 : low >= 1 check was replaced with ((low > 1) || (fabs(low-1) < 1e-6))
	// Basically  every place where there was comparison check with >=, = was removed. if(a==1) doesn't work in c++.

	// low is always guaranteed to be greater than 0. High is always guaranteed to be greater than or equal to low.


	/* bug2 : else if ((low <= 1) && (high <= 1))
			return (tri_poly -> PAED_integrate_piece_value(high,piecenum) - tri_poly -> PAED_integrate_piece_value(low,piecenum));
		else if ((low <= 1) && (high >= 1))		
			return (tri_poly -> PAED_integrate_piece_value(1,piecenum) - tri_poly -> PAED_integrate_piece_value(low,piecenum));	

		This code snippet was replaced with

		else if ((low < 1) && (high < 1))
			return (tri_poly -> PAED_integrate_piece_value(high,piecenum) - tri_poly -> PAED_integrate_piece_value(low,piecenum));
		else if ((low <= 1) && (high >= 1))		
			return (tri_poly -> PAED_integrate_piece_value(1,piecenum) - tri_poly -> PAED_integrate_piece_value(low,piecenum));	

*/
	if(e_ge(high,1))
	{
		high = 1;
	}
	
	if(e_le(low,0))
	{
		low = 0;
	}

	if (e_ge(low,1))
		return 0;
	else if (high - low < e_eps)
		return 0;
		//	return (tri_poly -> approx_PAED_integrate_piece_value(high,piecenum) - tri_poly -> approx_PAED_integrate_piece_value(low,piecenum));			
	else		
		return (epan_poly -> PQRS_integrate_piece_value(high,piecenum) - epan_poly -> PQRS_integrate_piece_value(low,piecenum));
}

double e_expected_value(double low, double high, int piecenum)
{	

	if(e_ge(high,1))
	{
		high = 1;
	}

	if(e_le(low,0))
	{
		low = 0;
	}
	
	// assuming low is always greater than 0 and  low<=high
	if (e_ge(low,1))
		return 0;
	else if (high - low < e_eps)
		return 0;
	//	return tri_poly -> approx_PAED_expected_piece_value(high,piecenum) - tri_poly -> approx_PAED_expected_piece_value(low,piecenum);
	else		
		return epan_poly -> PQRS_expected_piece_value(high,piecenum) - epan_poly -> PQRS_expected_piece_value(low,piecenum);
}

double e_second_moment(double low, double high, int piecenum)
{	

	if(e_ge(high,1))
	{
		high = 1;
	}
	
	if(e_le(low,0))
	{
		low = 0;
	}

	// assuming low is always greater than 0 and  low<=high

	if (e_ge(low,1))
		return 0;
	else if (high - low < e_eps)
		return 0;		
		//return tri_poly -> PAED_e_second_moment_piece_value(high,piecenum) - tri_poly -> PAED_e_second_moment_piece_value(low,piecenum);
	else		
		return epan_poly -> PQRS_second_moment_piece_value(high,piecenum) - epan_poly -> PQRS_second_moment_piece_value(low,piecenum);
}




//double e_integrate_pdf(double low, double high, int piecenum)
//{

	// assuming low is always greater than 0 and  low<=high

	// bug1 : low >= 1 check was replaced with ((low > 1) || (fabs(low-1) < 1e-6))
	// Basically  every place where there was comparison check with >=, = was removed. if(a==1) doesn't work in c++.

	// low is always guaranteed to be greater than 0. High is always guaranteed to be greater than or equal to low.


	/* bug2 : else if ((low <= 1) && (high <= 1))
			return (tri_poly -> PAED_integrate_piece_value(high,piecenum) - tri_poly -> PAED_integrate_piece_value(low,piecenum));
		else if ((low <= 1) && (high >= 1))		
			return (tri_poly -> PAED_integrate_piece_value(1,piecenum) - tri_poly -> PAED_integrate_piece_value(low,piecenum));	

		This code snippet was replaced with

		else if ((low < 1) && (high < 1))
			return (tri_poly -> PAED_integrate_piece_value(high,piecenum) - tri_poly -> PAED_integrate_piece_value(low,piecenum));
		else if ((low <= 1) && (high >= 1))		
			return (tri_poly -> PAED_integrate_piece_value(1,piecenum) - tri_poly -> PAED_integrate_piece_value(low,piecenum));	

*/
/*	if(ge(high,1))
	{
		high = 1;
	}

	if (sub_parallelogram_num == 0)
	{
		if (ge(low,1))
			return 0;
		else if (high - low < eps)
			return 0;
		//	return (tri_poly -> approx_PAED_integrate_piece_value(high,piecenum) - tri_poly -> approx_PAED_integrate_piece_value(low,piecenum));			
		else		
			return (tri_poly -> PAED_integrate_piece_value(high,piecenum) - tri_poly -> PAED_integrate_piece_value(low,piecenum));
	}
	else if (sub_parallelogram_num == 1)
	{
		if (ge(low,1))
			return 0;
		else if (high - low < eps)
			return 0;		
		//	return (tri_poly -> approx_AQBE_integrate_piece_value(high,piecenum) - tri_poly -> approx_AQBE_integrate_piece_value(low,piecenum));		
		else		
			return (tri_poly -> AQBE_integrate_piece_value(high,piecenum) - tri_poly -> AQBE_integrate_piece_value(low,piecenum));
	}
	else if (sub_parallelogram_num == 2)
	{
		if (ge(low,1))
			return 0;
		else if (high - low < eps)
			return 0;		
		//	return (tri_poly -> approx_EBRC_integrate_piece_value(high,piecenum) - tri_poly -> approx_EBRC_integrate_piece_value(low,piecenum));		
		else		
			return (tri_poly -> EBRC_integrate_piece_value(high,piecenum) - tri_poly -> EBRC_integrate_piece_value(low,piecenum));
	}
	else if (sub_parallelogram_num == 3)
	{
		if (ge(low,1)) 
			return 0;	
		else if (high - low < eps)
			return 0;		
		//	return (tri_poly -> approx_DECS_integrate_piece_value(high,piecenum) - tri_poly -> approx_DECS_integrate_piece_value(low,piecenum));
		else		
			return (tri_poly -> DECS_integrate_piece_value(high,piecenum) - tri_poly -> DECS_integrate_piece_value(low,piecenum));
	}
}*/

/*double e_expected_value(double low, double high, int piecenum)
{	

	if(ge(high,1))
	{
		high = 1;
	}

	// assuming low is always greater than 0 and  low<=high
	if (sub_parallelogram_num == 0)
	{
		if (ge(low,1))
			return 0;
		else if (high - low < eps)
			return 0;
		//	return tri_poly -> approx_PAED_expected_piece_value(high,piecenum) - tri_poly -> approx_PAED_expected_piece_value(low,piecenum);
		else		
			return tri_poly -> PAED_expected_piece_value(high,piecenum) - tri_poly -> PAED_expected_piece_value(low,piecenum);
	}
	else if (sub_parallelogram_num == 1)
	{
		if (ge(low,1))
			return 0;
		else if (high - low < eps)
			return 0;
		//	return tri_poly -> approx_AQBE_expected_piece_value(high,piecenum) - tri_poly -> approx_AQBE_expected_piece_value(low,piecenum);
		else 		
			return tri_poly -> AQBE_expected_piece_value(high,piecenum) - tri_poly -> AQBE_expected_piece_value(low,piecenum);
	}
	else if (sub_parallelogram_num == 2)
	{
		if (ge(low,1))
			return 0;
		else if (high - low < eps)
			return 0;
		//	tri_poly -> approx_EBRC_expected_piece_value(high,piecenum) - tri_poly -> approx_EBRC_expected_piece_value(low,piecenum);
		else		
			return tri_poly -> EBRC_expected_piece_value(high,piecenum) - tri_poly -> EBRC_expected_piece_value(low,piecenum);
	}
	else if (sub_parallelogram_num == 3)
	{
		if (ge(low,1)) 
			return 0;
		else if (high - low < eps)
			return 0;		
		//	return tri_poly -> approx_DECS_expected_piece_value(high,piecenum) - tri_poly -> approx_DECS_expected_piece_value(low,piecenum);			
		else			
			return tri_poly -> DECS_expected_piece_value(high,piecenum) - tri_poly -> DECS_expected_piece_value(low,piecenum);
	}	
}*/

/*double e_second_moment(double low, double high, int piecenum)
{	

	if(ge(high,1))
	{
		high = 1;
	}

	// assuming low is always greater than 0 and  low<=high

	if (sub_parallelogram_num == 0)
	{
		if (ge(low,1))
			return 0;
		else if (high - low < eps)
			return 0;		
			//return tri_poly -> PAED_e_second_moment_piece_value(high,piecenum) - tri_poly -> PAED_e_second_moment_piece_value(low,piecenum);
		else		
			return tri_poly -> PAED_e_second_moment_piece_value(high,piecenum) - tri_poly -> PAED_e_second_moment_piece_value(low,piecenum);
	}
	else if (sub_parallelogram_num == 1)
	{
		if (ge(low,1))
			return 0;
		else if (high - low < eps)
			return 0;
			//return tri_poly -> AQBE_e_second_moment_piece_value(high,piecenum) - tri_poly -> AQBE_e_second_moment_piece_value(low,piecenum);
		else		
			return tri_poly -> AQBE_e_second_moment_piece_value(high,piecenum) - tri_poly -> AQBE_e_second_moment_piece_value(low,piecenum);
	}
	else if (sub_parallelogram_num == 2)
	{		
		if (ge(low,1))
			return 0;
		else if (high - low < eps)
			return 0;
			//return tri_poly -> EBRC_e_second_moment_piece_value(high,piecenum) - tri_poly -> EBRC_e_second_moment_piece_value(low,piecenum);
		else 		
			return tri_poly -> EBRC_e_second_moment_piece_value(high,piecenum) - tri_poly -> EBRC_e_second_moment_piece_value(low,piecenum);
	}
	else if (sub_parallelogram_num == 3)
	{
		if (ge(low,1)) 
			return 0;
		else if (high - low < eps)
			return 0;
			//return tri_poly -> DECS_e_second_moment_piece_value(high,piecenum) - tri_poly -> DECS_e_second_moment_piece_value(low,piecenum);
		else		
			return tri_poly -> DECS_e_second_moment_piece_value(high,piecenum) - tri_poly -> DECS_e_second_moment_piece_value(low,piecenum);
	}
}*/

double e_path_1(double p,double a,double e,double d, double* exp, double* var, double *crossprob)
{    
    
	*exp = 0;
	*var = 0;
	*crossprob = 0;	

	// Vertex Order PADE 
	if (e_le(p,a) && e_le(a,d) && e_le(d,e))
	{   				
	*exp = e_expected_value(p,a,1) + e_expected_value(a,d,5) - e_expected_value(d,e,3);
	*var = e_second_moment(p,a,1) + e_second_moment(a,d,5) - e_second_moment(d,e,3);
	*crossprob = e_integrate_pdf(p,a,1) + e_integrate_pdf(a,d,5) - e_integrate_pdf(d,e,3);
	}
	
	// Vertex Order PAED 
	else if (e_le(p,a) && e_le(a,e) && e_le(e,d))
 	{
	*exp = e_expected_value(p,a,1) + e_expected_value(a,e,5) - e_expected_value(e,d,4);
	*var = e_second_moment(p,a,1) + e_second_moment(a,e,5) - e_second_moment(e,d,4);
	*crossprob = e_integrate_pdf(p,a,1) + e_integrate_pdf(a,e,5) - e_integrate_pdf(e,d,4);	
	}

	// Vertex Order PDAE 
	else if (e_le(p,d) && e_le(d,a) && e_le(a,e))
	{    	    		
	*exp = e_expected_value(p,d,1) + e_expected_value(d,a,6) - e_expected_value(a,e,3);
	*var = e_second_moment(p,d,1) + e_second_moment(d,a,6) - e_second_moment(a,e,3);
	*crossprob = e_integrate_pdf(p,d,1) + e_integrate_pdf(d,a,6) - e_integrate_pdf(a,e,3);
	}
	
}

double e_path_2(double p,double a,double e,double d, double* exp, double* var, double *crossprob)
{
	*exp = 0;
	*var = 0;
	*crossprob = 0;

	// Vertex Order EDAP (polynomial 1)
	if (e_le(e,d) && e_le(d,a) && e_le(a,p))
    	{   	 	
	*exp = e_expected_value(e,d,3) - e_expected_value(d,a,5) - e_expected_value(a,p,1);
	*var = e_second_moment(e,d,3) - e_second_moment(d,a,5) - e_second_moment(a,p,1);
	*crossprob = e_integrate_pdf(e,d,3) - e_integrate_pdf(d,a,5) - e_integrate_pdf(a,p,1);	
	}
    
    	// Vertex Order EDPA (polynomial 1)
	else if (e_le(e,d) && e_le(d,p) && e_le(p,a))
    	{      
	*exp = e_expected_value(e,d,3) - e_expected_value(d,p,5) - e_expected_value(p,a,2);
	*var = e_second_moment(e,d,3) - e_second_moment(d,p,5) - e_second_moment(p,a,2);
	*crossprob = e_integrate_pdf(e,d,3) - e_integrate_pdf(d,p,5) - e_integrate_pdf(p,a,2);
     	}

    	// Vertex Order EADP (polynomial 1)
	else if (e_le(e,a) && e_le(a,d) && e_le(d,p))
    	{    
	*exp = e_expected_value(e,a,3) - e_expected_value(a,d,6) - e_expected_value(d,p,1);
	*var = e_second_moment(e,a,3) - e_second_moment(a,d,6) - e_second_moment(d,p,1);
	*crossprob = e_integrate_pdf(e,a,3) - e_integrate_pdf(a,d,6) - e_integrate_pdf(d,p,1);		
	}	
}

double e_path_3(double p,double a,double e,double d, double* exp, double* var, double *crossprob)
{
	*exp = 0;
	*var = 0;
	*crossprob = 0;

	/* Only 1 possible ordering DPAE
 All equations derived from PAED_positive_1
 first and the fourth part should be identical*/
	*exp = -e_expected_value(0,d,6) - e_expected_value(d,p,1) + e_expected_value(a,e,2) - e_expected_value(e,INFINITY,6);
	*var = -e_second_moment(0,d,6) - e_second_moment(d,p,1) + e_second_moment(a,e,2) - e_second_moment(e,INFINITY,6);
	*crossprob = -e_integrate_pdf(0,d,6) - e_integrate_pdf(d,p,1) + e_integrate_pdf(a,e,2) - e_integrate_pdf(e,INFINITY,6);	
}

double e_path_4(double p,double a,double e,double d, double* exp, double* var, double *crossprob)
{
	*exp = 0;
	*var = 0;
	*crossprob = 0;
	
	/* Only 1 possible ordering AEDP
% All equations derived from PAED_positive_1
% first and the fourth part should be identical*/	
	*exp = e_expected_value(0,a,6) - e_expected_value(a,e,3) + e_expected_value(d,p,4) + e_expected_value(p,INFINITY,6);	
	*var = e_second_moment(0,a,6) - e_second_moment(a,e,3) + e_second_moment(d,p,4) + e_second_moment(p,INFINITY,6);	
	*crossprob = e_integrate_pdf(0,a,6) - e_integrate_pdf(a,e,3) + e_integrate_pdf(d,p,4) + e_integrate_pdf(p,INFINITY,6);	
}

double e_path_5(double p,double a,double e,double d, double* exp, double* var, double *crossprob)
{
	*exp = 0;
	*var = 0;
	*crossprob = 0;
	
	// Vertex Order DAPE 
	if (e_le(d,a) && e_le(a,p) && e_le(p,e))
	{  
	*exp = e_expected_value(0,d,7) + e_expected_value(0,d,8) + e_expected_value(d,a,7) + e_expected_value(d,a,10) + e_expected_value(a,p,9)+ e_expected_value(a,p,10) + e_expected_value(p,e,9)+ e_expected_value(p,e,12) + e_expected_value(e,INFINITY,11)+ e_expected_value(e,INFINITY,12);	
	*var = e_second_moment(0,d,7) + e_second_moment(0,d,8) + e_second_moment(d,a,7) + e_second_moment(d,a,10) + e_second_moment(a,p,9)+ e_second_moment(a,p,10) + e_second_moment(p,e,9)+ e_second_moment(p,e,12) + e_second_moment(e,INFINITY,11)+ e_second_moment(e,INFINITY,12);	
	*crossprob = e_integrate_pdf(0,d,7) + e_integrate_pdf(0,d,8) + e_integrate_pdf(d,a,7) + e_integrate_pdf(d,a,10) + e_integrate_pdf(a,p,9)+ e_integrate_pdf(a,p,10) + e_integrate_pdf(p,e,9)+ e_integrate_pdf(p,e,12) + e_integrate_pdf(e,INFINITY,11)+ e_integrate_pdf(e,INFINITY,12);
	}
    
    	// Vertex Order ADPE (polynomial 1)
	else if (e_le(a,d) && e_le(d,p) && e_le(p,e))
	{
	*exp = e_expected_value(0,a,7) + e_expected_value(0,a,8) + e_expected_value(a,d,9) + e_expected_value(a,d,8) + e_expected_value(d,p,9)+ e_expected_value(d,p,10) + e_expected_value(p,e,9)+ e_expected_value(p,e,12) + e_expected_value(e,INFINITY,11)+ e_expected_value(e,INFINITY,12);
	*var = e_second_moment(0,a,7) + e_second_moment(0,a,8) + e_second_moment(a,d,9) + e_second_moment(a,d,8) + e_second_moment(d,p,9)+ e_second_moment(d,p,10) + e_second_moment(p,e,9)+ e_second_moment(p,e,12) + e_second_moment(e,INFINITY,11)+ e_second_moment(e,INFINITY,12);
	*crossprob = e_integrate_pdf(0,a,7) + e_integrate_pdf(0,a,8) + e_integrate_pdf(a,d,9) + e_integrate_pdf(a,d,8) + e_integrate_pdf(d,p,9)+ e_integrate_pdf(d,p,10) + e_integrate_pdf(p,e,9)+ e_integrate_pdf(p,e,12) + e_integrate_pdf(e,INFINITY,11)+ e_integrate_pdf(e,INFINITY,12);
	}
    
	// Vertex Order DAEP (polynomial 1)
	else if (e_le(d,a) && e_le(a,e) && e_le(e,p))
	{  
	*exp = e_expected_value(0,d,7) + e_expected_value(0,d,8) + e_expected_value(d,a,7) + e_expected_value(d,a,10) + e_expected_value(a,e,9)+ e_expected_value(a,e,10) + e_expected_value(e,p,11)+ e_expected_value(e,p,10) + e_expected_value(p,INFINITY,11)+ e_expected_value(p,INFINITY,12);
	*var = e_second_moment(0,d,7) + e_second_moment(0,d,8) + e_second_moment(d,a,7) + e_second_moment(d,a,10) + e_second_moment(a,e,9)+ e_second_moment(a,e,10) + e_second_moment(e,p,11)+ e_second_moment(e,p,10) + e_second_moment(p,INFINITY,11)+ e_second_moment(p,INFINITY,12);
	*crossprob = e_integrate_pdf(0,d,7) + e_integrate_pdf(0,d,8) + e_integrate_pdf(d,a,7) + e_integrate_pdf(d,a,10) + e_integrate_pdf(a,e,9)+ e_integrate_pdf(a,e,10) + e_integrate_pdf(e,p,11)+ e_integrate_pdf(e,p,10) + e_integrate_pdf(p,INFINITY,11)+ e_integrate_pdf(p,INFINITY,12);
	}
   
   	// Vertex Order ADEP (polynomial 1)
	else if (e_le(a,d) && e_le(d,e) && e_le(e,p))
	{       
	*exp = e_expected_value(0,a,7) + e_expected_value(0,a,8) + e_expected_value(a,d,9) + e_expected_value(a,d,8) + e_expected_value(d,e,9)+ e_expected_value(d,e,10) + e_expected_value(e,p,11)+ e_expected_value(e,p,10) + e_expected_value(p,INFINITY,11)+ e_expected_value(p,INFINITY,12);
	*var = e_second_moment(0,a,7) + e_second_moment(0,a,8) + e_second_moment(a,d,9) + e_second_moment(a,d,8) + e_second_moment(d,e,9)+ e_second_moment(d,e,10) + e_second_moment(e,p,11)+ e_second_moment(e,p,10) + e_second_moment(p,INFINITY,11)+ e_second_moment(p,INFINITY,12);
	*crossprob = e_integrate_pdf(0,a,7) + e_integrate_pdf(0,a,8) + e_integrate_pdf(a,d,9) + e_integrate_pdf(a,d,8) + e_integrate_pdf(d,e,9)+ e_integrate_pdf(d,e,10) + e_integrate_pdf(e,p,11)+ e_integrate_pdf(e,p,10) + e_integrate_pdf(p,INFINITY,11)+ e_integrate_pdf(p,INFINITY,12);
	}	
}

double e_path_6(double p,double a,double e,double d, double* exp, double* var, double *crossprob)
{
	*exp = 0;
	*var = 0;
	*crossprob = 0;

//% only one possible vertex ordering DPA (Ignore -ve vertices)
	*exp = -e_expected_value(0,d,6) - e_expected_value(d,p,1) + e_expected_value(a,INFINITY,2);	
	*var = -e_second_moment(0,d,6) - e_second_moment(d,p,1) + e_second_moment(a,INFINITY,2);	
	*crossprob = -e_integrate_pdf(0,d,6) - e_integrate_pdf(d,p,1) + e_integrate_pdf(a,INFINITY,2);		
}

double e_path_7(double p,double a,double e,double d, double* exp, double* var, double *crossprob)
{
	*exp = 0;
	*var = 0;
	*crossprob = 0;

	// vertex order DAP
	if (e_le(d,a) && e_le(a,p))
	{  	
	*exp = e_expected_value(0,d,7) + e_expected_value(0,d,8) + e_expected_value(d,a,7) + e_expected_value(d,a,10) + e_expected_value(a,p,9)+ e_expected_value(a,p,10) + e_expected_value(p,INFINITY,9)+ e_expected_value(p,INFINITY,12);
	*var = e_second_moment(0,d,7) + e_second_moment(0,d,8) + e_second_moment(d,a,7) + e_second_moment(d,a,10) + e_second_moment(a,p,9)+ e_second_moment(a,p,10) + e_second_moment(p,INFINITY,9)+ e_second_moment(p,INFINITY,12);
	*crossprob = e_integrate_pdf(0,d,7) + e_integrate_pdf(0,d,8) + e_integrate_pdf(d,a,7) + e_integrate_pdf(d,a,10) + e_integrate_pdf(a,p,9)+ e_integrate_pdf(a,p,10) + e_integrate_pdf(p,INFINITY,9)+ e_integrate_pdf(p,INFINITY,12);
	}
    
    	// vertex order ADP
	else if (e_le(a,d) && e_le(d,p))
	{       
	*exp = e_expected_value(0,a,7) + e_expected_value(0,a,8) + e_expected_value(a,d,9) + e_expected_value(a,d,8) + e_expected_value(d,p,9)+ e_expected_value(d,p,10) + e_expected_value(p,INFINITY,9)+ e_expected_value(p,INFINITY,12);
	*var = e_second_moment(0,a,7) + e_second_moment(0,a,8) + e_second_moment(a,d,9) + e_second_moment(a,d,8) + e_second_moment(d,p,9)+ e_second_moment(d,p,10) + e_second_moment(p,INFINITY,9)+ e_second_moment(p,INFINITY,12);
	*crossprob = e_integrate_pdf(0,a,7) + e_integrate_pdf(0,a,8) + e_integrate_pdf(a,d,9) + e_integrate_pdf(a,d,8) + e_integrate_pdf(d,p,9)+ e_integrate_pdf(d,p,10) + e_integrate_pdf(p,INFINITY,9)+ e_integrate_pdf(p,INFINITY,12);        
	}	
}

double e_path_8(double p,double a,double e,double d, double* exp, double* var, double *crossprob)
{
	*exp = 0;
	*var = 0;
	*crossprob = 0;
	
	// only one possible vertex ordering AED (Ignore -ve vertices)
	*exp = e_expected_value(0,a,6) - e_expected_value(a,e,3) + e_expected_value(d,INFINITY,4);	
	*var = e_second_moment(0,a,6) - e_second_moment(a,e,3) + e_second_moment(d,INFINITY,4);	
	*crossprob = e_integrate_pdf(0,a,6) - e_integrate_pdf(a,e,3) + e_integrate_pdf(d,INFINITY,4);		
}

double e_path_9(double p,double a,double e,double d, double* exp, double* var, double *crossprob)
{
	*exp = 0;
	*var = 0;
	*crossprob = 0;
	
	// vertex order ADE
	if (e_le(a,d) && e_le(d,e))
	{       
	*exp = e_expected_value(0,a,7) + e_expected_value(0,a,8) + e_expected_value(a,d,9) + e_expected_value(a,d,8) + e_expected_value(d,e,9)+ e_expected_value(d,e,10) + e_expected_value(e,INFINITY,11)+ e_expected_value(e,INFINITY,10);
	*var = e_second_moment(0,a,7) + e_second_moment(0,a,8) + e_second_moment(a,d,9) + e_second_moment(a,d,8) + e_second_moment(d,e,9)+ e_second_moment(d,e,10) + e_second_moment(e,INFINITY,11)+ e_second_moment(e,INFINITY,10);
	*crossprob = e_integrate_pdf(0,a,7) + e_integrate_pdf(0,a,8) + e_integrate_pdf(a,d,9) + e_integrate_pdf(a,d,8) + e_integrate_pdf(d,e,9)+ e_integrate_pdf(d,e,10) + e_integrate_pdf(e,INFINITY,11)+ e_integrate_pdf(e,INFINITY,10);
	}
    
	// vertex order DAE
	else if (e_le(d,a) && e_le(a,e))
	{    
	*exp = e_expected_value(0,d,7) + e_expected_value(0,d,8) + e_expected_value(d,a,7) + e_expected_value(d,a,10) + e_expected_value(a,e,9)+ e_expected_value(a,e,10) + e_expected_value(e,INFINITY,11)+ e_expected_value(e,INFINITY,10);
	*var = e_second_moment(0,d,7) + e_second_moment(0,d,8) + e_second_moment(d,a,7) + e_second_moment(d,a,10) + e_second_moment(a,e,9)+ e_second_moment(a,e,10) + e_second_moment(e,INFINITY,11)+ e_second_moment(e,INFINITY,10);
	*crossprob = e_integrate_pdf(0,d,7) + e_integrate_pdf(0,d,8) + e_integrate_pdf(d,a,7) + e_integrate_pdf(d,a,10) + e_integrate_pdf(a,e,9)+ e_integrate_pdf(a,e,10) + e_integrate_pdf(e,INFINITY,11)+ e_integrate_pdf(e,INFINITY,10);
	}	
}

double e_path_10(double p,double a,double e,double d, double* exp, double* var, double *crossprob)
{
	*exp = 0;
	*var = 0;
	*crossprob = 0;

	// vertex order PAE
	if (e_le(p,a) && e_le(a,e))
	{ 
	*exp = e_expected_value(p,a,1) + e_expected_value(a,e,5) - e_expected_value(e,INFINITY,4);	
	*var = e_second_moment(p,a,1) + e_second_moment(a,e,5) - e_second_moment(e,INFINITY,4);	
	*crossprob = e_integrate_pdf(p,a,1) + e_integrate_pdf(a,e,5) - e_integrate_pdf(e,INFINITY,4);	
	}
    
	// vertex order APE    
	else if (e_le(a,p) && e_le(p,e)) 
	{    
	*exp = e_expected_value(a,p,2) + e_expected_value(p,e,5) - e_expected_value(e,INFINITY,4);	
	*var = e_second_moment(a,p,2) + e_second_moment(p,e,5) - e_second_moment(e,INFINITY,4);	
	*crossprob = e_integrate_pdf(a,p,2) + e_integrate_pdf(p,e,5) - e_integrate_pdf(e,INFINITY,4);	
	}    

	// vertex order AEP      
	else if (e_le(a,e) && e_le(e,p))  
	{    
	*exp = e_expected_value(a,e,2) - e_expected_value(e,p,6) - e_expected_value(p,INFINITY,4);	
	*var = e_second_moment(a,e,2) - e_second_moment(e,p,6) - e_second_moment(p,INFINITY,4);	
	*crossprob = e_integrate_pdf(a,e,2) - e_integrate_pdf(e,p,6) - e_integrate_pdf(p,INFINITY,4);		
	}	
}

double e_path_11(double p,double a,double e,double d, double* exp, double* var, double *crossprob)
{
	*exp = 0;
	*var = 0;
	*crossprob = 0;
	
	// vertex order EDP
	if (e_le(e,d) && e_le(d,p))
	{  
	*exp = e_expected_value(e,d,3) - e_expected_value(d,p,5) - e_expected_value(p,INFINITY,2);	
	*var = e_second_moment(e,d,3) - e_second_moment(d,p,5) - e_second_moment(p,INFINITY,2);		
	*crossprob = e_integrate_pdf(e,d,3) - e_integrate_pdf(d,p,5) - e_integrate_pdf(p,INFINITY,2);	 		
	}    

	// vertex order DEP    
	else if (e_le(d,e) && e_le(e,p))  
	{     
	*exp = e_expected_value(d,e,4) - e_expected_value(e,p,5) - e_expected_value(p,INFINITY,2);	
	*var = e_second_moment(d,e,4) - e_second_moment(e,p,5) - e_second_moment(p,INFINITY,2);	
	*crossprob = e_integrate_pdf(d,e,4) - e_integrate_pdf(e,p,5) - e_integrate_pdf(p,INFINITY,2);	
	}
    
	// vertex order DPE    
	else if (e_le(d,p) && e_le(p,e))  
	{	         
	*exp = e_expected_value(d,p,4) + e_expected_value(p,e,6) - e_expected_value(e,INFINITY,2);	
	*var = e_second_moment(d,p,4) + e_second_moment(p,e,6) - e_second_moment(e,INFINITY,2);	
	*crossprob = e_integrate_pdf(d,p,4) + e_integrate_pdf(p,e,6) - e_integrate_pdf(e,INFINITY,2);	
	}	
}

double e_path_12(double p,double a,double e,double d, double* exp, double* var, double *crossprob)
{
	*exp = 0;
	*var = 0;
	*crossprob = 0;

	// Vertex order PA
	if (e_le(p,a))
	{  
	*exp = e_expected_value(p,a,1) + e_expected_value(a,INFINITY,5);
	*var = e_second_moment(p,a,1) + e_second_moment(a,INFINITY,5);
	*crossprob = e_integrate_pdf(p,a,1) + e_integrate_pdf(a,INFINITY,5);
	}
    
	// Vertex order AP
	else if (e_le(a,p))
	{      
	*exp = e_expected_value(a,p,2) + e_expected_value(p,INFINITY,5);
	*var = e_second_moment(a,p,2) + e_second_moment(p,INFINITY,5);
	*crossprob = e_integrate_pdf(a,p,2) + e_integrate_pdf(p,INFINITY,5);
	}	
}

double e_path_13(double p,double a,double e,double d, double* exp, double* var, double *crossprob)
{
	*exp = 0;
	*var = 0;
	*crossprob = 0;
	
	// Vertex order ED
	if (e_le(e,d))
	{  
	*exp = e_expected_value(e,d,3) - e_expected_value(d,INFINITY,5);
	*var = e_second_moment(e,d,3) - e_second_moment(d,INFINITY,5);
	*crossprob = e_integrate_pdf(e,d,3) - e_integrate_pdf(d,INFINITY,5);
	}
    
	// Vertex order DE
	else if (e_le(d,e))
	{     
	*exp = e_expected_value(d,e,4) - e_expected_value(e,INFINITY,5);
	*var = e_second_moment(d,e,4) - e_second_moment(e,INFINITY,5);
	*crossprob = e_integrate_pdf(d,e,4) - e_integrate_pdf(e,INFINITY,5);    
	}	
}

double e_path_14(double p,double a,double e,double d, double* exp, double* var, double *crossprob)
{
	*exp = 0;
	*var = 0;
	*crossprob = 0;
	
	// vertex order AE (only one possible)  
	*exp = e_expected_value(0,a,6) - e_expected_value(a,e,3);
	*var = e_second_moment(0,a,6) - e_second_moment(a,e,3);
	*crossprob = e_integrate_pdf(0,a,6) - e_integrate_pdf(a,e,3);	
}

double e_path_15(double p,double a,double e,double d, double* exp, double* var, double *crossprob)
{

	*exp = 0;
	*var = 0;
	*crossprob = 0;
	
 	// vertex order DP (only one possible)
	*exp = -e_expected_value(0,d,6) - e_expected_value(d,p,1);	
	*var = -e_second_moment(0,d,6) - e_second_moment(d,p,1);	
	*crossprob = -e_integrate_pdf(0,d,6) - e_integrate_pdf(d,p,1);	
	
}

double e_path_16(double p,double a,double e,double d, double* exp, double* var, double *crossprob)
{
	*exp = 0;
	*var = 0;
	*crossprob = 0;
	
	// Vertex order DA
	if (e_le(d,a))
	{   
	*exp = e_expected_value(0,d,7) + e_expected_value(0,d,8) + e_expected_value(d,a,7) + e_expected_value(d,a,10) + e_expected_value(a,INFINITY,9) + e_expected_value(a,INFINITY,10);
	*var = e_second_moment(0,d,7) + e_second_moment(0,d,8) + e_second_moment(d,a,7) + e_second_moment(d,a,10) + e_second_moment(a,INFINITY,9) + e_second_moment(a,INFINITY,10);
	*crossprob = e_integrate_pdf(0,d,7) + e_integrate_pdf(0,d,8) + e_integrate_pdf(d,a,7) + e_integrate_pdf(d,a,10) + e_integrate_pdf(a,INFINITY,9) + e_integrate_pdf(a,INFINITY,10);  
	}
  
	// Vertex order AD
	else if (e_le(a,d))
	{     
	*exp = e_expected_value(0,a,7) + e_expected_value(0,a,8) + e_expected_value(a,d,9) + e_expected_value(a,d,8) + e_expected_value(d,INFINITY,9) + e_expected_value(d,INFINITY,10);
	*var = e_second_moment(0,a,7) + e_second_moment(0,a,8) + e_second_moment(a,d,9) + e_second_moment(a,d,8) + e_second_moment(d,INFINITY,9) + e_second_moment(d,INFINITY,10);
	*crossprob = e_integrate_pdf(0,a,7) + e_integrate_pdf(0,a,8) + e_integrate_pdf(a,d,9) + e_integrate_pdf(a,d,8) + e_integrate_pdf(d,INFINITY,9) + e_integrate_pdf(d,INFINITY,10);    
	}	
}

double e_path_17(double p,double a,double e,double d, double* exp, double* var, double *crossprob)
{

	*exp = 0;
	*var = 0;
	*crossprob = 0;
	
	*exp = e_expected_value(a,INFINITY,2);	
	*var = e_second_moment(a,INFINITY,2);	
	*crossprob = e_integrate_pdf(a,INFINITY,2);	
	
}

double e_path_18(double p,double a,double e,double d, double* exp, double* var, double *crossprob)
{

	*exp = 0;
	*var = 0;
	*crossprob = 0;
	
	*exp = e_expected_value(d,INFINITY,4);	
	*var = e_second_moment(d,INFINITY,4);	
	*crossprob = e_integrate_pdf(d,INFINITY,4);	
	
}

double e_path_19(double p,double a,double e,double d, double* exp, double* var, double *crossprob)
{

	*exp = 0;
	*var = 0;
	*crossprob = 0;

	*exp = e_expected_value(a,INFINITY,2);	
	*var = e_second_moment(a,INFINITY,2);	
	*crossprob = e_integrate_pdf(a,INFINITY,2);		

}

double e_path_20(double p,double a,double e,double d, double* exp, double* var, double *crossprob)
{

	*exp = 0;
	*var = 0;
	*crossprob = 0;
	
	*exp = e_expected_value(a,e,2) - e_expected_value(e,INFINITY,6);	
	*var = e_second_moment(a,e,2) - e_second_moment(e,INFINITY,6);	
	*crossprob = e_integrate_pdf(a,e,2) - e_integrate_pdf(e,INFINITY,6);		
}

double e_path_21(double p,double a,double e,double d, double* exp, double* var, double *crossprob)
{

	*exp = 0;
	*var = 0;
	*crossprob = 0;

	// vertex order DP
	*exp = e_expected_value(0,d,3) - e_expected_value(d,p,1);
	*var = e_second_moment(0,d,3) - e_second_moment(d,p,1);
	*crossprob = e_integrate_pdf(0,d,3) - e_integrate_pdf(d,p,1);	
}

double e_path_22(double p,double a,double e,double d, double* exp, double* var, double *crossprob)
{

	*exp = 0;
	*var = 0;
	*crossprob = 0;

	*exp = e_expected_value(d,INFINITY,4);	
	*var = e_second_moment(d,INFINITY,4);	
	*crossprob = e_integrate_pdf(d,INFINITY,4);		
}

double e_path_23(double p,double a,double e,double d, double* exp, double* var, double *crossprob)
{
	*exp = 0;
	*var = 0;
	*crossprob = 0;
	
	// vertex order DP

	*exp = e_expected_value(d,p,4) + e_expected_value(p,INFINITY,6);
	*var = e_second_moment(d,p,4) + e_second_moment(p,INFINITY,6);
	*crossprob = e_integrate_pdf(d,p,4) + e_integrate_pdf(p,INFINITY,6);	
}

double e_path_24(double p,double a,double e,double d, double* exp, double* var, double *crossprob)
{

	*exp = 0;
	*var = 0;
	*crossprob = 0;
	
	// vertex order AE
	*exp = e_expected_value(0,a,1) - e_expected_value(a,e,3);
	*var = e_second_moment(0,a,1) - e_second_moment(a,e,3);
	*crossprob = e_integrate_pdf(0,a,1) - e_integrate_pdf(a,e,3);	
}


void e_pdf_piece(double p,double a,double e, double d, double* exp, double* var, double* crossprob)
{
	double expected_pdf_piece = 0;

// Provide the coordinates of the the  the vertices, and coordinates of
// leftmost and the rightmost intersection of the parallelogram with 
//the horizontal axis.

double Px, Py, Ax, Ay, Ex, Ey, Dx, Dy, Z2right, Z2left;

// PQRS

// Following is the mapping

// P -> P
// Q -> A
// R -> E
// S -> D
	Px = e_mu2-e_mu1+e_delta2-e_delta1;
	Py = e_c-e_mu1-e_delta1;
	Ax = e_mu2-e_mu1+e_delta1+e_delta2;
	Ay = e_c-e_mu1+e_delta1;
	Ex = e_mu2-e_mu1+e_delta1-e_delta2;
	Ey = e_c-e_mu1+e_delta1;
	Dx = e_mu2-e_mu1-e_delta1-e_delta2;
	Dy = e_c-e_mu1-e_delta1;
	// intersection of AP with with horizontal axis (Z2). 
	Z2right = e_mu2-e_c+e_delta2;
	// intersection of ED with with horizontal axis (Z2). 
	Z2left = e_mu2-e_c-e_delta2;




//        e -  -  - a
//        /         /
//       /   P1    /
//      /         /
//     d- - - - - p

//% P1 is a polynomial.


//-------------------------------------------------------------------------
 // sets of 4 positive vertices 
 
// Slopes of P, A, E, D all are greater than 0. There are 5 possible parallelogram configurations for
// when all slopes are are greater than 0.

//  If whole parallelogram lies in the first quadrant   
if (e_ge(p,0) && e_ge(Px,0) && e_ge(Py,0) && e_ge(a,0) && e_ge(Ax,0) && e_ge(Ay,0) && e_ge(e,0) && e_ge(Ex,0) && e_ge(Ey,0) && e_ge(d,0) && e_ge(Dx,0) && e_ge(Dy,0))

    e_path_1(p,a,e,d,exp,var,crossprob);

//  If whole parallelogram lies in the third quadrant
else if (e_ge(p,0) && e_le(Px,0) && e_le(Py,0) && e_ge(a,0) && e_le(Ax,0) && e_le(Ay,0) && e_ge(e,0) && e_le(Ex,0) && e_le(Ey,0) && e_ge(d,0) && e_le(Dx,0) && e_le(Dy,0))
    
    e_path_2(p,a,e,d,exp,var,crossprob);
    
//  If E,A lie in the first quadrant and D,P lie in the third quadrant assuming edge AP crosses second quadrant;
else if (e_ge(p,0) && e_le(Px,0)  && e_le(Py,0) && e_ge(a,0) && e_ge(Ax,0) && e_ge(Ay,0) && e_ge(e,0) && e_ge(Ex,0) && e_ge(Ey,0) && e_ge(d,0) && e_le(Dx,0) && e_le(Dy,0) && e_le(Z2right,0))
    
    e_path_3(p,a,e,d,exp,var,crossprob);
    
//  If E,A lie in the first quadrant and D,P lie in the third quadrant assuming edge ED crosses fourh quadrant;
else if (e_ge(p,0) && e_le(Px,0) && e_le(Py,0) && e_ge(a,0) && e_ge(Ax,0) && e_ge(Ay,0) && e_ge(e,0) && e_ge(Ex,0) && e_ge(Ey,0) && e_ge(d,0) && e_le(Dx,0) && e_le(Dy,0) && e_ge(Z2left,0))

    e_path_4(p,a,e,d,exp,var,crossprob);	
   
//  If E,A lie in the first quadrant and D,P lie in the third quadrant assuming edge ED crosses second quadrant
//  and AP crosses the fourth quadrant;
else if (e_ge(p,0) && e_le(Px,0) && e_le(Py,0) && e_ge(a,0) && e_ge(Ax,0) && e_ge(Ay,0) && e_ge(e,0) && e_ge(Ex,0) && e_ge(Ey,0) && e_ge(d,0) && e_le(Dx,0) && e_le(Dy,0) && e_le(Z2left,0) && e_ge(Z2right,0))

    e_path_5(p,a,e,d,exp,var,crossprob);
    
//------------------------------------------------------------------------
    // sets of 3 positive vertices 
    
    // Slopes of P, A, D are greater than 0 and slope of E is less than 0. There are 2 possible parallelogram configurations
    // in this case.
    
    // If D, P lie in the third quadrant,A in first and E in second assuming edge ED
    // and AP cross the second quadrant;
else if (e_ge(p,0) && e_le(Px,0) && e_le(Py,0) && e_ge(a,0) && e_ge(Ax,0) && e_ge(Ay,0) && e_le(e,0) && e_le(Ex,0) && e_ge(Ey,0) && e_ge(d,0) && e_le(Dx,0)  && e_le(Dy,0) && e_le(Z2right,0))
    
    e_path_6(p,a,e,d,exp,var,crossprob);
    
    // If D, P lie in the third quadrant,A in first and E in second assuming edge ED crosses second quadrant
    // and AP crosses the fourth quadrant;
else if (e_ge(p,0) && e_le(Px,0) && e_le(Py,0) && e_ge(a,0) && e_ge(Ax,0) && e_ge(Ay,0) && e_le(e,0) && e_le(Ex,0) && e_ge(Ey,0) && e_ge(d,0) && e_le(Dx,0) && e_le(Dy,0) && e_le(Z2left,0) && e_ge(Z2right,0))
    
    e_path_7(p,a,e,d,exp,var,crossprob);
    
    // Slopes of A, E, D are greater than 0 and slope of P is less than 0. There are 2 possible parallelogram configurations
    //in this case.
    
    // If E, A lie in the first quadrant,D in third and P in fourth assuming edge ED
    //and AP cross the fourth quadrant;
else if (e_le(p,0) && e_ge(Px,0) && e_le(Py,0) && e_ge(a,0) && e_ge(Ax,0) && e_ge(Ay,0) && e_ge(e,0) && e_ge(Ex,0) && e_ge(Ey,0) && e_ge(d,0) && e_le(Dx,0) && e_le(Dy,0) && e_ge(Z2left,0))
    
    e_path_8(p,a,e,d,exp,var,crossprob);
    
    // If E, A lie in the first quadrant,D in third and P in fourth assuming edge ED crosses second quadrant
    //and AP crosses the fourth quadrant;
else if (e_le(p,0) && e_ge(Px,0) && e_le(Py,0) && e_ge(a,0) && e_ge(Ax,0) && e_ge(Ay,0) && e_ge(e,0) && e_ge(Ex,0) && e_ge(Ey,0) && e_ge(d,0) && e_le(Dx,0)  && e_le(Dy,0) && e_le(Z2left,0) && e_ge(Z2right,0))
    
    e_path_9(p,a,e,d,exp,var,crossprob);
    
    // Newly added after debugging
    // Slopes of P, A, E are greater than 0 and slope of D is less than 0.
    // P, A, E in the first quadrant and D in the second quadrant 
else if (e_ge(p,0) && e_ge(Px,0) && e_ge(Py,0) && e_ge(a,0) && e_ge(Ax,0) && e_ge(Ay,0) && e_ge(e,0) && e_ge(Ex,0) && e_ge(Ey,0) && e_le(d,0) && e_le(Dx,0) && e_ge(Dy,0))
    
    e_path_10(p,a,e,d,exp,var,crossprob);
    
    // Slopes of E, D, P are greater than 0 and slope of A is less than 0.
    // E, D, P in the third quadrant and A in the fourth quadrant     
else if (e_ge(p,0) && e_le(Px,0) && e_le(Py,0) && e_le(a,0) && e_ge(Ax,0) && e_le(Ay,0) && e_ge(e,0) && e_le(Ex,0) && e_le(Ey,0) && e_ge(d,0) && e_le(Dx,0) && e_le(Dy,0))
    
    e_path_11(p,a,e,d,exp,var,crossprob);
    
   //------------------------------------------------------------------------
    // sets of 2 positive vertices 
    // There are 5 possible parallelogram configurations in this case.
    
// E,D in second quadrant and A,P in the first quadrant    
else if (e_ge(p,0) && e_ge(Px,0) && e_ge(Py,0) && e_ge(a,0) && e_ge(Ax,0) && e_ge(Ay,0) && e_le(e,0) && e_le(Ex,0)  && e_ge(Ey,0) && e_le(d,0) && e_le(Dx,0)  && e_ge(Dy,0))
   
    e_path_12(p,a,e,d,exp,var,crossprob);
    
// E,D in third quadrant and A,P in the fourth quadrant        
else if (e_le(p,0) && e_ge(Px,0) && e_le(Py,0) && e_le(a,0) && e_ge(Ax,0) && e_le(Ay,0) && e_ge(e,0) && e_le(Ex,0) && e_le(Ey,0) && e_ge(d,0) && e_le(Dx,0) && e_le(Dy,0))
   
    e_path_13(p,a,e,d,exp,var,crossprob);
    
// E,A in first quadrant and D,P in the fourth quadrant        
else if (e_le(p,0) && e_ge(Px,0) && e_le(Py,0) && e_ge(a,0) && e_ge(Ax,0) && e_ge(Ay,0) && e_ge(e,0) && e_ge(Ex,0)  && e_ge(Ey,0) && e_le(d,0) && e_ge(Dx,0) && e_le(Dy,0))
   
   e_path_14(p,a,e,d,exp,var,crossprob);
    
// E,A in second quadrant and D,P in the third quadrant           
else if (e_ge(p,0) && e_le(Px,0) && e_le(Py,0) && e_le(a,0) && e_le(Ax,0) && e_ge(Ay,0) && e_le(e,0) && e_le(Ex,0) && e_ge(Ey,0) && e_ge(d,0) && e_le(Dx,0) && e_le(Dy,0))
   
    e_path_15(p,a,e,d,exp,var,crossprob);
    
// E in second quadrant, P in fourth quadrant, A in the first quadrant and D in the third quadrant           
else if (e_le(p,0) && e_ge(Px,0) && e_le(Py,0) && e_ge(a,0) && e_ge(Ax,0) && e_ge(Ay,0) && e_le(e,0) && e_le(Ex,0)  && e_ge(Ey,0) && e_ge(d,0) && e_le(Dx,0) && e_le(Dy,0))
   
    e_path_16(p,a,e,d,exp,var,crossprob);
   
  //------------------------------------------------------------------------
    //  1 positive vertex possibilities 
    // There are 2 possible parallelogram configurations in this case.

// A in first quadrant, E,D,P in the second quadrant           
else if (e_le(p,0) && e_le(Px,0) && e_ge(Py,0) && e_ge(a,0) && e_ge(Ax,0) && e_ge(Ay,0) && e_le(e,0) && e_le(Ex,0) && e_ge(Ey,0) && e_le(d,0) && e_le(Dx,0)  && e_ge(Dy,0))
    
    e_path_17(p,a,e,d,exp,var,crossprob);
    
// D in third quadrant, E,A,P in the fourth quadrant           
else if (e_le(p,0) && e_ge(Px,0) && e_le(Py,0) && e_le(a,0) && e_ge(Ax,0) && e_le(Ay,0) && e_le(e,0) && e_ge(Ex,0) && e_le(Ey,0) && e_ge(d,0) && e_le(Dx,0)  && e_le(Dy,0))

    e_path_18(p,a,e,d,exp,var,crossprob);
    
    //------------------------------------------------------------------------
    //  special case where one of the vertex is at origin(0,0) i.e. when
    //  slope of one vertex is NAN
    //  There are 6 possible parallelogram configurations in this case.

    // P at origin, A in first quadrant, E in the second quadrant 
else if (e_eq(Px,0) && e_eq(Py,0) && e_ge(a,0) && e_ge(Ax,0) && e_ge(Ay,0) && e_le(e,0) && e_le(Ex,0) && e_ge(Ey,0))
    
     e_path_19(p,a,e,d,exp,var,crossprob);
    
   // P at origin, E,A in first quadrant
else if (e_eq(Px,0) && e_eq(Py,0) && e_ge(a,0) && e_ge(Ax,0) && e_ge(Ay,0) && e_ge(e,0) && e_ge(Ex,0) && e_ge(Ey,0)) 
    
     e_path_20(p,a,e,d,exp,var,crossprob);
    
   // A at origin, D,P in the third quadrant 
else if (e_ge(p,0) && e_le(Px,0) && e_le(Py,0) && e_eq(Ax,0) && e_eq(Ay,0) && e_ge(d,0) && e_le(Dx,0) && e_le(Dy,0))
    
     e_path_21(p,a,e,d,exp,var,crossprob);
    
   // E at origin, D in the third quadrant and P in the fourth quadrant 
else if (e_le(p,0) && e_ge(Px,0) && e_le(Py,0) && e_eq(Ex,0) && e_eq(Ey,0) && e_ge(d,0) && e_le(Dx,0) && e_le(Dy,0))
    
     e_path_22(p,a,e,d,exp,var,crossprob);
    
   // E at origin, D,P in the third quadrant 
else if (e_ge(p,0) && e_le(Px,0) && e_le(Py,0) && e_eq(Ex,0) && e_eq(Ey,0) && e_ge(d,0) && e_le(Dx,0) && e_le(Dy,0)) 
    
     e_path_23(p,a,e,d,exp,var,crossprob);
    
   // D at origin, E,A in the first quadrant  
else if (e_ge(a,0) && e_ge(Ax,0) && e_ge(Ay,0) && e_ge(e,0) && e_ge(Ex,0) && e_ge(Ey,0) && e_eq(Dx,0) && e_eq(Dy,0))  

     e_path_24(p,a,e,d,exp,var,crossprob);

// Added during debugging
//If one of the e_paths is not taken means subparallelogram lies in the second/fourth quadrant where all slope values
// are negative		
else
{
	// no contribution to expected value, variance and crossing probability from the subparallelogram
	*exp = 0;
	*var = 0;
	*crossprob = 0;
}     	

}

double e_determine_slope_value(double z1, double z2)
{

	double slope_temp;

	if(fabs(z2) < DBL_EPSILON)
	{
		// Assign NAN value if both z1 and z2 are 0
		if(fabs(z1) < DBL_EPSILON)
		{
			slope_temp = std::numeric_limits<double>::quiet_NaN();
		}
		else if (z1 > 0)
		{
			slope_temp = INFINITY;
		}
		else
		{
			slope_temp = -INFINITY;
		}
	}
	else
	{
		slope_temp = (double)((z1) / (z2));
	}
	
	return slope_temp;

}

/*void pdf_piece(double p,double a,double e, double d, double* exp, double* var, double* crossprob)
{
	double expected_pdf_piece = 0;

// Provide the coordinates of the the  the vertices, and coordinates of
// leftmost and the rightmost intersection of the parallelogram with 
//the horizontal axis.

double Px, Py, Ax, Ay, Ex, Ey, Dx, Dy, Z2right, Z2left;

// PAED
if(sub_parallelogram_num == 0)
{
	Px = mu2-mu1+delta2-delta1;
	Py = c-mu1-delta1;
	Ax = mu2-mu1+delta2;
	Ay = c-mu1;
	Ex = mu2-mu1;
	Ey = c-mu1;
	Dx = mu2-mu1-delta1;
	Dy = c-mu1-delta1;
	// intersection of AP with with horizontal axis (Z2). 
	Z2right = mu2-c+delta2;
	// intersection of ED with with horizontal axis (Z2). 
	Z2left = mu2-c;
}

// AQBE

// Following is the mapping

// A -> P
// Q -> A
// B -> E
// E -> D

else if(sub_parallelogram_num == 1)
{
	Px = mu2-mu1+delta2;
	Py = c-mu1;
	Ax = mu2-mu1+delta1+delta2;
	Ay = c-mu1+delta1;
	Ex = mu2-mu1+delta1;
	Ey = c-mu1+delta1;
	Dx = mu2-mu1;
	Dy = c-mu1;
	// intersection of AP with with horizontal axis (Z2). 
	Z2right = mu2-c+delta2;
	// intersection of ED with with horizontal axis (Z2). 
	Z2left = mu2-c;
}


// EBRC
// Following is the mapping

// E -> P
// B -> A
// R -> E
// C -> D

else if(sub_parallelogram_num == 2)
{	
	Px = mu2-mu1;
	Py = c-mu1;
	Ax = mu2-mu1+delta1;
	Ay = c-mu1+delta1;
	Ex = mu2-mu1+delta1-delta2;
	Ey = c-mu1+delta1;
	Dx = mu2-mu1-delta2;
	Dy = c-mu1;
	// intersection of AP with with horizontal axis (Z2). 
	Z2right = mu2-c;
	// intersection of ED with with horizontal axis (Z2). 
	Z2left = mu2-c-delta2;
}

// DECS
// Following is the mapping

// D -> P
// E -> A
// C -> E
// S -> D

else if(sub_parallelogram_num == 3)
{
	Px = mu2-mu1-delta1;
	Py = c-mu1-delta1;
	Ax = mu2-mu1;
	Ay = c-mu1;
	Ex = mu2-mu1-delta2;
	Ey = c-mu1;
	Dx = mu2-mu1-delta1-delta2;
	Dy = c-mu1-delta1;
	// intersection of AP with with horizontal axis (Z2). 
	Z2right = mu2-c;
	// intersection of ED with with horizontal axis (Z2). 
	Z2left = mu2-c-delta2;	
}

//        e -  -  - a
//        /         /
//       /   P1    /
//      /         /
//     d- - - - - p

//% P1 is a polynomial.


//-------------------------------------------------------------------------
 // sets of 4 positive vertices 
 
// Slopes of P, A, E, D all are greater than 0. There are 5 possible parallelogram configurations for
// when all slopes are are greater than 0.

//  If whole parallelogram lies in the first quadrant   
if (ge(p,0) && ge(Px,0) && ge(Py,0) && ge(a,0) && ge(Ax,0) && ge(Ay,0) && ge(e,0) && ge(Ex,0) && ge(Ey,0) && ge(d,0) && ge(Dx,0) && ge(Dy,0))

    e_path_1(p,a,e,d,exp,var,crossprob);

//  If whole parallelogram lies in the third quadrant
else if (ge(p,0) && le(Px,0) && le(Py,0) && ge(a,0) && le(Ax,0) && le(Ay,0) && ge(e,0) && le(Ex,0) && le(Ey,0) && ge(d,0) && le(Dx,0) && le(Dy,0))
    
    e_path_2(p,a,e,d,exp,var,crossprob);
    
//  If E,A lie in the first quadrant and D,P lie in the third quadrant assuming edge AP crosses second quadrant;
else if (ge(p,0) && le(Px,0)  && le(Py,0) && ge(a,0) && ge(Ax,0) && ge(Ay,0) && ge(e,0) && ge(Ex,0) && ge(Ey,0) && ge(d,0) && le(Dx,0) && le(Dy,0) && le(Z2right,0))
    
    e_path_3(p,a,e,d,exp,var,crossprob);
    
//  If E,A lie in the first quadrant and D,P lie in the third quadrant assuming edge ED crosses fourh quadrant;
else if (ge(p,0) && le(Px,0) && le(Py,0) && ge(a,0) && ge(Ax,0) && ge(Ay,0) && ge(e,0) && ge(Ex,0) && ge(Ey,0) && ge(d,0) && le(Dx,0) && le(Dy,0) && ge(Z2left,0))

    e_path_4(p,a,e,d,exp,var,crossprob);	
   
//  If E,A lie in the first quadrant and D,P lie in the third quadrant assuming edge ED crosses second quadrant
//  and AP crosses the fourth quadrant;
else if (ge(p,0) && le(Px,0) && le(Py,0) && ge(a,0) && ge(Ax,0) && ge(Ay,0) && ge(e,0) && ge(Ex,0) && ge(Ey,0) && ge(d,0) && le(Dx,0) && le(Dy,0) && le(Z2left,0) && ge(Z2right,0))

    e_path_5(p,a,e,d,exp,var,crossprob);
    
//------------------------------------------------------------------------
    // sets of 3 positive vertices 
    
    // Slopes of P, A, D are greater than 0 and slope of E is less than 0. There are 2 possible parallelogram configurations
    // in this case.
    
    // If D, P lie in the third quadrant,A in first and E in second assuming edge ED
    // and AP cross the second quadrant;
else if (ge(p,0) && le(Px,0) && le(Py,0) && ge(a,0) && ge(Ax,0) && ge(Ay,0) && le(e,0) && le(Ex,0) && ge(Ey,0) && ge(d,0) && le(Dx,0)  && le(Dy,0) && le(Z2right,0))
    
    e_path_6(p,a,e,d,exp,var,crossprob);
    
    // If D, P lie in the third quadrant,A in first and E in second assuming edge ED crosses second quadrant
    // and AP crosses the fourth quadrant;
else if (ge(p,0) && le(Px,0) && le(Py,0) && ge(a,0) && ge(Ax,0) && ge(Ay,0) && le(e,0) && le(Ex,0) && ge(Ey,0) && ge(d,0) && le(Dx,0) && le(Dy,0) && le(Z2left,0) && ge(Z2right,0))
    
    e_path_7(p,a,e,d,exp,var,crossprob);
    
    // Slopes of A, E, D are greater than 0 and slope of P is less than 0. There are 2 possible parallelogram configurations
    //in this case.
    
    // If E, A lie in the first quadrant,D in third and P in fourth assuming edge ED
    //and AP cross the fourth quadrant;
else if (le(p,0) && ge(Px,0) && le(Py,0) && ge(a,0) && ge(Ax,0) && ge(Ay,0) && ge(e,0) && ge(Ex,0) && ge(Ey,0) && ge(d,0) && le(Dx,0) && le(Dy,0) && ge(Z2left,0))
    
    e_path_8(p,a,e,d,exp,var,crossprob);
    
    // If E, A lie in the first quadrant,D in third and P in fourth assuming edge ED crosses second quadrant
    //and AP crosses the fourth quadrant;
else if (le(p,0) && ge(Px,0) && le(Py,0) && ge(a,0) && ge(Ax,0) && ge(Ay,0) && ge(e,0) && ge(Ex,0) && ge(Ey,0) && ge(d,0) && le(Dx,0)  && le(Dy,0) && le(Z2left,0) && ge(Z2right,0))
    
    e_path_9(p,a,e,d,exp,var,crossprob);
    
    // Newly added after debugging
    // Slopes of P, A, E are greater than 0 and slope of D is less than 0.
    // P, A, E in the first quadrant and D in the second quadrant 
else if (ge(p,0) && ge(Px,0) && ge(Py,0) && ge(a,0) && ge(Ax,0) && ge(Ay,0) && ge(e,0) && ge(Ex,0) && ge(Ey,0) && le(d,0) && le(Dx,0) && ge(Dy,0))
    
    e_path_10(p,a,e,d,exp,var,crossprob);
    
    // Slopes of E, D, P are greater than 0 and slope of A is less than 0.
    // E, D, P in the third quadrant and A in the fourth quadrant     
else if (ge(p,0) && le(Px,0) && le(Py,0) && le(a,0) && ge(Ax,0) && le(Ay,0) && ge(e,0) && le(Ex,0) && le(Ey,0) && ge(d,0) && le(Dx,0) && le(Dy,0))
    
    e_path_11(p,a,e,d,exp,var,crossprob);
    
   //------------------------------------------------------------------------
    // sets of 2 positive vertices 
    // There are 5 possible parallelogram configurations in this case.
    
// E,D in second quadrant and A,P in the first quadrant    
else if (ge(p,0) && ge(Px,0) && ge(Py,0) && ge(a,0) && ge(Ax,0) && ge(Ay,0) && le(e,0) && le(Ex,0)  && ge(Ey,0) && le(d,0) && le(Dx,0)  && ge(Dy,0))
   
    e_path_12(p,a,e,d,exp,var,crossprob);
    
// E,D in third quadrant and A,P in the fourth quadrant        
else if (le(p,0) && ge(Px,0) && le(Py,0) && le(a,0) && ge(Ax,0) && le(Ay,0) && ge(e,0) && le(Ex,0) && le(Ey,0) && ge(d,0) && le(Dx,0) && le(Dy,0))
   
    e_path_13(p,a,e,d,exp,var,crossprob);
    
// E,A in first quadrant and D,P in the fourth quadrant        
else if (le(p,0) && ge(Px,0) && le(Py,0) && ge(a,0) && ge(Ax,0) && ge(Ay,0) && ge(e,0) && ge(Ex,0)  && ge(Ey,0) && le(d,0) && ge(Dx,0) && le(Dy,0))
   
   e_path_14(p,a,e,d,exp,var,crossprob);
    
// E,A in second quadrant and D,P in the third quadrant           
else if (ge(p,0) && le(Px,0) && le(Py,0) && le(a,0) && le(Ax,0) && ge(Ay,0) && le(e,0) && le(Ex,0) && ge(Ey,0) && ge(d,0) && le(Dx,0) && le(Dy,0))
   
    e_path_15(p,a,e,d,exp,var,crossprob);
    
// E in second quadrant, P in fourth quadrant, A in the first quadrant and D in the third quadrant           
else if (le(p,0) && ge(Px,0) && le(Py,0) && ge(a,0) && ge(Ax,0) && ge(Ay,0) && le(e,0) && le(Ex,0)  && ge(Ey,0) && ge(d,0) && le(Dx,0) && le(Dy,0))
   
    e_path_16(p,a,e,d,exp,var,crossprob);
   
  //------------------------------------------------------------------------
    //  1 positive vertex possibilities 
    // There are 2 possible parallelogram configurations in this case.

// A in first quadrant, E,D,P in the second quadrant           
else if (le(p,0) && le(Px,0) && ge(Py,0) && ge(a,0) && ge(Ax,0) && ge(Ay,0) && le(e,0) && le(Ex,0) && ge(Ey,0) && le(d,0) && le(Dx,0)  && ge(Dy,0))
    
    e_path_17(p,a,e,d,exp,var,crossprob);
    
// D in third quadrant, E,A,P in the fourth quadrant           
else if (le(p,0) && ge(Px,0) && le(Py,0) && le(a,0) && ge(Ax,0) && le(Ay,0) && le(e,0) && ge(Ex,0) && le(Ey,0) && ge(d,0) && le(Dx,0)  && le(Dy,0))

    e_path_18(p,a,e,d,exp,var,crossprob);
    
    //------------------------------------------------------------------------
    //  special case where one of the vertex is at origin(0,0) i.e. when
    //  slope of one vertex is NAN
    //  There are 6 possible parallelogram configurations in this case.

    // P at origin, A in first quadrant, E in the second quadrant 
else if (eq(Px,0) && eq(Py,0) && ge(a,0) && ge(Ax,0) && ge(Ay,0) && le(e,0) && le(Ex,0) && ge(Ey,0))
    
     e_path_19(p,a,e,d,exp,var,crossprob);
    
   // P at origin, E,A in first quadrant
else if (eq(Px,0) && eq(Py,0) && ge(a,0) && ge(Ax,0) && ge(Ay,0) && ge(e,0) && ge(Ex,0) && ge(Ey,0)) 
    
     e_path_20(p,a,e,d,exp,var,crossprob);
    
   // A at origin, D,P in the third quadrant 
else if (ge(p,0) && le(Px,0) && le(Py,0) && eq(Ax,0) && eq(Ay,0) && ge(d,0) && le(Dx,0) && le(Dy,0))
    
     e_path_21(p,a,e,d,exp,var,crossprob);
    
   // E at origin, D in the third quadrant and P in the fourth quadrant 
else if (le(p,0) && ge(Px,0) && le(Py,0) && eq(Ex,0) && eq(Ey,0) && ge(d,0) && le(Dx,0) && le(Dy,0))
    
     e_path_22(p,a,e,d,exp,var,crossprob);
    
   // E at origin, D,P in the third quadrant 
else if (ge(p,0) && le(Px,0) && le(Py,0) && eq(Ex,0) && eq(Ey,0) && ge(d,0) && le(Dx,0) && le(Dy,0)) 
    
     e_path_23(p,a,e,d,exp,var,crossprob);
    
   // D at origin, E,A in the first quadrant  
else if (ge(a,0) && ge(Ax,0) && ge(Ay,0) && ge(e,0) && ge(Ex,0) && ge(Ey,0) && eq(Dx,0) && eq(Dy,0))  

     e_path_24(p,a,e,d,exp,var,crossprob);

// Added during debugging
//If one of the e_paths is not taken means subparallelogram lies in the second/fourth quadrant where all slope values
// are negative		
else
{
	// no contribution to expected value, variance and crossing probability from the subparallelogram
	*exp = 0;
	*var = 0;
	*crossprob = 0;
}     	

}

double determine_slope_value(double z1, double z2)
{

	double slope_temp;

	if(fabs(z2) < DBL_EPSILON)
	{
		// Assign NAN value if both z1 and z2 are 0
		if(fabs(z1) < DBL_EPSILON)
		{
			slope_temp = std::numeric_limits<double>::quiet_NaN();
		}
		else if (z1 > 0)
		{
			slope_temp = INFINITY;
		}
		else
		{
			slope_temp = -INFINITY;
		}
	}
	else
	{
		slope_temp = (double)((z1) / (z2));
	}
	
	return slope_temp;

}*/

double z_density_epanechnikov :: alpha_density_Epanechnikov(double m1, double d1, double  m2, double  d2, double isovalue, double*  exp, double*  var, double*  crossprob)
{

	e_mu1 = m1;
	e_delta1 = d1;
	e_mu2 = m2;
	e_delta2 = d2;
	e_c = isovalue;

	//cout<<"\ne_mu1:"<<e_mu1<<" e_delta1:"<<e_delta1<<" e_mu2:"<<e_mu2<<" e_delta2:"<<e_delta2<<" isovalue:"<<c;
	
	double slope_OP, slope_OQ, slope_OR, slope_OS, slope_OA, slope_OB, slope_OC, slope_OD, slope_OE, mean_temp, delta_temp;

	// Assume e_mu2 is always greater than e_mu1
	if (e_mu1 > e_mu2)
	{
		
		mean_temp = e_mu1;
		e_mu1 = e_mu2;
		e_mu2 = mean_temp;

		delta_temp = e_delta1;
		e_delta1 = e_delta2;
		e_delta2 = delta_temp;
	}   

	// If e_c is out of range [e_mu1-e_delta1,e_mu2+e_delta2] return default expected value i.e. 0.5 (Removed this part. Latest!)
	/*double rangemin, rangemax;
	if((e_mu1-e_delta1) < (e_mu2-e_delta2))
		rangemin = (e_mu1-e_delta1);
	else
		rangemin = (e_mu2-e_delta2);
	
	if((e_mu1+e_delta1) > (e_mu2+e_delta2))
		rangemax = (e_mu1+e_delta1);
	else
		rangemax = (e_mu2+e_delta2);
	
	if((e_c < rangemin) || (e_c > rangemax))
			return 0.5;*/


	// Precomputation of slopes and variables which are part of final PDF    

	// Clean up stuff : add else part

//	if ((e_mu2 - e_mu1 + e_delta2 - e_delta1) != 0)
//		slope_OP = (double)((e_c - e_mu1 - e_delta1) / (e_mu2 - e_mu1 + e_delta2 - e_delta1));

	slope_OP = e_determine_slope_value(e_c - e_mu1 - e_delta1, e_mu2 - e_mu1 + e_delta2 - e_delta1);
	
//	if ((e_mu2 - e_mu1 + e_delta2 + e_delta1) != 0) 
//		slope_OQ = (double)((e_c - e_mu1 + e_delta1) / (e_mu2 - e_mu1 + e_delta2 + e_delta1)); 
	
	slope_OQ = e_determine_slope_value(e_c - e_mu1 + e_delta1, e_mu2 - e_mu1 + e_delta2 + e_delta1);

//	if ((e_mu2 - e_mu1 - e_delta2 - e_delta1) != 0)
//		slope_OS = (double)((e_c - e_mu1 - e_delta1) / (e_mu2 - e_mu1 - e_delta2 - e_delta1));

	slope_OS = e_determine_slope_value(e_c - e_mu1 - e_delta1, e_mu2 - e_mu1 - e_delta2 - e_delta1);		

//	if ((e_mu2 - e_mu1 + e_delta1 - e_delta2) != 0)
//		slope_OR = (double)((e_c - e_mu1 + e_delta1) / (e_mu2 - e_mu1 + e_delta1 - e_delta2));    

	slope_OR = e_determine_slope_value(e_c - e_mu1 + e_delta1, e_mu2 - e_mu1 + e_delta1 - e_delta2);
	

	double expected = 0;
	double crossP = 0;
	double secondM = 0; 
	

	// set parameters for computing the polynomials. Use follwing global object in whole program.
	//tri_poly = new triangular_kernel_polynomial(e_mu1,e_delta1,e_mu2,e_delta2,c);
	epan_poly = new epanechnikov_kernel_polynomial(e_mu1,e_delta1,e_mu2,e_delta2,e_c);


	// Parallelogram PQRS (only one polynomial over whole parallelogram)
        e_pdf_piece(slope_OP, slope_OQ, slope_OR, slope_OS, exp, var, crossprob);
	expected += *exp; 
	crossP += *crossprob;
	secondM += *var;   

	return crossP;
}

/*double z_density_triangular :: alpha_density_triangular(double m1, double d1, double  m2, double  d2, double isovalue, double*  exp, double*  var, double*  crossprob)
{

	mu1 = m1;
	delta1 = d1;
	mu2 = m2;
	delta2 = d2;
	c = isovalue;

	//cout<<"\nmu1:"<<mu1<<" delta1:"<<delta1<<" mu2:"<<mu2<<" delta2:"<<delta2<<" isovalue:"<<c;
	
	double slope_OP, slope_OQ, slope_OR, slope_OS, slope_OA, slope_OB, slope_OC, slope_OD, slope_OE, mean_temp, delta_temp;

	// Assume mu2 is always greater than mu1
	if (mu1 > mu2)
	{
		
		mean_temp = mu1;
		mu1 = mu2;
		mu2 = mean_temp;

		delta_temp = delta1;
		delta1 = delta2;
		delta2 = delta_temp;
	}   */

	// If c is out of range [mu1-delta1,mu2+delta2] return default expected value i.e. 0.5 (Removed this part. Latest!)
	/*double rangemin, rangemax;
	if((mu1-delta1) < (mu2-delta2))
		rangemin = (mu1-delta1);
	else
		rangemin = (mu2-delta2);
	
	if((mu1+delta1) > (mu2+delta2))
		rangemax = (mu1+delta1);
	else
		rangemax = (mu2+delta2);
	
	if((c < rangemin) || (c > rangemax))
			return 0.5;*/


	// Precomputation of slopes and variables which are part of final PDF    

	/*// Clean up stuff : add else part

//	if ((mu2 - mu1 + delta2 - delta1) != 0)
//		slope_OP = (double)((c - mu1 - delta1) / (mu2 - mu1 + delta2 - delta1));

	slope_OP = determine_slope_value(c - mu1 - delta1, mu2 - mu1 + delta2 - delta1);
	
//	if ((mu2 - mu1 + delta2 + delta1) != 0) 
//		slope_OQ = (double)((c - mu1 + delta1) / (mu2 - mu1 + delta2 + delta1)); 
	
	slope_OQ = determine_slope_value(c - mu1 + delta1, mu2 - mu1 + delta2 + delta1);

//	if ((mu2 - mu1 - delta2 - delta1) != 0)
//		slope_OS = (double)((c - mu1 - delta1) / (mu2 - mu1 - delta2 - delta1));

	slope_OS = determine_slope_value(c - mu1 - delta1, mu2 - mu1 - delta2 - delta1);		

//	if ((mu2 - mu1 + delta1 - delta2) != 0)
//		slope_OR = (double)((c - mu1 + delta1) / (mu2 - mu1 + delta1 - delta2));    

	slope_OR = determine_slope_value(c - mu1 + delta1, mu2 - mu1 + delta1 - delta2);
	
//	if ((mu2 - mu1 + delta2) != 0)
//		slope_OA = (double)(c-mu1)/(mu2-mu1+delta2);

	slope_OA = determine_slope_value(c-mu1, mu2-mu1+delta2);

//	if ((mu2 - mu1 + delta1) != 0)
//		slope_OB = (double)(c-mu1+delta1)/(mu2-mu1+delta1);

	slope_OB = determine_slope_value(c-mu1+delta1, mu2-mu1+delta1);

//	if ((mu2 - mu1 - delta2) != 0)
//		slope_OC = (double)(c-mu1)/(mu2-mu1-delta2);

	slope_OC = determine_slope_value(c-mu1, mu2-mu1-delta2);

//	if ((mu2 - mu1 - delta1) != 0)
//		slope_OD = (double)(c-mu1-delta1)/(mu2-mu1-delta1);

	slope_OD = determine_slope_value(c-mu1-delta1, mu2-mu1-delta1);

//	if ((mu2 - mu1) != 0)
//		slope_OE = (double)(c-mu1)/(mu2-mu1);

	slope_OE = determine_slope_value(c-mu1, mu2-mu1);

	double expected = 0;
	double crossP = 0;
	double secondM = 0; 
	

	// set parameters for computing the polynomials. Use follwing global object in whole program.
	tri_poly = new triangular_kernel_polynomial(mu1,delta1,mu2,delta2,c);*/
	/*cout<<"\n\nSubparallelogram PAED:"; 
  	cout<< "slope_OP:" << slope_OP << "\n";  
    cout<< "slope_OA:" << slope_OA << "\n";  
    cout<< "slope_OE:" << slope_OE << "\n";  
    cout<< "slope_OD:" << slope_OD << "\n"; */

	// Parallelogram PAED (polynomial 1)
       /* pdf_piece(slope_OP, slope_OA, slope_OE, slope_OD, exp, var, crossprob);
	expected += *exp; 
	crossP += *crossprob;
	secondM += *var;   */
	
	/*cout<< "\n Expected Value:"<< *exp << " ";
	cout<< "\n Crossing Probability:"<<*crossprob << " ";
	cout<< "\n Second Moment:"<<*var << " ";
	cout<< "\n";*/
	
	// Switch to next subparallelogram
	//sub_parallelogram_num++;
		/*cout<<"\n\nSubparallelogram AQBE:"; 
		cout<< "slope_OA:" << slope_OA << "\n";  
    cout<< "slope_OQ:" << slope_OQ << "\n";  
    cout<< "slope_OB:" << slope_OB << "\n";  
    cout<< "slope_OE:" << slope_OE << "\n";  */

    	// Parallelogram AQBE (polynomial 2)
    	/*pdf_piece(slope_OA, slope_OQ, slope_OB, slope_OE, exp, var, crossprob);
	expected += *exp;    
	crossP += *crossprob;
	secondM += *var;  */
	
	/*cout<< "\n Expected Value:"<<*exp << " ";
	cout<< "\n Crossing Probability:"<<*crossprob << " ";
	cout<<  "\n Second Moment:"<<*var << " ";
	cout<< "\n";*/

	// Switch to next subparallelogram
	//sub_parallelogram_num++;
   /* cout<<"\n\nSubparallelogram EBRC:"; 
    cout<< "slope_OE:" << slope_OE << "\n";  
    cout<< "slope_OB:" << slope_OB << "\n";  
    cout<< "slope_OR:" << slope_OR << "\n";  
    cout<< "slope_OC:" << slope_OC << "\n";  */
    	// Parallelogram EBRC (polynomial 3)
    	/*pdf_piece(slope_OE, slope_OB, slope_OR, slope_OC, exp, var, crossprob);
	expected += *exp;  
	crossP += *crossprob;
	secondM += *var;  */
	
	/*cout<< "\n Expected Value:"<<*exp << " ";  
	cout<< "\n Crossing Probability:"<<*crossprob << " ";
	cout<<  "\n Second Moment:"<<*var << " ";
	cout<< "\n";*/

	// Switch to next subparallelogram
	//sub_parallelogram_num++;
   /* cout<<"\n\nSubparallelogram DECS:"; 
    cout<< "slope_OD:" << slope_OD << "\n";  
    cout<< "slope_OE:" << slope_OE << "\n";  
    cout<< "slope_OC:" << slope_OC << "\n";  
    cout<< "slope_OS:" << slope_OS << "\n";  */
    	// Parallelogram DECS (polynomial 4)
    	/*pdf_piece(slope_OD, slope_OE, slope_OC, slope_OS, exp, var, crossprob);
	expected += *exp;
	crossP += *crossprob;
	secondM += *var;*/  
	
	/*cout<< "\n Expected Value:"<<*exp << " ";
	cout<< "\n Crossing Probability:"<<*crossprob << " ";
	cout<<  "\n Second Moment:"<<*var << " ";
	cout<< "\n";*/
	
	/**exp = expected;
	*crossprob = crossP;
	*var = secondM;*/
	
	// Major bug found!
	// Switch back to the first subparallelogram (PAED) in case the function is called again.
	//sub_parallelogram_num = 0;

	//return crossP;
//}

// set number of uniform distributions in kde1 
void z_density_epanechnikov :: setNumEpanechnikovInKde1(int a)
{
	numEpanechnikovInKde1 = a;
}

// get number of uniform distributions in kde2 
int z_density_epanechnikov :: getNumEpanechnikovInKde1()
{
	return numEpanechnikovInKde1;
}

// set number of uniform distributions in kde1 
void z_density_epanechnikov :: setNumEpanechnikovInKde2(int a)
{
	numEpanechnikovInKde2 = a;
}

// get number of uniform distributions in kde2 
int z_density_epanechnikov :: getNumEpanechnikovInKde2()
{
	return numEpanechnikovInKde2;
}

double z_density_epanechnikov :: kde_z_pdf_expected(float* mu1, double h1, float* mu2, double h2, double c)
{
	int numKde1 = getNumEpanechnikovInKde1();
	int numKde2 = getNumEpanechnikovInKde2();

	//cout<<"\nnumKde1:"<< numKde1; 
	//cout<<"\nnumKde2:"<< numKde2; 

	double expected, c_prob, sig, second, first;

	double kde_expected = 0;

	double total_crossing_prob = 0;
	
	double rangemin,rangemax;
	
	//double EPSILON = 0.000000001;

	for(int i=0; i<numKde1; i++)
	{
		for(int j=0; j<numKde2; j++)
		{
			//temp = alpha_pdf(mu1[i], h1, mu2[j], h2, c);
	
			//Compute0To1(temp,&expected,&c_prob,&sig,&second,&first);

			//cout<<"c_prob:"<<c_prob<<"\n";
			
			// If c is in the total range then only call alpha_density_triangular (Latest change)
			if((mu1[i]-h1) < (mu2[j]-h2))
				rangemin = (mu1[i]-h1);
			else
				rangemin = (mu2[j]-h2);
	
			if((mu1[i]+h1) > (mu2[j]+h2))
				rangemax = (mu1[i]+h1);
			else
				rangemax = (mu2[j]+h2);
	
	 		if((c < rangemin) || (c > rangemax))
		 	{
		 		expected = 0;
		 		c_prob =0; 
		 		sig = 0;		 	
			}
			else
		  	{
				//cout<<mu1[i]<<" "<<h1<<" "<<mu2[j] <<" "<< h2<< " "<<c<<"\n";
				alpha_density_Epanechnikov(mu1[i], h1, mu2[j], h2, c,&expected,&sig,&c_prob);
			}

			//cout<< "i:"<<i<<" j:"<<j<<" Expected:"<< expected<<"\n";
			// What happens when crossing probability is 0 skip
			if(mu1[i] > mu2[j])
			{
				if(fabs(c_prob)<DBL_EPSILON)
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
				if(fabs(c_prob)<DBL_EPSILON)
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


double z_density_epanechnikov :: kde_z_pdf_variance(float* mu1, double h1, float* mu2, double h2, double c, double expected_crossing)
{

	int numKde1 = getNumEpanechnikovInKde1();
	int numKde2 = getNumEpanechnikovInKde2();	

	//piecewise temp;

	double expected, c_prob, second;

	double kde_e_second_moment = 0;
	double total_crossing_prob = 0;	

	//double kde_first_moment = kde_z_pdf_expected(mu1, h1, mu2, h2, c);
	double kde_first_moment = expected_crossing;	
	
	double kde_variance = 0;
	
	double rangemin,rangemax;

	// kde second moment can be obtained in following two ways when mu1 > mu2.

	// reverse the pdf using function adjustPieceLimits and compute second moment in the forward direction (may be closed form compution is not possible because of things like log of negative number etc.)

        // or


        // use the following formula with current pdf instead of reversing the pdf
	// second moment in reverse direction = second moment in forward direction + 1 - 2*first moment in forward direction.

	for(int i=0; i<numKde1; i++)
	{
		for(int j=0; j<numKde2; j++)
		{
			//temp = alpha_pdf(mu1[i], h1, mu2[j], h2, c);

			//Compute0To1(temp,&expected,&c_prob,&sig,&second,&first);			
			
			// If c is in the total range then only call alpha_density_triangular (Latest change)
			if((mu1[i]-h1) < (mu2[j]-h2))
				rangemin = (mu1[i]-h1);
			else
				rangemin = (mu2[j]-h2);
	
			if((mu1[i]+h1) > (mu2[j]+h2))
				rangemax = (mu1[i]+h1);
			else
				rangemax = (mu2[j]+h2);
	
	 		if((c < rangemin) || (c > rangemax))
		 	{
		 		expected = 0;
		 		c_prob =0; 
		 		second = 0;		 	
			}
			else
		        {		
				alpha_density_Epanechnikov(mu1[i], h1, mu2[j], h2, c,&expected,&second,&c_prob);
			}

			if(mu1[i] > mu2[j])
			{
					if(fabs(c_prob)<DBL_EPSILON)
					{
						// kde_e_second_moment = kde_e_second_moment + 0.1;
						// Means given pair of kernels doesn't determine variance
					}
					else
					{
						//If pdf is not normalized 1 should be replace with crossing probability	
						kde_e_second_moment = kde_e_second_moment + second + c_prob - 2*expected;
						total_crossing_prob = total_crossing_prob + c_prob;
					}
			}	
			else
			{
					if(fabs(c_prob)<DBL_EPSILON)
					{
						// kde_e_second_moment = kde_e_second_moment + 0.1;
						// Means given pair of kernels doesn't determine variance
					}
					else
					{
						kde_e_second_moment =  kde_e_second_moment + second;
						total_crossing_prob = total_crossing_prob + c_prob;
					}
			}
			
		}
	}

	if(total_crossing_prob > 0)
	
		kde_variance  = (double)(kde_e_second_moment/total_crossing_prob) - kde_first_moment*kde_first_moment;

	else
		
		kde_variance = 0;	
	//cout<<"kde_variance:"<<kde_variance<<"\n";

	return kde_variance;
}

//int main()
//{

//z_density_epanechnikov z;

//float m1 [] = {189.465, 189.683, 189.576, 189.349, 189.62, 189.576, 188.958, 189.482, 189.774, 189.468, 189.424, 189.448, 189.516, 189.25, 189.399, 189.365, 189.463, 189.546, 189.682, 189.026, 188.932, 189.384, 189.52, 189.067, 189.307, 189.762, 189.182, 190.609, 190.932, 190.223, 190.246, 190.018, 190.104, 190.417, 190.189, 190.523, 189.222, 188.983, 189.443, 189.089, 188.83, 188.917, 189.24, 189.376, 189.348, 189.885, 189.48, 190.379, 189.371, 189.882, 189.666, 188.346, 188.8, 190.155, 189.496, 189.337, 189.298, 189.295, 189.402, 189.529, 189.275, 189.295, 189.141};


//float m2 [] = {188.874, 189.052, 189.01, 188.905, 189.019, 189.338, 188.526, 188.999, 188.964, 188.954, 188.914, 188.752, 188.916, 188.931, 188.923, 188.74, 188.963, 189.015, 189.242, 188.632, 188.196, 189.049, 188.69, 188.651, 188.833, 189.203, 188.915, 190.217, 190.518, 189.968, 189.777, 189.531, 189.553, 189.963, 189.957, 190.037, 188.807, 188.188, 188.749, 188.51, 188.572, 188.435, 189.154, 188.896, 188.922, 189.73, 189.405, 189.527, 189.266, 189.591, 189.315, 188.22, 187.866, 189.337, 189.259, 189.052, 189.138, 188.919, 189.171, 189.188, 189.026, 188.905, 188.842};

/*
int num1 = sizeof m1/sizeof (float);
int num2 = sizeof m2/sizeof (float);

double h1 = 0.219641;

double h2 = 0.230373;

double isoval = 189.433;
*/

/*double mu1 = 180, delta1 = 2, mu2 = 190, delta2 = 2, c = 191.5;
double *exp, *var, *crossprob; 

exp = new double;
var = new double;
crossprob = new double;

double cp = z.alpha_density_Epanechnikov(mu1,delta1,mu2,delta2,c,exp,var,crossprob);
cout<<"\nCrossing probability is:" << *crossprob <<"\n";
cout<<"\nExpected crossing location is:" << (*exp)/(*crossprob) <<"\n";*/
/*z.setNumTriangularInKde1(num1);
z.setNumTriangularInKde2(num2);

double e_expected_value = z.kde_z_pdf_expected(m1,h1,m2,h2,isoval);
cout<<"\nExpected crossing location is:" << e_expected_value <<"\n";
double variance = z.kde_z_pdf_variance(m1,h1,m2,h2,isoval,e_expected_value);
cout<<"\nVariance is:" << variance <<"\n";*/

//return 0;
//}

