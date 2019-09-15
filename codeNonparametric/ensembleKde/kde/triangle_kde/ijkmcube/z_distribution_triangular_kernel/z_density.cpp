#include "stdio.h"
#include <iostream>
#include "stdlib.h"
#include "math.h"
//#include "cmath.h"
#include <algorithm>
#include "z_density.h"
#include <cfloat>
#include <limits> 
#include "triangular_kernel_polynomial.h"

//# define infinity 10000
//# define minus_infinity -10000
#define _USE_MATH_DEFINES

using namespace std;

double mu1, delta1, mu2, delta2, c;

// 0 : PAED, 1: AQBE, 2: EBRC, 3: DECS
// set initially to PAED
int sub_parallelogram_num = 0;

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

triangular_kernel_polynomial* tri_poly;

double integrate_pdf(double low, double high, int piecenum)
{
	
	// assuming low is always greater than 0 and  low<=high

	if (sub_parallelogram_num == 0)
	{

		if (high - low < DBL_EPSILON)
			return 0;
		else if (low >= 1)
			return 0;
		else if ((low <= 1) && (high <= 1))
			return (tri_poly -> PAED_integrate_piece_value(high,piecenum) - tri_poly -> PAED_integrate_piece_value(low,piecenum));
		else if ((low <= 1) && (high >= 1))		
			return (tri_poly -> PAED_integrate_piece_value(1,piecenum) - tri_poly -> PAED_integrate_piece_value(low,piecenum));
	}
	else if (sub_parallelogram_num == 1)
	{
		if (high - low < DBL_EPSILON)
			return 0;
		else if (low >= 1)
			return 0;
		else if ((low <= 1) && (high <= 1))
			return (tri_poly -> AQBE_integrate_piece_value(high,piecenum) - tri_poly -> AQBE_integrate_piece_value(low,piecenum));
		else if ((low <= 1) && (high >= 1))		
			return (tri_poly -> AQBE_integrate_piece_value(1,piecenum) - tri_poly -> AQBE_integrate_piece_value(low,piecenum));
	}
	else if (sub_parallelogram_num == 2)
	{
		if (high - low < DBL_EPSILON)
			return 0;
		else if (low >= 1)
			return 0;
		else if ((low <= 1) && (high <= 1))
			return (tri_poly -> EBRC_integrate_piece_value(high,piecenum) - tri_poly -> EBRC_integrate_piece_value(low,piecenum));
		else if ((low <= 1) && (high >= 1))		
			return (tri_poly -> EBRC_integrate_piece_value(1,piecenum) - tri_poly -> EBRC_integrate_piece_value(low,piecenum));
	}
	else if (sub_parallelogram_num == 3)
	{
		if (high - low < DBL_EPSILON)
			return 0;
		else if (low >= 1)
			return 0;
		else if ((low <= 1) && (high <= 1))
			return (tri_poly -> DECS_integrate_piece_value(high,piecenum) - tri_poly -> DECS_integrate_piece_value(low,piecenum));
		else if ((low <= 1) && (high >= 1))		
			return (tri_poly -> DECS_integrate_piece_value(1,piecenum) - tri_poly -> DECS_integrate_piece_value(low,piecenum));
	}

}

double expected_value(double low, double high, int piecenum)
{	
	// assuming low is always greater than 0 and  low<=high
	if (sub_parallelogram_num == 0)
	{
		if (high - low < DBL_EPSILON)
			return 0;
		else if (low >= 1)
			return 0;
		else if ((low <= 1) && (high <= 1))
			return tri_poly -> PAED_expected_piece_value(high,piecenum) - tri_poly -> PAED_expected_piece_value(low,piecenum);
		else if ((low <= 1) && (high >= 1))		
			return tri_poly -> PAED_expected_piece_value(1,piecenum) - tri_poly -> PAED_expected_piece_value(low,piecenum);
	}
	else if (sub_parallelogram_num == 1)
	{
		if (high - low < DBL_EPSILON)
			return 0;
		else if (low >= 1)
			return 0;
		else if ((low <= 1) && (high <= 1))
			return tri_poly -> AQBE_expected_piece_value(high,piecenum) - tri_poly -> AQBE_expected_piece_value(low,piecenum);
		else if ((low <= 1) && (high >= 1))		
			return tri_poly -> AQBE_expected_piece_value(1,piecenum) - tri_poly -> AQBE_expected_piece_value(low,piecenum);
	}
	else if (sub_parallelogram_num == 2)
	{
		if (high - low < DBL_EPSILON)
			return 0;
		else if (low >= 1)
			return 0;
		else if ((low <= 1) && (high <= 1))
			return tri_poly -> EBRC_expected_piece_value(high,piecenum) - tri_poly -> EBRC_expected_piece_value(low,piecenum);
		else if ((low <= 1) && (high >= 1))		
			return tri_poly -> EBRC_expected_piece_value(1,piecenum) - tri_poly -> EBRC_expected_piece_value(low,piecenum);
	}
	else if (sub_parallelogram_num == 3)
	{
		if (high - low < DBL_EPSILON)
			return 0;
		else if (low >= 1)
			return 0;
		else if ((low <= 1) && (high <= 1))
			return tri_poly -> DECS_expected_piece_value(high,piecenum) - tri_poly -> DECS_expected_piece_value(low,piecenum);
		else if ((low <= 1) && (high >= 1))		
			return tri_poly -> DECS_expected_piece_value(1,piecenum) - tri_poly -> DECS_expected_piece_value(low,piecenum);
	}

	
}

double second_moment(double low, double high, int piecenum)
{	
	// assuming low is always greater than 0 and  low<=high

	if (sub_parallelogram_num == 0)
	{
		if (high - low < DBL_EPSILON)
			return 0;
		else if (low >= 1)
			return 0;
		else if ((low <= 1) && (high <= 1))
			return tri_poly -> PAED_second_moment_piece_value(high,piecenum) - tri_poly -> PAED_second_moment_piece_value(low,piecenum);
		else if ((low <= 1) && (high >= 1))		
			return tri_poly -> PAED_second_moment_piece_value(1,piecenum) - tri_poly -> PAED_second_moment_piece_value(low,piecenum);
	}
	else if (sub_parallelogram_num == 1)
	{
		if (high - low < DBL_EPSILON)
			return 0;
		else if (low >= 1)
			return 0;
		else if ((low <= 1) && (high <= 1))
			return tri_poly -> AQBE_second_moment_piece_value(high,piecenum) - tri_poly -> AQBE_second_moment_piece_value(low,piecenum);
		else if ((low <= 1) && (high >= 1))		
			return tri_poly -> AQBE_second_moment_piece_value(1,piecenum) - tri_poly -> AQBE_second_moment_piece_value(low,piecenum);
	}
	else if (sub_parallelogram_num == 2)
	{
		if (high - low < DBL_EPSILON)
			return 0;
		else if (low >= 1)
			return 0;
		else if ((low <= 1) && (high <= 1))
			return tri_poly -> EBRC_second_moment_piece_value(high,piecenum) - tri_poly -> EBRC_second_moment_piece_value(low,piecenum);
		else if ((low <= 1) && (high >= 1))		
			return tri_poly -> EBRC_second_moment_piece_value(1,piecenum) - tri_poly -> EBRC_second_moment_piece_value(low,piecenum);
	}
	else if (sub_parallelogram_num == 3)
	{
		if (high - low < DBL_EPSILON)
			return 0;
		else if (low >= 1)
			return 0;
		else if ((low <= 1) && (high <= 1))
			return tri_poly -> DECS_second_moment_piece_value(high,piecenum) - tri_poly -> DECS_second_moment_piece_value(low,piecenum);
		else if ((low <= 1) && (high >= 1))		
			return tri_poly -> DECS_second_moment_piece_value(1,piecenum) - tri_poly -> DECS_second_moment_piece_value(low,piecenum);
	}

}

double path_1(double p,double a,double e,double d, double* exp, double* var, double *crossprob)
{    
    
	*exp = 0;
	*var = 0;
	*crossprob = 0;
	

	// Vertex Order PADE 
	if ((p <= a) && (a <= d) && (d <= e))
	{   		
		
	*exp = expected_value(p,a,1) + expected_value(a,d,5) - expected_value(d,e,3);
	*var = second_moment(p,a,1) + second_moment(a,d,5) - second_moment(d,e,3);
	*crossprob = integrate_pdf(p,a,1) + integrate_pdf(a,d,5) - integrate_pdf(d,e,3);

	}
	
	// Vertex Order PAED 
	else if ((p <= a) && (a <= e) && (e <= d))
 	{

	*exp = expected_value(p,a,1) + expected_value(a,e,5) - expected_value(e,d,4);
	*var = second_moment(p,a,1) + second_moment(a,e,5) - second_moment(e,d,4);
	*crossprob = integrate_pdf(p,a,1) + integrate_pdf(a,e,5) - integrate_pdf(e,d,4);
	
	}

	// Vertex Order PDAE 
	else if ((p <= d) && (d <= a) && (a <= e))
	{    	
    		
	*exp = expected_value(p,d,1) + expected_value(d,a,6) - expected_value(a,e,3);
	*var = second_moment(p,d,1) + second_moment(d,a,6) - second_moment(a,e,3);
	*crossprob = integrate_pdf(p,d,1) + integrate_pdf(d,a,6) - integrate_pdf(a,e,3);

	}

	
}

double path_2(double p,double a,double e,double d, double* exp, double* var, double *crossprob)
{
	*exp = 0;
	*var = 0;
	*crossprob = 0;

	// Vertex Order EDAP (polynomial 1)
	if ((e <= d) && (d <= a) && (a <= p))
    	{
   	 	
	*exp = expected_value(e,d,3) - expected_value(d,a,5) - expected_value(a,p,1);
	*var = second_moment(e,d,3) - second_moment(d,a,5) - second_moment(a,p,1);
	*crossprob = integrate_pdf(e,d,3) - integrate_pdf(d,a,5) - integrate_pdf(a,p,1);
	
	}
    
    // Vertex Order EDPA (polynomial 1)
	else if ((e <= d) && (d <= p) && (p <= a))
    	{
      
	*exp = expected_value(e,d,3) - expected_value(d,p,5) - expected_value(p,a,2);
	*var = second_moment(e,d,3) - second_moment(d,p,5) - second_moment(p,a,2);
	*crossprob = integrate_pdf(e,d,3) - integrate_pdf(d,p,5) - integrate_pdf(p,a,2);

     	}
    // Vertex Order EADP (polynomial 1)
	else if ((e <= a) && (a <= d) && (d <= p))
    	{
    
	*exp = expected_value(e,a,3) - expected_value(a,d,6) - expected_value(d,p,1);
	*var = second_moment(e,a,3) - second_moment(a,d,6) - second_moment(d,p,1);
	*crossprob = integrate_pdf(e,a,3) - integrate_pdf(a,d,6) - integrate_pdf(d,p,1);	
	
	}	
}

double path_3(double p,double a,double e,double d, double* exp, double* var, double *crossprob)
{
	*exp = 0;
	*var = 0;
	*crossprob = 0;

	/* Only 1 possible ordering DPAE
 All equations derived from PAED_positive_1
 first and the fourth part should be identical*/

	*exp = -expected_value(0,d,6) - expected_value(d,p,1) + expected_value(a,e,2) - expected_value(e,INFINITY,6);
	*var = -second_moment(0,d,6) - second_moment(d,p,1) + second_moment(a,e,2) - second_moment(e,INFINITY,6);
	*crossprob = -integrate_pdf(0,d,6) - integrate_pdf(d,p,1) + integrate_pdf(a,e,2) - integrate_pdf(e,INFINITY,6);	
}

double path_4(double p,double a,double e,double d, double* exp, double* var, double *crossprob)
{
	*exp = 0;
	*var = 0;
	*crossprob = 0;
	
	/* Only 1 possible ordering AEDP
% All equations derived from PAED_positive_1
% first and the fourth part should be identical*/
	
	*exp = expected_value(0,a,6) - expected_value(a,e,3) + expected_value(d,p,4) + expected_value(p,INFINITY,6);	
	*var = second_moment(0,a,6) - second_moment(a,e,3) + second_moment(d,p,4) + second_moment(p,INFINITY,6);	
	*crossprob = integrate_pdf(0,a,6) - integrate_pdf(a,e,3) + integrate_pdf(d,p,4) + integrate_pdf(p,INFINITY,6);		
}

double path_5(double p,double a,double e,double d, double* exp, double* var, double *crossprob)
{
	*exp = 0;
	*var = 0;
	*crossprob = 0;
	
	// Vertex Order DAPE 
	if ((d <= a) && (a <= p) && (p <= e))
	{  

	*exp = expected_value(0,d,7) + expected_value(0,d,8) + expected_value(d,a,7) + expected_value(d,a,10) + expected_value(a,p,9)+ expected_value(a,p,10) + expected_value(p,e,9)+ expected_value(p,e,12) + expected_value(e,INFINITY,11)+ expected_value(e,INFINITY,12);	
	*var = second_moment(0,d,7) + second_moment(0,d,8) + second_moment(d,a,7) + second_moment(d,a,10) + second_moment(a,p,9)+ second_moment(a,p,10) + second_moment(p,e,9)+ second_moment(p,e,12) + second_moment(e,INFINITY,11)+ second_moment(e,INFINITY,12);	
	*crossprob = integrate_pdf(0,d,7) + integrate_pdf(0,d,8) + integrate_pdf(d,a,7) + integrate_pdf(d,a,10) + integrate_pdf(a,p,9)+ integrate_pdf(a,p,10) + integrate_pdf(p,e,9)+ integrate_pdf(p,e,12) + integrate_pdf(e,INFINITY,11)+ integrate_pdf(e,INFINITY,12);	

	}
    
    	// Vertex Order ADPE (polynomial 1)
	else if ((a <= d) && (d <= p) && (p <= e))
	{
	*exp = expected_value(0,a,7) + expected_value(0,a,8) + expected_value(a,d,9) + expected_value(a,d,8) + expected_value(d,p,9)+ expected_value(d,p,10) + expected_value(p,e,9)+ expected_value(p,e,12) + expected_value(e,INFINITY,11)+ expected_value(e,INFINITY,12);
	*var = second_moment(0,a,7) + second_moment(0,a,8) + second_moment(a,d,9) + second_moment(a,d,8) + second_moment(d,p,9)+ second_moment(d,p,10) + second_moment(p,e,9)+ second_moment(p,e,12) + second_moment(e,INFINITY,11)+ second_moment(e,INFINITY,12);
	*crossprob = integrate_pdf(0,a,7) + integrate_pdf(0,a,8) + integrate_pdf(a,d,9) + integrate_pdf(a,d,8) + integrate_pdf(d,p,9)+ integrate_pdf(d,p,10) + integrate_pdf(p,e,9)+ integrate_pdf(p,e,12) + integrate_pdf(e,INFINITY,11)+ integrate_pdf(e,INFINITY,12);

	}
    
	// Vertex Order DAEP (polynomial 1)
	else if ((d <= a) && (a <= e) && (e <= p))
	{  

	*exp = expected_value(0,d,7) + expected_value(0,d,8) + expected_value(d,a,7) + expected_value(d,a,10) + expected_value(a,e,9)+ expected_value(a,e,10) + expected_value(e,p,11)+ expected_value(e,p,10) + expected_value(p,INFINITY,11)+ expected_value(p,INFINITY,12);
	*var = second_moment(0,d,7) + second_moment(0,d,8) + second_moment(d,a,7) + second_moment(d,a,10) + second_moment(a,e,9)+ second_moment(a,e,10) + second_moment(e,p,11)+ second_moment(e,p,10) + second_moment(p,INFINITY,11)+ second_moment(p,INFINITY,12);
	*crossprob = integrate_pdf(0,d,7) + integrate_pdf(0,d,8) + integrate_pdf(d,a,7) + integrate_pdf(d,a,10) + integrate_pdf(a,e,9)+ integrate_pdf(a,e,10) + integrate_pdf(e,p,11)+ integrate_pdf(e,p,10) + integrate_pdf(p,INFINITY,11)+ integrate_pdf(p,INFINITY,12);

	}
   
   	// Vertex Order ADEP (polynomial 1)
	else if ((a <= d) && (d <= e) && (e <= p))
	{   
    
	*exp = expected_value(0,a,7) + expected_value(0,a,8) + expected_value(a,d,9) + expected_value(a,d,8) + expected_value(d,e,9)+ expected_value(d,e,10) + expected_value(e,p,11)+ expected_value(e,p,10) + expected_value(p,INFINITY,11)+ expected_value(p,INFINITY,12);
	*var = second_moment(0,a,7) + second_moment(0,a,8) + second_moment(a,d,9) + second_moment(a,d,8) + second_moment(d,e,9)+ second_moment(d,e,10) + second_moment(e,p,11)+ second_moment(e,p,10) + second_moment(p,INFINITY,11)+ second_moment(p,INFINITY,12);
	*crossprob = integrate_pdf(0,a,7) + integrate_pdf(0,a,8) + integrate_pdf(a,d,9) + integrate_pdf(a,d,8) + integrate_pdf(d,e,9)+ integrate_pdf(d,e,10) + integrate_pdf(e,p,11)+ integrate_pdf(e,p,10) + integrate_pdf(p,INFINITY,11)+ integrate_pdf(p,INFINITY,12);

	}	
}

double path_6(double p,double a,double e,double d, double* exp, double* var, double *crossprob)
{
	*exp = 0;
	*var = 0;
	*crossprob = 0;

//% only one possible vertex ordering DPA (Ignore -ve vertices)

	*exp = -expected_value(0,d,6) - expected_value(d,p,1) + expected_value(a,INFINITY,2);	
	*var = -second_moment(0,d,6) - second_moment(d,p,1) + second_moment(a,INFINITY,2);	
	*crossprob = -integrate_pdf(0,d,6) - integrate_pdf(d,p,1) + integrate_pdf(a,INFINITY,2);		
}

double path_7(double p,double a,double e,double d, double* exp, double* var, double *crossprob)
{
	*exp = 0;
	*var = 0;
	*crossprob = 0;

	// vertex order DAP
	if ((d <= a) && (a <= p))
	{  
	
	*exp = expected_value(0,d,7) + expected_value(0,d,8) + expected_value(d,a,7) + expected_value(d,a,10) + expected_value(a,p,9)+ expected_value(a,p,10) + expected_value(p,INFINITY,9)+ expected_value(p,INFINITY,12);
	*var = second_moment(0,d,7) + second_moment(0,d,8) + second_moment(d,a,7) + second_moment(d,a,10) + second_moment(a,p,9)+ second_moment(a,p,10) + second_moment(p,INFINITY,9)+ second_moment(p,INFINITY,12);
	*crossprob = integrate_pdf(0,d,7) + integrate_pdf(0,d,8) + integrate_pdf(d,a,7) + integrate_pdf(d,a,10) + integrate_pdf(a,p,9)+ integrate_pdf(a,p,10) + integrate_pdf(p,INFINITY,9)+ integrate_pdf(p,INFINITY,12);

	}
    
    	// vertex order ADP
	else if ((a <= d) && (d <= p))
	{    
   
	*exp = expected_value(0,a,7) + expected_value(0,a,8) + expected_value(a,d,9) + expected_value(a,d,8) + expected_value(d,p,9)+ expected_value(d,p,10) + expected_value(p,INFINITY,9)+ expected_value(p,INFINITY,12);
	*var = second_moment(0,a,7) + second_moment(0,a,8) + second_moment(a,d,9) + second_moment(a,d,8) + second_moment(d,p,9)+ second_moment(d,p,10) + second_moment(p,INFINITY,9)+ second_moment(p,INFINITY,12);
	*crossprob = integrate_pdf(0,a,7) + integrate_pdf(0,a,8) + integrate_pdf(a,d,9) + integrate_pdf(a,d,8) + integrate_pdf(d,p,9)+ integrate_pdf(d,p,10) + integrate_pdf(p,INFINITY,9)+ integrate_pdf(p,INFINITY,12);
        
	}
	
}

double path_8(double p,double a,double e,double d, double* exp, double* var, double *crossprob)
{
	*exp = 0;
	*var = 0;
	*crossprob = 0;
	
	// only one possible vertex ordering AED (Ignore -ve vertices)
	*exp = expected_value(0,a,6) - expected_value(a,e,3) + expected_value(d,INFINITY,4);	
	*var = second_moment(0,a,6) - second_moment(a,e,3) + second_moment(d,INFINITY,4);	
	*crossprob = integrate_pdf(0,a,6) - integrate_pdf(a,e,3) + integrate_pdf(d,INFINITY,4);		
}

double path_9(double p,double a,double e,double d, double* exp, double* var, double *crossprob)
{
	*exp = 0;
	*var = 0;
	*crossprob = 0;
	
	// vertex order ADE
	if ((a<=d) && (d<=e))
	{    
   
	*exp = expected_value(0,a,7) + expected_value(0,a,8) + expected_value(a,d,9) + expected_value(a,d,8) + expected_value(d,e,9)+ expected_value(d,e,10) + expected_value(e,INFINITY,11)+ expected_value(e,INFINITY,10);
	*var = second_moment(0,a,7) + second_moment(0,a,8) + second_moment(a,d,9) + second_moment(a,d,8) + second_moment(d,e,9)+ second_moment(d,e,10) + second_moment(e,INFINITY,11)+ second_moment(e,INFINITY,10);
	*crossprob = integrate_pdf(0,a,7) + integrate_pdf(0,a,8) + integrate_pdf(a,d,9) + integrate_pdf(a,d,8) + integrate_pdf(d,e,9)+ integrate_pdf(d,e,10) + integrate_pdf(e,INFINITY,11)+ integrate_pdf(e,INFINITY,10);

	}
    
	// vertex order DAE
	else if ((d<=a) && (a<=e))
	{
    
	*exp = expected_value(0,d,7) + expected_value(0,d,8) + expected_value(d,a,7) + expected_value(d,a,10) + expected_value(a,e,9)+ expected_value(a,e,10) + expected_value(e,INFINITY,11)+ expected_value(e,INFINITY,10);
	*var = second_moment(0,d,7) + second_moment(0,d,8) + second_moment(d,a,7) + second_moment(d,a,10) + second_moment(a,e,9)+ second_moment(a,e,10) + second_moment(e,INFINITY,11)+ second_moment(e,INFINITY,10);
	*crossprob = integrate_pdf(0,d,7) + integrate_pdf(0,d,8) + integrate_pdf(d,a,7) + integrate_pdf(d,a,10) + integrate_pdf(a,e,9)+ integrate_pdf(a,e,10) + integrate_pdf(e,INFINITY,11)+ integrate_pdf(e,INFINITY,10);
	
	}	
}

double path_10(double p,double a,double e,double d, double* exp, double* var, double *crossprob)
{
	*exp = 0;
	*var = 0;
	*crossprob = 0;

	// vertex order PAE
	if ((p <= a) && (a <= e))
	{
 
	*exp = expected_value(p,a,1) + expected_value(a,e,5) - expected_value(e,INFINITY,4);	
	*var = second_moment(p,a,1) + second_moment(a,e,5) - second_moment(e,INFINITY,4);	
	*crossprob = integrate_pdf(p,a,1) + integrate_pdf(a,e,5) - integrate_pdf(e,INFINITY,4);	

	}
    
	// vertex order APE    
	else if ((a <= p) && (p <= e)) 
	{
    
	*exp = expected_value(a,p,2) + expected_value(p,e,5) - expected_value(e,INFINITY,4);	
	*var = second_moment(a,p,2) + second_moment(p,e,5) - second_moment(e,INFINITY,4);	
	*crossprob = integrate_pdf(a,p,2) + integrate_pdf(p,e,5) - integrate_pdf(e,INFINITY,4);	
	}    

	// vertex order AEP    
	else if ((a <= e) && (e <= p))  
	{    
	*exp = expected_value(a,e,2) - expected_value(e,p,6) - expected_value(p,INFINITY,4);	
	*var = second_moment(a,e,2) - second_moment(e,p,6) - second_moment(p,INFINITY,4);	
	*crossprob = integrate_pdf(a,e,2) - integrate_pdf(e,p,6) - integrate_pdf(p,INFINITY,4);		
	}	
}

double path_11(double p,double a,double e,double d, double* exp, double* var, double *crossprob)
{
	*exp = 0;
	*var = 0;
	*crossprob = 0;
	
	// vertex order EDP
	if ((e <= d) && (d <= p))
	{
  
	*exp = expected_value(e,d,3) - expected_value(d,p,5) - expected_value(p,INFINITY,2);	
	*var = second_moment(e,d,3) - second_moment(d,p,5) - second_moment(p,INFINITY,2);		
	*crossprob = integrate_pdf(e,d,3) - integrate_pdf(d,p,5) - integrate_pdf(p,INFINITY,2);			
 		
	}    

	// vertex order DEP    
	else if ((d <= e) && (e <= p))  
	{  
   
	*exp = expected_value(d,e,4) - expected_value(e,p,5) - expected_value(p,INFINITY,2);	
	*var = second_moment(d,e,4) - second_moment(e,p,5) - second_moment(p,INFINITY,2);	
	*crossprob = integrate_pdf(d,e,4) - integrate_pdf(e,p,5) - integrate_pdf(p,INFINITY,2);	

	}
    
	// vertex order DPE    
	else if ((d <= p) && (p <= e))  
	{	     
    
	*exp = expected_value(d,p,4) + expected_value(p,e,6) - expected_value(e,INFINITY,2);	
	*var = second_moment(d,p,4) + second_moment(p,e,6) - second_moment(e,INFINITY,2);	
	*crossprob = integrate_pdf(d,p,4) + integrate_pdf(p,e,6) - integrate_pdf(e,INFINITY,2);	
	}	
}

double path_12(double p,double a,double e,double d, double* exp, double* var, double *crossprob)
{

	*exp = 0;
	*var = 0;
	*crossprob = 0;

	// Vertex order PA
	if (p <= a)
	{
  
	*exp = expected_value(p,a,1) + expected_value(a,INFINITY,5);
	*var = second_moment(p,a,1) + second_moment(a,INFINITY,5);
	*crossprob = integrate_pdf(p,a,1) + integrate_pdf(a,INFINITY,5);

	}
    
	// Vertex order AP
	else if (a <= p)
	{    
  
	*exp = expected_value(a,p,2) + expected_value(p,INFINITY,5);
	*var = second_moment(a,p,2) + second_moment(p,INFINITY,5);
	*crossprob = integrate_pdf(a,p,2) + integrate_pdf(p,INFINITY,5);

	} 
	
}

double path_13(double p,double a,double e,double d, double* exp, double* var, double *crossprob)
{
	*exp = 0;
	*var = 0;
	*crossprob = 0;
	
	// Vertex order ED
	if (e <= d)
	{
  
	*exp = expected_value(e,d,3) - expected_value(d,INFINITY,5);
	*var = second_moment(e,d,3) - second_moment(d,INFINITY,5);
	*crossprob = integrate_pdf(e,d,3) - integrate_pdf(d,INFINITY,5);

	}
    
	// Vertex order DE
	else if (d <= e)
	{   
  
	*exp = expected_value(d,e,4) - expected_value(e,INFINITY,5);
	*var = second_moment(d,e,4) - second_moment(e,INFINITY,5);
	*crossprob = integrate_pdf(d,e,4) - integrate_pdf(e,INFINITY,5);
    
	}	
}

double path_14(double p,double a,double e,double d, double* exp, double* var, double *crossprob)
{
	*exp = 0;
	*var = 0;
	*crossprob = 0;
	
	// vertex order AE (only one possible)  
	*exp = expected_value(0,a,6) - expected_value(a,e,3);
	*var = second_moment(0,a,6) - second_moment(a,e,3);
	*crossprob = integrate_pdf(0,a,6) - integrate_pdf(a,e,3);	
}

double path_15(double p,double a,double e,double d, double* exp, double* var, double *crossprob)
{

	*exp = 0;
	*var = 0;
	*crossprob = 0;
	
 	// vertex order DP (only one possible)
	*exp = -expected_value(0,d,6) - expected_value(d,p,1);	
	*var = -second_moment(0,d,6) - second_moment(d,p,1);	
	*crossprob = -integrate_pdf(0,d,6) - integrate_pdf(d,p,1);	
	
}

double path_16(double p,double a,double e,double d, double* exp, double* var, double *crossprob)
{
	*exp = 0;
	*var = 0;
	*crossprob = 0;
	
	// Vertex order DA
	if (d <= a)
	{
   
	*exp = expected_value(0,d,7) + expected_value(0,d,8) + expected_value(d,a,7) + expected_value(d,a,10) + expected_value(a,INFINITY,9) + expected_value(a,INFINITY,10);
	*var = second_moment(0,d,7) + second_moment(0,d,8) + second_moment(d,a,7) + second_moment(d,a,10) + second_moment(a,INFINITY,9) + second_moment(a,INFINITY,10);
	*crossprob = integrate_pdf(0,d,7) + integrate_pdf(0,d,8) + integrate_pdf(d,a,7) + integrate_pdf(d,a,10) + integrate_pdf(a,INFINITY,9) + integrate_pdf(a,INFINITY,10);
  
	}
  
	// Vertex order AD
	else if (a <= d)
	{    
 
	*exp = expected_value(0,a,7) + expected_value(0,a,8) + expected_value(a,d,9) + expected_value(a,d,8) + expected_value(d,INFINITY,9) + expected_value(d,INFINITY,10);
	*var = second_moment(0,a,7) + second_moment(0,a,8) + second_moment(a,d,9) + second_moment(a,d,8) + second_moment(d,INFINITY,9) + second_moment(d,INFINITY,10);
	*crossprob = integrate_pdf(0,a,7) + integrate_pdf(0,a,8) + integrate_pdf(a,d,9) + integrate_pdf(a,d,8) + integrate_pdf(d,INFINITY,9) + integrate_pdf(d,INFINITY,10);
    
	}	
}

double path_17(double p,double a,double e,double d, double* exp, double* var, double *crossprob)
{

	*exp = 0;
	*var = 0;
	*crossprob = 0;
	
	*exp = expected_value(a,INFINITY,2);	
	*var = second_moment(a,INFINITY,2);	
	*crossprob = integrate_pdf(a,INFINITY,2);	
	
}

double path_18(double p,double a,double e,double d, double* exp, double* var, double *crossprob)
{

	*exp = 0;
	*var = 0;
	*crossprob = 0;
	
	*exp = expected_value(d,INFINITY,4);	
	*var = second_moment(d,INFINITY,4);	
	*crossprob = integrate_pdf(d,INFINITY,4);	
	
}

double path_19(double p,double a,double e,double d, double* exp, double* var, double *crossprob)
{

	*exp = 0;
	*var = 0;
	*crossprob = 0;

	*exp = expected_value(a,INFINITY,2);	
	*var = second_moment(a,INFINITY,2);	
	*crossprob = integrate_pdf(a,INFINITY,2);		

}

double path_20(double p,double a,double e,double d, double* exp, double* var, double *crossprob)
{

	*exp = 0;
	*var = 0;
	*crossprob = 0;
	
	*exp = expected_value(a,e,2) - expected_value(e,INFINITY,6);	
	*var = second_moment(a,e,2) - second_moment(e,INFINITY,6);	
	*crossprob = integrate_pdf(a,e,2) - integrate_pdf(e,INFINITY,6);		
}

double path_21(double p,double a,double e,double d, double* exp, double* var, double *crossprob)
{

	*exp = 0;
	*var = 0;
	*crossprob = 0;

	// vertex order DP
	*exp = expected_value(0,d,3) - expected_value(d,p,1);
	*var = second_moment(0,d,3) - second_moment(d,p,1);
	*crossprob = integrate_pdf(0,d,3) - integrate_pdf(d,p,1);	
}

double path_22(double p,double a,double e,double d, double* exp, double* var, double *crossprob)
{

	*exp = 0;
	*var = 0;
	*crossprob = 0;

	*exp = expected_value(d,INFINITY,4);	
	*var = second_moment(d,INFINITY,4);	
	*crossprob = integrate_pdf(d,INFINITY,4);		
}

double path_23(double p,double a,double e,double d, double* exp, double* var, double *crossprob)
{
	*exp = 0;
	*var = 0;
	*crossprob = 0;
	
	// vertex order DP

	*exp = expected_value(d,p,4) + expected_value(p,INFINITY,6);
	*var = second_moment(d,p,4) + second_moment(p,INFINITY,6);
	*crossprob = integrate_pdf(d,p,4) + integrate_pdf(p,INFINITY,6);	
}

double path_24(double p,double a,double e,double d, double* exp, double* var, double *crossprob)
{

	*exp = 0;
	*var = 0;
	*crossprob = 0;
	
	// vertex order AE
	*exp = expected_value(0,a,1) - expected_value(a,e,3);
	*var = second_moment(0,a,1) - second_moment(a,e,3);
	*crossprob = integrate_pdf(0,a,1) - integrate_pdf(a,e,3);	
}

void pdf_piece(double p,double a,double e, double d, double* exp, double* var, double* crossprob)
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

/* Following is the mapping

 A -> P
 Q -> A
 B -> E
 E -> D*/

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
/* Following is the mapping

 E -> P
 B -> A
 R -> E
 C -> D*/

else if(sub_parallelogram_num == 2)
{	
	//cout<< "p:" << p << "\n";  
    //cout<< "a:" << a << "\n";  
    //cout<< "e:" << e << "\n";  
    //cout<< "d:" << d << "\n";  

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
/* Following is the mapping

 D -> P
 E -> A
 C -> E
 S -> D*/

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




/*        e -  -  - a
        /         /
       /   P1    /
      /         /
     d- - - - - p

% P1 is a polynomial.*/


//-------------------------------------------------------------------------
 // sets of 4 positive vertices 
 
// Slopes of P, A, E, D all are greater than 0. There are 5 possible parallelogram configurations for
// when all slopes are are greater than 0.

//  If whole parallelogram lies in the first quadrant


if ((p >= 0) && (Px >= 0) && (Py >= 0) && (a >= 0) && (Ax >= 0) && (Ay >= 0) && (e >= 0) && (Ex >= 0) && (Ey >= 0) && (d >= 0) && (Dx >= 0) && (Dy >= 0)) 
{
    
    path_1(p,a,e,d,exp,var,crossprob);

}
//  If whole parallelogram lies in the third quadrant
else if ((p >= 0) && (Px <= 0) && (Py <= 0) && (a >= 0) && (Ax <= 0) && (Ay <= 0) && (e >= 0) && (Ex <= 0) && (Ey <= 0)&& (d >= 0) && (Dx <= 0) && (Dy <= 0))
    
    path_2(p,a,e,d,exp,var,crossprob);
    
//  If E,A lie in the first quadrant and D,P lie in the third quadrant assuming edge AP crosses second quadrant;
else if ((p >= 0) && (Px <= 0)  && (Py <= 0) && (a >= 0) && (Ax >= 0)  && (Ay >= 0) && (e >= 0) && (Ex >= 0)  && (Ey >= 0) && (d >= 0) && (Dx <= 0) && (Dy <= 0) && (Z2right <= 0))
    
    path_3(p,a,e,d,exp,var,crossprob);
    
//  If E,A lie in the first quadrant and D,P lie in the third quadrant assuming edge ED crosses fourh quadrant;
else if ((p >= 0) && (Px <= 0)  && (Py <= 0) && (a >= 0) && (Ax >= 0)  && (Ay >= 0) && (e >= 0) && (Ex >= 0)  && (Ey >= 0) && (d >= 0) && (Dx <= 0) && (Dy <= 0) && (Z2left >= 0))
    
    path_4(p,a,e,d,exp,var,crossprob);
    
//  If E,A lie in the first quadrant and D,P lie in the third quadrant assuming edge ED crosses second quadrant
//  and AP crosses the fourth quadrant;
else if ((p >= 0) && (Px <= 0)  && (Py <= 0) && (a >= 0) && (Ax >= 0)  && (Ay >= 0) && (e >= 0) && (Ex >= 0)  && (Ey >= 0) && (d >= 0) && (Dx <= 0) && (Dy <= 0) && (Z2left <= 0) && (Z2right >= 0))
    
    path_5(p,a,e,d,exp,var,crossprob);
    
//------------------------------------------------------------------------
    // sets of 3 positive vertices 
    
    // Slopes of P, A, D are greater than 0 and slope of E is less than 0. There are 2 possible parallelogram configurations
    // in this case.
    
    // If D, P lie in the third quadrant,A in first and E in second assuming edge ED
    // and AP cross the second quadrant;
    
else if ((p >= 0) && (Px <= 0) && (Py <= 0) && (a >= 0) && (Ax >= 0) && (Ay >= 0) && (e <= 0) && (Ex <= 0)  && (Ey >= 0) && (d >= 0) && (Dx <= 0)  && (Dy <= 0) && (Z2right <= 0))
    
    path_6(p,a,e,d,exp,var,crossprob);
    
    // If D, P lie in the third quadrant,A in first and E in second assuming edge ED crosses second quadrant
    // and AP crosses the fourth quadrant;
    
else if ((p >= 0) && (Px <= 0) && (Py <= 0) && (a >= 0) && (Ax >= 0) && (Ay >= 0) && (e <= 0) && (Ex <= 0)  && (Ey >= 0) && (d >= 0) && (Dx <= 0)  && (Dy <= 0) && (Z2left <= 0) && (Z2right >= 0))
    
    path_7(p,a,e,d,exp,var,crossprob);
    
    // Slopes of A, E, D are greater than 0 and slope of P is less than 0. There are 2 possible parallelogram configurations
    //in this case.
    
    // If E, A lie in the first quadrant,D in third and P in fourth assuming edge ED
    //and AP cross the fourth quadrant;
    
else if ((p <= 0) && (Px >= 0) && (Py <= 0) && (a >= 0) && (Ax >= 0) && (Ay >= 0) && (e >= 0) && (Ex >= 0)  && (Ey >= 0) && (d >= 0) && (Dx <= 0)  && (Dy <= 0) && (Z2left >= 0))
    
    path_8(p,a,e,d,exp,var,crossprob);
    
    // If E, A lie in the first quadrant,D in third and P in fourth assuming edge ED crosses second quadrant
    //and AP crosses the fourth quadrant;
    
else if ((p <= 0) && (Px >= 0) && (Py <= 0) && (a >= 0) && (Ax >= 0) && (Ay >= 0) && (e >= 0) && (Ex >= 0)  && (Ey >= 0) && (d >= 0) && (Dx <= 0)  && (Dy <= 0) && (Z2left <= 0) && (Z2right >= 0))
    
    path_9(p,a,e,d,exp,var,crossprob);
    
    // Newly added after debugging
    // Slopes of P, A, E are greater than 0 and slope of D is less than 0.
    // P, A, E in the first quadrant and D in the second quadrant 
    
else if ((p >= 0) && (Px >= 0) && (Py >= 0) && (a >= 0) && (Ax >= 0) && (Ay >= 0) && (e >= 0) && (Ex >= 0)  && (Ey >= 0) && (d <= 0) && (Dx <= 0)  && (Dy >= 0))
    
    path_10(p,a,e,d,exp,var,crossprob);
    
    // Slopes of E, D, P are greater than 0 and slope of A is less than 0.
    // E, D, P in the third quadrant and A in the fourth quadrant     
else if ((p >= 0) && (Px <= 0) && (Py <= 0) && (a <= 0) && (Ax >= 0) && (Ay <= 0) && (e >= 0) && (Ex <= 0)  && (Ey <= 0) && (d >= 0) && (Dx <= 0)  && (Dy <= 0))
    
    path_11(p,a,e,d,exp,var,crossprob);
    
   //------------------------------------------------------------------------
    // sets of 2 positive vertices 
    // There are 5 possible parallelogram configurations in this case.
    
// E,D in second quadrant and A,P in the first quadrant    
else if ((p >= 0) && (Px >= 0) && (Py >= 0) && (a >= 0) && (Ax >= 0) && (Ay >= 0) && (e <= 0) && (Ex <= 0)  && (Ey >= 0) && (d <= 0) && (Dx <= 0)  && (Dy >= 0))
   
    path_12(p,a,e,d,exp,var,crossprob);
    
// E,D in third quadrant and A,P in the fourth quadrant        
else if ((p <= 0) && (Px >= 0) && (Py <= 0) && (a <= 0) && (Ax >= 0) && (Ay <= 0) && (e >= 0) && (Ex <= 0)  && (Ey <= 0) && (d >= 0) && (Dx <= 0)  && (Dy <= 0))
   
    path_13(p,a,e,d,exp,var,crossprob);
    
// E,A in first quadrant and D,P in the fourth quadrant        
else if ((p <= 0) && (Px >= 0) && (Py <= 0) && (a >= 0) && (Ax >= 0) && (Ay >= 0) && (e >= 0) && (Ex >= 0)  && (Ey >= 0) && (d <= 0) && (Dx >= 0)  && (Dy <= 0))
   
   path_14(p,a,e,d,exp,var,crossprob);
    
// E,A in second quadrant and D,P in the third quadrant           
else if ((p >= 0) && (Px <= 0) && (Py <= 0) && (a <= 0) && (Ax <= 0) && (Ay >= 0) && (e <= 0) && (Ex <= 0)  && (Ey >= 0) && (d >= 0) && (Dx <= 0)  && (Dy <= 0))
   
    path_15(p,a,e,d,exp,var,crossprob);
    
// E in second quadrant, P in fourth quadrant, A in the first quadrant and D in the third quadrant           
else if ((p <= 0) && (Px >= 0) && (Py <= 0) && (a >= 0) && (Ax >= 0) && (Ay >= 0) && (e <= 0) && (Ex <= 0)  && (Ey >= 0) && (d >= 0) && (Dx <= 0)  && (Dy <= 0))
   
    path_16(p,a,e,d,exp,var,crossprob);
   
  //------------------------------------------------------------------------
    //  1 positive vertex possibilities 
    // There are 2 possible parallelogram configurations in this case.

// A in first quadrant, E,D,P in the second quadrant           
else if ((p <= 0) && (Px <= 0) && (Py >= 0) && (a >= 0) && (Ax >= 0) && (Ay >= 0) && (e <= 0) && (Ex <= 0)  && (Ey >= 0) && (d <= 0) && (Dx <= 0)  && (Dy >= 0))
    
    path_17(p,a,e,d,exp,var,crossprob);
    
// D in third quadrant, E,A,P in the fourth quadrant           
else if ((p <= 0) && (Px >= 0) && (Py <= 0) && (a <= 0) && (Ax >= 0) && (Ay <= 0) && (e <= 0) && (Ex >= 0)  && (Ey <= 0) && (d >= 0) && (Dx <= 0)  && (Dy <= 0))

    path_18(p,a,e,d,exp,var,crossprob);
    
    //------------------------------------------------------------------------
    //  special case where one of the vertex is at origin(0,0) i.e. when
    //  slope of one vertex is NAN
    //  There are 6 possible parallelogram configurations in this case.

    // P at origin, A in first quadrant, E in the second quadrant 
else if ((Px == 0) && (Py == 0) && (a >= 0) && (Ax >= 0) && (Ay >= 0) && (e <= 0) && (Ex <= 0)  && (Ey >= 0))
    
     path_19(p,a,e,d,exp,var,crossprob);
    
   // P at origin, E,A in first quadrant
else if ((Px == 0) && (Py == 0) && (a >= 0) && (Ax >= 0) && (Ay >= 0) && (e >= 0) && (Ex >= 0)  && (Ey >= 0)) 
    
     path_20(p,a,e,d,exp,var,crossprob);
    
   // A at origin, D,P in the third quadrant 
else if ((p >= 0) && (Px <= 0) && (Py <= 0) && (Ax == 0) && (Ay == 0) && (d >= 0) && (Dx <= 0)  && (Dy <= 0))
    
     path_21(p,a,e,d,exp,var,crossprob);
    
   // E at origin, D in the third quadrant and P in the fourth quadrant 
else if ((p <= 0) && (Px >= 0) && (Py <= 0) && (Ex == 0) && (Ey == 0) && (d >= 0) && (Dx <= 0)  && (Dy <= 0))
    
     path_22(p,a,e,d,exp,var,crossprob);
    
   // E at origin, D,P in the third quadrant 
else if ((p >= 0) && (Px <= 0) && (Py <= 0) && (Ex == 0) && (Ey == 0) && (d >= 0) && (Dx <= 0)  && (Dy <= 0)) 
    
     path_23(p,a,e,d,exp,var,crossprob);
    
   // D at origin, E,A in the first quadrant 
else if ((a >= 0) && (Ax >= 0) && (Ay >= 0) && (e >= 0) && (Ex >= 0)  && (Ey >= 0) && (Dx == 0) && (Dy == 0))  

     path_24(p,a,e,d,exp,var,crossprob);

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

}

double alpha_density_triangular(double m1, double d1, double  m2, double  d2, double isovalue, double*  exp, double*  var, double*  crossprob)
{

	mu1 = m1;
	delta1 = d1;
	mu2 = m2;
	delta2 = d2;
	c = isovalue;
	
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
	}   

	// If c is out of range [mu1-delta1,mu2+delta2] return default expected value i.e. 0.5
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
			return 0.5;


	// Precomputation of slopes and variables which are part of final PDF    

	// Clean up stuff : add else part

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
	tri_poly = new triangular_kernel_polynomial(mu1,delta1,mu2,delta2,c);
	cout<<"\n\nSubparallelogram PAED:"; 
  	cout<< "slope_OP:" << slope_OP << "\n";  
    cout<< "slope_OA:" << slope_OA << "\n";  
    cout<< "slope_OE:" << slope_OE << "\n";  
    cout<< "slope_OD:" << slope_OD << "\n"; 

	// Parallelogram PAED (polynomial 1)
        pdf_piece(slope_OP, slope_OA, slope_OE, slope_OD, exp, var, crossprob);
	expected += *exp; 
	crossP += *crossprob;
	secondM += *var;   
	
	cout<< "\n Expected Value:"<< *exp << " ";
	cout<< "\n Crossing Probability:"<<*crossprob << " ";
	cout<< "\n Second Moment:"<<*var << " ";
	cout<< "\n";
	
	// Switch to next subparallelogram
	sub_parallelogram_num++;
		cout<<"\n\nSubparallelogram AQBE:"; 
		cout<< "slope_OA:" << slope_OA << "\n";  
    cout<< "slope_OQ:" << slope_OQ << "\n";  
    cout<< "slope_OB:" << slope_OB << "\n";  
    cout<< "slope_OE:" << slope_OE << "\n";  

    	// Parallelogram AQBE (polynomial 2)
    	pdf_piece(slope_OA, slope_OQ, slope_OB, slope_OE, exp, var, crossprob);
	expected += *exp;    
	crossP += *crossprob;
	secondM += *var;  
	
	cout<< "\n Expected Value:"<<*exp << " ";
	cout<< "\n Crossing Probability:"<<*crossprob << " ";
	cout<<  "\n Second Moment:"<<*var << " ";
	cout<< "\n";

	// Switch to next subparallelogram
	sub_parallelogram_num++;
    cout<<"\n\nSubparallelogram EBRC:"; 
    cout<< "slope_OE:" << slope_OE << "\n";  
    cout<< "slope_OB:" << slope_OB << "\n";  
    cout<< "slope_OR:" << slope_OR << "\n";  
    cout<< "slope_OC:" << slope_OC << "\n";  
    	// Parallelogram EBRC (polynomial 3)
    	pdf_piece(slope_OE, slope_OB, slope_OR, slope_OC, exp, var, crossprob);
	expected += *exp;  
	crossP += *crossprob;
	secondM += *var;  
	
	cout<< "\n Expected Value:"<<*exp << " ";  
	cout<< "\n Crossing Probability:"<<*crossprob << " ";
	cout<<  "\n Second Moment:"<<*var << " ";
	cout<< "\n";

	// Switch to next subparallelogram
	sub_parallelogram_num++;
    cout<<"\n\nSubparallelogram DECS:"; 
    cout<< "slope_OD:" << slope_OD << "\n";  
    cout<< "slope_OE:" << slope_OE << "\n";  
    cout<< "slope_OC:" << slope_OC << "\n";  
    cout<< "slope_OS:" << slope_OS << "\n";  
    	// Parallelogram DECS (polynomial 4)
    	pdf_piece(slope_OD, slope_OE, slope_OC, slope_OS, exp, var, crossprob);
	expected += *exp;
	crossP += *crossprob;
	secondM += *var;  
	
	cout<< "\n Expected Value:"<<*exp << " ";
	cout<< "\n Crossing Probability:"<<*crossprob << " ";
	cout<<  "\n Second Moment:"<<*var << " ";
	cout<< "\n";
	
	*exp = expected;
	*crossprob = crossP;
	*var = secondM;
	
	return crossP;
}

/*int main()
{

double mu1 = 3, delta1 = 1, mu2 = 10, delta2 = 2, c = 6;
double *exp, *var, *crossprob; 

exp = new double;
var = new double;
crossprob = new double;

double cp = alpha_density_triangular(mu1,delta1,mu2,delta2,c,exp,var,crossprob);
cout<<"\nCrossing probability is:" << *crossprob <<"\n";
cout<<"\nExpected crossing location is:" << (*exp)/(*crossprob) <<"\n";
return 0;
}*/
