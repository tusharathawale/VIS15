/// \file ijkinterpolate.txx
/// ijk templates for linear and multilinear interpolation
/// Version 0.1.0

/*
  IJK: Isosurface Jeneration Kode
  Copyright (C) 2009 Rephael Wenger

  This library is free software; you can redistribute it and/or
  modify it under the terms of the GNU Lesser General Public License
  (LGPL) as published by the Free Software Foundation; either
  version 2.1 of the License, or (at your option) any later version.

  This library is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
  Lesser General Public License for more details.

  You should have received a copy of the GNU Lesser General Public
  License along with this library; if not, write to the Free Software
  Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
*/

#ifndef _IJKINTERPOLATE_
#define _IJKINTERPOLATE_

#include "ijkcoord.txx"
#include "kdeWithUniformKernel.h"
#include "kdeWithTriangleKernel.h"
#include "kdeWithEpanechnikovKernel.h"

//#include "z_density.h"
//#include "even_z_density.h"

namespace IJK {

// **********************************************************
// Linear and multilinear interpolation functions.
// **********************************************************

  /// Compute coordinates using linear interpolation.
  template <class DTYPE, class WTYPE,
	    class CTYPE0, class CTYPE1, class CTYPE2>
  inline void linear_interpolate_coord
  (const DTYPE dimension, const WTYPE w0, const CTYPE0 coord0,
   const CTYPE1 coord1, CTYPE2 coord2)
  {
    const WTYPE w1 = 1-w0;

    for (int d = 0; d < dimension; d++)
      { coord2[d] = w0*coord0[d] + w1*coord1[d]; }
  }  

  /// Newly added. Compute coordinates using linear interpolation.
  template <class DTYPE, class WTYPE,
	    class CTYPE0, class CTYPE1, class CTYPE2>
  inline void linear_interpolate_coord_second
  (const DTYPE dimension, const WTYPE w0, const CTYPE0 coord0,
   const CTYPE1 coord1, CTYPE2 coord2)
  {
    const WTYPE w1 = 1-w0;

    for (int d = 0; d < dimension; d++)
      { coord2[d] = w0*coord1[d] + w1*coord0[d]; }
  }  
	

  /// Compute coordinates using linear interpolation.
  /// @param PTYPE = Precision type to be used in calculations. Should be float of double.
  template <class PTYPE, class DTYPE, class STYPE, class VTYPE,
	    class CTYPE0, class CTYPE1, class CTYPE2>
  inline void linear_interpolate_coord
  (const DTYPE dimension, const STYPE s0, const CTYPE0 coord0,
   const STYPE s1, const CTYPE1 coord1,
   const VTYPE isovalue, CTYPE2 coord2)
  {
    PTYPE w0;
    STYPE s_diff = s1 - s0;
    const PTYPE EPSILON = 0.00001;
    if (s_diff > EPSILON || s_diff < -EPSILON) { 
      w0 = (s1 - isovalue) / s_diff;
    }
    else {
      // arbitrarily set weights to 0.5
      w0 = 0.5;
    };

    linear_interpolate_coord
      (dimension, w0, coord0, coord1, coord2);
  }


   /// Compute coordinates using linear interpolation.
  template <class DTYPE, class STYPE, class VTYPE,
	    class CTYPE0, class CTYPE1, class CTYPE2>
  inline void linear_interpolate_coord
  (const DTYPE dimension, const STYPE s0, const CTYPE0 coord0,
   const STYPE s1, const CTYPE1 coord1,
   const VTYPE isovalue, CTYPE2 coord2)
  {
    linear_interpolate_coord<double>
      (dimension, s0, coord0, s1, coord1, isovalue, coord2);
  } 	

   /// Compute coordinates using linear interpolation.
  template <class DTYPE, class STYPE, class VTYPE,
	    class CTYPE0, class CTYPE1, class CTYPE2>
  inline void alpha_uncertainty_coord
  (const DTYPE dimension, const STYPE s0, const CTYPE0 coord0,
   const STYPE s1, const CTYPE1 coord1,
   const VTYPE isovalue, CTYPE2 coord2, const STYPE s2, const STYPE s3, float* max_variance, float* color)
  {
    alpha_uncertainty_coord<double>
      (dimension, s0, coord0, s1, coord1, isovalue, coord2, s2, s3,max_variance,color);
  } 	

  /// change!
/// Compute coordinates using alpha uncertainty.
  /// @param PTYPE = Precision type to be used in calculations. Should be float of double.
  template <class PTYPE, class DTYPE, class STYPE, class VTYPE,
	    class CTYPE0, class CTYPE1, class CTYPE2>
  inline void kde_alpha_uncertainty_coord
  (const DTYPE dimension, STYPE* s0, const CTYPE0 coord0,
   STYPE* s1, const CTYPE1 coord1,
   const VTYPE isovalue, CTYPE2 coord2, const STYPE s2, const STYPE s3, float* max_variance, float* color, int numKernels)
  {
  
  		PTYPE w0;
  		PTYPE s2_new, s3_new;
    	double lin;
    	double m0=0,m1=0;   	
    	STYPE m_diff;
      const PTYPE EPSILON = 0.00001;
       
    			z_density_uniform ad;      	
			z_density_triangular tri;
			z_density_epanechnikov ep;
    	piecewise pic;  
    	double expected, c_prob, sig;
       		
    	//cout<<"\nnumber of kernels:"<<numKernels;
			tri.setNumTriangularInKde1(numKernels); 
			tri.setNumTriangularInKde2(numKernels); 	
			
			// use uniform kernels for the expected value and variance computation when triangular returns -ve value for the expected crossing.
			// This can happen because when s2 and s3 are small (below 0.6 in general from observations), the greens theorem computation becomes unstable.
			
			// ****************** IMP obseravation (for code triangle.cxx and traingular_kernel_polynomial.cxx of vis sub 1 following snippet is necessary) ***********************
		/*	// For tangle dataset many times variance is less than 0.6. So in vis14 submission 1 the kde result obtained is most of the times with uniform kernel 
			ad.setNumUniformsInKde1(numKernels);
      ad.setNumUniformsInKde2(numKernels);
				
			// use uniform kernels for the expected value and variance computation when triangular returns -ve value for the expected crossing.
			if(s2 < 0.6 || s3 < 0.6)
      {
      	//	cout<<"\ns2:"<<s2<<" s3:"<<s3;
      	//  w0 = 0.5;
     		w0 = ad.kde_alpha_pdf_expected(s0,s2,s1,s3,isovalue);    	 		
     	 	sig = ad.kde_alpha_pdf_variance(s0,s2,s1,s3,isovalue);      
      }*/
      
      //************************************************************************************************************************************************************

			// get expected crossing assuming triangular base.
			//else
			//{
			
			
				   // triangular kernel.
				   //Bandwidth was stored assuming the Gaussian kernel in MATLAB script. Thus, rescale bandwidth to triangular kernel (Refer to Table 1 in Appendix).
				   s2_new = (double)(s2*2.4320); 
				   s3_new = (double)(s3*2.4320); 
                               
                                   // uniform kernel.
				   //Bandwidth was stored assuming the Gaussian kernel in MATLAB script. Thus, rescale bandwidth to triangular kernel.
				   //s2_new = (double)(s2*1.7401); 
				   //s3_new = (double)(s3*1.7401); 

				   // epanechnikov kernel.
				   //Bandwidth was stored assuming the Gaussian kernel in MATLAB script. Thus, rescale bandwidth to triangular kernel.
				   //s2_new = (double)(s2*2.2138); 
				   //s3_new = (double)(s3*2.2138); 	

			
				w0 = tri.kde_z_pdf_expected(s0,s2_new,s1,s3_new,isovalue); 
				sig = tri.kde_z_pdf_variance(s0,s2_new,s1,s3_new,isovalue,w0);		

				// If w0 value is 0.5 and sig is 0 means we explicitly set these values since crossing probability, expected value and variance values are 0.
				// In this case we find the w0 with inverse linear interpolation of means and set sig to 0.
				// Code for linear interpolation in a case where crossing probability, expected value and variance values are 0 (This is when parallelogram size is very small).				
				if ((fabs(w0-0.5)<0.00000001) && (fabs(sig)<0.00000001))
				{
    	 			for(int i=0; i<numKernels; i++)
    	 			{
    	 				m0+=s0[i];
    	    		m1+=s1[i];
    	 			}
    	 			m0 = m0/(double)numKernels;
    	 			m1 = m1/(double)numKernels;
						
						m_diff = m1 - m0;
					  if (m_diff > EPSILON || m_diff < -EPSILON) { 
      				w0 = (m1 - isovalue) / m_diff;
    				}
    				else {
      			// arbitrarily set weights to 0.5
      			w0 = 0.5;
    				};
				}


				/*if (s2<0.6 || s3<0.6)
				{
	                		cout<<"\ns2:"<<s2<<" s3:"<<s3<< " w0:"<<w0;
				}*/

			//}
      if ((w0 < 0.00000001)||(w0 > 0.99999999))
      {
	    
		 cout<<"\ns0:";
     for(int i=0; i<numKernels; i++)
     	cout<<s0[i]<<" ";     	
     	
     	cout<<"\ns2:"<<s2;
     	
     cout<<"\ns1:";
     for(int i=0; i<numKernels; i++)
     	cout<<s1[i]<<" ";	
     	
     	cout<<"\ns3:"<<s3;

	cout<<"\nw0:"<<w0; 
    }

    /* cout<<"\nsig:"<<sig;*/

		// find maximum variance
		
    if(sig > *max_variance)
			*max_variance = sig; 		
		    
    // store unnormalized variance	
      *color = sig;    
    	  
    linear_interpolate_coord_second(dimension, w0, coord0, coord1, coord2);
}  
  	  
  
   /// Compute coordinates using linear interpolation.
  template <class DTYPE, class STYPE, class VTYPE,
	    class CTYPE0, class CTYPE1, class CTYPE2>
  inline void kde_alpha_uncertainty_coord
  (const DTYPE dimension, STYPE* s0, const CTYPE0 coord0,
   STYPE* s1, const CTYPE1 coord1,
   const VTYPE isovalue, CTYPE2 coord2, const STYPE s2, const STYPE s3, float* max_variance, float* color, int numKernels)
  {
    kde_alpha_uncertainty_coord<double>
      (dimension, s0, coord0, s1, coord1, isovalue, coord2, s2, s3,max_variance,color,numKernels);
  } 	

 

  /// Multilinear interpolate scalar values at unit cube vertices.
  /// @tparam ITYPE = Type of interpolated scalar value.
  /// @tparam WTYPE = Weight type.
  template <class ITYPE, class WTYPE,
	    class DTYPE, class CTYPE, class NTYPE, class STYPE>
  ITYPE multilinear_interpolate_scalar
  (const DTYPE dimension, const CTYPE coord,
   const NTYPE num_cube_vertices, const STYPE * cube_scalar)
  {
    ITYPE s = 0;
    for (NTYPE i = 0; i < num_cube_vertices; i++) {
      WTYPE weight = 1.0;
      NTYPE mask = 1;
      for (DTYPE d = 0; d < dimension; d++) {
	if ((i & mask) == 0) { weight = weight * (1-coord[d]); }
	else { weight = weight*coord[d]; };

	mask = (mask << 1L);
      }
      s += weight * cube_scalar[i];
    }
    return(s);
  }

  /// Multilinear interpolate scalar values at unit cube vertices.
  template <class DTYPE, class CTYPE, class NTYPE, class STYPE>
  double multilinear_interpolate_scalar
  (const DTYPE dimension, const CTYPE coord,
   const NTYPE num_cube_vertices, const STYPE * cube_scalar)
  {
    multilinear_interpolate_scalar<double, double>
      (dimension, coord, num_cube_vertices, cube_scalar);
  }

  namespace {

    template <class STYPE, class VTYPE>
    inline bool get_vertex_sign
    (const STYPE s, const VTYPE isovalue)
    {
      if (s < isovalue) { return(false); }
      else { return(true); };
    }

  };


  /// Compute location of isosurface vertex on line segment [A,B].
  /// Use multilinear interpolation on unit cube vertices to calculate scalar values.
  /// @tparam ITYPE = Type of interpolated scalar value.
  /// @tparam WTYPE = Weight type.
  template <class ITYPE, class WTYPE, class DTYPE, 
	    class CTYPE, class VTYPE, class NTYPE, class STYPE>
  void multilinear_interpolate_coord
  (const DTYPE dimension, const CTYPE * coordA, const CTYPE * coordB,
   const VTYPE isovalue, const NTYPE num_cube_vertices,
   const STYPE * cube_scalar, const NTYPE num_iter, CTYPE * coordC)
  {

    ITYPE s0 = multilinear_interpolate_scalar<ITYPE,WTYPE>
      (dimension, coordA, num_cube_vertices, cube_scalar);
    ITYPE s1 = multilinear_interpolate_scalar<ITYPE,WTYPE>
      (dimension, coordB, num_cube_vertices, cube_scalar);

    bool sign0 = get_vertex_sign(s0, isovalue);
    bool sign1 = get_vertex_sign(s1, isovalue);

    if (sign0 == sign1) {
      linear_interpolate_coord(dimension, 0.5, coordA, coordB, coordC);
      return;
    }

    if (s0 == isovalue) {
      IJK::copy_coord(dimension, coordA, coordC);
      return;
    }

    if (s1 == isovalue) {
      IJK::copy_coord(dimension, coordB, coordC);
      return;
    }

    WTYPE t0 = 0.0;
    WTYPE t1 = 1.0;
    WTYPE t2;
    ITYPE s2;
    for (NTYPE i = 0; i  < num_iter; i++) {
      t2 = (t0+t1)/2.0;
      linear_interpolate_coord(dimension, 1.0-t2, coordA, coordB, coordC);
      s2 = multilinear_interpolate_scalar<ITYPE,WTYPE>
	(dimension, coordC, num_cube_vertices, cube_scalar);
      bool sign2 = get_vertex_sign(s2, isovalue);

      if (sign0 == sign2) { t0 = t2; }
      else { t1 = t2; }
    }

  }


  // Compute location of isosurface vertex on line segment [A,B].
  // Use multilinear interpolation on unit cube vertices to calculate scalar values.
  template <class DTYPE, 
	    class CTYPE, class VTYPE, class NTYPE, class STYPE>
  void multilinear_interpolate_coord
  (const DTYPE dimension, const CTYPE * coordA, const CTYPE * coordB,
   const VTYPE isovalue, const NTYPE num_cube_vertices,
   const STYPE * cube_scalar, const NTYPE num_iter, CTYPE * coordC)
  {
    multilinear_interpolate_coord<double,double>
      (dimension, coordA, coordB, isovalue, num_cube_vertices,
       cube_scalar, num_iter, coordC);
  }

// **********************************************************
// Calculation of saddle points of multilinear interpolants
// **********************************************************

  /// Compute parity of integer i.  Parity is 0 or 1.
  template <class DTYPE, class ITYPE>
  ITYPE compute_parity(const DTYPE dimension, const ITYPE i)
  {
    ITYPE i2 = i;
    ITYPE parity = 0;
    for (DTYPE d = 0; d < dimension; d++) {
      parity = (parity+i2)%2;
      i2 = (i2 >> 1);
    }
    return(parity);
  }

  /// Convert parity to sign. sign = (-1)^parity.
  template <class ITYPE>
  ITYPE convert_parity_to_sign(const ITYPE parity)
  {
    if (parity == 0) { return(1); }
    else { return(-1); }
  }

  /// Compute coefficient of highest degree term in multilinear representation.
  template <class CTYPE, class DTYPE, class STYPE, class ITYPE>
  CTYPE compute_multilinear_a
  (const DTYPE dimension, const STYPE * scalar, const ITYPE * parity_sign)
  {
    // number of cube vertices
    const long numv = (1L << dimension);

    CTYPE a = 0;
    for (long i = 0; i < numv; i++) {
      ITYPE sign = convert_parity_to_sign(dimension%2); 
      sign = sign * parity_sign[i];
      a += scalar[i] * sign;
    }

    return(a);
  }

  /// Compute coefficient of second highest degree terms in multilinear representation.
  /// @param j = Index of variable missing in second highest degree term.
  template <class CTYPE, class DTYPE, class STYPE, class ITYPE>
  CTYPE compute_multilinear_b
  (const DTYPE dimension, const STYPE * scalar, const ITYPE * parity_sign,
   const DTYPE j)
  {
    // number of cube vertices
    const long numv = (1L << dimension);

    CTYPE b = 0;
    long mask = (1L << j);
    for (long i = 0; i < numv; i++) {

      if ((i & mask) == 0) {
	ITYPE sign = convert_parity_to_sign((dimension+1)%2); 
	sign = sign * parity_sign[i];
	b += scalar[i] * sign;
      }
    }

    return(b);
  }

  /// Compute coefficient of third highest degree terms in multilinear representation.
  /// @param j1 = Index of first variable missing in second highest degree term.
  /// @param j2 = Index of second variable missing in second highest degree term.
  template <class CTYPE, class DTYPE, class STYPE, class ITYPE>
  CTYPE compute_multilinear_c
  (const DTYPE dimension, const STYPE * scalar, const ITYPE * parity_sign,
   const DTYPE j1, const DTYPE j2)
  {
    // number of cube vertices
    const long numv = (1L << dimension);

    CTYPE c = 0;
    long mask = ((1L << j1) | (1L << j2));
    for (long i = 0; i < numv; i++) {

      if ((i & mask) == 0) {
	ITYPE sign = convert_parity_to_sign(dimension%2); 
	sign = sign * parity_sign[i];
	c += scalar[i] * sign;
      }
    }

    return(c);
  }

  /// Compute saddle points of 3D unit cube.
  /// Saddle points are trilinear_mipdoint[] +/- saddle_offset[].
  template <class STYPE, class ITYPE, class ATYPE, class CTYPE,
	    class NTYPE, class ETYPE>
  void compute_saddle_3D
  (const STYPE * scalar, const ITYPE * parity_sign,
   ATYPE & a, CTYPE trilinear_midpoint[3],
   CTYPE saddle_offset[3], NTYPE & num_saddle,
   const ETYPE EPSILON)
  {
    const int dimension = 3;
    ATYPE b[3];
    ATYPE c[3];
    CTYPE q[3];
    CTYPE v[3];

    // initialize
    num_saddle = 0;
    a = 0;
    for (int d = 0; d < dimension; d++) {
      trilinear_midpoint[d] = 0;
      saddle_offset[d] = 0;
    }

    a = compute_multilinear_a<STYPE>(dimension, scalar, parity_sign);

    if (-EPSILON < a && a < EPSILON) { 
      a = 0;
      return; 
    }

    for (int i = 0; i < 3; i++) {
      b[i] = compute_multilinear_b<STYPE>
	(dimension, scalar, parity_sign, i);

      int j1 = (i+1)%dimension;
      int j2 = (i+2)%dimension;
      c[i] = compute_multilinear_c<STYPE>
	(dimension, scalar, parity_sign, j1, j2);
    }

    for (int d = 0; d < 3; d++) 
      { trilinear_midpoint[d] = -b[d]/a; }

    for (int i = 0; i < 3; i++) {
      int j1 = (i+1)%dimension;
      int j2 = (i+2)%dimension;

      q[i] = b[j1]*b[j2] - a * c[i];

      // Round q[i] to zero if it is very close to 0.
      if (-EPSILON < q[i] && q[i] < EPSILON)
	{ q[i] = 0; };
    }

    CTYPE t = q[0]*q[1]*q[2];

    if (q[0] == 0 && q[1] == 0 && q[2] == 0) {
      // One saddle.
      num_saddle = 1;
      return;
    }

    if (t <= 0) { 
      // Zero saddles.
      num_saddle = 0;
      return; 
    }

    CTYPE r = 0;
    r = sqrt(t);
    for (int i = 0; i < 3; i++) { saddle_offset[i] = r / (a*q[i]); }
    num_saddle = 2;
  }

}

#endif