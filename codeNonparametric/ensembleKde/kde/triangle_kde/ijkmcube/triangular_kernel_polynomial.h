#include "stdio.h"
#include <iostream>
#include "stdlib.h"
#include "math.h"

using namespace std;

class triangular_kernel_polynomial{

	private: double mu1, delta1, mu2, delta2, c;

	public:
		triangular_kernel_polynomial(double m1,double d1,double m2,double d2,double isoval)
		{
			mu1 = m1;
			delta1 = d1;
			mu2 = m2;
			delta2 = d2;
			c = isoval;
		}

		double PAED_integrate_piece_value(double z, int piecenum);
		double PAED_expected_piece_value(double z, int piecenum);
		double PAED_second_moment_piece_value(double z, int piecenum);

		double AQBE_integrate_piece_value(double z, int piecenum);
		double AQBE_expected_piece_value(double z, int piecenum);
		double AQBE_second_moment_piece_value(double z, int piecenum);

		double EBRC_integrate_piece_value(double z, int piecenum);
		double EBRC_expected_piece_value(double z, int piecenum);
		double EBRC_second_moment_piece_value(double z, int piecenum);

		double DECS_integrate_piece_value(double z, int piecenum);
		double DECS_expected_piece_value(double z, int piecenum);
		double DECS_second_moment_piece_value(double z, int piecenum);		


		double approx_PAED_integrate_piece_value(double z, int piecenum);
		double approx_PAED_expected_piece_value(double z, int piecenum);
		double approx_PAED_second_moment_piece_value(double z, int piecenum);

		double approx_AQBE_integrate_piece_value(double z, int piecenum);
		double approx_AQBE_expected_piece_value(double z, int piecenum);
		double approx_AQBE_second_moment_piece_value(double z, int piecenum);

		double approx_EBRC_integrate_piece_value(double z, int piecenum);
		double approx_EBRC_expected_piece_value(double z, int piecenum);
		double approx_EBRC_second_moment_piece_value(double z, int piecenum);

		double approx_DECS_integrate_piece_value(double z, int piecenum);
		double approx_DECS_expected_piece_value(double z, int piecenum);
		double approx_DECS_second_moment_piece_value(double z, int piecenum);		
};

