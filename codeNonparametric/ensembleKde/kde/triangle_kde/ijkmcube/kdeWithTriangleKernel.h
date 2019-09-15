using namespace std;

class z_density_triangular{

	private :
		double numTriangularInKde1, numTriangularInKde2;
		
	public :

        double alpha_density_triangular(double m1, double d1, double  m2, double  d2, double isovalue, double*  exp, double*  var, double*  crossprob);

	double kde_z_pdf_expected(float* mu1, double h1, float* mu2, double h2, double c);

	double kde_z_pdf_variance(float* mu1, double h1, float* mu2, double h2, double c, double expected_crossing);
	
	// set number of uniform distributions in kde_1 
	void setNumTriangularInKde1(int a);

	// get number of uniform distributions in kde_1 
	int getNumTriangularInKde1();

	// set number of uniform distributions in kde_2 
	void setNumTriangularInKde2(int a);

	// get number of uniform distributions in kde_2
	int getNumTriangularInKde2();
		
};

