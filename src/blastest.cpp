#include "blas_wrap.h"
#include <assert.h>
#include <iostream>

using namespace std;
int main( int argc, char **argv )
{
	// test case 1
	{
		Matrix a(2, 3);
		Matrix b(2, 4);
		Matrix c(3, 4);

		a.set(0,0, -3.0);
		a.set(0,1, 1.0);
		a.set(0,2, 1.0);
		a.set(1,0, 1.0);
		a.set(1,1, -1.0);
		a.set(1,2, 4.0);

		a.transpose();

		b.set(0,0, 2.0);
		b.set(0,1, 1.0);
		b.set(0,2, 9.0);
		b.set(0,3, 4.1);
		b.set(1,0, -3.2);
		b.set(1,1, 1.4);
		b.set(1,2, 0.0);
		b.set(1,3, 0.0);

		c.zero_fill();

		// multiply matrices
		matrix_multiply(a, b, c, 0.1, 2.0);

		// check result
		assert(almost_equal( c.get(0,0), -0.92 ));
		assert(almost_equal( c.get(0,1), -0.16 ));
		assert(almost_equal( c.get(0,2), -2.7  ));
		assert(almost_equal( c.get(0,3), -1.23 ));
		assert(almost_equal( c.get(1,0), 0.52  ));
		assert(almost_equal( c.get(1,1), -0.04 ));
		assert(almost_equal( c.get(1,2), 0.9   ));
		assert(almost_equal( c.get(1,3), 0.41  ));
		assert(almost_equal( c.get(2,0), -1.08 ));
		assert(almost_equal( c.get(2,1), 0.66  ));
		assert(almost_equal( c.get(2,2), 0.9   ));
		assert(almost_equal( c.get(2,3), 0.41  ));
	}

	// test case 2
	{
		Matrix a(3,1);
		Matrix b(1,3);
		Matrix c(3,3);
		a.one_fill();
		b.one_fill();
		c.zero_fill();
		

		// multiply matrices
		matrix_multiply(a, b, c, 0.3, 2.0);

		// check result
		for (int i = 0; i < c.first_dim(); ++i)
		{
			for (int j = 0; j < c.second_dim(); ++j)
			{
				assert(almost_equal( c.get(i,j), 0.3 ));
			}
		}
	}

	// test case 3: external memory management
	{
		int m = 10;
		int n = 30;
		double *data = (double*) malloc( m*n * sizeof(double) );
		Matrix a(m,n,data);
		a.one_fill();
		for (int i = 0; i < m*n; ++i)
		{
			assert( data[i] == 1.0 );
		}
		free(data);
	}

	// test case 4: matrix copy
	{
		int m = 11;
		int n = 14;
		Matrix a(m,n);
		Matrix b(m,n);
		a.one_fill();
		matrix_copy(a, b);
		for (int i = 0; i < m; ++i)
		{
			for (int j = 0; j < n; ++j)
			{
				assert( b.get(i, j) == 1.0 );
			}
		}
	}

	std::cout << "All tests passed." << std::endl;

	return 0;
}