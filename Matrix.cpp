#include <Matrix.h>

void Print(const float* A, int row, int col, String label)
{
    Serial.println();
    Serial.println(label);
    for(int i=0; i<row; i++)
    {
        for(int j=0; j<col; j++)
        {
            Serial.print(A[i*col+j]);
            Serial.print("\t");
        }
        Serial.print("\n");
    }
}

void Print(const float* A, int dim, String label)
{
	Print(A, dim, dim, label);
}

void Copy(const float* A, int row, int col, float* B)
{
    for(int i=0; i < row; i++)
    {
        for(int j=0; j<col; j++)
        {
            B[i*col+j]=A[i*col+j];
        }
    }
}

void Copy(const float* A, int dim, float* B)
{
    Copy(A, dim, dim, B);
}

void Add(const float* A, const float* B, const int& row, const int& col, float* C)
{
    for(int i=0; i<row; i++)
    {
        for(int j=0; j<col; j++)
        {
            C[i*col+j]=A[i*col+j]+B[i*col+j];
        }
    }
}

void Add(const float* A, const float* B, const int& dim, float* C)
{
    Add(A, B, dim, dim, C);
}

void Sub(const float* A, const float* B, const int& row, const int& col, float* C)
{
    for(int i=0; i<row; i++)
    {
        for(int j=0; j<col; j++)
        {
            C[i*col+j]=A[i*col+j]-B[i*col+j];
        }
    }
}

void Sub(const float* A, const float* B, const int& dim, float* C)
{
    Sub(A, B, dim, dim, C);
}

void Multiply(const float* A, const int& A_row, const int& A_col,
              const float* B, const int& B_row, const int& B_col,
              float* C)
{
    if(A_col != B_row) return;

    for(int i=0; i<A_row; i++)
    {
        for(int j=0; j<B_col; j++)
        {
            C[i*B_col+j]=0;
            for(int k=0; k<A_col; k++)
            {
                C[i*B_col+j]+=A[i*A_col+k]*B[k*B_col+j];
            }
        }
    }
}

void Multiply(const float* A, const float* B, const int& dim, float* C)
{
    Multiply(A, dim, dim, B, dim, dim, C);
}

void Transpose(const float* A, const int& A_row, const int& A_col, float* B)
{
    for(int i=0; i<A_row; i++)
    {
        for(int j=0; j<A_col; j++)
        {
            B[j*A_row+i]=A[i*A_col+j];
        }
    }
}

void Transpose(const float* A, const int& dim, float* B)
{
    Transpose(A, dim, dim, B);
}

void Scale(float* A, const int& A_row, const int& A_col, float scale)
{
    for(int i=0; i<A_row; i++)
    {
        for(int j=0; j<A_col; j++)
        {
            A[i*A_col+j]*=scale;
        }
    }
}

void Scale(float* A, const int& dim, float scale)
{
    Scale(A, dim, dim, scale);
}

//Matrix Inversion Routine
// * This function inverts a matrix based on the Gauss Jordan method.
// * Specifically, it uses partial pivoting to improve numeric stability.
// * The algorithm is drawn from those presented in
//	 NUMERICAL RECIPES: The Art of Scientific Computing.
// * The function returns 1 on success, 0 on failure.
// * NOTE: The argument is ALSO the result matrix, meaning the input matrix is REPLACED
int Invert(float* A, int dim)
{
	int n=dim;
	// A = input matrix AND result matrix
	// n = number of rows = number of columns in A (n x n)
	int pivrow;		// keeps track of current pivot row
	int k, i, j;		// k: overall index along diagonal; i: row index; j: col index
	int pivrows[n]; // keeps track of rows swaps to undo at end
	float tmp;		// used for finding max value and making column swaps

	for (k = 0; k < n; k++)
	{
		// find pivot row, the row with biggest entry in current column
		tmp = 0;
		for (i = k; i < n; i++)
		{
			if (abs(A[i * n + k]) >= tmp)	// 'Avoid using other functions inside abs()?'
			{
				tmp = abs(A[i * n + k]);
				pivrow = i;
			}
		}

		// check for singular matrix
		if (A[pivrow * n + k] == 0.0f)
		{
			Serial.println("Inversion failed due to singular matrix");
			return 0;
		}

		// Execute pivot (row swap) if needed
		if (pivrow != k)
		{
			// swap row k with pivrow
			for (j = 0; j < n; j++)
			{
				tmp = A[k * n + j];
				A[k * n + j] = A[pivrow * n + j];
				A[pivrow * n + j] = tmp;
			}
		}
		pivrows[k] = pivrow;	// record row swap (even if no swap happened)

		tmp = 1.0f / A[k * n + k];	// invert pivot element
		A[k * n + k] = 1.0f;		// This element of input matrix becomes result matrix

		// Perform row reduction (divide every element by pivot)
		for (j = 0; j < n; j++)
		{
			A[k * n + j] = A[k * n + j] * tmp;
		}

		// Now eliminate all other entries in this column
		for (i = 0; i < n; i++)
		{
			if (i != k)
			{
				tmp = A[i * n + k];
				A[i * n + k] = 0.0f; // The other place where in matrix becomes result mat
				for (j = 0; j < n; j++)
				{
					A[i * n + j] = A[i * n + j] - A[k * n + j] * tmp;
				}
			}
		}
	}

	// Done, now need to undo pivot row swaps by doing column swaps in reverse order
	for (k = n - 1; k >= 0; k--)
	{
		if (pivrows[k] != k)
		{
			for (i = 0; i < n; i++)
			{
				tmp = A[i * n + k];
				A[i * n + k] = A[i * n + pivrows[k]];
				A[i * n + pivrows[k]] = tmp;
			}
		}
	}
	return 1;
}

Matrix::Matrix(const int& _row, const int& _col)
{
	float* p;
	p = (float*)malloc(_row*_col*sizeof(float));
	mat = p;
	row = _row;
	col = _col;

	for(int i=0; i<row; i++)
	{
		for(int j=0; j<col; j++)
		{
			((float*)mat)[i*col+j]=0.f;
		}
	}
}

Matrix::Matrix(const int& _dim)
{
	float* p;
	p = (float*)malloc(_dim*_dim*sizeof(float));
	mat = p;
	row = _dim;
	col = _dim;

	for(int i=0; i<row; i++)
	{
		for(int j=0; j<col; j++)
		{
			((float*)mat)[i*col+j]=0.f;
		}
	}
}

void Matrix::zero()
{
	for(int i=0; i<row; i++)
	{
		for(int j=0; j<col; j++)
		{
			((float*)mat)[i*col+j]=0.f;
		}
	}
}

void Matrix::eye()
{
	for(int i=0; i<row; i++)
	{
		for(int j=0; j<col; j++)
		{
			if(i==j)
				((float*)mat)[i*col+j]=1.f;
			else
				((float*)mat)[i*col+j]=0.f;
		}
	}
}

void Matrix::one()
{
	for(int i=0; i<row; i++)
	{
		for(int j=0; j<col; j++)
		{
			((float*)mat)[i*col+j]=1.f;
		}
	}
}

Matrix Matrix::add(const Matrix& _mat)
{
	Matrix res(row,col);
	Add((float*)mat, (float*)(_mat.mat), row, col, (float*)(res.mat));

	return res;
}

Matrix Matrix::sub(const Matrix& _mat)
{
	Matrix res(row,col);
	Sub((float*)mat, (float*)(_mat.mat), row, col, (float*)(res.mat));

	return res;
}

Matrix Matrix::multiply(const Matrix& _mat)
{
	Matrix res(row, _mat.col);
	
	if(col!=_mat.row)
		return res;

	Multiply((float*)mat, row, col, (float*)(_mat.mat), _mat.row, _mat.col, (float*)(res.mat));

	return res;
}

Matrix Matrix::invert()
{
	Matrix res(row, col);

	if(row!=col)
		return res;
	
	for(int i=0; i<row; i++)
	{
		for(int j=0; j<col; j++)
		{
			res.mat[i*col+j]=mat[i*col+j];
		}
	}

	Invert((float*)res.mat, row);

	return res;
}

Matrix Matrix::scale(const int& k)
{
	Matrix res(row,col);
	
	Copy((float*)mat, row, col, (float*)(res.mat));
	Scale((float*)(res.mat), row, col, k);

	return res;
}

Matrix Matrix::transpose()
{
	Matrix res(row,col);
	Transpose((float*)mat, row, col, (float*)(res.mat));

	return res;
}

Matrix Matrix::copy()
{
	Matrix res(row, col);
	Copy((float*)mat, row, col, (float*)(res.mat));

	return res;
}

void Matrix::print(String label)
{
	Print((float*)mat, row, col, label);
}

float Matrix::at(const unsigned int& _row, const unsigned int& _col)
{
	return ((float*)mat)[_row*col+_col];
}

float Matrix::at(const unsigned int& _row, const unsigned int& _col, const float value)
{
	((float*)mat)[_row*col+_col] = value;
	return value;
}
