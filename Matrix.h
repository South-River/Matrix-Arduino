#ifndef MATRIX_H
#define MATRIX_H

#include <Arduino.h>

void Print(const float* A, int row, int col, String label="");
void Print(const float* A, int dim, String label="");

// from A to B
void Copy(const float* A, int row, int col, float* B);
void Copy(const float* A, int dim, float* B);

// C=A+B
void Add(const float* A, const float* B, const int& row, const int& col, float* C);
void Add(const float* A, const float* B, const int& dim, float* C);

// C=A-B
void Sub(const float* A, const float* B, const int& row, const int& col, float* C);
void Sub(const float* A, const float* B, const int& dim, float* C);

// C=A*B
void Multiply(const float* A, const int& A_row, const int& A_col,
              const float* B, const int& B_row, const int& B_col,
              float* C);
void Multiply(const float* A, const float* B, const int& dim, float* C);

//C=A^T
void Transpose(const float* A, const int& A_row, const int& A_col, float* B);
void Transpose(const float* A, const int& dim, float* B);

//A=kA
void Scale(float* A, const int& row, const int& col, float scale);
void Scale(float* A, const int& dim, float scale);

//A=A^(-1)
int Invert(float* A, int dim);

namespace Matrix
{
    class Matrix
    {
    public: 
        Matrix(const int& _row, const int& _col);   
        Matrix(const int& _dim);

        void zero();
        void eye();
        void one();
    
    public:
        Matrix copy();

    public:
        Matrix add(const Matrix& _mat);
        Matrix sub(const Matrix& _mat);
        Matrix multiply(const Matrix& _mat);
        Matrix invert();

    public:
        Matrix scale(const int& k);
        Matrix transpose();        

    public:
        Matrix Block(const int& _row, const int& _col, const int& height, const int& width);
        void Block(const int& _row, const int& _col, const Matrix& _mat);

        void swapRow(const int& _row1, const int& _row2);
        void swapCol(const int& _col1, const int& _col2);

        float at(const unsigned int& _row, const unsigned int& _col);
        float at(const unsigned int& _row, const unsigned int& _col, const float value);

    public:
        void resize(const int& _row, const int& _col);
        void resize(const int& _dim);
        
    public:
        void print(String label="");

    public:
        Matrix operator+(const Matrix& _mat);
        Matrix operator-(const Matrix& _mat);
        Matrix operator*(const Matrix& _mat);

        bool operator==(const Matrix& _mat);
        bool operator!=(const Matrix& _mat);

    public:
        float* mat;
        int row, col;
    };

    Matrix zero(const int& _row, const int& _col);
    Matrix zero(const int& _dim);
    Matrix eye(const int& _row, const int& _col);
    Matrix eye(const int& _dim);
    Matrix one(const int& _row, const int& _col);
    Matrix one(const int& _dim);
}

#endif