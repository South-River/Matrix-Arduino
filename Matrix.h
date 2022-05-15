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
void Scale(float* A, const int& A_row, const int& A_col, float scale);
void Scale(float* A, const int& dim, float scale);

//A=A^(-1)
int Invert(float* A, int dim);

class Matrix
{
public:    
    Matrix(const int& _dim);
    Matrix(const int& _row, const int& _col);

    void zero();
    void eye();
    void one();

    Matrix add(const Matrix& _mat);
    Matrix sub(const Matrix& _mat);
    Matrix multiply(const Matrix& _mat);
    Matrix scale(const int& k);
    Matrix transpose();
    Matrix copy();
    Matrix invert();

    void print(String label="");
    
    float at(const unsigned int& _row, const unsigned int& _col);
    float at(const unsigned int& _row, const unsigned int& _col, const float value);

public:
    float* mat;
    int row, col;
};

#endif