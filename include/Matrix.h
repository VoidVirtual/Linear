#pragma once
#include "Constants.h"
#include <vector>
#include <iostream>
#include <cmath>
#include "Vector.h"
#include <ctime>
const unsigned char RANDOM = 1;
const unsigned char DIAG = 2;
const unsigned char SYMMETRIC = 4;
const unsigned char ANTISYMMETRIC = 8;
class Matrix
{
public:
/**Constructors, destructor**/
    Matrix();
    Matrix(int N);
    Matrix(int N, unsigned char flags);
    Matrix(const Vector& v);
    Matrix(const Matrix& other);
    Matrix(const std::initializer_list<double>& diag);
    Matrix(std::initializer_list<std::initializer_list<double>> _lists);
    Matrix(const char* filename);
    virtual ~Matrix();
/**Operators**/
    Vector operator*(const Vector& v)const;
    Vector& operator[](const int index);
    Vector& operator[](const int index)const;
    Matrix operator*(const Matrix& other)const;
    Matrix operator*(const double& coef)const;
    Matrix operator-(const Matrix& other)const;
    Matrix operator+(const Matrix& other)const;
    Matrix& operator=(const Matrix& other);
    Matrix& operator*=(const double coef);
/**Main methods**/
    Matrix getOpposite()const;
    Matrix _T()const;
    int getDim()const {return dim;}
    double powerMethod(const Vector& X);
    double getMinSelfValue(const Vector& X);
    double getMaxSelfValue(const Vector& X);
    double getMinSelfValue();
    double getMaxSelfValue();
    double Det();
    Vector Jakobi();
/**Friend Operators**/
    friend Matrix operator*(const double coef, const Matrix &m);
    friend std::ostream& operator<<(std::ostream &out, const Matrix &m);
    friend std::istream& operator>>(std::istream &in, Matrix &m);
protected:
    Vector* m_elements;
    int i_max;
    int j_max;
    int dim;
/**Elemenary**/
    void Init(int N);
    void Clean();
    void Copy(const Matrix& other);
    void SwapStrings(int index1, int index2);
    void AddString(double coef, int s1, int s2);
    void AddRow(double coef, int r1, int r2);
    double maxNonDiag();
};
