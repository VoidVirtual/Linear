#pragma once
#include "Matrix.h"
#include "Constants.h"
class ExtendedMatrix:public Matrix
{
public:
/**Constructors, destructor:**/
    ExtendedMatrix();
    ExtendedMatrix(int N);
    ExtendedMatrix(int N, const char flags):Matrix(N, flags), m_column(N){}
    ExtendedMatrix(ExtendedMatrix& other);
    ExtendedMatrix& operator=(ExtendedMatrix& other);
    ExtendedMatrix(const Matrix& m, const Vector& v);
    virtual ~ExtendedMatrix();
/**Elementary**/
    friend std::ostream& operator<<(std::ostream &out, const ExtendedMatrix &m);
    friend std::istream& operator>>(std::istream &in, ExtendedMatrix &m);
    void SwapStrings(int index1, int index2);
    void ScalarString(double coef, int index);
    void AddString(double coef, int s1, int s2);
    void PrintSolution();
    Vector& getSolution();
    Vector& getColumn();
    Matrix getMatrix();
/**Solving linear systems:**/
    void Jakobi(Vector x);
    void Gauss();
    void SolveGauss();
    void SolveRightTriangle();
    void SolveLeftTriangle();
    void TransposedGauss();
    void Diag();
    void Seidel(Vector X);
    void Seidel();
    void Square();
protected:
    Vector m_column;
    Vector m_solution;
    double det;
    bool m_solvable = false;
};

