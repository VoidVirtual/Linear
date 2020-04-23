#include "Matrix.h"
#include "ExtendedMatrix.h"
#include <fstream>
#include <cmath>
#include <cassert>
void Matrix::Init(int N)
{
    dim = N;
    m_elements = new Vector[N];
    for(int i = 0; i < dim; i++)
    {
        m_elements[i] = Vector(N);
    }
}
Matrix::Matrix(int N)
{
    Init(N);
}
Matrix::Matrix(int N, unsigned char flags):Matrix(N)
{
    if(!(flags&RANDOM))
            return;
    for(int i = 0; i < dim; i++)
        for(int j = 0; j < dim; j++)
        {
            if(i==j || !(flags&DIAG))
                m_elements[i][j] = (unsigned char)rand();
        }
    if(flags&SYMMETRIC)
    {
        for(int i = 0; i < dim; i++)
            for(int j = 0; j < i; j++)
                m_elements[i][j] = m_elements[j][i];
    }
    if(flags&ANTISYMMETRIC)
    {
        for(int i = 0; i < dim; i++)
            for(int j = 0; j < i; j++)
                m_elements[i][j] = -m_elements[j][i];
    }
}
Matrix::Matrix(const std::initializer_list<double>& diag):Matrix(diag.size())
 {
     int i = 0;
     for(auto x: diag)
     {
         (*this)[i][i] = x;
         i++;
     }
 }
Matrix::Matrix(const Vector& v):Matrix(v.size())
{
    for(unsigned int i = 0; i < dim; i++)
        (*this)[i][i] = v[i];
}
Matrix::Matrix(const char* filename)
{
    std::ifstream in;
    in.exceptions(std::ifstream::failbit|std::ifstream::badbit);
    try
    {
        in.open(filename);
    }
    catch(std::exception& e)
    {
        std::cerr<<"Wrong filename\n"<<e.what();
        exit(1);
    }
    in>>dim;
    Init(dim);
    in>>(*this);
    in.close();
}
Matrix::Matrix():Matrix(3)
{

}
void Matrix::Clean()
{
     delete[] m_elements;
}
Matrix::~Matrix()
{
    Clean();
}
void Matrix::SwapStrings(int index1, int index2)
{
    std::swap((*this)[index1], (*this)[index2]);
}
void Matrix::AddString(double coef, int s1, int s2)
{
    m_elements[s2]+=coef*m_elements[s1];
}
void Matrix::AddRow(double coef, int r1, int r2)
{
    for(unsigned int i = 0; i < dim; i++)
    {
        m_elements[i][r2] += coef*m_elements[i][r1];
    }
}
Matrix Matrix::getOpposite()const
{
    Matrix res(dim);
    for(int j = 0; j < dim; j++)
    {
        Vector v(dim);
        v[j] = 1;
        ExtendedMatrix ext(*this, v);
        try
        {
            ext.Diag();
        }
        catch(const char* exception)
        {
            throw "The opposite matrix doesn't exist.\n";
        }
        Vector tempColumn = ext.getColumn();
        for(int i = 0; i < dim; i++)
        {
            res[i][j] = tempColumn[i];
        }
    }
    return res;
}
Matrix Matrix::operator*(const Matrix& other)const
{
    if(dim!=other.dim)
        throw "Different dimentions\n";
    Matrix res(dim);
    for(int i = 0; i < dim; i++)
        for(int j = 0; j < dim; j++)
            for(int k = 0; k < dim; k++)
                res[i][j] += (*this)[i][k]*other[k][j];
    return res;
}
Vector Matrix::operator*(const Vector& v)const
{
    if(dim!=v.size())
        throw "Different dimentions\n";
    Vector res(dim);
    for(int i = 0; i < dim; i++)
        for(int k = 0; k < dim; k++)
            res[i] += (*this)[i][k]*v[k];
    return res;
}
double Matrix::powerMethod(const Vector& X)
{
    int iterations = 0;
    Vector Xn(dim);
	Vector Xn1(dim);
	Xn = (*this)*X;
	double _eps = 1;
	while (_eps > 0.00001)
	{
		Xn1 = (*this)*Xn;
		for (int i = 0; i < dim; i++)
        {
			Xn1[i] = Xn1[i] / Xn.getNorm();
		}
		_eps = fabs(Xn.getNorm() - Xn1.getNorm());
        Xn = Xn1;
        iterations++;
	}
	std::cout<<"Iterations: "<<iterations<<std::endl;
	return Xn1.getNorm();
}
double Matrix::getMaxSelfValue(const Vector& X)
{
   return powerMethod(X);
}
double Matrix::getMinSelfValue(const Vector& X)
{
    double maxSelfValue = getMaxSelfValue(X);
    Matrix B(dim, maxSelfValue);
    B = B - (*this);
    double res = maxSelfValue - B.getMaxSelfValue(X);
    return res;
}
double Matrix::getMaxSelfValue()
{
    Vector v(dim);
    v[0] = 1;
    return getMaxSelfValue(v);
}
double Matrix::getMinSelfValue()
{
    Vector v(dim);
    v[0] = 1;
    return getMinSelfValue(v);
}
void Matrix::Copy(const Matrix& other)
{
    for(int i = 0; i < dim; i++)
        for(int j = 0; j < dim; j++)
            (*this)[i][j] = other[i][j];
}

Matrix::Matrix(const Matrix& other)
{
    //Clear();
    Init(other.getDim());
    Copy(other);
}
Matrix& Matrix::operator=(const Matrix& other)
{
    if(this == &other)
        return *this;
    Clean();
    Init(other.getDim());
    Copy(other);
    return *this;
}
Matrix Matrix::operator*(const double& scalar)const
{
    Matrix res(dim);
    for(int i = 0; i < dim; i++)
        for(int j = 0; j < dim; j++)
            res[i][j] = m_elements[i][j]*scalar;
    return res;
}
Matrix Matrix::operator-(const Matrix& other)const
{
    if(dim!=other.dim)
        throw "Different dimentions\n";
    Matrix res(dim);
    for(int i = 0; i < dim; i++)
        for(int j = 0; j < dim; j++)
            res[i][j] = (*this)[i][j] - other[i][j];
    return res;
}
Matrix Matrix::operator+(const Matrix& other)const
{
    if(dim!=other.dim)
        throw "Different dimentions\n";
    Matrix res(dim);
    for(int i = 0; i < dim; i++)
        for(int j = 0; j < dim; j++)
            res[i][j] = (*this)[i][j] + other[i][j];
    return res;
}

Matrix::Matrix(std::initializer_list<std::initializer_list<double>> _lists):Matrix(_lists.size())
{
    int i = 0;
    for(auto _list: _lists)
    {
        m_elements[i] = _list;
        i++;
    }
}
double Matrix::Det()
{
    Vector v(dim);
    v[0] = 1;
    ExtendedMatrix Ext((*this),v);
    Ext.Gauss();
    double res = 1;
    for(int i = 0; i < dim; i++)
        res*=Ext[i][i];
    return res;
}
Matrix Matrix::_T()const
{
    Matrix res(dim);
    for(int i = 0; i < dim; i++)
    {
        for(int j = 0; j < dim; j++)
            res[i][j] = (*this)[j][i];
    }
    return res;
}
Vector& Matrix::operator[](const int index)
{
    return m_elements[index];
}
Vector& Matrix::operator[](const int index) const
{
    return m_elements[index];
}
double Matrix::maxNonDiag()
{
    double curMax = 0;
    for(int i = 0; i < dim; i++)
        for(int j = 0; j < dim; j++)
            if(i!=j && fabs((*this)[i][j]) > fabs(curMax))
            {
                curMax = (*this)[i][j];
                i_max = i;
                j_max = j;
            }
    return(*this)[i_max][j_max];
}
Vector Matrix::Jakobi()
{
    Vector Ai(dim);
    Vector Aj(dim);
    double m;
    do
    {
        m = maxNonDiag();
        double cos;
        double sin;
        if((*this)[i_max][i_max]!=(*this)[j_max][j_max])
        {
            double d = (*this)[i_max][i_max]- (*this)[j_max][j_max];
            double tau = 1.0/(2*m/d);
            double t = -tau + (tau/fabs(tau))*sqrt(1.0 + tau*tau);
            cos = 1.0/sqrt(1.0 + t*t);
            sin = t*cos;
        }
        else
        {
            cos = 1.0/sqrt(2.0);
            sin = 1.0/sqrt(2.0);
        }
        for(int k = 0; k < dim; k++)
        {
            Ai[k] = cos*(*this)[i_max][k] + sin*(*this)[j_max][k];
            Aj[k] = cos*(*this)[j_max][k] - sin*(*this)[i_max][k];
        }
        for(int k = 0; k < dim; k++)
        {
            (*this)[i_max][k] = Ai[k];
            (*this)[j_max][k] = Aj[k];
        }
        for(int k = 0; k < dim; k++)
        {
            Ai[k] = cos*(*this)[k][i_max] + sin*(*this)[k][j_max];
            Aj[k] = cos*(*this)[k][j_max] - sin*(*this)[k][i_max];
        }
        for(int k = 0; k < dim; k++)
        {
            (*this)[k][i_max] = Ai[k];
            (*this)[k][j_max] = Aj[k];
        }
    }
    while(fabs(m)>eps);
    Vector res(dim);
    for(unsigned int i = 0; i < dim; i++)
        res[i] = (*this)[i][i];
    return res;

}
Matrix operator*(const double coef, const Matrix &m)
{
    Matrix res(m);
    for(int i = 0; i < res.dim; i++)
        for(int j = 0; j < res.dim; j++)
            res[i][j]*=coef;
    return res;
}
Matrix& Matrix::operator*=(const double coef)
{
    for(int i = 0; i < dim; i++)
        for(int j = 0; j < dim; j++)
            (*this)[i][j] *= coef;
    return *this;
}
std::ostream& operator<<(std::ostream &out, const Matrix &m)
{
    for(int i = 0; i < m.dim; i++)
    {
        for(auto x: m[i])
            out<<x<<" ";
        out<<std::endl;
    }
    return out;
}
std::istream& operator >> (std::istream &in, Matrix &m)
{
    for(int i = 0; i < m.dim; i++)
        for(auto &x: m[i])
            in>>x;
    return in;
}


