#include "ExtendedMatrix.h"
#include "Vector.h"
#include <cassert>
#include <cmath>
double solve_linear(double k, double b)
{
    assert(k!=0);
    return -b/k;
}
ExtendedMatrix::ExtendedMatrix():ExtendedMatrix(3){}
ExtendedMatrix::ExtendedMatrix(int N):Matrix(N),m_column(N),m_solution(N)
{

}
Matrix ExtendedMatrix::getMatrix()
{
    Matrix res(dim);
    for(int i = 0; i < dim; i++)
        res[i] = (*this)[i];
    return res;
}
void ExtendedMatrix::SwapStrings(int s1, int s2)
{
    Matrix::SwapStrings(s1,s2);
    std::swap(m_column[s1], m_column[s2]);
}
void ExtendedMatrix::AddString(double coef, int s1, int s2)
{
    m_elements[s2]+=coef*m_elements[s1];
    m_column[s2]+=coef*m_column[s1];
}
void ExtendedMatrix::ScalarString(double coef, int index)
{
    m_elements[index]*=coef;
    m_column[index]*=coef;
}
void ExtendedMatrix::Gauss()
{
    m_solvable = true;
    for (int key = 0; key < dim-1; key++)
    {
        bool nonZeroColumn = true;
        if(m_elements[key][key]==0)
        {
            nonZeroColumn = false;
            for(int i = key+1; i < dim; i++)
            {
                if(m_elements[i][key]!=0)
                {
                    nonZeroColumn = true;
                    SwapStrings(key, i);
                }
            }
        }
        if(nonZeroColumn)
        {
            for(int i = key +1; i < dim; i++)
            {
                double coef = -1*m_elements[i][key]/m_elements[key][key];
                AddString(coef,key,i);
            }
        }
        else
        {
            m_solvable = false;
            break;
        }
    }
}
void ExtendedMatrix::SolveGauss()
{
    Gauss();
    if(!m_solvable)
        return;
    SolveRightTriangle();
}
void ExtendedMatrix::SolveRightTriangle()
{
    Vector x(dim);
    x[dim-1] = solve_linear(m_elements[dim-1][dim-1],-m_elements[dim-1][dim]);
    for(int key = dim-2; key >=0; key--)
    {
        double b = -m_elements[key][dim];
        for(int j = dim -1; j > key; j--)
        {
            b+=(m_elements[key][j]*x[j]);
        }
        x[key] = solve_linear(m_elements[key][key],b);
    }
    m_solvable = true;
    m_solution = x;
}
void ExtendedMatrix::SolveLeftTriangle()
{
    Vector x(dim);
    x[0] = solve_linear(m_elements[0][0], -m_column[0]);
    for(int key = 1; key < dim; key++)
    {
        double b = -m_column[key];
        for(int j = 0; j <key; j++)
            b+=m_elements[key][j]*x[j];
        x[key] = solve_linear(m_elements[key][key],b);
    }
    m_solution = x;
    m_solvable = true;
}
void ExtendedMatrix::TransposedGauss()
{
    m_solvable = true;
    for (int key = dim-1; key > 0; key--)
    {
        bool nonZeroColumn = true;
        if(m_elements[key][key]==0)
        {
            nonZeroColumn = false;
            for(int i = key-1; i >=0; i--)
            {
                if(m_elements[i][key]!=0)
                {
                    nonZeroColumn = true;
                    SwapStrings(key, i);
                }
            }
        }
        if(nonZeroColumn)
        {
            for(int i = key -1; i >= 0; i--)
            {
                double coef = -1*m_elements[i][key]/m_elements[key][key];
               AddString(coef,key,i);
            }
        }
        else
        {
            m_solvable = false;
            break;
        }
    }
}
void ExtendedMatrix::Seidel(Vector x)
{
    double norm;
    int iterations = 0;
    do
    {
        Vector previous = x;
        for(int i = 0; i < dim; i++)
        {
            double a = 0;
            for(int j = 0; j < i; j++)
                a+=((*this)[i][j]*x[j]);
            for(int j = i+1; j < dim; j++)
                a+=((*this)[i][j]*previous[j]);
            x[i] = (m_column[i] - a)/(*this)[i][i];
        }
        norm = (x-previous).getNorm();
        iterations++;
    } while(norm > eps);
    std::cout<<"Iterations: "<<iterations<<"\n";
    m_solution = x;
    m_solvable = true;
}
void ExtendedMatrix::Seidel()
{
    Vector v(dim);
    v[0] = 1;
    Seidel(v);
}
void ExtendedMatrix::Square()
{
    Matrix G(dim);
    Matrix A = getMatrix();
    for(int i = 0; i < dim; i++)
    {
        double sum1 = 0;
        for(int k = 0; k < i; k++)
            sum1+=G[k][i]*G[k][i];
        G[i][i] = sqrt((*this)[i][i] - sum1);
        for(int j = i+1; j < dim; j++)
        {
            double sum2 = 0;
            for(int k = 0; k < i; k++)
                sum2+=G[k][i]*G[k][j];
            G[i][j] = ((*this)[i][j] - sum2)/G[i][i];
        }
    }
    Matrix Gt = G._T(); // Left Triangle
    ExtendedMatrix Gt_ext(G._T(),m_column);
    Gt_ext.SolveLeftTriangle();

    Vector y = Gt_ext.getSolution();
    ExtendedMatrix G_ext(G,y); //Right Triangle
    G_ext.SolveRightTriangle();

    m_solution = G_ext.getSolution();
    m_solvable = true;
}
void ExtendedMatrix::Jakobi(Vector X)
{
    double norm;
    int it = 0;
    Vector tempX(dim);
    do {
        tempX = m_column;
		for (int i = 0; i < dim; i++)
        {
			for (int j= 0; j < dim; j++)
			{
				if (i != j)
					tempX[i] -= (*this)[i][j] * X[j];
			}
			tempX[i] /= (*this)[i][i];
		}
        norm = fabs(X[0] - tempX[0]);
		for (int i = 0; i < dim; i++)
		 {
			if (fabs(X[i] - tempX[i]) > norm)
				norm = fabs(X[i] - tempX[i]);
		 }
		 X = tempX;
		 it++;
	} while (norm > eps);
	std::cout<<"Iterations: "<<it<<std::endl;
    m_solution = X;
    m_solvable = true;
}
void ExtendedMatrix::Diag()
{
    Gauss();
    TransposedGauss();
    if(!m_solvable)
        throw "Can't be diagonalized\n";
    for(int i = 0; i < dim; i++)
    {
        double coef = 1.0/m_elements[i][i];
        ScalarString(coef,i);
    }

}
void ExtendedMatrix::PrintSolution()
{
    if(!m_solvable)
        std::cout <<"No solution"<<std::endl;
    else
    {
        for(int i = 0; i < dim; i++)
            std::cout <<m_solution[i]<<" ";
        std::cout<<std::endl;
    }
}
ExtendedMatrix::ExtendedMatrix(ExtendedMatrix& other)
{
   Copy(other);
   m_column = other.m_column;
}
ExtendedMatrix& ExtendedMatrix::operator=(ExtendedMatrix& other)
{
    if(this==&other)
        return *this;
    Matrix::Clean();
    Init(other.getDim());
    Copy(other);
    return *this;
}
ExtendedMatrix::ExtendedMatrix(const Matrix& m, const Vector& v):Matrix(m)
{
    m_column = v;
}
Vector& ExtendedMatrix::getSolution()
{
    if(!m_solvable)
        throw "No solution";
    return m_solution;
}
Vector& ExtendedMatrix::getColumn()
{
    return m_column;
}
ExtendedMatrix::~ExtendedMatrix()
{

}
std::ostream& operator<<(std::ostream &out, const ExtendedMatrix &m)
{
    for(int i = 0; i < m.dim; i++)
    {
        out<<m[i];
        out<<m.m_column[i];
        out<<std::endl;
    }
    return out;
}
std::istream& operator >> (std::istream &in, ExtendedMatrix &m)
{
    for(int i = 0; i < m.dim; i++)
    {
        in>>m[i];
        in>>m.m_column[i];
    }
    return in;
}
