#include "Vector.h"
double Vector::getNorm()const
{
    double res = 0;
    for(auto x: *this)
        res += x*x;
    return sqrt(res);
}
void Vector::Normalize()
{
    double norm = getNorm();
    for(auto &x: *this)
        x/=norm;
}
std::ostream& operator<<(std::ostream& out, const Vector& v)
{
    for(auto x:v)
       out<<x<<" ";
    return out;
}
std::istream& operator>>(std::istream& in, Vector& v)
{
    for(auto& x:v)
        in>>x;
    return in;
}
Vector& Vector::operator*=(const double coef)
{
    for(int i = 0; i < size(); i++)
        (*this)[i]*=coef;
    return *this;
}
Vector& Vector::operator+=(const Vector& other)
{
    for(int i = 0; i < size(); i++)
        (*this)[i]+=other[i];
    return *this;
}
Vector operator*(const double coef, const Vector& v)
{
    Vector res(v.size());
    for(int i = 0; i < v.size(); i++)
        res[i] = coef*v[i];
    return res;
}
Vector operator*(const Vector& v, const double coef)
{
   return coef*v;
}
Vector Vector::operator+(const Vector& other)const
{
    Vector res(size());
    for(int i = 0; i < size(); i++)
    {
        res[i] = (*this)[i] + other[i];
    }
    return res;
}
Vector Vector::operator-(const Vector& other)const
{
    Vector res(size());
    for(int i = 0; i < size(); i++)
    {
        res[i] = (*this)[i] - other[i];
    }
    return res;
}
double Vector::Dot(const Vector& other)const
{
    double res = 0;
    for(int i = 0; i < size(); i++)
        res+=(*this)[i]*other[i];
    return res;
}
