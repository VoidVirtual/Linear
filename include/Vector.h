#pragma once
#include <vector>
#include <array>
#include <cmath>
#include <iostream>
class Vector:public::std::vector<double>
{
public:
    Vector(){}
    Vector(std::initializer_list<double> _list):std::vector<double>(_list){}
    Vector(int N):std::vector<double>(N){}
    virtual ~Vector(){}
    Vector& operator+=(const Vector& other);
    Vector& operator*=(const double coef);
    Vector operator+(const Vector& other)const;
    Vector operator-(const Vector& other)const;
    double Dot(const Vector& other)const;
    Vector operator*(const Vector& other)const;
    double getNorm()const;
    void Normalize();
    friend std::ostream& operator<<(std::ostream& out, const Vector& v);
    friend std::istream& operator>>(std::istream& in, Vector& v);
    friend Vector operator*(const double coef, const Vector& v);
    friend Vector operator*(const Vector& v,const double coef);
};
