#include <iostream>
#include <fstream>
#include "Matrix.h"
#include "ExtendedMatrix.h"
#include "Vector.h"
#include <ctime>
int main(int argc, char** argv)
{
    srand(time(0));
    Matrix A("Matrix.txt");
    Matrix B = {
                    {3,0,0},
                    {0,3,0},
                    {0,0,3}
               };
    ExtendedMatrix C(3,RANDOM|SYMMETRIC);
    std::cout<<A*B<<"\n";
    std::cout<<C<<"\n";
}
