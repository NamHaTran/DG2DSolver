#ifndef CONSTDECLARATION_H_INCLUDED
#define CONSTDECLARATION_H_INCLUDED
#include <string>

//Maximum size of Points Array & Elements1D Array
int const pointsArrSize(50000);

//Maximum size of Elements2D Array
int const elements2DArrSize(200000);

//Maximum size of boundary conditions Array
int const bcSize(10);

//Maximum Gauss Points
int const maxGauss(3);

//Maximum Order of accuracy
int const maxOrder(3);

//Maximum size for Gauss-Seidel function
int const maxGS(5);

#endif // CONSTDECLARATION_H_INCLUDED