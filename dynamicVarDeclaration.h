#ifndef DYNAMICVARDECLARATION_H_INCLUDED
#define DYNAMICVARDECLARATION_H_INCLUDED
#include "VarDeclaration.h"
#include <vector>

namespace meshVar
{
	/*Gauss points on edges*/
	extern std::vector<std::vector<double>> edgeGaussPoints_a, edgeGaussPoints_b;
	
	/*Vector contents BC edges name and location of them on BC values arrays*/
	extern std::vector<int>adressOfBCVals;
}

/*Conservative variables declaration*/
extern std::vector<std::vector<double>> rho, rhou, rhov, rhoE;
extern std::vector<std::vector<double>> rhoN, rhouN, rhovN, rhoEN;

/*Primary variables declaration*/
extern std::vector<std::vector<double>>u, v, e, p, T, mu;

/*Auxilary variables
//X direction*/
extern std::vector<std::vector<double>> rhoX, rhouX, rhovX, rhoEX;

/*Y direction*/
extern std::vector<std::vector<double>> rhoY, rhouY, rhovY, rhoEY;

//time step
extern double dt, runTime;

//Limiting coefficients
extern std::vector<double>
theta1Arr,
theta2Arr;

namespace SurfaceBCFields
{
	extern std::vector<std::vector<double>> rhoBc, rhouBc, rhovBc, rhoEBc;
}
#endif // DYNAMICVARDECLARATION_H_INCLUDED
