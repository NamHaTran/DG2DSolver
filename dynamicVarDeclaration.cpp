#include "VarDeclaration.h"
#include "dynamicVarDeclaration.h"
#include <vector>

namespace meshVar
{
	/*Gauss points on edges*/
	std::vector<std::vector<double>>
		edgeGaussPoints_a(1, std::vector<double>(1, 0.0)),
		edgeGaussPoints_b(1, std::vector<double>(1, 0.0));

	/*Vector contents BC edges name and location of them on BC values arrays*/
	std::vector<int>adressOfBCVals(0, 0);
}

/*Conservative variables declaration*/
std::vector<std::vector<double>>
rho(1, std::vector<double>(1, 0.0)),
rhou(1, std::vector<double>(1, 0.0)),
rhov(1, std::vector<double>(1, 0.0)),
rhoE(1, std::vector<double>(1, 0.0));
std::vector<std::vector<double>>
rhoN(1, std::vector<double>(1, 0.0)),
rhouN(1, std::vector<double>(1, 0.0)),
rhovN(1, std::vector<double>(1, 0.0)),
rhoEN(1, std::vector<double>(1, 0.0));

/*Primary variables declaration*/
std::vector<std::vector<double>>
u(1, std::vector<double>(1, 0.0)),
v(1, std::vector<double>(1, 0.0)),
e(1, std::vector<double>(1, 0.0)),
p(1, std::vector<double>(1, 0.0)),
T(1, std::vector<double>(1, 0.0)),
mu(1, std::vector<double>(1, 0.0));

/*Auxilary variables
//X direction*/
std::vector<std::vector<double>>
rhoX(1, std::vector<double>(1, 0.0)),
rhouX(1, std::vector<double>(1, 0.0)),
rhovX(1, std::vector<double>(1, 0.0)),
rhoEX(1, std::vector<double>(1, 0.0));

/*Y direction*/
std::vector<std::vector<double>>
rhoY(1, std::vector<double>(1, 0.0)),
rhouY(1, std::vector<double>(1, 0.0)),
rhovY(1, std::vector<double>(1, 0.0)),
rhoEY(1, std::vector<double>(1, 0.0));

//Limiting coefficients
std::vector<double>
theta1Arr(1, 1.0),
theta2Arr(1, 1.0);

namespace SurfaceBCFields
{
	std::vector<std::vector<double>> rhoBc(1, std::vector<double>(1, 0.0)),
		rhouBc(1, std::vector<double>(1, 0.0)),
		rhovBc(1, std::vector<double>(1, 0.0)),
		rhoEBc(1, std::vector<double>(1, 0.0));
}