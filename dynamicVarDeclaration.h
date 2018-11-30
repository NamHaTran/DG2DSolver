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

/*Primary variables declaration
extern std::vector<std::vector<double>>u, v, e, p, T, mu;*/

/*Auxilary variables
//X direction*/
extern std::vector<std::vector<double>> rhoX, rhouX, rhovX, rhoEX;

/*Y direction*/
extern std::vector<std::vector<double>> rhoY, rhouY, rhovY, rhoEY;

/*Interface values*/
//conservative variable
extern std::vector<std::vector<double>> interface_rho, interface_rhou, interface_rhov, interface_rhoE;

//auxilary equaiton
extern std::vector<std::vector<double>> aux_interface_rho, aux_interface_rhou, aux_interface_rhov, aux_interface_rhoE;

//X direction*/
extern std::vector<std::vector<double>> invis_interface_rhoX, invis_interface_rhouX, invis_interface_rhovX, invis_interface_rhoEX,
Vis_interface_rhoX, Vis_interface_rhouX, Vis_interface_rhovX, Vis_interface_rhoEX;

/*Y direction*/
extern std::vector<std::vector<double>> invis_interface_rhoY, invis_interface_rhouY, invis_interface_rhovY, invis_interface_rhoEY,
Vis_interface_rhoY, Vis_interface_rhouY, Vis_interface_rhovY, Vis_interface_rhoEY;

//Lax-Friedrich constant
extern std::vector<double> LxFConst;

//time step
extern double dt, runTime;

//Limiting coefficients
extern std::vector<double>
theta1Arr,
theta2Arr;

namespace SurfaceBCFields
{
	extern std::vector<std::vector<double>> rhoBc, rhouBc, rhovBc, rhoEBc;
	extern std::vector<int>BCPoints;
}
#endif // DYNAMICVARDECLARATION_H_INCLUDED
