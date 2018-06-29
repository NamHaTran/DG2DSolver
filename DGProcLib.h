#ifndef DGPROCLIB_H_INCLUDED
#define DGPROCLIB_H_INCLUDED
#include <vector>

namespace meshParam
{
	/*Function calculate Gaussian constants*/
	void GaussParam();

	/*Function calculate Jacobian*/
	void JacobianParam();

	/*Function calculate basis function*/
	void basisFcParam();

	/*Function saves coordinates derivatives to array*/
	void derivCoordinates();
}

namespace process
{
	/*Function sets initial value to all elements*/
	void setIniValues();

	/*Function distributes initial value to each order of accuracy of element*/
	std::vector<double> calcIniValues(double iniVal);

	namespace auxEq
	{
		/*Function solves auxilary equation at all elements for auxilary variables*/
		void solveAuxEquation();

		/*Function calculates right hand side term of auxilary equations at ONLY one order*/
		std::vector<double> CalcRHSTerm(int ielem, int order, int dir);

		/*Function returns matrix content Gauss values of <valType> conservative variable at all Gauss points in element volume,
		use only for auxilary equation*/
		std::vector<std::vector<double>> getGaussMatrixOfConserVar(int element, int valType);

		/*Function returns vector content Gauss values of flux of <valType> conservative variable at all Gauss points on the edge,
		use only for auxilary equation*/
		std::vector<double> getGaussVectorOfConserVarFluxesAtInternal(int edge, int element, int nG, double nVectorComp);

		/*Function calculates Auxilary Stiff matrix*/
		std::vector<std::vector<double>> calculateStiffMatrix(int element);

		/*Function calculates value of element of Auxilary Stiff matrix*/
		double calculateStiffMatrixElement(int element, int order1, int order2);
	}

	namespace NSFEq
	{
		/*Function calculates right hand side terms of all conservative variables at ONLY one order*/
		std::vector<double> CalcRHSTerm(int element, int order);

		/*Function calculates Inviscid terms at Gauss point (a, b) and returns 2D matrix
		InviscidTerm 2D array has 4 row 2 column:
		- column 1: Ox direction
		- column 2: Oy direction*/
		std::vector<std::vector<double>> calcGaussInviscidTerm(int element, double a, double b);

		/*Function calculates Viscous terms at Gauss point (a, b) and returns 2D matrix
		ViscousTerm 2D array has 4 row 2 column:
		- column 1: Ox direction
		- column 2: Oy direction*/
		std::vector<std::vector<double>> calcGaussViscousTerm(int element, double a, double b);

		/*Function calculates volume integral terms in NSF equation at ONLY ONE ORDER*/
		std::vector<double> calcVolumeIntegralTerms(int element, int order);

		/*Function calculates surface integral terms in NSF equation at ONLY ONE ORDER*/
		std::vector<double> calcSurfaceIntegralTerms(int element, int order);

		/*Function calculates flux at nGauss point of all conservative variables at internal egde*/
		std::vector<std::vector<double>> getGaussVectorOfConserVarFluxesAtInternal(int edgeName, int element, int nGauss);
	}

	/*Function calculates volume integral terms of auxilary equations of ONLY one order
	User's guide:
	elem: element index
	Ui: array of Gauss values
	direction: index of direction
	1: x direction
	2: y direction*/
	double volumeInte(int elem, std::vector< std::vector<double> > &Ui, int order, int direction);

	/*Function calculates surface integral at 1 surface of auxilary equations of ONLY one order
	User's guide:
	direction: index of direction
	1: x direction
	2: y direction*/
	double surfaceInte(int elem, int edge, std::vector<double> &FluxVector, int order);
}

#endif // DGPROCLIB_H_INCLUDED