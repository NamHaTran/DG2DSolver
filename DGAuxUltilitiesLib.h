#ifndef DGAUXULTILITIESLIB_H_INCLUDED
#define DGAUXULTILITIESLIB_H_INCLUDED
#include <tuple>
#include <vector>
namespace auxUlti
{
	/*Funtion finds order of edge respective to element*/
	int findEdgeOrder(int element, int edge);

	/*Function checks element's type, return 3 if element is tri, 4 if element is quad*/
	int checkType(int element);

	/*Function gets coordinates of index'th vertex of element*/
	std::tuple<double, double> getElemCornerCoord(int elem, int index);

	/*Function calculates primary variables from conservative variables*/
	void ConserToPri();

	/*Function gets working directory*/
	std::string workingdir();

	/*Function return true if considering element is master of considering edge, otherwise it return false*/
	bool checkMaster(int elem, int edge);

	/*Function gets Jacobi 1D value from J1D array*/
	double getJ1D(int elem, int edge);

	/*Function gets normal vector component (nx or ny) from normalVector array*/
	double getNormVectorComp(int elem, int edge, int dir);

	/*Function gets (a,b) coordinates of Gauss point on surface*/
	std::tuple<double, double> getGaussSurfCoor(int edge, int elem, int nG);

	/*Function executes file .exe*/
	void openFileEXE(std::string location);

	/*Function gets value of all order of accuracy of conservative variables at inputted element and returns output as a vector
	Type 1: rho
	2: rhou
	3: rhov
	4: rhoE*/
	std::vector<double> getElementConserValuesOfOrder(int element, int type);

	/*Function gets value of all order of accuracy of primary variables at inputted element and returns output as a vector
	Type 1: rho
	2: u
	3: v
	4: e
	5: p
	6: T
	7: mu*/
	std::vector<double> getElementPriValuesOfOrder(int element, int type);

	/*Function gets value of all order of accuracy of auxilary variables at inputted element and returns output as a vector
	Type 1: drho
	2: drhou
	3: drhov
	4: drhoE
	
	dir 1: Ox
	dir 2: Oy*/
	std::vector<double> getElementAuxValuesOfOrder(int element, int type, int dir);

	/*Function gets (a,b) coordinates of Gauss point at inside element*/
	std::tuple<double, double> getGaussCoor(int na, int nb);

	/*Function gets group index of edge*/
	int getGrpOfEdge(int edge);

	/*Function gets boundary type of edge*/
	int getBCType(int edge);

	/*Function returns true if problem is subsonic*/
	bool checkSubSonic(double TInf, double uInf, double vInf);

	/*Function returns master element and servant element of edge*/
	std::tuple<int, int> getMasterServantOfEdge(int edge);
}
#endif // DGAUXULTILITIESLIB_H_INCLUDED