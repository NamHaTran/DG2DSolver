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

	/*Function gets (a,b) coordinates of Gauss point on surface of cell master*/
	std::tuple<double, double> getGaussSurfCoorMaster(int edge, int elem, int nG);

	/*Function executes file .exe*/
	void openFileEXE(std::string location);

	/*Function gets value of all order of accuracy of conservative variables at inputted element and returns output as a vector
	Type 1: rho
	2: rhou
	3: rhov
	4: rhoE*/
	std::vector<double> getElementConserValuesOfOrder(int element, int type);

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
	bool checkSubSonic();

	/*Function checks subsonic flow locally*/
	bool checkSubSonicLocally(double TVal, double uVal, double vVal);

	/*Function returns master element and servant element of edge*/
	std::tuple<int, int> getMasterServantOfEdge(int edge);

	//Function returns cell centroid coordinates and size (cell area)
	std::tuple<double, double, double> getCellMetrics(int element);

	/*Function resize 2D array*/
	void resize2DArray(std::vector<std::vector<double>> &Array, int row, int column);

	/*Function resizes all dynamic arrays, it helps to reduce amount of consumed RAM*/
	void resizeDGArrays();

	/*Function computes coordinates of Gauss point on all edges*/
	void mappingEdges();

	//This function supports for inverse coodinates mapping
	std::vector<std::vector<double>> getVectorGaussSurfCoor(int edge, int elem);

	//Function returns location of input edge on BC values array
	int getAdressOfBCEdgesOnBCValsArray(int edge);

	//Function gets centroid coordinates of inputted cell
	std::tuple<double, double> getCellCentroid(int element);
}
#endif // DGAUXULTILITIESLIB_H_INCLUDED
