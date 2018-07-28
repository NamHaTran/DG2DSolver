#include "VarDeclaration.h"
#include <string>
#include "ConstDeclaration.h"
#include <vector>

namespace systemVar
{
	std::string wD;
	std::string caseName;
	std::string pwd;
	std::string cmd;
	bool endKey(false);

	double CFL(0.5); //Courant number
	double Ttime(0.0); //Total time
	int wrtI(0); //write interval
	bool wrtLog("true");
}

namespace meshVar
{
	double Points[pointsArrSize][3] = {};
	int Elements1D[pointsArrSize][3] = {}, Elements2D[elements2DArrSize][4] = {}, BoundaryType[bcSize][3] = {};
	int nelem1D(0), nelem2D(0), npoin(0), nBc(0);
	double normalVector[2][2 * elements2DArrSize] = {};  //row 1 contents nx, row 2 contents ny
	int MasterElemOfEdge[2 * elements2DArrSize] = {};  //array content master element of iedge, use it with normalVector to get information of normal vector of edge

	/*Default value*/
	//number of nodes per element
	int const nnode(4);
	//number of edges per element (default is 4)
	int const nedel(4);  //change nedel value at file .h

	/*Elements surrounding point*/
	int esup1[4 * pointsArrSize] = {}, esup2[pointsArrSize + 1] = {}, inpoel[5][elements2DArrSize] = {};

	/*Points surrounding point*/
	int psup1[5 * pointsArrSize] = {}, psup2[pointsArrSize + 1] = {};

	/*Elements surrounding element*/
	int esuel[4][elements2DArrSize] = {};  //Default is 4 faces

	/*Edges informations*/
	int inpoed[4][2 * elements2DArrSize] = {};
	/*column 3 contents group which edge belongs to (group 0 is internal group),
	column 4 contents type of boundary (type 0 is internal edge)*/

	/*Edges of element*/
	int inedel[4][elements2DArrSize] = {}, //column index is element index, each row in column contents index of edge belong to element, number of row is 4 because of default quad element
		ineled[3][2 * elements2DArrSize] = {}; //column index is edge index, each row in column contents index of element which edge is belong to, row 3 contents pointer

	/*Variables help to save mesh data*/
	int inpoedCount(0);  //can be used for normalVector, MasterElemOfEdge, ineled

	/*Jacobian*/
	double J2D[elements2DArrSize][maxGauss][maxGauss] = {}, J1D[pointsArrSize][2] = {};

	/*derivatives dx/da, dx/db, dy/da, dy/db*/
	double dxa[elements2DArrSize][maxGauss][maxGauss] = {},
		dxb[elements2DArrSize][maxGauss][maxGauss] = {},
		dya[elements2DArrSize][maxGauss][maxGauss] = {},
		dyb[elements2DArrSize][maxGauss][maxGauss] = {};

	std::vector<std::vector<double>> geoCenter(elements2DArrSize, std::vector<double>(2, 0.0));
	std::vector<double> cellSize(elements2DArrSize, 0.0);
}

namespace mathVar
{
	int nGauss(0), orderElem(0);
	double wGauss[maxGauss] = {}, xGauss[maxGauss] = {}, wGaussLobatto[maxGauss] = {}, xGaussLobatto[maxGauss] = {};
	double B[maxOrder] = {}, dBa[maxOrder] = {}, dBb[maxOrder] = {};
	double BPts[maxOrder][maxGauss][maxGauss] = {}, dBaPts[maxOrder][maxGauss][maxGauss] = {}, dBbPts[maxOrder][maxGauss][maxGauss] = {};
	double GaussPts[maxGauss][maxGauss][2] = {}, //coordinate a is array (..,..,0), coordinate b is array (..,..,1)
		wGaussPts[maxGauss][maxGauss][2] = {}, //weights on a direction (w1) is array (..,..,1), weights on b direction (w2) is array (..,..,2)
		GaussLobattoPts[maxGauss][maxGauss][2] = {},
		wGaussLobattoPts[maxGauss][maxGauss][2] = {};
}

namespace material
{
	double gamma(1.4), R(287.0), Pr(0.72), As(0.001), Ts(110.4), Cp(0.0), Cv(0.0);
}

namespace iniValues
{
	double uIni(0.0), vIni(0.0), wIni(0.0), pIni(0.0), TIni(0.0), muIni(0.0), rhoIni(0.0), eIni(0.0);
}

namespace bcValues
{
	double uBC[bcSize] = { 0.0 },
		vBC[bcSize] = { 0.0 },
		wBC[bcSize] = { 0.0 },
		pBC[bcSize] = { 0.0 },
		TBC[bcSize] = { 0.0 };

	/*Values for Maxwell-Smoluchovsky condition*/
	double TWall[bcSize] = {}, uWall[bcSize] = {}, vWall[bcSize] = {}, wWall[bcSize] = {};

	/*Variables help identify boundary condition at each group*/
	int UBcType[bcSize] = {};
	int pBcType[bcSize] = {};
	int TBcType[bcSize] = {};
}

namespace refValues
{
	double Ma(0.0);
	bool subsonic(false);
}

/*Conservative variables declaration
double rho[elements2DArrSize][maxOrder] = {},
rhou[elements2DArrSize][maxOrder] = {},
rhov[elements2DArrSize][maxOrder] = {},
rhoE[elements2DArrSize][maxOrder] = {};*/
std::vector<std::vector<double>> 
rho(elements2DArrSize, std::vector<double>(maxOrder, 0.0)),
rhou(elements2DArrSize, std::vector<double>(maxOrder, 0.0)),
rhov(elements2DArrSize, std::vector<double>(maxOrder, 0.0)),
rhoE(elements2DArrSize, std::vector<double>(maxOrder, 0.0));

/*Primary variables declaration
double u[elements2DArrSize][maxOrder] = {},
v[elements2DArrSize][maxOrder] = {},
e[elements2DArrSize][maxOrder] = {},
p[elements2DArrSize][maxOrder] = {},
T[elements2DArrSize][maxOrder] = {},
mu[elements2DArrSize][maxOrder] = {};*/
std::vector<std::vector<double>>
u(elements2DArrSize, std::vector<double>(maxOrder, 0.0)),
v(elements2DArrSize, std::vector<double>(maxOrder, 0.0)),
e(elements2DArrSize, std::vector<double>(maxOrder, 0.0)),
p(elements2DArrSize, std::vector<double>(maxOrder, 0.0)),
T(elements2DArrSize, std::vector<double>(maxOrder, 0.0)),
mu(elements2DArrSize, std::vector<double>(maxOrder, 0.0));

/*Auxilary variables
//X direction
double rhoX[elements2DArrSize][maxOrder] = {},
rhouX[elements2DArrSize][maxOrder] = {},
rhovX[elements2DArrSize][maxOrder] = {},
rhoEX[elements2DArrSize][maxOrder] = {};*/
std::vector<std::vector<double>>
rhoX(elements2DArrSize, std::vector<double>(maxOrder, 0.0)),
rhouX(elements2DArrSize, std::vector<double>(maxOrder, 0.0)),
rhovX(elements2DArrSize, std::vector<double>(maxOrder, 0.0)),
rhoEX(elements2DArrSize, std::vector<double>(maxOrder, 0.0));

/*Y direction
double rhoY[elements2DArrSize][maxOrder] = {},
rhouY[elements2DArrSize][maxOrder] = {},
rhovY[elements2DArrSize][maxOrder] = {},
rhoEY[elements2DArrSize][maxOrder] = {};*/
std::vector<std::vector<double>>
rhoY(elements2DArrSize, std::vector<double>(maxOrder, 0.0)),
rhouY(elements2DArrSize, std::vector<double>(maxOrder, 0.0)),
rhovY(elements2DArrSize, std::vector<double>(maxOrder, 0.0)),
rhoEY(elements2DArrSize, std::vector<double>(maxOrder, 0.0));

//time step
double dt(1e-5);

//Limiting coefficients
std::vector<double>
theta1Arr(elements2DArrSize, 0.0),
theta2Arr(elements2DArrSize, 0.0);

/*Mean values
row1: mean rho
row2: mean rhou
row3: mean rhov
row4: mean rhoE*/
std::vector<std::vector<double>> meanVals(elements2DArrSize, std::vector<double>(4, 0.0));

//system settings
namespace sysSetting
{
	int ddtScheme(1);
	int limiter(1);
	double epsilon(1e-13);
}