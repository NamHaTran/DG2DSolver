#include "DGPostProcessLib.h"
#include "VarDeclaration.h"
#include "dynamicVarDeclaration.h"
#include "DGAuxUltilitiesLib.h"
#include "DGMath.h"
#include "DGMessagesLib.h"
#include <string>
#include <iostream>
#include <fstream>
#include <direct.h>

namespace debugTool
{
	void checkElemSurPt(int ipoin)
	{
		std::cout << "In SALOME, point ";
		std::cout << ipoin << " is surrounded by elements \n";
		ipoin--;
		for (int ii = meshVar::esup2[ipoin]+1; ii <= meshVar::esup2[ipoin+1]; ii++)
		{
			std::cout << meshVar::esup1[ii]+ meshVar::nelem1D + 1 << std::endl;
		}
	}

	void checkPtsSurPt(int ipoin)
	{
		std::cout << "In SALOME, point ";
		std::cout << ipoin << " is surrounded by points \n";
		ipoin--;
		for (int ii = meshVar::psup2[ipoin] + 1; ii <= meshVar::psup2[ipoin + 1]; ii++)
		{
			std::cout << meshVar::psup1[ii] + 1 << std::endl;
		}
	}

	void checkElemsSurElem(int ielem)
	{
		std::cout << "In SALOME, element ";
		std::cout << ielem << " is surrounded by elements \n";
		ielem = ielem - 1 - meshVar::nelem1D;
		for (int i = 0; i < 4; i++)
		{
			if (meshVar::esuel[i][ielem]>=0)
			{
				std::cout << meshVar::esuel[i][ielem] + 1 + meshVar::nelem1D << std::endl;
			}
		}
	}

	void checkElemInfor(int elem)
	{
		std::cout << "In SALOME, element ";
		std::cout << elem << " has the following edges \n";
		elem = elem - 1 - meshVar::nelem1D;
		int edgeName(0), point1(0), point2(0), grp(0), bcType(0), edgeOrder(0), elemType(0);
		double nx(0.0), ny(0.0);
		bool master(true);
		elemType = auxUlti::checkType(elem);
		
		if (elem>0)
		{
			for (int iedge = 0; iedge < elemType; iedge++)
			{
				edgeName = meshVar::inedel[iedge][elem];
				point1 = meshVar::inpoed[0][edgeName] + 1;
				point2 = meshVar::inpoed[1][edgeName] + 1;
				grp = auxUlti::getGrpOfEdge(edgeName);
				bcType = auxUlti::getBCType(edgeName);
				nx = (auxUlti::getNormVectorComp(elem, edgeName, 1));
				ny = (auxUlti::getNormVectorComp(elem, edgeName, 2));
				edgeOrder = auxUlti::findEdgeOrder(elem, edgeName);
				master = auxUlti::checkMaster(elem, edgeName);
				std::cout << "+ Edge " << edgeName << " created by two points: " << point1 << " and " << point2 << std::endl
					<< "	Coordinates of normal vector are: " << nx << ", " << ny << std::endl
					<< "	Edge belongs to group " << grp << " and has bc type is " << bcType << std::endl
					<< "	Order of edge is " << edgeOrder << std::endl
					<< "	Edge is shared by two elements that are " << meshVar::ineled[0][edgeName] + 1 + meshVar::nelem1D << " and " << meshVar::ineled[1][edgeName] + 1 + meshVar::nelem1D << std::endl;
				if (master)
				{
					std::cout << "	Considering element is a master of edge\n";
				}
				else
				{
					std::cout << "	Considering element is not a master of edge\n";
				}
			}
		}
		else
		{
			std::cout << "	Considering element is not a polygonal element!\n";
		}

	}

	void checkPointValue(int element)
	{
		double rhoVal(math::pointValue(element, 0.0, 0.0, 1, 1)),
			uVal(math::pointValue(element, 0.0, 0.0, 2, 1)),
			vVal(math::pointValue(element, 0.0, 0.0, 3, 1)),
			//eVal(math::pointValue(element, 0.0, 0.0, 4, 1)),
			pVal(math::pointValue(element, 0.0, 0.0, 5, 1)),
			TVal(math::pointValue(element, 0.0, 0.0, 6, 1));
			//muVal(math::pointValue(element, 0.0, 0.0, 7, 1));
		std::cout << "- Element: " << element + meshVar::nelem1D + 1 << std::endl
			<< "- rho: " << rhoVal << std::endl
			<< "- u: " << uVal << std::endl
			<< "- v: " << vVal << std::endl
			//<< "- e: " << eVal << std::endl
			<< "- p: " << pVal << std::endl
			<< "- T: " << TVal << std::endl << std::endl;
			//<< "- mu: " << muVal << std::endl << std::endl;
	}
}

namespace DG2Matlab
{
	void createMatlabCode()
	{
		std::string fileName("MeshPlot.m"), Loc(systemVar::wD + "\\CASES\\" + systemVar::caseName + "\\matlabFile");
		std::string meshPlotLoc(Loc + "\\" + fileName);
		std::ofstream meshPlotFlux(meshPlotLoc.c_str());
		std::string code(" ");
		code = R"(
clear;
clc;
inpoed=dlmread('inpoed.txt',' ',13,1);
points=dlmread('Points.txt');
nedge=size(inpoed,2);
figure;
hold on;
for iedge=1:nedge
    point1=inpoed(1,iedge);
    point2=inpoed(2,iedge);
    XCoor=[points(point1+1,2) points(point2+1,2)];
    YCoor=[points(point1+1,3) points(point2+1,3)];
    plot(XCoor,YCoor,'-b');
end
axis equal;
)";
		if (meshPlotFlux)
		{
			meshPlotFlux << code;
		}
		else
		{
			message::writeLog(systemVar::pwd, systemVar::caseName, message::opFError("MeshPlot.m", meshPlotLoc));
		}
	}

	void exportData(int iter)
	{
		double xC(0.0), yC(0.0), aC(0.0), bC(0.0);
		std::string iter_str = std::to_string(iter);
		std::string fileName("rho"+ iter_str +".txt"), Loc(systemVar::wD + "\\CASES\\" + systemVar::caseName + "\\matlabFile");
		std::string rhoLoc(Loc + "\\" + fileName);
		std::ofstream rhoFlux(rhoLoc.c_str());
		if (rhoFlux)
		{
			for (int i = 0; i < meshVar::nelem2D; i++)
			{
				for (int j = 0; j <= mathVar::orderElem; j++)
				{
					rhoFlux << rho[i][j] << " ";
				}
				std::tie(xC, yC) = auxUlti::getCellCentroid(i);
				std::tie(aC, bC) = math::inverseMapping(i, xC, yC);
				rhoFlux << math::pointValue(i, aC, bC, 1, 2) << "\n";
			}
		}
		else
		{
			message::writeLog(systemVar::pwd, systemVar::caseName, message::opFError("rho.txt", rhoLoc));
		}
		fileName = "rhou" + iter_str + ".txt";
		std::string rhouLoc(Loc + "\\" + fileName);
		std::ofstream rhouFlux(rhouLoc.c_str());
		if (rhouFlux)
		{
			for (int i = 0; i < meshVar::nelem2D; i++)
			{
				for (int j = 0; j <= mathVar::orderElem; j++)
				{
					rhouFlux << rhou[i][j] << " ";
				}
				std::tie(xC, yC) = auxUlti::getCellCentroid(i);
				std::tie(aC, bC) = math::inverseMapping(i, xC, yC);
				rhouFlux << math::pointValue(i, aC, bC, 2, 2) << "\n";
			}
		}
		else
		{
			message::writeLog(systemVar::pwd, systemVar::caseName, message::opFError("rhou.txt", rhoLoc));
		}
		fileName = "rhov" + iter_str + ".txt";
		std::string rhovLoc(Loc + "\\" + fileName);
		std::ofstream rhovFlux(rhovLoc.c_str());
		if (rhovFlux)
		{
			for (int i = 0; i < meshVar::nelem2D; i++)
			{
				for (int j = 0; j <= mathVar::orderElem; j++)
				{
					rhovFlux << rhov[i][j] << " ";
				}
				std::tie(xC, yC) = auxUlti::getCellCentroid(i);
				std::tie(aC, bC) = math::inverseMapping(i, xC, yC);
				rhovFlux << math::pointValue(i, aC, bC, 3, 2) << "\n";
			}
		}
		else
		{
			message::writeLog(systemVar::pwd, systemVar::caseName, message::opFError("rhov.txt", rhoLoc));
		}
		fileName = "rhoE" + iter_str + ".txt";
		std::string rhoELoc(Loc + "\\" + fileName);
		std::ofstream rhoEFlux(rhoELoc.c_str());
		if (rhoEFlux)
		{
			for (int i = 0; i < meshVar::nelem2D; i++)
			{
				for (int j = 0; j <= mathVar::orderElem; j++)
				{
					rhoEFlux << rhoE[i][j] << " ";
				}
				std::tie(xC, yC) = auxUlti::getCellCentroid(i);
				std::tie(aC, bC) = math::inverseMapping(i, xC, yC);
				rhoEFlux << math::pointValue(i, aC, bC, 4, 2) << "\n";
			}
		}
		else
		{
			message::writeLog(systemVar::pwd, systemVar::caseName, message::opFError("rhoE.txt", rhoLoc));
		}

		fileName = "cellCentroid.txt"; 
		std::string cellCentroidLoc(Loc + "\\" + fileName);
		std::ofstream cellCentroidFlux(cellCentroidLoc.c_str());
		if (cellCentroidFlux)
		{
			for (int i = 0; i < meshVar::nelem2D; i++)
			{
				for (int j = 0; j < 2; j++)
				{
					cellCentroidFlux << meshVar::geoCenter[i][j] << " ";
				}
				cellCentroidFlux << "\n";
			}
		}
		else
		{
			message::writeLog(systemVar::pwd, systemVar::caseName, message::opFError("cellCentroid.txt", cellCentroidLoc));
		}
	}
}

namespace DG2Tecplot
{
	std::vector<double> calcNodeValues(int valType)
	{
		//Only for primary variables
		std::vector<int>elemSurPt;
		std::vector<double>nodeValues(meshVar::npoin,0.0);
		int numElemSurPt(0), elemId(-1);
		double aP(0.0), bP(0.0), localValue(0.0);
		for (int ipoint = 0; ipoint < meshVar::npoin; ipoint++)
		{
			localValue = 0.0;
			elemSurPt = auxUlti::postProcess::getElementsSurroundingPoint(ipoint);
			numElemSurPt = elemSurPt.size();
			for (int nelem = 0; nelem < numElemSurPt; nelem++)
			{
				elemId = elemSurPt[nelem];
				std::tie(aP, bP) = auxUlti::postProcess::findPointCoorInStandardSpace(ipoint, elemId);
				localValue += math::pointValue(elemId, aP, bP, valType, 1);
			}
			nodeValues[ipoint] = (localValue / numElemSurPt);
		}
		return nodeValues;
	}

	std::vector<double> calcCellCenteredValues(int valType)
	{
		std::vector<double>cellCenteredValues(meshVar::nelem2D, 0.0);
		double xC(0.0), yC(0.0), aC(0.0), bC(0.0);
		for (int ielem = 0; ielem < meshVar::nelem2D; ielem++)
		{
			std::tie(xC, yC) = auxUlti::getCellCentroid(ielem);
			std::tie(aC, bC) = math::inverseMapping(ielem, xC, yC);
			cellCenteredValues[ielem] = math::pointValue(ielem, aC, bC, valType, 1);
		}
		return cellCenteredValues;
	}

	void exportNodeData(int iter)
	{
		std::vector<double>nodeRho, node_u, node_v, node_p, node_T, node_uMag(meshVar::npoin, 0.0);
		nodeRho = DG2Tecplot::calcNodeValues(1);
		node_u = DG2Tecplot::calcNodeValues(2);
		node_v = DG2Tecplot::calcNodeValues(3);
		node_p = DG2Tecplot::calcNodeValues(5);
		node_T = DG2Tecplot::calcNodeValues(6);
		for (int i = 0; i < meshVar::npoin; i++)
		{
			node_uMag[i] = sqrt(pow(node_u[i], 2) + pow(node_v[i], 2));
		}

		std::string iter_str = std::to_string(iter);
		std::string fileName(systemVar::caseName + ".dat"), Loc(systemVar::wD + "\\CASES\\" + systemVar::caseName + "\\TecplotFile\\" + iter_str), code;
		_mkdir(Loc.c_str());
		
		std::string fileLoc(Loc + "\\" + fileName);
		std::ofstream fileFlux(fileLoc.c_str());

		code = R"(
TITLE     = "DG2D to Tecplot"
VARIABLES = "X", "Y", "RHO", "U" ,"V", "VELOCITY_MAG", "P", "T"
ZONE T="ZONE 1"
ZONETYPE=FEQUADRILATERAL
DATAPACKING=BLOCK)";

		if (fileFlux)
		{
			fileFlux << code << std::endl << "NODES=" << std::to_string(meshVar::npoin) << ", " << "ELEMENTS=" << std::to_string(meshVar::nelem2D) << std::endl;
			//X
			int counter(0);
			for (int i = 0; i < meshVar::npoin; i++)
			{
				counter++;
				if (counter == 1000)
				{
					fileFlux << std::endl;
					counter = 0;
				}
				fileFlux << meshVar::Points[i][0] << " ";
			}
			fileFlux << std::endl;

			//Y
			counter = 0;
			for (int i = 0; i < meshVar::npoin; i++)
			{
				counter++;
				if (counter == 1000)
				{
					fileFlux << std::endl;
					counter = 0;
				}
				fileFlux << meshVar::Points[i][1] << " ";
			}
			fileFlux << std::endl;

			//RHO
			counter = 0;
			for (int i = 0; i < meshVar::npoin; i++)
			{
				counter++;
				if (counter == 1000)
				{
					fileFlux << std::endl;
					counter = 0;
				}
				fileFlux << nodeRho[i] << " ";
			}
			fileFlux << std::endl;

			//U
			counter = 0;
			for (int i = 0; i < meshVar::npoin; i++)
			{
				counter++;
				if (counter == 1000)
				{
					fileFlux << std::endl;
					counter = 0;
				}
				fileFlux << node_u[i] << " ";
			}
			fileFlux << std::endl;

			//U
			counter = 0;
			for (int i = 0; i < meshVar::npoin; i++)
			{
				counter++;
				if (counter == 1000)
				{
					fileFlux << std::endl;
					counter = 0;
				}
				fileFlux << node_v[i] << " ";
			}
			fileFlux << std::endl;

			//VELOCITY_MAG
			counter = 0;
			for (int i = 0; i < meshVar::npoin; i++)
			{
				counter++;
				if (counter == 1000)
				{
					fileFlux << std::endl;
					counter = 0;
				}
				fileFlux << node_uMag[i] << " ";
			}
			fileFlux << std::endl;

			//P
			counter = 0;
			for (int i = 0; i < meshVar::npoin; i++)
			{
				counter++;
				if (counter == 1000)
				{
					fileFlux << std::endl;
					counter = 0;
				}
				fileFlux << node_p[i] << " ";
			}
			fileFlux << std::endl;

			//T
			counter = 0;
			for (int i = 0; i < meshVar::npoin; i++)
			{
				counter++;
				if (counter == 1000)
				{
					fileFlux << std::endl;
					counter = 0;
				}
				fileFlux << node_T[i] << " ";
			}
			fileFlux << std::endl;

			//CONNECTIVITY
			for (int ielem = 0; ielem < meshVar::nelem2D; ielem++)
			{
				int elemType(auxUlti::checkType(ielem));
				for (int ipoin = 0; ipoin < 3; ipoin++)
				{
					fileFlux << meshVar::Elements2D[ielem][ipoin] + 1 << " ";
				}
				switch (elemType)
				{
				case 3:
				{
					fileFlux << meshVar::Elements2D[ielem][0] + 1 << std::endl;
					break;
				}
				case 4:
				{
					fileFlux << meshVar::Elements2D[ielem][3] + 1 << std::endl;
					break;
				}
				default:
					break;
				}
			}
		}
		else
		{
			message::writeLog(systemVar::pwd, systemVar::caseName, message::opFError(systemVar::caseName + ".dat", fileLoc));
		}
	}

	void exportCellCenteredData(int iter)
	{
		std::vector<double>nodeRho, node_u, node_v, node_p, node_T, node_uMag(meshVar::npoin, 0.0);
		nodeRho = DG2Tecplot::calcCellCenteredValues(1);
		node_u = DG2Tecplot::calcCellCenteredValues(2);
		node_v = DG2Tecplot::calcCellCenteredValues(3);
		node_p = DG2Tecplot::calcCellCenteredValues(5);
		node_T = DG2Tecplot::calcCellCenteredValues(6);
		double xC(0.0), yC(0.0);
		for (int i = 0; i < meshVar::npoin; i++)
		{
			node_uMag[i] = sqrt(pow(node_u[i], 2) + pow(node_v[i], 2));
		}

		std::string iter_str = std::to_string(iter);
		std::string fileName(systemVar::caseName + ".dat"), Loc(systemVar::wD + "\\CASES\\" + systemVar::caseName + "\\TecplotFile\\" + iter_str), code;
		_mkdir(Loc.c_str());

		std::string fileLoc(Loc + "\\" + fileName);
		std::ofstream fileFlux(fileLoc.c_str());

		code = R"(
TITLE     = "DG2D to Tecplot"
VARIABLES = "X", "Y", "RHO", "U" ,"V", "VELOCITY_MAG", "P", "T"
ZONE T="ZONE 1"
ZONETYPE=FEQUADRILATERAL
DATAPACKING=BLOCK
VARLOCATION=CELLCENTERED)";

		if (fileFlux)
		{
			fileFlux << code << std::endl << "NODES=" << std::to_string(meshVar::npoin) << ", " << "ELEMENTS=" << std::to_string(meshVar::nelem2D) << std::endl;
			//X
			int counter(0);
			for (int i = 0; i < meshVar::nelem2D; i++)
			{
				counter++;
				std::tie(xC, yC) = auxUlti::getCellCentroid(i);
				if (counter==1000)
				{
					fileFlux << std::endl;
					counter = 0;
				}
				fileFlux << xC << " ";
			}
			fileFlux << std::endl;

			//Y
			counter = 0;
			for (int i = 0; i < meshVar::nelem2D; i++)
			{
				counter++;
				if (counter == 1000)
				{
					fileFlux << std::endl;
					counter = 0;
				}
				std::tie(xC, yC) = auxUlti::getCellCentroid(i);
				fileFlux << yC << " ";
			}
			fileFlux << std::endl;

			//RHO
			counter = 0;
			for (int i = 0; i < meshVar::nelem2D; i++)
			{
				counter++;
				if (counter == 1000)
				{
					fileFlux << std::endl;
					counter = 0;
				}
				fileFlux << nodeRho[i] << " ";
			}
			fileFlux << std::endl;

			//U
			counter = 0;
			for (int i = 0; i < meshVar::nelem2D; i++)
			{
				counter++;
				if (counter == 1000)
				{
					fileFlux << std::endl;
					counter = 0;
				}
				fileFlux << node_u[i] << " ";
			}
			fileFlux << std::endl;

			//U
			counter = 0;
			for (int i = 0; i < meshVar::nelem2D; i++)
			{
				counter++;
				if (counter == 1000)
				{
					fileFlux << std::endl;
					counter = 0;
				}
				fileFlux << node_v[i] << " ";
			}
			fileFlux << std::endl;

			//VELOCITY_MAG
			counter = 0;
			for (int i = 0; i < meshVar::nelem2D; i++)
			{
				counter++;
				if (counter == 1000)
				{
					fileFlux << std::endl;
					counter = 0;
				}
				fileFlux << node_uMag[i] << " ";
			}
			fileFlux << std::endl;

			//P
			counter = 0;
			for (int i = 0; i < meshVar::nelem2D; i++)
			{
				counter++;
				if (counter == 1000)
				{
					fileFlux << std::endl;
					counter = 0;
				}
				fileFlux << node_p[i] << " ";
			}
			fileFlux << std::endl;

			//T
			counter = 0;
			for (int i = 0; i < meshVar::nelem2D; i++)
			{
				counter++;
				if (counter == 1000)
				{
					fileFlux << std::endl;
					counter = 0;
				}
				fileFlux << node_T[i] << " ";
			}
			fileFlux << std::endl;

			//CONNECTIVITY
			for (int ielem = 0; ielem < meshVar::nelem2D; ielem++)
			{
				int elemType(auxUlti::checkType(ielem));
				for (int ipoin = 0; ipoin < 3; ipoin++)
				{
					fileFlux << meshVar::Elements2D[ielem][ipoin] + 1 << " ";
				}
				switch (elemType)
				{
				case 3:
				{
					fileFlux << meshVar::Elements2D[ielem][0] + 1 << std::endl;
					break;
				}
				case 4:
				{
					fileFlux << meshVar::Elements2D[ielem][3] + 1 << std::endl;
					break;
				}
				default:
					break;
				}
			}
		}
		else
		{
			message::writeLog(systemVar::pwd, systemVar::caseName, message::opFError(systemVar::caseName + ".dat", fileLoc));
		}
	}
}