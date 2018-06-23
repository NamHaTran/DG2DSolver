#include "DGAuxUltilitiesLib.h"
#include "VarDeclaration.h"
#include "DGMath.h"
#include <vector>
#include <math.h>
#include <windows.h>
#include <tuple>  //Include this for returning multiple values in function

namespace auxUlti
{
	int findEdgeOrder(int element, int edge)
	{
		int order(0);
		int pt1(meshVar::inpoed[0][edge]), pt2(meshVar::inpoed[1][edge]);
		int typeElem(checkType(element));
		int maxEdge(0), ABpt1(0), ABpt2(0), BCpt1(0), BCpt2(0), CDpt1(0), CDpt2(0), DApt1(0), DApt2(0), CApt1(0), CApt2(0);

		if (typeElem==4)  //Quad element
		{
			maxEdge = 4;
			ABpt1 = meshVar::Elements2D[element][0];
			ABpt2 = meshVar::Elements2D[element][1];

			BCpt1 = ABpt2;
			BCpt2 = meshVar::Elements2D[element][2];

			CDpt1 = BCpt2;
			CDpt2 = meshVar::Elements2D[element][3];

			DApt1 = CDpt2;
			DApt2 = ABpt1;

			for (int iedge = 0; iedge < maxEdge; iedge++)
			{
				if ((pt1 == ABpt1 && pt2 == ABpt2) || (pt1 == ABpt2 && pt2 == ABpt1))
				{
					order = 0;
					break;
				}
				else if ((pt1 == BCpt1 && pt2 == BCpt2) || (pt1 == BCpt2 && pt2 == BCpt1))
				{
					order = 1;
					break;
				}
				else if ((pt1 == CDpt1 && pt2 == CDpt2) || (pt1 == CDpt2 && pt2 == CDpt1))
				{
					order = 2;
					break;
				}
				else if ((pt1 == DApt1 && pt2 == DApt2) || (pt1 == DApt2 && pt2 == DApt1))
				{
					order = 3;
					break;
				}
			}
		}
		else if (typeElem == 3)  //Tri element
		{
			maxEdge = 3;
			ABpt1 = meshVar::Elements2D[element][0];
			ABpt2 = meshVar::Elements2D[element][1];

			BCpt1 = ABpt2;
			BCpt2 = meshVar::Elements2D[element][2];

			CApt1 = BCpt2;
			CApt2 = ABpt1;

			for (int iedge = 0; iedge < maxEdge; iedge++)
			{
				if ((pt1 == ABpt1 && pt2 == ABpt2) || (pt1 == ABpt2 && pt2 == ABpt1))
				{
					order = 0;
					break;
				}
				else if ((pt1 == BCpt1 && pt2 == BCpt2) || (pt1 == BCpt2 && pt2 == BCpt1))
				{
					order = 1;
					break;
				}
				else if ((pt1 == CApt1 && pt2 == CApt2) || (pt1 == CApt2 && pt2 == CApt1))
				{
					order = 2;
					break;
				}
			}
		}
		return order;
	}

	int checkType(int element)
	{
		int typeElem(0);
		int typeFlag(meshVar::Elements2D[element][3]);
		if (typeFlag<0)
		{
			typeElem = 3;
		}
		else
		{
			typeElem = 4;
		}
		return typeElem;
	}

	std::tuple<double, double> getElemCornerCoord(int elem, int index)
	{
		/*Note: index starts from 0*/
		int pt(meshVar::Elements2D[elem][index]);
		double x(meshVar::Points[pt][0]), y(meshVar::Points[pt][1]);
		return std::make_tuple(x, y);
	}
	
	void ConserToPri()
	{
		for (int ielem = 0; ielem < meshVar::nelem2D; ielem++)
		{
			for (int i = 0; i <= mathVar::orderElem; i++)
			{
				u[ielem][i] = rhou[ielem][i] / rho[ielem][i];
				v[ielem][i] = rhov[ielem][i] / rho[ielem][i];
				T[ielem][i] = math::CalcTFromPriVar(rho[ielem][i], rhou[ielem][i], rhov[ielem][i], rhoE[ielem][i]);
				p[ielem][i] = math::CalcP(T[ielem][i], rho[ielem][i]);
				mu[ielem][i] = math::CalcVisCoef(T[ielem][i]);
			}
		}
	}

	std::string workingdir()
	{
		char buf[256];
		GetCurrentDirectoryA(256, buf);
		return std::string(buf);
	}

	bool checkMaster(int elem, int edge)
	{
		bool master(true);
		if (meshVar::MasterElemOfEdge[edge]==elem)
		{
			master = true;
		}
		else
		{
			master = false;
		}
		return master;
	}

	double getJ1D(int elem, int edge)
	{
		double J(0.0);
		bool master(auxUlti::checkMaster(elem,edge));
		if (master==true)
		{
			J = meshVar::J1D[edge][0];
		}
		else
		{
			J = meshVar::J1D[edge][1];
		}
		return J;
	}

	std::tuple<double, double> getGaussSurfCoor(int edge, int elem, int nG)
	{
		double a(0.0), b(0.0);
		int edgeOrder(auxUlti::findEdgeOrder(elem, edge));
		int elemType(auxUlti::checkType(elem));

		if (elemType == 4)  //quad element
		{
			if (edgeOrder == 0)
			{
				a = mathVar::xGauss[nG];
				b = -1.0;
			}
			else if (edgeOrder == 1)
			{
				a = 1.0;
				b = mathVar::xGauss[nG];
			}
			else if (edgeOrder == 2)
			{
				a = mathVar::xGauss[nG];
				b = 1.0;
			}
			else if (edgeOrder == 3)
			{
				a = -1.0;
				b = mathVar::xGauss[nG];
			}
		}
		else if (elemType == 3)  //tri element
		{
			if (edgeOrder == 0)
			{
				a = mathVar::xGauss[nG];
				b = -1.0;
			}
			else if (edgeOrder == 1)
			{
				a = 1.0;
				b = mathVar::xGauss[nG];
			}
			else if (edgeOrder == 2)
			{
				a = -1.0;
				b = mathVar::xGauss[nG];
			}
		}
		return std::make_tuple(a, b);
	}

	double getNormVectorComp(int elem, int edge, int dir)
	{
		double n(0.0);
		bool master(auxUlti::checkMaster(elem, edge));
		if (dir == 1)  //x direction
		{
			n = meshVar::normalVector[0][edge];
		}
		else if (dir == 2)  //y direction
		{
			n = meshVar::normalVector[1][edge];
		}
		if (master == false)
		{
			n = -n;
		}
		return n;
	}

	void openFileEXE(std::string location)
	{
		LPCSTR sw = location.c_str();
		ShellExecute(GetDesktopWindow(), "open", sw, NULL, NULL, SW_SHOWNORMAL);
	}

	std::vector<double> getElementConserValuesOfOrder(int element, int type)
	{
		std::vector<double> Out(mathVar::orderElem + 1, 0.0);
		if (type == 1)  //rho
		{
			for (int iorder = 0; iorder < mathVar::orderElem; iorder++)
			{
				Out[iorder] = rho[element][iorder];
			}
		}
		else if (type == 2)  //rhou
		{
			for (int iorder = 0; iorder < mathVar::orderElem; iorder++)
			{
				Out[iorder] = rhou[element][iorder];
			}
		}
		else if (type == 3)  //rhov
		{
			for (int iorder = 0; iorder < mathVar::orderElem; iorder++)
			{
				Out[iorder] = rhov[element][iorder];
			}
		}
		else if (type == 4)  //rhoE
		{
			for (int iorder = 0; iorder < mathVar::orderElem; iorder++)
			{
				Out[iorder] = rhoE[element][iorder];
			}
		}
		
		return Out;
	}

	std::vector<double> getElementPriValuesOfOrder(int element, int type)
	{
		std::vector<double> Out(mathVar::orderElem + 1, 0.0);
		if (type == 1)  //rho
		{
			for (int iorder = 0; iorder < mathVar::orderElem; iorder++)
			{
				Out[iorder] = rho[element][iorder];
			}
		}
		else if (type == 2)  //u
		{
			for (int iorder = 0; iorder < mathVar::orderElem; iorder++)
			{
				Out[iorder] = u[element][iorder];
			}
		}
		else if (type == 3)  //v
		{
			for (int iorder = 0; iorder < mathVar::orderElem; iorder++)
			{
				Out[iorder] = v[element][iorder];
			}
		}
		else if (type == 4)  //e
		{
			for (int iorder = 0; iorder < mathVar::orderElem; iorder++)
			{
				Out[iorder] = e[element][iorder];
			}
		}
		else if (type == 5)  //p
		{
			for (int iorder = 0; iorder < mathVar::orderElem; iorder++)
			{
				Out[iorder] = p[element][iorder];
			}
		}
		else if (type == 6)  //T
		{
			for (int iorder = 0; iorder < mathVar::orderElem; iorder++)
			{
				Out[iorder] = T[element][iorder];
			}
		}
		else if (type == 7)  //mu
		{
			for (int iorder = 0; iorder < mathVar::orderElem; iorder++)
			{
				Out[iorder] = mu[element][iorder];
			}
		}

		return Out;
	}

	std::vector<double> getElementAuxValuesOfOrder(int element, int type, int dir)
	{
		std::vector<double> Out(mathVar::orderElem + 1, 0.0);
		if (dir==1)  //Ox direction
		{
			if (type == 1)  //d(rho)x
			{
				for (int iorder = 0; iorder < mathVar::orderElem; iorder++)
				{
					Out[iorder] = rhoX[element][iorder];
				}
			}
			else if (type == 2)  //d(rhou)x
			{
				for (int iorder = 0; iorder < mathVar::orderElem; iorder++)
				{
					Out[iorder] = rhouX[element][iorder];
				}
			}
			else if (type == 3)  //d(rhov)x
			{
				for (int iorder = 0; iorder < mathVar::orderElem; iorder++)
				{
					Out[iorder] = rhovX[element][iorder];
				}
			}
			else if (type == 4)  //d(rhoE)x
			{
				for (int iorder = 0; iorder < mathVar::orderElem; iorder++)
				{
					Out[iorder] = rhoEX[element][iorder];
				}
			}
		}
		else if (dir==2)  //Oy direction
		{
			if (type == 1)  //d(rho)y
			{
				for (int iorder = 0; iorder < mathVar::orderElem; iorder++)
				{
					Out[iorder] = rhoY[element][iorder];
				}
			}
			else if (type == 2)  //d(rhou)y
			{
				for (int iorder = 0; iorder < mathVar::orderElem; iorder++)
				{
					Out[iorder] = rhouY[element][iorder];
				}
			}
			else if (type == 3)  //d(rhov)y
			{
				for (int iorder = 0; iorder < mathVar::orderElem; iorder++)
				{
					Out[iorder] = rhovY[element][iorder];
				}
			}
			else if (type == 4)  //d(rhoE)y
			{
				for (int iorder = 0; iorder < mathVar::orderElem; iorder++)
				{
					Out[iorder] = rhoEY[element][iorder];
				}
			}
		}

		return Out;
	}

	std::tuple<double, double> getGaussCoor(int na, int nb)
	{
		double a(0.0), b(0.0);
		a = mathVar::GaussPts[na][nb][0];
		b= mathVar::GaussPts[na][nb][1];
		return std::make_tuple(a, b);
	}

	int getGrpOfEdge(int edge)
	{
		int grp(meshVar::inpoed[2][edge]);
		return grp;
	}

	int getBCType(int edge)
	{
		int bcType(meshVar::inpoed[3][edge]);
		return bcType;
	}

	bool checkSubSonic(double TInf, double uInf, double vInf)
	{
		bool Out(true);
		double SpeedOfSound(sqrt(material::gamma*material::R*TInf)), Mach(0.0), Velocity(sqrt(uInf*uInf + vInf * vInf));
		Mach = Velocity / SpeedOfSound;
		if (Mach >= 1.0)
		{
			Out = false;
		}

		return Out;
	}

	std::tuple<int, int> getMasterServantOfEdge(int edge)
	{
		int master(meshVar::MasterElemOfEdge[edge]), servant(0), elem1(meshVar::ineled[0][edge]), elem2(meshVar::ineled[1][edge]);
		if (master==elem1)
		{
			servant = elem2;
		}
		else if (master==elem2)
		{
			servant = elem1;
		}
		return std::make_tuple(master, servant);
	}
}