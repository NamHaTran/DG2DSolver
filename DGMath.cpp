#include "DGMath.h"
#include "VarDeclaration.h"
#include "dynamicVarDeclaration.h"
#include "DGAuxUltilitiesLib.h"
#include <math.h>
#include <tuple>  //Include this for returning multiple values in function
#include <algorithm>

namespace math
{
	void Gauss(int nGauss)
	{
		if (nGauss==0)
		{
			mathVar::xGauss[nGauss] = 0.0;
			mathVar::wGauss[nGauss] = 2.0;
		}
		else if (nGauss==1)
		{
			mathVar::xGauss[nGauss - 1] = -0.577350269189625764509148780502;
			mathVar::xGauss[nGauss] = 0.577350269189625764509148780502;

			mathVar::wGauss[nGauss - 1] = 1.0;
			mathVar::wGauss[nGauss] = 1.0;
		}
		else if (nGauss==2)
		{
			mathVar::xGauss[nGauss - 2] = -0.774596669241483377035853079956;
			mathVar::xGauss[nGauss - 1] = 0.0;
			mathVar::xGauss[nGauss] = 0.774596669241483377035853079956;

			mathVar::wGauss[nGauss - 2] = 0.555555555555555555555;
			mathVar::wGauss[nGauss - 1] = 0.888888888888888888888;
			mathVar::wGauss[nGauss] = 0.55555555555555555555555;
		}
	}

	void GaussLobatto(int nGauss)
	{
		if (nGauss == 0)
		{
			mathVar::xGaussLobatto[nGauss] = 0.0;
			mathVar::wGaussLobatto[nGauss] = 2.0;
		}
		else if (nGauss == 1)
		{
			mathVar::xGaussLobatto[nGauss - 1] = -1.0;
			mathVar::xGaussLobatto[nGauss] = 1.0;

			mathVar::wGaussLobatto[nGauss - 1] = 1.0;
			mathVar::wGaussLobatto[nGauss] = 1.0;
		}
		else if (nGauss == 2)
		{
			mathVar::xGaussLobatto[nGauss - 2] = -1.0;
			mathVar::xGaussLobatto[nGauss - 1] = 0.0;
			mathVar::xGaussLobatto[nGauss] = 1.0;

			mathVar::wGaussLobatto[nGauss - 2] = 1.0/3.0;
			mathVar::wGaussLobatto[nGauss - 1] = 4.0/3.0;
			mathVar::wGaussLobatto[nGauss] = 1.0/3.0;
		}
	}

	void basisFc(double a, double b)
	{
		for (int i = 0; i <= mathVar::orderElem; i++)
		{
			if (i==0)
			{
				mathVar::B[i] = 1.0;
			}
			else if (i==2)
			{
				mathVar::B[i] = (3.0*b + 1.0) / 2.0;
			}
			else if (i==1)
			{
				mathVar::B[i] = a * (1 - b);
			}
			else if (i==3)
			{
				mathVar::B[i] = -0.5 + b + 5.0*pow(b, 2) / 2.0;
			}
		}
	}

	void dBasisFc(double a, double b)
	{
		for (int i = 0; i <= mathVar::orderElem; i++)
		{
			if (i == 0)
			{
				mathVar::dBa[i] = 0;
				mathVar::dBb[i] = 0;
			}
			else if (i == 2)
			{
				mathVar::dBa[i] = 0;
				mathVar::dBb[i] = 3.0/2.0;
			}
			else if (i == 1)
			{
				mathVar::dBa[i] = 1 - b;
				mathVar::dBb[i] = -a;
			}
			else if (i == 3)
			{
				
			}
		}
	}

	double J2DCal(int elem, double a, double b)
	{
		double Jacobi(0.0), xA(0.0), xB(0.0), xC(0.0), xD(0.0), yA(0.0), yB(0.0), yC(0.0), yD(0.0);
		int elemType(auxUlti::checkType(elem));
		if (elemType == 4)  //Quad
		{
			std::tie(xA, yA) = auxUlti::getElemCornerCoord(elem, 0);
			std::tie(xB, yB) = auxUlti::getElemCornerCoord(elem, 1);
			std::tie(xC, yC) = auxUlti::getElemCornerCoord(elem, 2);
			std::tie(xD, yD) = auxUlti::getElemCornerCoord(elem, 3);

			Jacobi = math::jacobianQuad(xA, xB, xC, xD, yA, yB, yC, yD, a, b);
		}
		else if (elemType == 3)  //Tri
		{
			std::tie(xA, yA) = auxUlti::getElemCornerCoord(elem, 0);
			std::tie(xB, yB) = auxUlti::getElemCornerCoord(elem, 1);
			std::tie(xC, yC) = auxUlti::getElemCornerCoord(elem, 2);

			Jacobi = math::jacobianTri(xA, xB, xC, yA, yB, yC, a, b);
		}

		return Jacobi;
	}

	std::tuple<double, double> J1DCal(int edge)
	{
		int master(0), servant(0), masterIndex(0), servantIndex(0), masterType(0), servantType(0);
		double JMaster(0.0), JServant(0.0), xA(0.0), xB(0.0), xC(0.0), xD(0.0), yA(0.0), yB(0.0), yC(0.0), yD(0.0);

		if (meshVar::ineled[0][edge]>meshVar::ineled[1][edge])
		{
			master = meshVar::ineled[0][edge];
			servant = meshVar::ineled[1][edge];
		}
		else if (meshVar::ineled[0][edge]<meshVar::ineled[1][edge])
		{
			master = meshVar::ineled[1][edge];
			servant = meshVar::ineled[0][edge];
		}

		masterIndex = auxUlti::findEdgeOrder(master, edge);
		servantIndex = auxUlti::findEdgeOrder(servant, edge);
		masterType = auxUlti::checkType(master);
		servantType = auxUlti::checkType(servant);

		if (masterType==4)
		{
			std::tie(xA, yA) = auxUlti::getElemCornerCoord(master, 0);
			std::tie(xB, yB) = auxUlti::getElemCornerCoord(master, 1);
			std::tie(xC, yC) = auxUlti::getElemCornerCoord(master, 2);
			std::tie(xD, yD) = auxUlti::getElemCornerCoord(master, 3);

			JMaster = math::jacobian1DQuad(masterIndex, xA, xB, xC, xD, yA, yB, yC, yD);
		}
		else if (masterType==3)
		{
			std::tie(xA, yA) = auxUlti::getElemCornerCoord(master, 0);
			std::tie(xB, yB) = auxUlti::getElemCornerCoord(master, 1);
			std::tie(xC, yC) = auxUlti::getElemCornerCoord(master, 2);

			JMaster = math::jacobian1DTri(masterIndex, xA, xB, xC, yA, yB, yC);
		}

		if (servantType == 4)
		{
			std::tie(xA, yA) = auxUlti::getElemCornerCoord(servant, 0);
			std::tie(xB, yB) = auxUlti::getElemCornerCoord(servant, 1);
			std::tie(xC, yC) = auxUlti::getElemCornerCoord(servant, 2);
			std::tie(xD, yD) = auxUlti::getElemCornerCoord(servant, 3);

			JServant = math::jacobian1DQuad(servantIndex, xA, xB, xC, xD, yA, yB, yC, yD);
		}
		else if (servantType == 3)
		{
			std::tie(xA, yA) = auxUlti::getElemCornerCoord(servant, 0);
			std::tie(xB, yB) = auxUlti::getElemCornerCoord(servant, 1);
			std::tie(xC, yC) = auxUlti::getElemCornerCoord(servant, 2);

			JServant = math::jacobian1DTri(servantIndex, xA, xB, xC, yA, yB, yC);
		}

		return std::make_tuple(JMaster, JServant);
	}

	double iniIntegral(int order)
	{
		double w1(0.0), w2(0.0), integral(0.0);
		for (int na = 0; na <= mathVar::nGauss; na++)
		{
			for (int nb = 0; nb <= mathVar::nGauss; nb++)
			{
				w1 = mathVar::wGaussPts[na][nb][0];
				w2 = mathVar::wGaussPts[na][nb][1];
				integral += w1 * w2* mathVar::BPts[order][na][nb];
			}
		}
		return integral;
	}

	double volumeInte(std::vector< std::vector<double> > &Fvalue, int elem)
	{
		double J2D(0.0), w1(0.0), w2(0.0), integral(0.0);
		for (int na = 0; na <= mathVar::nGauss; na++)
		{
			for (int nb = 0; nb <= mathVar::nGauss; nb++)
			{
				w1 = mathVar::wGaussPts[na][nb][0];
				w2 = mathVar::wGaussPts[na][nb][1];
				J2D = meshVar::J2D[elem][na][nb];
				integral += w1 * w2*(J2D)*Fvalue[na][nb];
			}
		}
		return integral;
	}

	double jacobianQuad(double xA, double xB, double xC, double xD, double yA, double yB, double yC, double yD, double a, double b)
	{
		double dxa(0.0), dxb(0.0), dya(0.0), dyb(0.0);
		double jQuad(0.0);

		dxa = (1.0 / 4.0)*(-xA + xB - xD + xC) + (1.0 / 4.0)*(xA - xB - xD + xC)*b;
		dxb = (1.0 / 4.0)*(-xA - xB + xD + xC) + (1.0 / 4.0)*(xA - xB - xD + xC)*a;
		dya = (1.0 / 4.0)*(-yA + yB - yD + yC) + (1.0 / 4.0)*(yA - yB - yD + yC)*b;
		dyb = (1.0 / 4.0)*(-yA - yB + yD + yC) + (1.0 / 4.0)*(yA - yB - yD + yC)*a;
		jQuad = dxa * dyb - dxb * dya;
		return jQuad;
	}

	double jacobianTri(double xA, double xB, double xC, double yA, double yB, double yC, double a, double b)
	{
		double dxa(0.0), dxb(0.0), dya(0.0), dyb(0.0);
		double jTri(0.0);

		dxa = (1 - b)*(xB - xA) / 4.0;
		dxb = a * (xA - xB) / 4.0 + (-xA - xB + 2 * xC) / 4.0;
		dya = (1 - b)*(yB - yA) / 4.0;
		dyb = a * (yA - yB) / 4.0 + (-yA - yB + 2 * yC) / 4.0;
		jTri = dxa * dyb - dxb * dya;
		return jTri;
	}

	double jacobian1DQuad(int edgeIndex, double xA, double xB, double xC, double xD, double yA, double yB, double yC, double yD)
	{
		double dx(0.0), dy(0.0), C(0.0);
		double jacobi(0.0);
		if ((edgeIndex == 0) || (edgeIndex == 2))  //AB or CD
		{
			if (edgeIndex == 0)
			{
				C = -1.0;
			}
			else if (edgeIndex == 2)
			{
				C = 1.0;
			}
			dx = (1.0 / 4.0)*(-xA + xB - xD + xC) + (1.0 / 4.0)*(xA - xB - xD + xC)*C;
			dy = (1.0 / 4.0)*(-yA + yB - yD + yC) + (1.0 / 4.0)*(yA - yB - yD + yC)*C;
		}
		else if ((edgeIndex == 1) || (edgeIndex == 3))  //BC or DA
		{
			if (edgeIndex == 1)
			{
				C = 1.0;
			}
			else if (edgeIndex == 3)
			{
				C = -1.0;
			}
			dx = (1.0 / 4.0)*(-xA - xB + xD + xC) + (1.0 / 4.0)*(xA - xB - xD + xC)*C;
			dy = (1.0 / 4.0)*(-yA - yB + yD + yC) + (1.0 / 4.0)*(yA - yB - yD + yC)*C;
		}

		jacobi = std::sqrt(pow(dx, 2) + pow(dy, 2));
		return jacobi;
	}

	double jacobian1DTri(int edgeIndex, double xA, double xB, double xC, double yA, double yB, double yC)
	{
		double dx(0.0), dy(0.0), C(0.0);
		double jacobi(0.0);
		if ((edgeIndex == 1) || (edgeIndex == 2))  //BC or CA
		{
			if (edgeIndex == 1)
			{
				C = 1.0;
			}
			else if (edgeIndex == 2)
			{
				C = -1.0;
			}
			dx = C * (xA - xB) / 4.0 + (-xA - xB + 2 * xC) / 4.0;
			dy = C * (yA - yB) / 4.0 + (-yA - yB + 2 * yC) / 4.0;
		}
		else if ((edgeIndex == 0))  //AB
		{
			C = -1.0;
			dx = (1 - C)*(xB - xA) / 4.0;
			dy = (1 - C)*(yB - yA) / 4.0;
		}

		jacobi = std::sqrt(pow(dx, 2) + pow(dy, 2));
		return jacobi;
	}

	std::vector<double> SolveSysEqs(std::vector< std::vector<double> > &a, std::vector<double> &b)
	{
		int n = static_cast<int>(b.size());
		double eMax(1e-9), e(1.0), sum(0.0), xi(0.0);
		std::vector<double> results(n, 1.0);
		int counter(0);
		//eMax = fabs(*std::min_element(b.begin(), b.end())) / 1e3;

		while (e>eMax && counter<=15)
		{
			for (int i = 0; i < n; i++)
			{
				sum = b[i];
				for (int j = 0; j < n; j++)
				{
					if (i != j)
					{
						sum -= a[i][j] * results[j];
					}
				}
				xi = sum / a[i][i];
				results[i] = xi;
			}
			e = math::errorGS(a, b, results);
			counter++;
		}
		return results;
	}

	double errorGS(std::vector< std::vector<double> > &a, std::vector<double> &b, std::vector<double> &No)
	{
		/* Calculate a*No-b */
		int n = static_cast<int>(b.size());
		double rVal(0.0), error(1.0);
		std::vector<double> R(n, 0.0);

		for (int i = 0; i < n; i++)
		{
			for (int j = 0; j < n; j++)
			{
				R[i] += a[i][j] * No[j];
			}
			R[i] = fabs(R[i] - b[i]);
		}
		error = *std::max_element(R.begin(), R.end());  //find max value of vector
		return error;
	}

	double CalcVisCoef(double T)
	{
		double mu(0.0);
		mu = material::As*pow(T, 1.5) / (T + material::Ts);
		return mu;
	}

	double CalcT(int elem, double a, double b)
	{
		double T(0.0), rhoGs(0.0), rhouGs(0.0), rhovGs(0.0), rhoEGs(0.0);
		rhoGs = math::pointValue(elem, a, b, 1, 2);
		rhouGs = math::pointValue(elem, a, b, 2, 2);
		rhovGs = math::pointValue(elem, a, b, 3, 2);
		rhoEGs = math::pointValue(elem, a, b, 4, 2);
		T = (material::gamma - 1)*(rhoEGs - 0.5*(pow(rhouGs, 2) + pow(rhovGs, 2)) / rhoGs) / (material::R*rhoGs);
		return T;
	}

	double CalcTFromConsvVar(double rho, double rhou, double rhov, double rhoE)
	{
		double T(0.0);
		return (T = (material::gamma - 1)*(rhoE - 0.5*(pow(rhou, 2) + pow(rhov, 2)) / rho) / (material::R*rho));
	}

	double CalcP(double T, double rho)
	{
		double p = T * (material::R*rho);
		return p;
	}

	std::tuple<double, double, double, double> Calc_dxydab(int elem, double a, double b)
	{
		double dxa(0.0), dxb(0.0), dya(0.0), dyb(0.0);
		double xA(0.0), xB(0.0), xC(0.0), xD(0.0), yA(0.0), yB(0.0), yC(0.0), yD(0.0);
		int elemType(0);

		elemType = auxUlti::checkType(elem);
		std::tie(xA, yA) = auxUlti::getElemCornerCoord(elem, 0);
		std::tie(xB, yB) = auxUlti::getElemCornerCoord(elem, 1);
		std::tie(xC, yC) = auxUlti::getElemCornerCoord(elem, 2);
		
		if (elemType==4)  //Quad
		{
			std::tie(xD, yD) = auxUlti::getElemCornerCoord(elem, 3);
			dxa = 0.25*(-xA + xB - xD + xC) + 0.25*(xA - xB - xD + xC)*b;
			dxb = 0.25*(-xA - xB + xD + xC) + 0.25*(xA - xB - xD + xC)*a;
			dya = 0.25*(-yA + yB - yD + yC) + 0.25*(yA - yB - yD + yC)*b;
			dyb = 0.25*(-yA - yB + yD + xC) + 0.25*(yA - yB - yD + yC)*a;
		}
		else if (elemType == 3)  //Tri
		{
			dxa = 0.25*(-xA + xB) + 0.25*(xA - xB)*b;
			dxb = 0.25*(-xA - xB + 2 * xC) + 0.25*(xA - xB)*a;
			dya = 0.25*(-yA + yB) + 0.25*(yA - yB)*b;
			dyb = 0.25*(-yA - yB + 2 * yC) + 0.25*(yA - yB)*a;
		}
		
		return std::make_tuple(dxa, dxb, dya, dyb);
	}

	double Calc_dBxdBy(int elem, int order, int na, int nb, int opt)
	{
		double dB(0.0);
		if (opt==1)  //x direction
		{
			dB = (1 / meshVar::J2D[elem][na][nb]) * (mathVar::dBaPts[order][na][nb] * meshVar::dyb[elem][na][nb] - mathVar::dBbPts[order][na][nb] * meshVar::dya[elem][na][nb]);
		}
		else if (opt==2)  //y direction
		{
			dB = (1 / meshVar::J2D[elem][na][nb]) * (mathVar::dBbPts[order][na][nb] * meshVar::dxa[elem][na][nb] - mathVar::dBaPts[order][na][nb] * meshVar::dxb[elem][na][nb]);
		}
		return dB;
	}

	double surfaceInte(std::vector<double> &Fvalue, int edge, int elem)
	{
		double inte(0.0), J(auxUlti::getJ1D(elem, edge)), w(0.0);
		for (int nG = 0; nG <= mathVar::nGauss; nG++)
		{
			w = mathVar::wGauss[nG];
			inte += w * Fvalue[nG] * J;
		}
		return inte;
	}

	double pointValue(int element, double a, double b, int valType, int valKind)
	{
		//Compute primary variables from conservative variables
		double out(0.0);
		if (valKind==1)  //primary variables
		{
			if (valType == 1)  //rho
			{
				out= limiter::calcConsvVarWthLimiter(element, a, b, valType);
			}
			else if (valType == 2)  //u
			{
				double rhoVal(limiter::calcConsvVarWthLimiter(element, a, b, 1)),
					rhouVal(limiter::calcConsvVarWthLimiter(element, a, b, 2));
				out = rhouVal / rhoVal;
			}
			else if (valType == 3)  //v
			{
				double rhoVal(limiter::calcConsvVarWthLimiter(element, a, b, 1)),
					rhovVal(limiter::calcConsvVarWthLimiter(element, a, b, 3));
				out = rhovVal / rhoVal;
			}
			else if (valType == 4)  //e
			{
				double rhoVal(limiter::calcConsvVarWthLimiter(element, a, b, 1)),
					rhouVal(limiter::calcConsvVarWthLimiter(element, a, b, 2)),
					rhovVal(limiter::calcConsvVarWthLimiter(element, a, b, 3)),
					rhoEVal(limiter::calcConsvVarWthLimiter(element, a, b, 4));
				out = material::Cv*math::CalcTFromConsvVar(rhoVal, rhouVal, rhovVal, rhoEVal);
			}
			else if (valType == 5)  //p
			{
				double rhoVal(limiter::calcConsvVarWthLimiter(element, a, b, 1)),
					rhouVal(limiter::calcConsvVarWthLimiter(element, a, b, 2)),
					rhovVal(limiter::calcConsvVarWthLimiter(element, a, b, 3)),
					rhoEVal(limiter::calcConsvVarWthLimiter(element, a, b, 4));
				double TVal(math::CalcTFromConsvVar(rhoVal, rhouVal, rhovVal, rhoEVal));
				out = math::CalcP(TVal, rhoVal);
			}
			else if (valType == 6)  //T
			{
				double rhoVal(limiter::calcConsvVarWthLimiter(element, a, b, 1)),
					rhouVal(limiter::calcConsvVarWthLimiter(element, a, b, 2)),
					rhovVal(limiter::calcConsvVarWthLimiter(element, a, b, 3)),
					rhoEVal(limiter::calcConsvVarWthLimiter(element, a, b, 4));
				out = math::CalcTFromConsvVar(rhoVal, rhouVal, rhovVal, rhoEVal);
			}
			else if (valType == 7)  //mu
			{
				double rhoVal(limiter::calcConsvVarWthLimiter(element, a, b, 1)),
					rhouVal(limiter::calcConsvVarWthLimiter(element, a, b, 2)),
					rhovVal(limiter::calcConsvVarWthLimiter(element, a, b, 3)),
					rhoEVal(limiter::calcConsvVarWthLimiter(element, a, b, 4));
				double TVal(math::CalcTFromConsvVar(rhoVal, rhouVal, rhovVal, rhoEVal));
				out = math::CalcVisCoef(TVal);
			}
		}
		else if (valKind==2)  //conservative variables
		{	
			out = limiter::calcConsvVarWthLimiter(element, a, b, valType);
		}
		return out;
	}

	double pointValueNoLimiter(int element, double a, double b, int valType, int valKind)
	{
		double out(0.0);
		std::vector<double> Value(mathVar::orderElem + 1, 0.0);

		if (valKind == 1)  //primary variables
		{
			Value = auxUlti::getElementPriValuesOfOrder(element, valType);
		}
		else if (valKind == 2)  //conservative variables
		{
			Value = auxUlti::getElementConserValuesOfOrder(element, valType);
		}

		math::basisFc(a, b);
		for (int order = 0; order <= mathVar::orderElem; order++)
		{
			out += Value[order] * mathVar::B[order];
		}
		return out;
	}

	double vectorDotProduct(std::vector<double> &a, std::vector<double> &b)
	{
		double out(0.0);
		out = a[0] * b[0] + a[1] * b[1];
		return out;
	}

	std::tuple <double, double> internalSurfaceValue(int edge, int element, int nG, int valType, int valKind)
	{
		int masterElem(0), servantElem(0);
		std::tie(masterElem, servantElem) = auxUlti::getMasterServantOfEdge(edge);
		double valPlus(0.0), valMinus(0.0), aMaster(0.0), bMaster(0.0), aServant(0.0), bServant(0.0);

		std::tie(aMaster, bMaster) = auxUlti::getGaussSurfCoor(edge, masterElem, nG);
		std::tie(aServant, bServant) = auxUlti::getGaussSurfCoor(edge, servantElem, nG);

		if (masterElem==element)  //considering element is master
		{
			valPlus = math::pointValue(masterElem, aMaster, bMaster, valType, valKind);
			valMinus = math::pointValue(servantElem, aServant, bServant, valType, valKind);
		}
		else
		{
			valMinus = math::pointValue(masterElem, aMaster, bMaster, valType, valKind);
			valPlus = math::pointValue(servantElem, aServant, bServant, valType, valKind);
		}

		return std::make_tuple(valPlus, valMinus);
	}

	std::tuple <double, double> internalSurfaceDerivativeValue(int edge, int element, int nG, int valType, int dir)
	{
		int masterElem(0), servantElem(0);
		std::tie(masterElem, servantElem) = auxUlti::getMasterServantOfEdge(edge);
		double valPlus(0.0), valMinus(0.0), aMaster(0.0), bMaster(0.0), aServant(0.0), bServant(0.0);

		std::tie(aMaster, bMaster) = auxUlti::getGaussSurfCoor(edge, masterElem, nG);
		std::tie(aServant, bServant) = auxUlti::getGaussSurfCoor(edge, servantElem, nG);

		if (masterElem == element)  //considering element is master
		{
			valPlus = math::pointAuxValue(masterElem, aMaster, bMaster, valType, dir);
			valMinus = math::pointAuxValue(servantElem, aServant, bServant, valType, dir);
		}
		else
		{
			valMinus = math::pointAuxValue(masterElem, aMaster, bMaster, valType, dir);
			valPlus = math::pointAuxValue(servantElem, aServant, bServant, valType, dir);
		}

		return std::make_tuple(valPlus, valMinus);
	}

	double SurfaceValueFromMaster(int edge, int nG, int valType, int valKind)
	{
		int masterElem(0), servantElem(0);
		std::tie(masterElem, servantElem) = auxUlti::getMasterServantOfEdge(edge);
		double valPlus(0.0), aMaster(0.0), bMaster(0.0);
		std::tie(aMaster, bMaster) = auxUlti::getGaussSurfCoor(edge, masterElem, nG);
		valPlus = math::pointValue(masterElem, aMaster, bMaster, valType, valKind);
		return valPlus;
	}

	double CalcSpeedOfSound(double T)
	{
		double C(sqrt(material::gamma*material::R*T));
		return C;
	}

	std::vector<double> vectorSum(std::vector<double> &a, std::vector<double> &b)
	{
		static int size(a.size());
		std::vector<double> out(size, 0.0);
		for (int i = 0; i < size; i++)
		{
			out[i] = a[i] + b[i];
		}
		return out;
	}

	double calcRhouvEDeriv(double dRhoUVal, double dRhoVal, double rhoUVal, double rhoVal)
	{
		/*NOTE! direction of derivatives based on input values, for example:
		d(u)/dx = [d(rhou)/dx - d(rho)/dx.(rhou)/rho]/(rho)
		so input values are
		d(rho)/dx	----> dRhoVal
		d(rhou)/dx	----> dRhoUVal
		rhou		----> rhoUVal
		rho			----> rhoVal*/
		double outVal(0.0);
		outVal = (dRhoUVal - dRhoVal * rhoUVal / rhoVal) / rhoVal;
		return outVal;
	}

	double calcTDeriv(double dE, double du, double dv, double u, double v)
	{
		/*NOTE! direction of derivatives based on input values, for example:
		d(T)/dx = f(dE/dx, du/dx, dv/dx, u, v)
		so input values are
		dE/dx	----> dE
		du/dx	----> du
		dv/dx	----> dv*/
		double outVal(0.0);
		outVal = (material::gamma - 1)*(dE - (u*du + v * dv)) / material::R;
		return outVal;
	}

	double pointAuxValue(int element, double a, double b, int valType, int dir)
	{
		double out(0.0);
		std::vector<double> Value(mathVar::orderElem + 1, 0.0);

		Value = auxUlti::getElementAuxValuesOfOrder(element, valType, dir);
		double muVal(math::pointValue(element, a, b, 7, 1));

		math::basisFc(a, b);
		for (int order = 0; order <= mathVar::orderElem; order++)
		{
			out += Value[order] * mathVar::B[order];
		}
		out = out / muVal;
		return out;
	}

	double calcThermalConductivity(double muVal)
	{
		double k(0.0);
		k = material::Cp*muVal / material::Pr;
		return k;
	}

	std::tuple<bool, double, double> solvQuadraticEq(double A, double B, double C)
	{
		double delta(B*B - 4 * A*C), out(0.0), root1(0.0), root2(0.0);
		bool realRoot(true);

		if (A != 0)
		{
			if (delta>0)
			{
				root1 = ((-B + sqrt(delta)) / (2 * A));
				root2 = ((-B - sqrt(delta)) / (2 * A));
			}
			else if (delta == 0)
			{
				root1 = (-B / (2 * A));
				root2 = root1;
			}
			else
			{
				realRoot = false;
			}
		}
		else if (A == 0.0 && B !=0)
		{
			root1 = -C / B;
			root2 = root1;
		}
		else
		{
			realRoot = false;
		}
		return std::make_tuple(realRoot, root1, root2);
	}

	double centerValue(int element, int valType, int valKind)
	{
		double xC(-1.0 / 3.0), yC(1.0 / 3.0), output(0.0);

		if (auxUlti::checkType(element) == 4)
		{
			xC = 0.0;
			yC = 0.0;
		}
		output = math::pointValue(element, xC, yC, valType, valKind);
		return output;
	}

	double centerAuxValue(int element, int valType, int dir)
	{
		double xC(-1.0 / 3.0), yC(1.0 / 3.0), output(0.0);

		if (auxUlti::checkType(element) == 4)
		{
			xC = 0.0;
			yC = 0.0;
		}
		output = math::pointAuxValue(element, xC, yC, valType, dir);
		return output;
	}

	std::tuple<double, double> mappingStdToReal(int element, double aCoor, double bCoor)
	{
		int elemType(auxUlti::checkType(element));
		double xCoor(0.0), yCoor(0.0), xA(0.0), xB(0.0), xC(0.0),
			yA(0.0), yB(0.0), yC(0.0);

		std::tie(xA, yA) = auxUlti::getElemCornerCoord(element, 0);
		std::tie(xB, yB) = auxUlti::getElemCornerCoord(element, 1);
		std::tie(xC, yC) = auxUlti::getElemCornerCoord(element, 2);

		if (elemType == 3) //Tri element
		{
			xCoor = 0.25*(1 - aCoor)*(1 - bCoor)*xA + 0.25*(1 + aCoor)*(1 - bCoor)*xB + 0.5*(1 + bCoor)*xC;
			yCoor = 0.25*(1 - aCoor)*(1 - bCoor)*yA + 0.25*(1 + aCoor)*(1 - bCoor)*yB + 0.5*(1 + bCoor)*yC;
		}
		else if (elemType == 4) //Quad element
		{
			double xD(0.0), yD(0.0);
			std::tie(xD, yD) = auxUlti::getElemCornerCoord(element, 3);

			xCoor = 0.25*(1 - aCoor)*(1 - bCoor)*xA + 0.25*(1 + aCoor)*(1 - bCoor)*xB + 0.25*(1 - aCoor)*(1 + bCoor)*xD + 0.25*(1 + aCoor)*(1 + bCoor)*xC;
			yCoor = 0.25*(1 - aCoor)*(1 - bCoor)*yA + 0.25*(1 + aCoor)*(1 - bCoor)*yB + 0.25*(1 - aCoor)*(1 + bCoor)*yD + 0.25*(1 + aCoor)*(1 + bCoor)*yC;
		}
		return std::make_tuple(xCoor, yCoor);
	}

	std::tuple<double, double> mappingRealToStd(int edge, int element, double xCoor, double yCoor)
	{
		double A1(0.0), B1(0.0), C1(0.0), D1(0.0),
			A2(0.0), B2(0.0), C2(0.0), D2(0.0),
			xA(0.0), xB(0.0), xC(0.0),
			yA(0.0), yB(0.0), yC(0.0),
			aCoor(0.0), bCoor(0.0),
			xCal(0.0), yCal(0.0), eX(100.0), eY(100.0);
		int elemType(auxUlti::checkType(element));
		std::vector<std::vector<double>> vectorGaussPoints(mathVar::nGauss + 1, std::vector<double>(2, 0.0));

		std::tie(xA, yA) = auxUlti::getElemCornerCoord(element, 0);
		std::tie(xB, yB) = auxUlti::getElemCornerCoord(element, 1);
		std::tie(xC, yC) = auxUlti::getElemCornerCoord(element, 2);

		vectorGaussPoints = auxUlti::getVectorGaussSurfCoor(edge, element);

		if (elemType==3) //tri
		{
			for (int nG = 0; nG <= mathVar::nGauss; nG++)
			{
				aCoor = vectorGaussPoints[nG][0];
				bCoor = vectorGaussPoints[nG][1];
				xCal = 0.25*(1 - aCoor)*(1 - bCoor)*xA + 0.25*(1 + aCoor)*(1 - bCoor)*xB + 0.5*(1 + bCoor)*xC;
				yCal = 0.25*(1 - aCoor)*(1 - bCoor)*yA + 0.25*(1 + aCoor)*(1 - bCoor)*yB + 0.5*(1 + bCoor)*yC;
				eX = fabs(xCal - xCoor) * 100 / xCoor;
				eY = fabs(yCal - yCoor) * 100 / yCoor;
				if ((eX < 0.05) && (eY < 0.05))
				{
					break;
				}
			}
		}
		else if (elemType==4) //quad
		{
			double xD(0.0), yD(0.0);
			std::tie(xD, yD) = auxUlti::getElemCornerCoord(element, 3);

			for (int nG = 0; nG <= mathVar::nGauss; nG++)
			{
				aCoor = vectorGaussPoints[nG][0];
				bCoor = vectorGaussPoints[nG][1];
				xCal = 0.25*(1 - aCoor)*(1 - bCoor)*xA + 0.25*(1 + aCoor)*(1 - bCoor)*xB + 0.25*(1 - aCoor)*(1 + bCoor)*xD + 0.25*(1 + aCoor)*(1 + bCoor)*xC;
				yCal = 0.25*(1 - aCoor)*(1 - bCoor)*yA + 0.25*(1 + aCoor)*(1 - bCoor)*yB + 0.25*(1 - aCoor)*(1 + bCoor)*yD + 0.25*(1 + aCoor)*(1 + bCoor)*yC;
				eX = fabs((xCal - xCoor) * 100 / xCoor);
				eY = fabs((yCal - yCoor) * 100 / yCoor);
				if ((eX < 0.05) && (eY < 0.05))
				{
					break;
				}
			}
		}

		return std::make_tuple(aCoor, bCoor);
	}

	//Function supports for math::mappingRealToStd
	/*
	double solve_abQuad(int option, double A, double B, double D, double C, double inVar)
	{
		double outVar(0.0);
		if (option==1) //a
		{
			outVar = (inVar - (A + B + D + C)) / (-A + B - D + C);
		}
		else if (option==2) //b
		{
			outVar = (inVar - (A + B + D + C)) / (-A - B + D + C);
		}
		return outVar;
	}

	//Function supports for math::mappingRealToStd
	double solve_abTri(int option, double A, double B, double C, double inVar)
	{
		double outVar(0.0);
		if (option == 1) //a
		{
			outVar = (inVar - (A + B + C)) / (-A + B);
		}
		else if (option == 2) //b
		{
			outVar = (inVar - (A + B + C)) / (-A - B + C);
		}
		return outVar;
	}*/

	namespace numericalFluxes
	{
		double auxFlux(double MinusVal, double PlusVar, double vectorComp)
		{
			/*use central numerical flux*/
			double flux(0.5*(MinusVal + PlusVar)*vectorComp);
			return flux;
		}

		double advectiveFlux(double FPlus, double FMinus, double UPlus, double UMinus, double C, double vectorComp)
		{
			/*use Lax - Friedrich numerical flux*/
			double flux(0.5*vectorComp*(FPlus + FMinus - C * (UPlus - UMinus)));
			return flux;
		}

		double diffusiveFlux(double MinusVal, double PlusVar, double vectorComp)
		{
			/*use central numerical flux*/
			double flux(0.5*(MinusVal + PlusVar)*vectorComp);
			return flux;
		}

		std::vector<std::vector<double>> NSFEqFluxFromConserVars(std::vector<double> &UPlus, std::vector<double> &UMinus, std::vector<double> &dUXPlus, std::vector<double> &dUXMinus, std::vector<double> &dUYPlus, std::vector<double> &dUYMinus, std::vector<double> &normVector)
		{
			/*Fluxes array has the following form:
			- column 0: advective fluxes
			- column 1: diffusive fluxes*/
			std::vector<std::vector<double>> Fluxes(4, std::vector<double>(2, 0.0));

			/*StressHeat matrix has form:
			[tauXx		tauXy		Qx]
			[tauYx		tauYy		Qy]
			*/
			std::vector<std::vector<double>> StressHeatP(2, std::vector<double>(3, 0.0));
			std::vector<std::vector<double>> StressHeatM(2, std::vector<double>(3, 0.0));

			double rhoPlus(UPlus[0]), rhouPlus(UPlus[1]), rhovPlus(UPlus[2]), rhoEPlus(UPlus[3]),
				rhoMinus(UMinus[0]), rhouMinus(UMinus[1]), rhovMinus(UMinus[2]), rhoEMinus(UMinus[3]),
				nx(normVector[0]), ny(normVector[1]);

			double
				uPlus(rhouPlus / rhoPlus),
				uMinus(rhouMinus / rhoMinus),

				vPlus(rhovPlus / rhoPlus),
				vMinus(rhovMinus / rhoMinus),

				totalEPlus(rhoEPlus / rhoPlus),
				totalEMinus(rhoEMinus / rhoMinus),

				TPlus(0.0),
				TMinus(0.0),

				pPlus(0.0),
				pMinus(0.0),
				
				muPlus(0.0),
				muMinus(0.0);

			double
				termX1P(0.0), termX1M(0.0),  //(rho*u)					or 0
				termX2P(0.0), termX2M(0.0),  //(rho*u^2 + p)			or tauxx
				termX3P(0.0), termX3M(0.0),  //(rho*u*v)				or tauxy
				termX4P(0.0), termX4M(0.0),  //(rho*totalE + p)*u		or tauxx*u + tauxy*v + Qx

				termY1P(0.0), termY1M(0.0),  //(rho*v)					or 0
				termY2P(0.0), termY2M(0.0),  //(rho*u*v)				or tauxy
				termY3P(0.0), termY3M(0.0),  //(rho*v^2 + p)			or tauyy
				termY4P(0.0), termY4M(0.0);  //(rho*totalE + p)*v		or tauxy*u + tauyy*v + Qy

			double C(0.0),
				uMagP(0.0),
				uMagM(0.0),
				aP(0.0),
				aM(0.0);

			/*INVISCID TERMS*/
			/*calculate velocity magnitude*/
			uMagP = sqrt(pow(uPlus, 2) + pow(vPlus, 2));
			uMagM = sqrt(pow(uMinus, 2) + pow(vMinus, 2));

			/*calculate T and P*/
			TPlus = math::CalcTFromConsvVar(rhoPlus, rhouPlus, rhovPlus, rhoEPlus);
			TMinus = math::CalcTFromConsvVar(rhoMinus, rhouMinus, rhovMinus, rhoEMinus);
			pPlus = math::CalcP(TPlus, rhoPlus);
			pMinus = math::CalcP(TMinus, rhoMinus);
			muPlus = math::CalcVisCoef(TPlus);
			muMinus = math::CalcVisCoef(TMinus);

			/*calculate speed of sound*/
			aP = math::CalcSpeedOfSound(TPlus);
			aM = math::CalcSpeedOfSound(TMinus);

			/*calculate constant for Lax-Friederich flux*/
			C = math::numericalFluxes::constantC(uMagP, uMagM, aP, aM);

			/*calculate inviscid terms*/
			std::tie(termX1P, termX2P, termX3P, termX4P) = math::inviscidTerms::calcInvisTermsFromPriVars(rhoPlus, uPlus, vPlus, totalEPlus, pPlus, 1);
			std::tie(termY1P, termY2P, termY3P, termY4P) = math::inviscidTerms::calcInvisTermsFromPriVars(rhoPlus, uPlus, vPlus, totalEPlus, pPlus, 2);

			std::tie(termX1M, termX2M, termX3M, termX4M) = math::inviscidTerms::calcInvisTermsFromPriVars(rhoMinus, uMinus, vMinus, totalEMinus, pMinus, 1);
			std::tie(termY1M, termY2M, termY3M, termY4M) = math::inviscidTerms::calcInvisTermsFromPriVars(rhoMinus, uMinus, vMinus, totalEMinus, pMinus, 2);

			/*Calculate fluxes*/
			Fluxes[0][0] = math::numericalFluxes::advectiveFlux(termX1P, termX1M, rhoPlus, rhoMinus, C, nx) + math::numericalFluxes::advectiveFlux(termY1P, termY1M, rhoPlus, rhoMinus, C, ny);
			Fluxes[1][0] = math::numericalFluxes::advectiveFlux(termX2P, termX2M, rhouPlus, rhouMinus, C, nx) + math::numericalFluxes::advectiveFlux(termY2P, termY2M, rhouPlus, rhouMinus, C, ny);
			Fluxes[2][0] = math::numericalFluxes::advectiveFlux(termX3P, termX3M, rhovPlus, rhovMinus, C, nx) + math::numericalFluxes::advectiveFlux(termY3P, termY3M, rhovPlus, rhovMinus, C, ny);
			Fluxes[3][0] = math::numericalFluxes::advectiveFlux(termX4P, termX4M, rhoEPlus, rhoEMinus, C, nx) + math::numericalFluxes::advectiveFlux(termY4P, termY4M, rhoEPlus, rhoEMinus, C, ny);
			
			/*VISCOUS TERMS*/
			/*calculate inviscid terms*/
			StressHeatP = math::viscousTerms::calcStressTensorAndHeatFlux(muPlus, UPlus, dUXPlus, dUYPlus);
			std::tie(termX1P, termX2P, termX3P, termX4P) = math::viscousTerms::calcViscousTermsFromStressHeatFluxMatrix(StressHeatP, uPlus, vPlus, 1);
			std::tie(termY1P, termY2P, termY3P, termY4P) = math::viscousTerms::calcViscousTermsFromStressHeatFluxMatrix(StressHeatP, uPlus, vPlus, 2);

			StressHeatM = math::viscousTerms::calcStressTensorAndHeatFlux(muMinus, UMinus, dUXMinus, dUYMinus);
			std::tie(termX1M, termX2M, termX3M, termX4M) = math::viscousTerms::calcViscousTermsFromStressHeatFluxMatrix(StressHeatM, uMinus, vMinus, 1);
			std::tie(termY1M, termY2M, termY3M, termY4M) = math::viscousTerms::calcViscousTermsFromStressHeatFluxMatrix(StressHeatM, uMinus, vMinus, 2);

			/*Calculate fluxes*/
			Fluxes[0][1] = math::numericalFluxes::diffusiveFlux(termX1M, termX1P, nx) + math::numericalFluxes::diffusiveFlux(termY1M, termY1P, ny);
			Fluxes[1][1] = math::numericalFluxes::diffusiveFlux(termX2M, termX2P, nx) + math::numericalFluxes::diffusiveFlux(termY2M, termY2P, ny);
			Fluxes[2][1] = math::numericalFluxes::diffusiveFlux(termX3M, termX3P, nx) + math::numericalFluxes::diffusiveFlux(termY3M, termY3P, ny);
			Fluxes[0][1] = math::numericalFluxes::diffusiveFlux(termX4M, termX4P, nx) + math::numericalFluxes::diffusiveFlux(termY4M, termY4P, ny);

			return Fluxes;
		}

		double constantC(double uMagP, double uMagM, double aP, double aM)
		{
			std::vector<double> CArray(2, 0.0);
			double C(0.0);
			CArray[0] = uMagP + aP / refValues::Ma;
			CArray[1] = uMagM + aM / refValues::Ma;
			C = *std::max_element(CArray.begin(), CArray.end());  //find max value of vector
			return C;
		}
	}//end of namespace numericalFluxes

	namespace inviscidTerms
	{
		std::tuple<double, double, double, double> calcInvisTermsFromPriVars(double rhoVal, double uVal, double vVal, double totalE, double pVal, int dir)
		{
			double term1(0.0), term2(0.0), term3(0.0), term4(0.0);

			if (dir==1)  //Ox direction
			{ 
				term1 = rhoVal * uVal;
				term2 = rhoVal * pow(uVal, 2) + pVal;
				term3 = rhoVal * uVal*vVal;
				term4 = (rhoVal*totalE + pVal)*uVal;
			}
			else if (dir==2)  //Oy direction
			{
				term1 = rhoVal * vVal;
				term2 = rhoVal * uVal*vVal;
				term3 = rhoVal * pow(vVal, 2) + pVal;
				term4 = (rhoVal*totalE + pVal)*vVal;
			}
			return std::make_tuple(term1, term2, term3, term4);
		}
	}//end of namespace invicidTerms

	namespace viscousTerms
	{
		std::vector<std::vector<double>> calcStressTensorAndHeatFlux(double muVal, std::vector<double> &U, std::vector<double> &dUx, std::vector<double> &dUy)
		{
			/*Output matrix has form:
			[tauXx		tauXy		Qx]
			[tauYx		tauYy		Qy]
			*/
			std::vector<std::vector<double>> OutputMatrix(2, std::vector<double>(3, 0.0));

			double dux(0.0), duy(0.0), dvx(0.0), dvy(0.0), dEx(0.0), dEy(0.0), rhouVal(0.0), rhovVal, rhoVal(0.0), rhoEVal(0.0);
			double drhox(0.0), drhoy(0.0),
				drhoux(0.0), drhouy(0.0),
				drhovx(0.0), drhovy(0.0),
				drhoEx(0.0), drhoEy(0.0),
				dTx(0.0), dTy(0.0);
			double uVal(0.0), vVal(0.0), Qx(0.0), Qy(0.0), k(0.0);
			int index(0);

			rhoVal = U[0];
			rhouVal = U[1];
			rhovVal = U[2];
			rhoEVal = U[3];

			uVal = rhouVal / rhoVal;
			vVal = rhovVal / rhoVal;

			drhox = dUx[0];
			drhoy = dUy[0];

			drhoux = dUx[1];
			drhouy = dUy[1];

			drhovx = dUx[2];
			drhovy = dUy[2];

			drhoEx = dUx[3];
			drhoEy = dUy[3];

			dux = math::calcRhouvEDeriv(drhoux, drhox, rhouVal, rhoVal);
			duy = math::calcRhouvEDeriv(drhouy, drhoy, rhouVal, rhoVal);

			dvx = math::calcRhouvEDeriv(drhovx, drhox, rhovVal, rhoVal);
			dvy = math::calcRhouvEDeriv(drhovy, drhoy, rhovVal, rhoVal);

			dEx = math::calcRhouvEDeriv(drhoEx, drhox, rhoEVal, rhoVal);
			dEy = math::calcRhouvEDeriv(drhoEy, drhoy, rhoEVal, rhoVal);

			/*calculate stresses*/
			for (int i = 0; i < 2; i++)
			{
				for (int j = 0; j < 2; j++)
				{
					index = 10 * (i + 1) + (j + 1);
					if (index == 11)  //x normal stress (tau_xx)
					{
						OutputMatrix[i][j] = math::viscousTerms::calcStressComponent(index, muVal, dux, dvy);
					}
					else if (index == 22)  //y normal stress (tau_yy)
					{
						OutputMatrix[i][j] = math::viscousTerms::calcStressComponent(index, muVal, dvy, dux);
					}
					else  //shear stress (tau_xy)
					{
						OutputMatrix[i][j] = math::viscousTerms::calcStressComponent(index, muVal, duy, dvx);
					}
				}
			}

			/*calculate heat flux*/
			dTx = math::calcTDeriv(dEx, dux, dvx, uVal, vVal);
			dTy = math::calcTDeriv(dEy, duy, dvy, uVal, vVal);
			k = math::calcThermalConductivity(muVal);
			std::tie(Qx, Qy) = math::viscousTerms::calcHeatFluxTerms(dTx, dTy, k);

			OutputMatrix[0][2] = Qx;
			OutputMatrix[1][2] = Qy;
			return OutputMatrix;
		}

		double calcStressComponent(int index, double muVal, double fstDeriv, double sndDeriv)
		{
			/*Formular of stress
			tau_xx = -mu((4/3)*du/dx - (2/3)*dv/dy)
			tau_yy = -mu((4/3)*dv/dy - (2/3)*du/dx)
			tau_xy = -mu(du/dy + dv/dx)

			input variables:
			index			fstDeriv		sndDeriv
			11 (tau_xx)		du/dx			dv/dy
			22 (tau_yy)		dv/dy			du/dx
			12 (tau_xy)		du/dy			dv/dx
			*/
			double tau(0.0);
			if (index == 11 || index == 22)
			{
				tau = -muVal * ((4.0 / 3.0)*fstDeriv - (2.0 / 3.0)*sndDeriv);
			}
			else if (index == 12 || index == 21)
			{
				tau = -muVal * (fstDeriv + sndDeriv);
			}
			return tau;
		}

		std::tuple<double, double> calcHeatFluxTerms(double dTx, double dTy, double k)
		{
			double Qx(0.0), Qy(0.0);
			Qx = -k * dTx;
			Qy = -k * dTy;
			return std::make_tuple(Qx, Qy);
		}

		std::tuple<double, double, double, double> calcViscousTermsFromStressHeatFluxMatrix(std::vector< std::vector<double> > &StressHeatFlux, double uVal, double vVal, int dir)
		{
			/*StressHeatFlux is 2D array contents stress and heat flux component, has the following form:
			[tauXx		tauXy		Qx]
			[tauYx		tauYy		Qy]*/

			double tauXy(0.0), tauXx(0.0), tauYy(0.0), Qx(0.0), Qy(0.0);
			double viscTerm1(0.0), viscTerm2(0.0), viscTerm3(0.0), viscTerm4(0.0);

			tauXx = StressHeatFlux[0][0];
			tauYy = StressHeatFlux[1][1];
			tauXy = StressHeatFlux[0][1];
			Qx = StressHeatFlux[0][2];
			Qy = StressHeatFlux[1][2];

			if (dir==1)
			{
				/*1. Ox direction*/
				viscTerm1 = 0.0;
				viscTerm2 = tauXx;
				viscTerm3 = tauXy;
				viscTerm4 = tauXx * uVal + tauXy * vVal + Qx;
			}
			else if (dir==2)
			{
				/*2. Oy direction*/
				viscTerm1 = 0.0;
				viscTerm2 = tauXy;
				viscTerm3 = tauYy;
				viscTerm4 = tauXy * uVal + tauYy * vVal + Qy;
			}
			return std::make_tuple(viscTerm1, viscTerm2, viscTerm3, viscTerm4);
		}
	}//end of namespace viscousTerms

	namespace limiter
	{
		//Function calculates mean value of conservative variables of quad element
		double calcMeanConsvVarQuad(int element, int valType)
		{
			double meanVal(0.0);

			meanVal += math::pointValueNoLimiter(element, -1.0, -1.0, valType, 2);
			meanVal += math::pointValueNoLimiter(element, -1.0, 1.0, valType, 2);
			meanVal += math::pointValueNoLimiter(element, 1.0, -1.0, valType, 2);
			meanVal += math::pointValueNoLimiter(element, 1.0, 1.0, valType, 2);
			meanVal = meanVal / 4.0;
			return meanVal;
		}

		//Function calculates mean value of conservative variables of tri element
		double calcMeanConsvVarTri(int element, int valType)
		{
			double meanVal(0.0);

			meanVal += math::pointValueNoLimiter(element, -1.0, -1.0, valType, 2);
			meanVal += math::pointValueNoLimiter(element, -1.0, 1.0, valType, 2);
			meanVal += math::pointValueNoLimiter(element, 1.0, -1.0, valType, 2);
			meanVal = meanVal / 3.0;
			return meanVal;
		}

		//Function calculates minimum value of rho of quad element
		double calcMinRhoQuad(int element)
		{
			std::vector<double> vectorRho(2 * (mathVar::nGauss + 1) * (mathVar::nGauss + 1), 0.0);
			double aG(0.0), bG(0.0), aGL(0.0), bGL(0.0), min(0.0);
			int index(0);
			for (int na = 0; na <= mathVar::nGauss; na++)
			{
				for (int nb = 0; nb <= mathVar::nGauss; nb++)
				{
					aG = mathVar::GaussPts[na][nb][0];
					bG = mathVar::GaussPts[na][nb][1];

					aGL = mathVar::GaussLobattoPts[na][nb][0];
					bGL = mathVar::GaussLobattoPts[na][nb][1];

					vectorRho[index] = math::pointValueNoLimiter(element, aG, bGL, 1, 1);
					index++;
					vectorRho[index] = math::pointValueNoLimiter(element, aGL, bG, 1, 1);
					index++;
				}
			}
			min = *std::min_element(vectorRho.begin(), vectorRho.end());  //find min value of vector
			return min;
		}

		//Function calculates minimum value of rho of tri element
		double calcMinRhoTri(int element)
		{
			std::vector<double> vectorRho(3, 0.0);
			double minVal(0.0);

			vectorRho[0] = math::pointValueNoLimiter(element, -1.0, -1.0, 1, 2);
			vectorRho[1] = math::pointValueNoLimiter(element, -1.0, 1.0, 1, 2);
			vectorRho[2] = math::pointValueNoLimiter(element, 1.0, -1.0, 1, 2);
			minVal = *std::min_element(vectorRho.begin(), vectorRho.end());
			return minVal;
		}

		//Function calculates minimum value of p of tri element
		double calcMinPTri(int element)
		{
			double rhoVal(0.0), rhouVal(0.0), rhovVal(0.0), rhoEVal(0.0), TVal(0.0);
			std::vector<double> vectorP(3, 0.0);
			double minVal(0.0);

			rhoVal = math::pointValueNoLimiter(element, -1.0, -1.0, 1, 2);
			rhouVal = math::pointValueNoLimiter(element, -1.0, -1.0, 2, 2);
			rhovVal = math::pointValueNoLimiter(element, -1.0, -1.0, 3, 2);
			rhoEVal = math::pointValueNoLimiter(element, -1.0, -1.0, 4, 2);
			TVal = math::CalcTFromConsvVar(rhoVal, rhouVal, rhovVal, rhoEVal);
			vectorP[0] = math::CalcP(TVal, rhoVal);

			rhoVal = math::pointValueNoLimiter(element, -1.0, 1.0, 1, 2);
			rhouVal = math::pointValueNoLimiter(element, -1.0, 1.0, 2, 2);
			rhovVal = math::pointValueNoLimiter(element, -1.0, 1.0, 3, 2);
			rhoEVal = math::pointValueNoLimiter(element, -1.0, 1.0, 4, 2);
			TVal = math::CalcTFromConsvVar(rhoVal, rhouVal, rhovVal, rhoEVal);
			vectorP[1] = math::CalcP(TVal, rhoVal);

			rhoVal = math::pointValueNoLimiter(element, 1.0, -1.0, 1, 2);
			rhouVal = math::pointValueNoLimiter(element, 1.0, -1.0, 2, 2);
			rhovVal = math::pointValueNoLimiter(element, 1.0, -1.0, 3, 2);
			rhoEVal = math::pointValueNoLimiter(element, 1.0, -1.0, 4, 2);
			TVal = math::CalcTFromConsvVar(rhoVal, rhouVal, rhovVal, rhoEVal);
			vectorP[2] = math::CalcP(TVal, rhoVal);

			minVal = *std::min_element(vectorP.begin(), vectorP.end());
			return minVal;
		}

		//Function calculates modified value of Rho at abitrary point (for calculating theta2)
		double calcRhoModified(int element, double a, double b, double theta1, double rhoMean)
		{
			double rhoOrigin(math::pointValueNoLimiter(element, a, b, 1, 1)), rhoMod(0.0);
			rhoMod = theta1*(rhoOrigin - rhoMean) + rhoMean;
			return rhoMod;
		}

		/*Function computes value of conservative variables at abitrary point with applying limiter
		valType:
		1: rho
		2: rhou
		3: rhov
		4: rhoE*/
		double calcConsvVarWthLimiter(int element, double a, double b, int valType)
		{
			double out(0.0), rhoVal(0.0);
			std::vector<double> Value(mathVar::orderElem + 1, 0.0);
			Value = auxUlti::getElementConserValuesOfOrder(element, valType);

			//Compute value at point (a, b) without limiter
			math::basisFc(a, b);
			for (int order = 0; order <= mathVar::orderElem; order++)
			{
				out += Value[order] * mathVar::B[order];
			}

			/*1st approaching:
			In this approaching, we follow instruction shown on Zhang's paper exactly, this means we need to check limiter condition at each time 
			we compute value of conservative variables. So this approaching leads to low computing performance. To improve performace, we use 2nd approaching.
			//Check limiter condition
			bool needLimiter(false);
			std::tie(needLimiter, rhoVal) = math::limiter::checkLimiterForQuad(element, a, b);

			if (needLimiter)
			{
				//Limit value at point (a, b)
				//rho is limited 2 times, other variables are limited 1 time
				if (valType == 1)  //rho
				{
					//Limit rho second time
					out = theta2Arr[element] * (rhoVal - meanVals[element][valType - 1]) + meanVals[element][valType - 1];
				}
				else
				{
					//Limit other value second time
					out = theta2Arr[element] * (out - meanVals[element][valType - 1]) + meanVals[element][valType - 1];
				}
				
			}
			*/

			/*2nd approaching:
			In this approaching, we compute values of theta1 and theta2 of every element, by limiting roots of quadratic equation of coefficient t, value of theta2
			always fall into 0-1 segment. Limiter is applied at all elements, whatever that element is needed to limited or not. because of that, limiting condition
			is not needed to be checked any more*/
			
			//Limit value at point (a, b)
			//rho is limited 2 times, other variables are limited 1 time
			if (valType == 1)  //rho
			{
				out = theta1Arr[element] * (out - meanVals[element][valType - 1]) + meanVals[element][valType - 1];
			}
			out = theta2Arr[element] * (out - meanVals[element][valType - 1]) + meanVals[element][valType - 1];

			return out;
		}

		//Function returns true if element is needed to limit
		std::tuple<bool, double> checkLimiterForQuad(int element, double a, double b)
		{
			double rhoVal(0.0), pVal(0.0), TVal(0.0), rhouVal(0.0), rhovVal(0.0), rhoEVal(0.0);
			bool needLimiter(false);

			rhoVal = math::pointValueNoLimiter(element, a, b, 1, 2);

			//Modify rho
			rhoVal = theta1Arr[element] * (rhoVal - meanVals[element][0]) + meanVals[element][0];

			rhouVal = math::pointValueNoLimiter(element, a, b, 2, 2);
			rhovVal = math::pointValueNoLimiter(element, a, b, 3, 2);
			rhoEVal = math::pointValueNoLimiter(element, a, b, 4, 2);

			TVal = math::CalcTFromConsvVar(rhoVal, rhouVal, rhovVal, rhoEVal);
			pVal = math::CalcP(TVal, rhoVal);

			if (pVal < systemVar::epsilon)
			{
				needLimiter = true;
			}
			return std::make_tuple(needLimiter, rhoVal);
		}

		std::tuple<double, double> calcTheta1Coeff(double meanRho, double minRho, double meanP)
		{
			double temp1(0.0), theta1(0.0);
			std::vector<double> vectorOmega(3, 0.0);
			vectorOmega[0] = systemVar::epsilon;
			vectorOmega[1] = meanP;
			vectorOmega[2] = meanRho;
			double omega(*std::min_element(vectorOmega.begin(), vectorOmega.end()));  //find min value of vector

			temp1 = (meanRho - omega) / (meanRho - minRho);
			if (temp1 < 1.0)
			{
				theta1 = temp1;
			}
			else
			{
				theta1 = 1.0;
			}
			return std::make_tuple(theta1, omega);
		}

		//Function computes theta2 at 1 Gauss point in input direction
		double calcTheta2Coeff(int element, int na, int nb, double theta1, double omega ,double meanRho, double meanRhou, double meanRhov, double meanRhoE, int dir)
		{
			double theta2(0.0);
			//coefficients of t equation
			double A1(0.0), A2(0.0), A3(0.0), A4(0.0), B1(0.0), B2(0.0), B3(0.0), B4(0.0), ACoef(0.0), BCoef(0.0), CCoef(0.0);
			bool realRoot(true);
			double root1(0.0), root2(0.0), rhouOrigin(0.0), rhovOrigin(0.0), rhoEOrigin(0.0), rhoMod(0.0),
				aG(0.0), bG(0.0), aGL(0.0), bGL(0.0), min(0.0);
			
			aG = mathVar::GaussPts[na][nb][0];
			bG = mathVar::GaussPts[na][nb][1];

			aGL = mathVar::GaussLobattoPts[na][nb][0];
			bGL = mathVar::GaussLobattoPts[na][nb][1];

			switch (dir)
			{
			case 1:
				rhoMod = math::limiter::calcRhoModified(element, aG, bGL, theta1, meanRho);
				rhouOrigin = math::pointValueNoLimiter(element, aG, bGL, 2, 2);
				rhovOrigin = math::pointValueNoLimiter(element, aG, bGL, 3, 2);
				rhoEOrigin = math::pointValueNoLimiter(element, aG, bGL, 4, 2);
				break;
			case 2:
				rhoMod = math::limiter::calcRhoModified(element, aGL, bG, theta1, meanRho);
				rhouOrigin = math::pointValueNoLimiter(element, aGL, bG, 2, 2);
				rhovOrigin = math::pointValueNoLimiter(element, aGL, bG, 3, 2);
				rhoEOrigin = math::pointValueNoLimiter(element, aGL, bG, 4, 2);
				break;
			}

			A1 = rhoMod - meanRho;
			A2 = rhouOrigin - meanRhou;
			A3 = rhovOrigin - meanRhov;
			A4 = rhoEOrigin - meanRhoE;

			ACoef = A4 * A1 - 0.5*(A2*A2 + A3 * A3);
			BCoef = A4 * A1 + A1 * rhoEOrigin - A2 * rhouOrigin - A3 * rhovOrigin - omega * A1 / (material::gamma - 1);
			CCoef = rhoEOrigin * rhoMod - 0.5*(pow(rhouOrigin, 2) + pow(rhovOrigin, 2)) - omega * rhoMod / (material::gamma - 1);

			std::tie(realRoot, root1, root2) = math::solvQuadraticEq(ACoef, BCoef, CCoef);
			if (realRoot)
			{
				if ((root1>0.0) & (root1<1.0))
				{
					theta2 = root1;
				}
				else if ((root2>0.0) & (root2<1.0))
				{
					theta2 = root2;
				}
				else
				{
					theta2 = 1.0;
				}
			}
			else
			{
				theta2 = 1.0;
			}
			return theta2;
		}
	}

	namespace geometricOp
	{
		std::tuple<double, double> calcGeoCenter(std::vector<double> &xCoor, std::vector<double> &yCoor, int type)
		{
			double x(0.0), y(0.0), xCG(0.0), yCG(0.0);
			for (int i = 0; i < type; i++)
			{
				x = xCoor[i];
				y = yCoor[i];
				xCG += x;
				yCG += y;
			}

			xCG = xCG / type;
			yCG = yCG / type;
			return std::make_tuple(xCG, yCG);
		}

		double calcPolygonArea(std::vector<double> &xCoor, std::vector<double> &yCoor, int type)
		{
			double x1(xCoor[0]), x2(xCoor[1]), x3(xCoor[2]), x4(0.0), y1(yCoor[0]), y2(yCoor[1]), y3(yCoor[2]), y4(0.0),Area(0.0);
			if (type == 3)
			{
				Area = fabs(0.5*(x1 * y2 - x2 * y1 + x2 * y3 - x3 * y2 + x3 * y1 - x1 * y3));
			}
			else if (type == 4)
			{
				x4 = xCoor[3];
				y4 = yCoor[3];

				Area = fabs(0.5*(x1 * y2 - x2 * y1 + x2 * y3 - x3 * y2 + x3 * y4 - x4 * y3 + x4 * y1 - x1 * y4));
			}
			return Area;
		}

		std::tuple<double, double> calcQuadCentroid(int element, double xCG, double yCG, double area)
		{
			std::vector<double> xSubTriCoor(3, 0.0), ySubTriCoor(3, 0.0);
			std::vector<double> xCGSubTri(4, 0.0), yCGSubTri(4, 0.0), subTriArea(4,0.0);
			double xC(0.0), yC(0.0);
			//1. point 0, 1
			std::tie(xSubTriCoor[0], ySubTriCoor[0]) = auxUlti::getElemCornerCoord(element, 0);
			std::tie(xSubTriCoor[1], ySubTriCoor[1]) = auxUlti::getElemCornerCoord(element, 1);
			xSubTriCoor[2] = xCG;
			ySubTriCoor[2] = yCG;
			std::tie(xCGSubTri[0], yCGSubTri[0]) = math::geometricOp::calcGeoCenter(xSubTriCoor, ySubTriCoor, 3);
			subTriArea[0] = math::geometricOp::calcPolygonArea(xSubTriCoor, ySubTriCoor, 3);

			//2. point 1, 2
			std::tie(xSubTriCoor[0], ySubTriCoor[0]) = auxUlti::getElemCornerCoord(element, 1);
			std::tie(xSubTriCoor[1], ySubTriCoor[1]) = auxUlti::getElemCornerCoord(element, 2);
			xSubTriCoor[2] = xCG;
			ySubTriCoor[2] = yCG;
			std::tie(xCGSubTri[1], yCGSubTri[1]) = math::geometricOp::calcGeoCenter(xSubTriCoor, ySubTriCoor, 3);
			subTriArea[1] = math::geometricOp::calcPolygonArea(xSubTriCoor, ySubTriCoor, 3);

			//3. point 2, 3
			std::tie(xSubTriCoor[0], ySubTriCoor[0]) = auxUlti::getElemCornerCoord(element, 2);
			std::tie(xSubTriCoor[1], ySubTriCoor[1]) = auxUlti::getElemCornerCoord(element, 3);
			xSubTriCoor[2] = xCG;
			ySubTriCoor[2] = yCG;
			std::tie(xCGSubTri[2], yCGSubTri[2]) = math::geometricOp::calcGeoCenter(xSubTriCoor, ySubTriCoor, 3);
			subTriArea[2] = math::geometricOp::calcPolygonArea(xSubTriCoor, ySubTriCoor, 3);

			//4. point 3, 0
			std::tie(xSubTriCoor[0], ySubTriCoor[0]) = auxUlti::getElemCornerCoord(element, 3);
			std::tie(xSubTriCoor[1], ySubTriCoor[1]) = auxUlti::getElemCornerCoord(element, 0);
			xSubTriCoor[2] = xCG;
			ySubTriCoor[2] = yCG;
			std::tie(xCGSubTri[3], yCGSubTri[3]) = math::geometricOp::calcGeoCenter(xSubTriCoor, ySubTriCoor, 3);
			subTriArea[3] = math::geometricOp::calcPolygonArea(xSubTriCoor, ySubTriCoor, 3);

			for (int i = 0; i < 4; i++)
			{
				xC += subTriArea[i] * xCGSubTri[i];
				yC += subTriArea[i] * yCGSubTri[i];
			}
			xC = xC / area;
			yC = yC / area;
			return std::make_tuple(xC, yC);
		}
	}
}//end of namespace math