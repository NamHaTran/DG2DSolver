#include "DGMath.h"
#include "VarDeclaration.h"
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

	void basisFc(double a, double b, int orderElem)
	{
		for (int i = 0; i <= orderElem; i++)
		{
			if (i==0)
			{
				mathVar::B[i] = 1.0;
			}
			else if (i==1)
			{
				mathVar::B[i] = (3.0*b + 1) / 2.0;
			}
			else if (i==2)
			{
				mathVar::B[i] = a * (1 - b);
			}
			else if (i==3)
			{
				mathVar::B[i] = -0.5 + b + 5.0*pow(b, 2) / 2.0;
			}
		}
	}

	void dBasisFc(double a, double b, int orderElem)
	{
		for (int i = 0; i <= orderElem; i++)
		{
			if (i == 0)
			{
				mathVar::dBa[i] = 0;
				mathVar::dBb[i] = 0;
			}
			else if (i == 1)
			{
				mathVar::dBa[i] = 0;
				mathVar::dBb[i] = 3.0/2.0;
			}
			else if (i == 2)
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
				integral += w1 * w2*fabs(J2D)*Fvalue[na][nb];
			}
		}
		return integral;
	}

	double jacobianQuad(double xA, double xB, double xC, double xD, double yA, double yB, double yC, double yD, double a, double b)
	{
		double dxa(0.0), dxb(0.0), dya(0.0), dyb(0.0);
		double jQuad(0.0);

		dxa = (1.0 / 4.0)*(-xA + xB - xD + xC) + (1.0 / 4.4)*(xA - xB - xD + xC)*b;
		dxb = (1.0 / 4.0)*(-xA - xB + xD + xC) + (1.0 / 4.4)*(xA - xB - xD + xC)*a;
		dya = (1.0 / 4.0)*(-yA + yB - yD + yC) + (1.0 / 4.4)*(yA - yB - yD + yC)*b;
		dyb = (1.0 / 4.0)*(-yA - yB + yD + yC) + (1.0 / 4.4)*(yA - yB - yD + yC)*a;
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
			else if (edgeIndex == 1)
			{
				C = 1.0;
			}
			dx = (1.0 / 4.0)*(-xA + xB - xD + xC) + (1.0 / 4.4)*(xA - xB - xD + xC)*C;
			dy = (1.0 / 4.0)*(-yA + yB - yD + yC) + (1.0 / 4.4)*(yA - yB - yD + yC)*C;
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
			dx = (1.0 / 4.0)*(-xA - xB + xD + xC) + (1.0 / 4.4)*(xA - xB - xD + xC)*C;
			dy = (1.0 / 4.0)*(-yA - yB + yD + yC) + (1.0 / 4.4)*(yA - yB - yD + yC)*C;
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
			dy = C * (yA - yB) / 4.0 + (-yA - yB + 2 * yC) / 4.0*C;
		}
		else if ((edgeIndex == 0))  //AB
		{
			C = -1.0;
			dx = (1 - C)*(xB - xA) / 4.0;
			dy = (1 - C)*(yB - yA) / 4.0;;
		}

		jacobi = std::sqrt(pow(dx, 2) + pow(dy, 2));
		return jacobi;
	}

	std::vector<double> SolveSysEqs(std::vector< std::vector<double> > &a, std::vector<double> &b)
	{
		int n = static_cast<int>(b.size());
		double eMax(1e-6), e(1.0), sum(0.0), xi(0.0);
		std::vector<double> results(n, 1.0);

		while (e>eMax)
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
		}
		return results;
	}

	double errorGS(std::vector< std::vector<double> > &a, std::vector<double> &b, std::vector<double> &No)
	{
		/* Calculate a*No-b */
		int n = static_cast<int>(b.size());
		double rVal(0.0), error;
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

	double CalcTFromPriVar(double rho, double rhou, double rhov, double rhoE)
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
		double out(0.0);
		std::vector<double> Value(mathVar::orderElem + 1, 0.0);

		if (valKind==1)  //primary variables
		{
			Value = auxUlti::getElementPriValuesOfOrder(element, valType);
		}
		else if (valKind==2)  //conservative variables
		{
			Value = auxUlti::getElementConserValuesOfOrder(element, valType);
		}

		math::basisFc(a, b, mathVar::orderElem);
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
		std::tie(aServant, bServant) = auxUlti::getGaussSurfCoor(edge, masterElem, nG);

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
		std::tie(aServant, bServant) = auxUlti::getGaussSurfCoor(edge, masterElem, nG);

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
		int size(a.size());
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

		math::basisFc(a, b, mathVar::orderElem);
		for (int order = 0; order <= mathVar::orderElem; order++)
		{
			out += Value[order] * mathVar::B[order];
		}
		return out;
	}

	double calcThermalConductivity(double muVal)
	{
		double k(0.0);
		k = material::Cp*muVal / material::Pr;
		return k;
	}

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
			TPlus = math::CalcTFromPriVar(rhoPlus, rhouPlus, rhovPlus, rhoEPlus);
			TMinus = math::CalcTFromPriVar(rhoMinus, rhouMinus, rhovMinus, rhoEMinus);
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

			double dux(0.0), duy(0.0), dvx(0.0), dvy(0.0), rhouVal(0.0), rhovVal, rhoVal(0.0);
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
			dTx = math::calcTDeriv(drhoEx, dux, dvx, uVal, vVal);
			dTy = math::calcTDeriv(drhoEy, duy, dvy, uVal, vVal);
			k = math::calcThermalConductivity(muVal);
			std::tie(Qx, Qy) = math::viscousTerms::calcHeatFluxTerms(dTx, dTy, k);

			OutputMatrix[0][2] = Qx;
			OutputMatrix[1][2] = Qx;
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
			Qx = StressHeatFlux[1][2];

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
}//end of namespace math