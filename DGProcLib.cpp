#include "DGProcLib.h"
#include "DGMath.h"
#include "varDeclaration.h"
#include "DGAuxUltilitiesLib.h"
#include "dynamicVarDeclaration.h"
#include <tuple>
#include "DGBCsLib.h"
#include <algorithm>
#include "DGIOLib.h"
#include <iostream>

namespace meshParam
{
	void GaussParam()
	{
		math::Gauss(mathVar::nGauss);
		math::GaussLobatto(mathVar::nGauss);  //run GaussLobatto for applying limiter
		for (int na = 0; na <= mathVar::nGauss; na++)  //nGauss is started from 0
		{
			for (int nb = 0; nb <= mathVar::nGauss; nb++)
			{
				mathVar::GaussPts[na][nb][0] = mathVar::xGauss[na];
				mathVar::GaussPts[na][nb][1] = mathVar::xGauss[nb];
				mathVar::GaussLobattoPts[na][nb][0] = mathVar::xGaussLobatto[na];
				mathVar::GaussLobattoPts[na][nb][1] = mathVar::xGaussLobatto[nb];

				mathVar::wGaussPts[na][nb][0] = mathVar::wGauss[na];
				mathVar::wGaussPts[na][nb][1] = mathVar::wGauss[nb];
				mathVar::wGaussLobattoPts[na][nb][0] = mathVar::wGaussLobatto[na];
				mathVar::wGaussLobattoPts[na][nb][1] = mathVar::wGaussLobatto[nb];
			}
		}
	}

	void JacobianParam()
	{
		/*2D Jacobi*/
		double a(0.0), b(0.0);
		for (int ielem = 0; ielem < meshVar::nelem2D; ielem++)
		{
			for (int na = 0; na <= mathVar::nGauss; na++)  //nGauss is started from 0
			{
				for (int nb = 0; nb <= mathVar::nGauss; nb++)
				{
					a = mathVar::GaussPts[na][nb][0];
					b = mathVar::GaussPts[na][nb][1];
					meshVar::J2D[ielem][na][nb] = math::J2DCal(ielem, a, b);
				}
			}
		}

		/*1D Jacobi*/
		for (int iedge = 0; iedge < meshVar::inpoedCount; iedge++)
		{
			std::tie(meshVar::J1D[iedge][0], meshVar::J1D[iedge][1]) = math::J1DCal(iedge);
		}
	}

	void basisFcParam()
	{
		double a(0.0), b(0.0);
		for (int na = 0; na <= mathVar::nGauss; na++)
		{
			for (int nb = 0; nb <= mathVar::nGauss; nb++)
			{
				a = mathVar::GaussPts[na][nb][0];
				b = mathVar::GaussPts[na][nb][1];
				math::basisFc(a, b);
				math::dBasisFc(a, b);
				for (int order = 0; order <= mathVar::orderElem; order++)
				{
					mathVar::BPts[order][na][nb] = mathVar::B[order];

					mathVar::dBaPts[order][na][nb] = mathVar::dBa[order];
					mathVar::dBbPts[order][na][nb] = mathVar::dBb[order];
				}
			}
		}
	}

	void derivCoordinates()
	{
		double dxa(0.0), dxb(0.0), dya(0.0), dyb(0.0);
		double a(0.0), b(0.0);
		for (int ielem = 0; ielem < meshVar::nelem2D; ielem++)
		{
			for (int na = 0; na <= mathVar::nGauss; na++)
			{
				for (int nb = 0; nb <= mathVar::nGauss; nb++)
				{
					a = mathVar::GaussPts[na][nb][0];
					b = mathVar::GaussPts[na][nb][1];
					std::tie(dxa, dxb, dya, dyb) = math::Calc_dxydab(ielem, a, b);
					meshVar::dxa[ielem][na][nb] = dxa;
					meshVar::dxb[ielem][na][nb] = dxb;
					meshVar::dya[ielem][na][nb] = dya;
					meshVar::dyb[ielem][na][nb] = dyb;
				}
			}
		}
	}

	void calcCellMetrics()
	{
		int elemType(0);
		double xCG(0.0), yCG(0.0);
		for (int nelem = 0; nelem < meshVar::nelem2D; nelem++)
		{
			elemType = auxUlti::checkType(nelem);
			std::vector<double> xCoor(elemType, 0.0),
				yCoor(elemType, 0.0),
				xSubTriCoor(elemType, 0.0),
				ySubTriCoor(elemType, 0.0);
			for (int i = 0; i < elemType; i++)
			{
				std::tie(xCoor[i], yCoor[i]) = auxUlti::getElemCornerCoord(nelem, i);
			}
			//1. Calculate cell area (cell size)
			meshVar::cellSize[nelem] = math::geometricOp::calcPolygonArea(xCoor, yCoor, elemType);

			//2. Compute geometric center of cell
			std::tie(xCG, yCG) = math::geometricOp::calcGeoCenter(xCoor, yCoor, elemType);

			//3. Compute centroid of cell
			if (elemType == 3)
			{
				meshVar::geoCenter[nelem][0] = xCG;
				meshVar::geoCenter[nelem][1] = yCG;
			}
			else  //quad element
			{
				//3. Compute geometric center of sub-triangles which creates polygon
				std::tie(meshVar::geoCenter[nelem][0], meshVar::geoCenter[nelem][1]) = math::geometricOp::calcQuadCentroid(nelem, xCG, yCG, meshVar::cellSize[nelem]);
			}
		}
	}
}

namespace process
{
	void setIniValues()
	{
		iniValues::rhoIni = iniValues::pIni / (material::R*iniValues::TIni);
		material::Cp = material::R*material::gamma / (material::gamma - 1);
		material::Cv = material::Cp - material::R;
		iniValues::eIni = material::Cv*iniValues::TIni; //+0.5*(pow(iniValues::uIni, 2) + pow(iniValues::vIni, 2) + pow(iniValues::wIni, 2));
		iniValues::muIni = math::CalcVisCoef(iniValues::TIni);
		std::vector<double> iniRho(mathVar::orderElem + 1),
			iniRhou(mathVar::orderElem + 1),
			iniRhov(mathVar::orderElem + 1),
			iniRhoE(mathVar::orderElem + 1);

		for (int nelement = 0; nelement < meshVar::nelem2D; nelement++)
		{
			iniRho = process::calcIniValues(iniValues::rhoIni, nelement);
			iniRhou = process::calcIniValues(iniValues::rhoIni*iniValues::uIni, nelement);
			iniRhov = process::calcIniValues(iniValues::rhoIni*iniValues::vIni, nelement);
			iniRhoE = process::calcIniValues(iniValues::rhoIni*(iniValues::eIni + 0.5*(pow(iniValues::uIni, 2) + pow(iniValues::vIni, 2))), nelement);

			for (int i = 0; i <= mathVar::orderElem; i++)
			{
				rho[nelement][i] = iniRho[i];
				rhou[nelement][i] = iniRhou[i];
				rhov[nelement][i] = iniRhov[i];
				rhoE[nelement][i] = iniRhoE[i];
			}
		}

		//Calculate limit of rhoE
		//limitVal::rhoEUp = material::Cv*limitVal::TUp*limitVal::rhoUp;
		//limitVal::rhoEDwn = material::Cv*limitVal::TDwn*limitVal::rhoDwn;
	}

	std::vector<double> calcIniValues(double iniVal, int element)
	{
		std::vector<std::vector<double>> matrix(mathVar::orderElem + 1, std::vector<double>(mathVar::orderElem + 1, 0.0));
		std::vector<double> RHS(mathVar::orderElem + 1, 0.0);
		std::vector<double> iniVector(mathVar::orderElem + 1, 0.0);

		matrix = process::calculateStiffMatrix(element);
		RHS = process::calcIniValuesRHS(element, iniVal);
		iniVector = math::SolveSysEqs(matrix, RHS);
		return iniVector;
	}

	std::vector<double> calcIniValuesRHS(int element, double iniVal)
	{
		std::vector<double> Out(mathVar::orderElem + 1, 0.0);
		std::vector<std::vector<double>> matrix(mathVar::nGauss + 1, std::vector<double>(mathVar::nGauss + 1, 0.0));
		double a(0.0), b(0.0);
		for (int order = 0; order <= mathVar::orderElem; order++)
		{
			for (int na = 0; na <= mathVar::nGauss; na++)
			{
				for (int nb = 0; nb <= mathVar::nGauss; nb++)
				{
					std::tie(a, b) = auxUlti::getGaussCoor(na, nb);
					math::basisFc(a, b);
					matrix[na][nb] = mathVar::B[order];
				}
			}
			Out[order] = math::volumeInte(matrix, element)*iniVal;
		}
		return Out;
	}

	namespace auxEq
	{
		void calcValuesAtInterface()
		{
			//mu included
			int masterCell(-1), slaveCell(-1), bcGrp(0);
			double muMaster(0.0), muSlave(0.0), tempMaster(0.0), tempSlave(0.0);
			for (int iedge = 0; iedge < meshVar::inpoedCount; iedge++)
			{
				bcGrp = auxUlti::getBCType(iedge);
				if (bcGrp == 0)
				{
					std::tie(masterCell, slaveCell) = auxUlti::getMasterServantOfEdge(iedge);
					for (int nG = 0; nG <= mathVar::nGauss; nG++)
					{
						std::tie(muMaster, muSlave) = math::internalSurfaceValue(iedge, masterCell, nG, 7, 1);

						//rho*mu
						std::tie(tempMaster, tempSlave) = math::internalSurfaceValue(iedge, masterCell, nG, 1, 2);
						aux_interface_rho[iedge][nG] = tempMaster * muMaster;
						aux_interface_rho[iedge][nG + mathVar::nGauss + 1] = tempSlave * muSlave;
						interface_rho[iedge][nG] = tempMaster;
						interface_rho[iedge][nG + mathVar::nGauss + 1] = tempSlave;

						//rhou*mu
						std::tie(tempMaster, tempSlave) = math::internalSurfaceValue(iedge, masterCell, nG, 2, 2);
						aux_interface_rhou[iedge][nG] = tempMaster * muMaster;
						aux_interface_rhou[iedge][nG + mathVar::nGauss + 1] = tempSlave * muSlave;
						interface_rhou[iedge][nG] = tempMaster;
						interface_rhou[iedge][nG + mathVar::nGauss + 1] = tempSlave;

						//rhov*mu
						std::tie(tempMaster, tempSlave) = math::internalSurfaceValue(iedge, masterCell, nG, 3, 2);
						aux_interface_rhov[iedge][nG] = tempMaster * muMaster;
						aux_interface_rhov[iedge][nG + mathVar::nGauss + 1] = tempSlave * muSlave;
						interface_rhov[iedge][nG] = tempMaster;
						interface_rhov[iedge][nG + mathVar::nGauss + 1] = tempSlave;

						//rhoE*mu
						std::tie(tempMaster, tempSlave) = math::internalSurfaceValue(iedge, masterCell, nG, 4, 2);
						aux_interface_rhoE[iedge][nG] = tempMaster * muMaster;
						aux_interface_rhoE[iedge][nG + mathVar::nGauss + 1] = tempSlave * muSlave;
						interface_rhoE[iedge][nG] = tempMaster;
						interface_rhoE[iedge][nG + mathVar::nGauss + 1] = tempSlave;
					}
				}
			}
		}

		std::tuple<double, double> getInternalValuesFromCalculatedArrays(int edge, int element, int nG, int valType)
		{
			//valType = 1(rho*mu), 2(rhou*mu), 3(rhov*mu), 4(rhoE*mu)
			bool isMaster(auxUlti::checkMaster(element, edge));
			int locationPlus(-1), locationMinus(-1);
			double valPlus(0.0), valMinus(0.0);
			if (isMaster)
			{
				locationPlus = nG;
				locationMinus = nG + mathVar::nGauss + 1;
			}
			else
			{
				locationPlus = nG + mathVar::nGauss + 1;
				locationMinus = nG;
			}

			switch (valType)
			{
			case 1:
			{
				valPlus = aux_interface_rho[edge][locationPlus];
				valMinus = aux_interface_rho[edge][locationMinus];
			}
			break;
			case 2:
			{
				valPlus = aux_interface_rhou[edge][locationPlus];
				valMinus = aux_interface_rhou[edge][locationMinus];
			}
			break;
			case 3:
			{
				valPlus = aux_interface_rhov[edge][locationPlus];
				valMinus = aux_interface_rhov[edge][locationMinus];
			}
			break;
			case 4:
			{
				valPlus = aux_interface_rhoE[edge][locationPlus];
				valMinus = aux_interface_rhoE[edge][locationMinus];
			}
			break;
			default:
				break;
			}
			return std::make_tuple(valPlus, valMinus);
		}

		void solveAuxEquation()
		{
			std::vector<std::vector<double>> StiffMatrix(mathVar::orderElem + 1, std::vector<double>(mathVar::orderElem + 1, 0.0));
			std::vector<double> rhoRHSTermOxDir(mathVar::orderElem + 1, 0.0),
				rhoRHSTermOyDir(mathVar::orderElem + 1, 0.0),
				rhouRHSTermOxDir(mathVar::orderElem + 1, 0.0),
				rhouRHSTermOyDir(mathVar::orderElem + 1, 0.0),
				rhovRHSTermOxDir(mathVar::orderElem + 1, 0.0),
				rhovRHSTermOyDir(mathVar::orderElem + 1, 0.0),
				rhoERHSTermOxDir(mathVar::orderElem + 1, 0.0),
				rhoERHSTermOyDir(mathVar::orderElem + 1, 0.0);

			std::vector<double> rhoxVector(mathVar::orderElem + 1, 0.0),
				rhoyVector(mathVar::orderElem + 1, 0.0),
				rhouxVector(mathVar::orderElem + 1, 0.0),
				rhouyVector(mathVar::orderElem + 1, 0.0),
				rhovxVector(mathVar::orderElem + 1, 0.0),
				rhovyVector(mathVar::orderElem + 1, 0.0),
				rhoExVector(mathVar::orderElem + 1, 0.0),
				rhoEyVector(mathVar::orderElem + 1, 0.0);

			for (int nelement = 0; nelement < meshVar::nelem2D; nelement++)
			{
				//1) Calculate Stiff matrix
				StiffMatrix = process::calculateStiffMatrix(nelement);
				
				//2) Calculate Right hand side terms
				process::auxEq::CalcRHSTerm(nelement, rhoRHSTermOxDir, rhoRHSTermOyDir, rhouRHSTermOxDir, rhouRHSTermOyDir, rhovRHSTermOxDir, rhovRHSTermOyDir, rhoERHSTermOxDir, rhoERHSTermOyDir);

				//3) Solve for auxilary variables
				//Ox direction
				rhoxVector = math::SolveSysEqs(StiffMatrix, rhoRHSTermOxDir);
				rhouxVector = math::SolveSysEqs(StiffMatrix, rhouRHSTermOxDir);
				rhovxVector = math::SolveSysEqs(StiffMatrix, rhovRHSTermOxDir);
				rhoExVector = math::SolveSysEqs(StiffMatrix, rhoERHSTermOxDir);

				//Oy direction
				rhoyVector = math::SolveSysEqs(StiffMatrix, rhoRHSTermOyDir);
				rhouyVector = math::SolveSysEqs(StiffMatrix, rhouRHSTermOyDir);
				rhovyVector = math::SolveSysEqs(StiffMatrix, rhovRHSTermOyDir);
				rhoEyVector = math::SolveSysEqs(StiffMatrix, rhoERHSTermOyDir);

				//4) Save results to auxilary variables array

				for (int order = 0; order <= mathVar::orderElem; order++)
				{
					rhoX[nelement][order] = rhoxVector[order];
					rhoY[nelement][order] = rhoyVector[order];
					rhouX[nelement][order] = rhouxVector[order];
					rhouY[nelement][order] = rhouyVector[order];
					rhovX[nelement][order] = rhovxVector[order];
					rhovY[nelement][order] = rhovyVector[order];
					rhoEX[nelement][order] = rhoExVector[order];
					rhoEY[nelement][order] = rhoEyVector[order];
				}
			}
		}

		/*Function calculates right hand side terms of all conservative variables at all order in all directions*/
		void CalcRHSTerm(int element, std::vector<double> &rhoRHSOx, std::vector<double> &rhoRHSOy, std::vector<double> &rhouRHSOx, std::vector<double> &rhouRHSOy, std::vector<double> &rhovRHSOx, std::vector<double> &rhovRHSOy, std::vector<double> &rhoERHSOx, std::vector<double> &rhoERHSOy)
		{
			std::vector<double>
				//vectors for volume integrals
				rhoVolIntOx(mathVar::orderElem + 1, 0.0),
				rhoVolIntOy(mathVar::orderElem + 1, 0.0),
				rhouVolIntOx(mathVar::orderElem + 1, 0.0),
				rhouVolIntOy(mathVar::orderElem + 1, 0.0),
				rhovVolIntOx(mathVar::orderElem + 1, 0.0),
				rhovVolIntOy(mathVar::orderElem + 1, 0.0),
				rhoEVolIntOx(mathVar::orderElem + 1, 0.0),
				rhoEVolIntOy(mathVar::orderElem + 1, 0.0),

				//vectors for surface integrals
				rhoSurfIntOx(mathVar::orderElem + 1, 0.0),
				rhoSurfIntOy(mathVar::orderElem + 1, 0.0),
				rhouSurfIntOx(mathVar::orderElem + 1, 0.0),
				rhouSurfIntOy(mathVar::orderElem + 1, 0.0),
				rhovSurfIntOx(mathVar::orderElem + 1, 0.0),
				rhovSurfIntOy(mathVar::orderElem + 1, 0.0),
				rhoESurfIntOx(mathVar::orderElem + 1, 0.0),
				rhoESurfIntOy(mathVar::orderElem + 1, 0.0);

			/*1. Calculate volume integral term*/
			process::auxEq::calcVolumeIntegralTerms(element, rhoVolIntOx, rhouVolIntOx, rhovVolIntOx, rhoEVolIntOx, rhoVolIntOy, rhouVolIntOy, rhovVolIntOy, rhoEVolIntOy);

			/*2. Calculate surface integral term*/
			process::auxEq::calcSurfaceIntegralTerms(element, rhoSurfIntOx, rhouSurfIntOx, rhovSurfIntOx, rhoESurfIntOx, rhoSurfIntOy, rhouSurfIntOy, rhovSurfIntOy, rhoESurfIntOy);

			for (int order = 0; order <= mathVar::orderElem; order++)
			{
				rhoRHSOx[order] = -rhoVolIntOx[order] + rhoSurfIntOx[order];
				rhoRHSOy[order] = -rhoVolIntOy[order] + rhoSurfIntOy[order];

				rhouRHSOx[order] = -rhouVolIntOx[order] + rhouSurfIntOx[order];
				rhouRHSOy[order] = -rhouVolIntOy[order] + rhouSurfIntOy[order];

				rhovRHSOx[order] = -rhovVolIntOx[order] + rhovSurfIntOx[order];
				rhovRHSOy[order] = -rhovVolIntOy[order] + rhovSurfIntOy[order];

				rhoERHSOx[order] = -rhoEVolIntOx[order] + rhoESurfIntOx[order];
				rhoERHSOy[order] = -rhoEVolIntOy[order] + rhoESurfIntOy[order];
			}
		}

		void calcVolumeIntegralTerms(int element, std::vector<double> &rhoVolIntX, std::vector<double> &rhouVolIntX, std::vector<double> &rhovVolIntX, std::vector<double> &rhoEVolIntX, std::vector<double> &rhoVolIntY, std::vector<double> &rhouVolIntY, std::vector<double> &rhovVolIntY, std::vector<double> &rhoEVolIntY)
		{
			std::vector<std::vector<double>> rhoGsVol(mathVar::nGauss + 1, std::vector<double>(mathVar::nGauss + 1, 0.0)),
				rhouGsVol(mathVar::nGauss + 1, std::vector<double>(mathVar::nGauss + 1, 0.0)),
				rhovGsVol(mathVar::nGauss + 1, std::vector<double>(mathVar::nGauss + 1, 0.0)),
				rhoEGsVol(mathVar::nGauss + 1, std::vector<double>(mathVar::nGauss + 1, 0.0));

			//Calculates Gauss matrix
			//rho -------------------------------------------------------------------------------------------
			rhoGsVol = process::auxEq::getGaussMatrixOfConserVar(element, 1);
			//rhou ------------------------------------------------------------------------------------------
			rhouGsVol = process::auxEq::getGaussMatrixOfConserVar(element, 2);
			//rhov ------------------------------------------------------------------------------------------
			rhovGsVol = process::auxEq::getGaussMatrixOfConserVar(element, 3);
			//rhou ------------------------------------------------------------------------------------------
			rhoEGsVol = process::auxEq::getGaussMatrixOfConserVar(element, 4);

			for (int order = 0; order <= mathVar::orderElem; order++)
			{
				rhoVolIntX[order] = process::volumeInte(element, rhoGsVol, order, 1);
				rhouVolIntX[order] = process::volumeInte(element, rhouGsVol, order, 1);
				rhovVolIntX[order] = process::volumeInte(element, rhovGsVol, order, 1);
				rhoEVolIntX[order] = process::volumeInte(element, rhoEGsVol, order, 1);

				rhoVolIntY[order] = process::volumeInte(element, rhoGsVol, order, 2);
				rhouVolIntY[order] = process::volumeInte(element, rhouGsVol, order, 2);
				rhovVolIntY[order] = process::volumeInte(element, rhovGsVol, order, 2);
				rhoEVolIntY[order] = process::volumeInte(element, rhoEGsVol, order, 2);
			}
		}

		void calcSurfaceIntegralTerms(int element, std::vector<double> &rhoSurfIntX, std::vector<double> &rhouSurfIntX, std::vector<double> &rhovSurfIntX, std::vector<double> &rhoESurfIntX, std::vector<double> &rhoSurfIntY, std::vector<double> &rhouSurfIntY, std::vector<double> &rhovSurfIntY, std::vector<double> &rhoESurfIntY)
		{
			/*User's guide:
			Input array rhoRHSTerm, rhouRHSTerm, rhovRHSTerm, rhoERHSTerm have following form:
			- number of row: orderElem + 1*/
			int elemType(auxUlti::checkType(element)), edgeName(0);
			std::vector<std::vector<double>> rhoFluxX(mathVar::nGauss + 1, std::vector<double>(elemType, 0.0)),
				rhouFluxX(mathVar::nGauss + 1, std::vector<double>(elemType, 0.0)),
				rhovFluxX(mathVar::nGauss + 1, std::vector<double>(elemType, 0.0)),
				rhoEFluxX(mathVar::nGauss + 1, std::vector<double>(elemType, 0.0)),
				rhoFluxY(mathVar::nGauss + 1, std::vector<double>(elemType, 0.0)),
				rhouFluxY(mathVar::nGauss + 1, std::vector<double>(elemType, 0.0)),
				rhovFluxY(mathVar::nGauss + 1, std::vector<double>(elemType, 0.0)),
				rhoEFluxY(mathVar::nGauss + 1, std::vector<double>(elemType, 0.0));

			std::vector<double> rhoFluxXTemp(mathVar::nGauss + 1, 0.0),
				rhouFluxXTemp(mathVar::nGauss + 1, 0.0),
				rhovFluxXTemp(mathVar::nGauss + 1, 0.0),
				rhoEFluxXTemp(mathVar::nGauss + 1, 0.0),
				rhoFluxYTemp(mathVar::nGauss + 1, 0.0),
				rhouFluxYTemp(mathVar::nGauss + 1, 0.0),
				rhovFluxYTemp(mathVar::nGauss + 1, 0.0),
				rhoEFluxYTemp(mathVar::nGauss + 1, 0.0);

			/*1. Calculate flux of conservative variables at all Gauss points on all faces of element*/
			process::auxEq::getGaussVectorOfConserVar(element, rhoFluxX, rhouFluxX, rhovFluxX, rhoEFluxX, rhoFluxY, rhouFluxY, rhovFluxY, rhoEFluxY);

			/*2. Calculates surface integrals of all conservative variables at all order*/
			
			for (int nface = 0; nface < elemType; nface++)
			{
				edgeName = meshVar::inedel[nface][element];
				for (int nG = 0; nG <= mathVar::nGauss; nG++)
				{
					rhoFluxXTemp[nG] = rhoFluxX[nG][nface];
					rhouFluxXTemp[nG] = rhouFluxX[nG][nface];
					rhovFluxXTemp[nG] = rhovFluxX[nG][nface];
					rhoEFluxXTemp[nG] = rhoEFluxX[nG][nface];

					rhoFluxYTemp[nG] = rhoFluxY[nG][nface];
					rhouFluxYTemp[nG] = rhouFluxY[nG][nface];
					rhovFluxYTemp[nG] = rhovFluxY[nG][nface];
					rhoEFluxYTemp[nG] = rhoEFluxY[nG][nface];
				}

				for (int order = 0; order <= mathVar::orderElem; order++)
				{
					rhoSurfIntX[order] += process::surfaceInte(element, edgeName, rhoFluxXTemp, order);
					rhouSurfIntX[order] += process::surfaceInte(element, edgeName, rhouFluxXTemp, order);
					rhovSurfIntX[order] += process::surfaceInte(element, edgeName, rhovFluxXTemp, order);
					rhoESurfIntX[order] += process::surfaceInte(element, edgeName, rhoEFluxXTemp, order);

					rhoSurfIntY[order] += process::surfaceInte(element, edgeName, rhoFluxYTemp, order);
					rhouSurfIntY[order] += process::surfaceInte(element, edgeName, rhouFluxYTemp, order);
					rhovSurfIntY[order] += process::surfaceInte(element, edgeName, rhovFluxYTemp, order);
					rhoESurfIntY[order] += process::surfaceInte(element, edgeName, rhoEFluxYTemp, order);
				}
			}
		}

		std::vector<std::vector<double>> getGaussMatrixOfConserVar(int element, int valType)
		{
			std::vector<std::vector<double>> GaussMatrix(mathVar::nGauss + 1, std::vector<double>(mathVar::nGauss + 1, 0.0));
			double muGs(0.0), a(0.0), b(0.0);
			for (int na = 0; na <= mathVar::nGauss; na++)
			{
				for (int nb = 0; nb <= mathVar::nGauss; nb++)
				{
					std::tie(a, b) = auxUlti::getGaussCoor(na, nb);
					muGs = math::pointValue(element, a, b, 7, 1);
					GaussMatrix[na][nb] = math::pointValue(element, a, b, valType, 2)*muGs;
				}
			}
			return GaussMatrix;
		}

		std::vector<std::vector<double>> getVectorOfConserVarFluxesAtInternal(int edge, int element, int nG, double nx, double ny)
		{
			std::vector<std::vector<double>> gaussVector(4, std::vector<double>(2, 0.0)); //columns 0, 1 are Ox, Oy values
			double tempP(0.0), tempM(0.0);
			for (int i = 0; i < 4; i++)
			{
				std::tie(tempP, tempM) = process::auxEq::getInternalValuesFromCalculatedArrays(edge, element, nG, i+1);
				gaussVector[i][0] = math::numericalFluxes::auxFlux(tempM, tempP, nx);
				gaussVector[i][1] = math::numericalFluxes::auxFlux(tempM, tempP, ny);
			}

			return gaussVector;
		}

		void getGaussVectorOfConserVar(int element, std::vector<std::vector<double>> &rhoFluxX, std::vector<std::vector<double>> &rhouFluxX, std::vector<std::vector<double>> &rhovFluxX, std::vector<std::vector<double>> &rhoEFluxX, std::vector<std::vector<double>> &rhoFluxY, std::vector<std::vector<double>> &rhouFluxY, std::vector<std::vector<double>> &rhovFluxY, std::vector<std::vector<double>> &rhoEFluxY)
		{
			/*User's guide:
			Input array rhoFlux, rhouFlux, rhovFlux, rhoEFlux have following form:
			- number of row: nGauss + 1*/

			int elemType(auxUlti::checkType(element)), edgeName(0);
			int faceBcType(0);
			double nVectorComp(0.0);
			std::vector<std::vector<double>> gaussVector(4, std::vector<double>(2, 0.0));

			for (int nface = 0; nface < elemType; nface++)
			{
				edgeName = meshVar::inedel[nface][element];
				faceBcType = auxUlti::getBCType(edgeName);
				double nx(auxUlti::getNormVectorComp(element, edgeName, 1)), ny(auxUlti::getNormVectorComp(element, edgeName, 2));
				if (faceBcType == 0)  //internal edge
				{
					for (int nGauss = 0; nGauss <= mathVar::nGauss; nGauss++)
					{
						gaussVector = process::auxEq::getVectorOfConserVarFluxesAtInternal(edgeName, element, nGauss, nx, ny);
						rhoFluxX[nGauss][nface] = gaussVector[0][0];
						rhouFluxX[nGauss][nface] = gaussVector[1][0];
						rhovFluxX[nGauss][nface] = gaussVector[2][0];
						rhoEFluxX[nGauss][nface] = gaussVector[3][0];

						rhoFluxY[nGauss][nface] = gaussVector[0][1];
						rhouFluxY[nGauss][nface] = gaussVector[1][1];
						rhovFluxY[nGauss][nface] = gaussVector[2][1];
						rhoEFluxY[nGauss][nface] = gaussVector[3][1];
					}
				}
				else  //boundary edge
				{
					for (int nGauss = 0; nGauss <= mathVar::nGauss; nGauss++)
					{
						gaussVector = auxEqBCsImplement(element, edgeName, nGauss, nx, ny);
						rhoFluxX[nGauss][nface] = gaussVector[0][0];
						rhouFluxX[nGauss][nface] = gaussVector[1][0];
						rhovFluxX[nGauss][nface] = gaussVector[2][0];
						rhoEFluxX[nGauss][nface] = gaussVector[3][0];

						rhoFluxY[nGauss][nface] = gaussVector[0][1];
						rhouFluxY[nGauss][nface] = gaussVector[1][1];
						rhovFluxY[nGauss][nface] = gaussVector[2][1];
						rhoEFluxY[nGauss][nface] = gaussVector[3][1];
					}
				}
			}
		}
	}//end namespace auxEq

	namespace NSFEq
	{
		void calcValuesAtInterface()
		{
			//mu included
			int masterCell(-1), slaveCell(-1), bcGrp(0);
			double uMaster(0.0), vMaster(0.0), totalEMaster(0.0), TMaster(0.0), pMaster(0.0),
				uSlave(0.0), vSlave(0.0), totalESlave(0.0), TSlave(0.0), pSlave(0.0),
				uMagM(0.0), uMagP(0.0), aM(0.0), aP(0.0),
				dRhoXMaster(0.0), dRhouXMaster(0.0), dRhovXMaster(0.0), dRhoEXMaster(0.0),
				dRhoYMaster(0.0), dRhouYMaster(0.0), dRhovYMaster(0.0), dRhoEYMaster(0.0),
				dRhoXSlave(0.0), dRhouXSlave(0.0), dRhovXSlave(0.0), dRhoEXSlave(0.0),
				dRhoYSlave(0.0), dRhouYSlave(0.0), dRhovYSlave(0.0), dRhoEYSlave(0.0);
			std::vector<double> UMaster(4, 0.0), dUXMaster(4, 0.0), dUYMaster(4, 0.0),
				USlave(4, 0.0), dUXSlave(4, 0.0), dUYSlave(4, 0.0);
			/*StressHeat matrix has form:
			[tauXx		tauXy		Qx]
			[tauYx		tauYy		Qy]
			*/
			std::vector<std::vector<double>> StressHeatP(2, std::vector<double>(3, 0.0));
			std::vector<std::vector<double>> StressHeatM(2, std::vector<double>(3, 0.0));
			for (int iedge = 0; iedge < meshVar::inpoedCount; iedge++)
			{
				bcGrp = auxUlti::getBCType(iedge);
				if (bcGrp == 0)
				{
					std::tie(masterCell, slaveCell) = auxUlti::getMasterServantOfEdge(iedge);
					for (int nG = 0; nG <= mathVar::nGauss; nG++)
					{
						for (int i = 0; i < 4; i++)
						{
							/*INVISCID TERMS*/
							std::tie(UMaster[i], USlave[i]) = math::internalSurfaceValue(iedge, masterCell, nG, i + 1, 2);
							/*VISCOUS TERMS*/
							std::tie(dUXMaster[i], dUXSlave[i]) = math::internalSurfaceDerivativeValue(iedge, masterCell, nG, i + 1, 1);  //dUx
							std::tie(dUYMaster[i], dUYSlave[i]) = math::internalSurfaceDerivativeValue(iedge, masterCell, nG, i + 1, 2);  //dUy
						}

						uMaster = (UMaster[1] / UMaster[0]);
						uSlave = (USlave[1] / USlave[0]);

						vMaster = (UMaster[2] / UMaster[0]);
						vSlave = (USlave[2] / USlave[0]);

						totalEMaster = (UMaster[3] / UMaster[0]);
						totalESlave = (USlave[3] / USlave[0]);

						/*calculate T and P*/
						TMaster = math::CalcTFromConsvVar(UMaster[0], UMaster[1], UMaster[2], UMaster[3]);
						TSlave = math::CalcTFromConsvVar(USlave[0], USlave[1], USlave[2], USlave[3]);
						pMaster = math::CalcP(TMaster, UMaster[0]);
						pSlave = math::CalcP(TSlave, USlave[0]);

						/*INVISCID TERMS*/
						/*calculate velocity magnitude*/
						uMagP = sqrt(pow(uMaster, 2) + pow(vMaster, 2));
						uMagM = sqrt(pow(uSlave, 2) + pow(vSlave, 2));

						/*calculate speed of sound*/
						aP = math::CalcSpeedOfSound(TMaster);
						aM = math::CalcSpeedOfSound(TSlave);

						/*calculate constant for Lax-Friederich flux*/
						LxFConst[iedge] = math::numericalFluxes::constantC(uMagP, uMagM, aP, aM);

						/*calculate inviscid terms*/
						
						std::tie(invis_interface_rhoX[iedge][nG], invis_interface_rhouX[iedge][nG], invis_interface_rhovX[iedge][nG], invis_interface_rhoEX[iedge][nG]) = math::inviscidTerms::calcInvisTermsFromPriVars(UMaster[0], uMaster, vMaster, totalEMaster, pMaster, 1);
						std::tie(invis_interface_rhoY[iedge][nG], invis_interface_rhouY[iedge][nG], invis_interface_rhovY[iedge][nG], invis_interface_rhoEY[iedge][nG]) = math::inviscidTerms::calcInvisTermsFromPriVars(UMaster[0], uMaster, vMaster, totalEMaster, pMaster, 2);

						std::tie(invis_interface_rhoX[iedge][nG + mathVar::nGauss + 1], invis_interface_rhouX[iedge][nG + mathVar::nGauss + 1], invis_interface_rhovX[iedge][nG + mathVar::nGauss + 1], invis_interface_rhoEX[iedge][nG + mathVar::nGauss + 1]) = math::inviscidTerms::calcInvisTermsFromPriVars(USlave[0], uSlave, vSlave, totalESlave, pSlave, 1);
						std::tie(invis_interface_rhoY[iedge][nG + mathVar::nGauss + 1], invis_interface_rhouY[iedge][nG + mathVar::nGauss + 1], invis_interface_rhovY[iedge][nG + mathVar::nGauss + 1], invis_interface_rhoEY[iedge][nG + mathVar::nGauss + 1]) = math::inviscidTerms::calcInvisTermsFromPriVars(USlave[0], uSlave, vSlave, totalESlave, pSlave, 2);

						/*calculate viscous terms*/
						StressHeatP = math::viscousTerms::calcStressTensorAndHeatFlux(UMaster, dUXMaster, dUYMaster);
						std::tie(Vis_interface_rhoX[iedge][nG], Vis_interface_rhouX[iedge][nG], Vis_interface_rhovX[iedge][nG], Vis_interface_rhoEX[iedge][nG]) = math::viscousTerms::calcViscousTermsFromStressHeatFluxMatrix(StressHeatP, uMaster, vMaster, 1);
						std::tie(Vis_interface_rhoY[iedge][nG], Vis_interface_rhouY[iedge][nG], Vis_interface_rhovY[iedge][nG], Vis_interface_rhoEY[iedge][nG]) = math::viscousTerms::calcViscousTermsFromStressHeatFluxMatrix(StressHeatP, uMaster, vMaster, 2);

						StressHeatM = math::viscousTerms::calcStressTensorAndHeatFlux(USlave, dUXSlave, dUYSlave);
						std::tie(Vis_interface_rhoX[iedge][nG + mathVar::nGauss + 1], Vis_interface_rhouX[iedge][nG + mathVar::nGauss + 1], Vis_interface_rhovX[iedge][nG + mathVar::nGauss + 1], Vis_interface_rhoEX[iedge][nG + mathVar::nGauss + 1]) = math::viscousTerms::calcViscousTermsFromStressHeatFluxMatrix(StressHeatM, uSlave, vSlave, 1);
						std::tie(Vis_interface_rhoY[iedge][nG + mathVar::nGauss + 1], Vis_interface_rhouY[iedge][nG + mathVar::nGauss + 1], Vis_interface_rhovY[iedge][nG + mathVar::nGauss + 1], Vis_interface_rhoEY[iedge][nG + mathVar::nGauss + 1]) = math::viscousTerms::calcViscousTermsFromStressHeatFluxMatrix(StressHeatM, uSlave, vSlave, 2);
					}
				}
			}
		}

		std::tuple<double, double> getInternalValuesFromCalculatedArrays(int edge, int element, int nG, int mod, int direction, int valType)
		{
			/*
			-mod: 1 for inviscid, 2 for viscous
			-direction: 1 for Ox, 2 for Oy
			*/
			bool isMaster(auxUlti::checkMaster(element, edge));
			int locationPlus(-1), locationMinus(-1);
			double valPlus(0.0), valMinus(0.0);
			if (isMaster)
			{
				locationPlus = nG;
				locationMinus = nG + mathVar::nGauss + 1;
			}
			else
			{
				locationPlus = nG + mathVar::nGauss + 1;
				locationMinus = nG;
			}

			switch (mod)
			{
			case 1: //inviscid
			{
				switch (direction)
				{
				case 1:
				{
					switch (valType)
					{
					case 1:
					{
						valPlus = invis_interface_rhoX[edge][locationPlus];
						valMinus = invis_interface_rhoX[edge][locationMinus];
					}
					break;
					case 2:
					{
						valPlus = invis_interface_rhouX[edge][locationPlus];
						valMinus = invis_interface_rhouX[edge][locationMinus];
					}
					break;
					case 3:
					{
						valPlus = invis_interface_rhovX[edge][locationPlus];
						valMinus = invis_interface_rhovX[edge][locationMinus];
					}
					break;
					case 4:
					{
						valPlus = invis_interface_rhoEX[edge][locationPlus];
						valMinus = invis_interface_rhoEX[edge][locationMinus];
					}
					break;
					default:
						break;
					}
				}
				break;
				case 2:
				{
					switch (valType)
					{
					case 1:
					{
						valPlus = invis_interface_rhoY[edge][locationPlus];
						valMinus = invis_interface_rhoY[edge][locationMinus];
					}
					break;
					case 2:
					{
						valPlus = invis_interface_rhouY[edge][locationPlus];
						valMinus = invis_interface_rhouY[edge][locationMinus];
					}
					break;
					case 3:
					{
						valPlus = invis_interface_rhovY[edge][locationPlus];
						valMinus = invis_interface_rhovY[edge][locationMinus];
					}
					break;
					case 4:
					{
						valPlus = invis_interface_rhoEY[edge][locationPlus];
						valMinus = invis_interface_rhoEY[edge][locationMinus];
					}
					break;
					default:
						break;
					}
				}
				break;
				default:
					break;
				}
			}
				break;
			case 2: //viscous
			{
				switch (direction)
				{
				case 1:
				{
					switch (valType)
					{
					case 1:
					{
						valPlus = Vis_interface_rhoX[edge][locationPlus];
						valMinus = Vis_interface_rhoX[edge][locationMinus];
					}
					break;
					case 2:
					{
						valPlus = Vis_interface_rhouX[edge][locationPlus];
						valMinus = Vis_interface_rhouX[edge][locationMinus];
					}
					break;
					case 3:
					{
						valPlus = Vis_interface_rhovX[edge][locationPlus];
						valMinus = Vis_interface_rhovX[edge][locationMinus];
					}
					break;
					case 4:
					{
						valPlus = Vis_interface_rhoEX[edge][locationPlus];
						valMinus = Vis_interface_rhoEX[edge][locationMinus];
					}
					break;
					default:
						break;
					}
				}
				break;
				case 2:
				{
					switch (valType)
					{
					case 1:
					{
						valPlus = Vis_interface_rhoY[edge][locationPlus];
						valMinus = Vis_interface_rhoY[edge][locationMinus];
					}
					break;
					case 2:
					{
						valPlus = Vis_interface_rhouY[edge][locationPlus];
						valMinus = Vis_interface_rhouY[edge][locationMinus];
					}
					break;
					case 3:
					{
						valPlus = Vis_interface_rhovY[edge][locationPlus];
						valMinus = Vis_interface_rhovY[edge][locationMinus];
					}
					break;
					case 4:
					{
						valPlus = Vis_interface_rhoEY[edge][locationPlus];
						valMinus = Vis_interface_rhoEY[edge][locationMinus];
					}
					break;
					default:
						break;
					}
				}
				break;
				default:
					break;
				}
			}
				break;
			default:
				break;
			}
			return std::make_tuple(valPlus, valMinus);
		}

		void solveNSFEquation()
		{
			std::vector<double> rhoError(meshVar::nelem2D, 1.0),
				rhouError(meshVar::nelem2D, 1.0),
				rhovError(meshVar::nelem2D, 1.0),
				rhoEError(meshVar::nelem2D, 1.0);

			std::vector<std::vector<double>> StiffMatrix(mathVar::orderElem + 1, std::vector<double>(mathVar::orderElem + 1, 0.0));
			std::vector<double>
				RHSTerm1(mathVar::orderElem + 1, 0.0),
				RHSTerm2(mathVar::orderElem + 1, 0.0),
				RHSTerm3(mathVar::orderElem + 1, 0.0),
				RHSTerm4(mathVar::orderElem + 1, 0.0),

				ddtRhoVector(mathVar::orderElem + 1, 0.0),
				ddtRhouVector(mathVar::orderElem + 1, 0.0),
				ddtRhovVector(mathVar::orderElem + 1, 0.0),
				ddtRhoEVector(mathVar::orderElem + 1, 0.0),

				rhoVectorN(mathVar::orderElem + 1, 0.0),
				rhouVectorN(mathVar::orderElem + 1, 0.0),
				rhovVectorN(mathVar::orderElem + 1, 0.0),
				rhoEVectorN(mathVar::orderElem + 1, 0.0),
				
				UnVector(mathVar::orderElem + 1, 0.0),
				
				timeStepArr(meshVar::nelem2D,1.0);

			double rhoRes(1.0), rhouRes(1.0), rhovRes(1.0), rhoERes(1.0);

			for (int nelement = 0; nelement < meshVar::nelem2D; nelement++)
			{
				//1) Calculate Stiff matrix
				StiffMatrix = process::calculateStiffMatrix(nelement);

				//2) Calculate Right hand side terms
				process::NSFEq::CalcRHSTerm(nelement, RHSTerm1, RHSTerm2, RHSTerm3, RHSTerm4);

				//3) Solve for derivartives of conservative variables
				ddtRhoVector = math::SolveSysEqs(StiffMatrix, RHSTerm1);
				ddtRhouVector = math::SolveSysEqs(StiffMatrix, RHSTerm2);
				ddtRhovVector = math::SolveSysEqs(StiffMatrix, RHSTerm3);
				ddtRhoEVector = math::SolveSysEqs(StiffMatrix, RHSTerm4);

				//4) Solve time marching
				//rho
				for (int order = 0; order <= mathVar::orderElem; order++)
				{
					UnVector[order] = rho[nelement][order];
				}
				rhoVectorN = process::NSFEq::solveTimeMarching(ddtRhoVector, UnVector);
				//rhou
				for (int order = 0; order <= mathVar::orderElem; order++)
				{
					UnVector[order] = rhou[nelement][order];
				}
				rhouVectorN = process::NSFEq::solveTimeMarching(ddtRhouVector, UnVector);
				//rhov
				for (int order = 0; order <= mathVar::orderElem; order++)
				{
					UnVector[order] = rhov[nelement][order];
				}
				rhovVectorN = process::NSFEq::solveTimeMarching(ddtRhovVector, UnVector);
				//rhoE
				for (int order = 0; order <= mathVar::orderElem; order++)
				{
					UnVector[order] = rhoE[nelement][order];
				}
				rhoEVectorN = process::NSFEq::solveTimeMarching(ddtRhoEVector, UnVector);

				//5) Save results to conservative variables array
				for (int order = 0; order <= mathVar::orderElem; order++)
				{
					rhoN[nelement][order] = rhoVectorN[order];
					rhouN[nelement][order] = rhouVectorN[order];
					rhovN[nelement][order] = rhovVectorN[order];
					rhoEN[nelement][order] = rhoEVectorN[order];
				}

				//6) Estimate Residuals
				rhoError[nelement] = (process::Euler::localErrorEstimate(nelement, ddtRhoVector));
				rhouError[nelement] = (process::Euler::localErrorEstimate(nelement, ddtRhouVector));
				rhovError[nelement] = (process::Euler::localErrorEstimate(nelement, ddtRhovVector));
				rhoEError[nelement] = (process::Euler::localErrorEstimate(nelement, ddtRhoEVector));
			}
			std::tie(rhoRes, rhouRes, rhovRes, rhoERes) = process::Euler::globalErrorEstimate(rhoError, rhouError, rhovError, rhoEError);
			IO::residualOutput(rhoRes, rhouRes, rhovRes, rhoERes);
		}

		void updateVariables()
		{
			for (int nelement = 0; nelement < meshVar::nelem2D; nelement++)
			{
				for (int order = 0; order <= mathVar::orderElem; order++)
				{
					rho[nelement][order] = rhoN[nelement][order];
					rhou[nelement][order] = rhouN[nelement][order];
					rhov[nelement][order] = rhovN[nelement][order];
					rhoE[nelement][order] = rhoEN[nelement][order];
				}
			}
		}

		/*Function calculates right hand side terms of all conservative variables at ONLY one order*/
		void CalcRHSTerm(int element, std::vector<double> &term1RHS, std::vector<double> &term2RHS, std::vector<double> &term3RHS, std::vector<double> &term4RHS)
		{
			std::vector<double>
				VolIntTerm1(mathVar::orderElem + 1, 0.0),
				VolIntTerm2(mathVar::orderElem + 1, 0.0),
				VolIntTerm3(mathVar::orderElem + 1, 0.0),
				VolIntTerm4(mathVar::orderElem + 1, 0.0),

				SurfIntTerm1(mathVar::orderElem + 1, 0.0),
				SurfIntTerm2(mathVar::orderElem + 1, 0.0),
				SurfIntTerm3(mathVar::orderElem + 1, 0.0),
				SurfIntTerm4(mathVar::orderElem + 1, 0.0);

			/*Volume integral term===========================================================================*/
			process::NSFEq::calcVolumeIntegralTerms(element, VolIntTerm1, VolIntTerm2, VolIntTerm3, VolIntTerm4);
			
			/*Surface integral term===========================================================================*/
			process::NSFEq::calcSurfaceIntegralTerms(element, SurfIntTerm1, SurfIntTerm2, SurfIntTerm3, SurfIntTerm4);
			
			for (int order = 0; order <= mathVar::orderElem; order++)
			{
				term1RHS[order] = VolIntTerm1[order] - SurfIntTerm1[order];
				term2RHS[order] = VolIntTerm2[order] - SurfIntTerm2[order];
				term3RHS[order] = VolIntTerm3[order] - SurfIntTerm3[order];
				term4RHS[order] = VolIntTerm4[order] - SurfIntTerm4[order];
			}
		}

		/*Function calculates Inviscid terms at Gauss point (a, b)*/
		std::vector<std::vector<double>> calcGaussInviscidTerm(int element, double a, double b)
		{
			/*InviscidTerm 2D array has 4 row 2 column:
			- column 1: Ox direction
			- column 2: Oy direction*/
			std::vector<std::vector<double>> InviscidTerm(4, std::vector<double>(2, 0.0));
			double
				rhoVal(math::pointValue(element, a, b, 1, 1)),
				uVal(math::pointValue(element, a, b, 2, 1)),
				vVal(math::pointValue(element, a, b, 3, 1)),
				eVal(math::pointValue(element, a, b, 4, 1)),
				pVal(math::pointValue(element, a, b, 5, 1));

			double totalE(eVal + 0.5*(uVal*uVal + vVal * vVal));

			/*1. Ox direction*/
			std::tie(InviscidTerm[0][0], InviscidTerm[1][0], InviscidTerm[2][0], InviscidTerm[3][0]) = math::inviscidTerms::calcInvisTermsFromPriVars(rhoVal, uVal, vVal, totalE, pVal, 1);

			/*2. Oy direction*/
			std::tie(InviscidTerm[0][1], InviscidTerm[1][1], InviscidTerm[2][1], InviscidTerm[3][1]) = math::inviscidTerms::calcInvisTermsFromPriVars(rhoVal, uVal, vVal, totalE, pVal, 2);

			return InviscidTerm;
		}

		/*Function calculates Viscous terms at Gauss point (a, b)*/
		std::vector<std::vector<double>> calcGaussViscousTerm(int element, double a, double b)
		{
			double uVal(math::pointValue(element, a, b, 2, 1)), vVal(math::pointValue(element, a, b, 3, 1));//, muVal(math::pointValue(element, a, b, 7, 1));
			/*ViscousTerm 2D array has 4 row 2 column:
			- column 1: Ox direction
			- column 2: Oy direction*/
			std::vector<std::vector<double>> ViscousTerm(4, std::vector<double>(2, 0.0));
			std::vector<std::vector<double>> StressHeatFlux(2, std::vector<double>(3, 0.0));
			std::vector<double> vectorU(4, 0.0);
			std::vector<double> vectordUx(4, 0.0);
			std::vector<double> vectordUy(4, 0.0);

			/*calculate conservative and derivative variables*/
			for (int i = 0; i < 4; i++)
			{
				vectorU[i] = math::pointValue(element, a, b, i + 1, 2);
				vectordUx[i] = math::pointAuxValue(element, a, b, i + 1, 1);
				vectordUy[i] = math::pointAuxValue(element, a, b, i + 1, 2);
			}

			/*calculate stresses and heat fluxes*/
			StressHeatFlux = math::viscousTerms::calcStressTensorAndHeatFlux(vectorU, vectordUx, vectordUy);
			std::tie(ViscousTerm[0][0], ViscousTerm[1][0], ViscousTerm[2][0], ViscousTerm[3][0]) = math::viscousTerms::calcViscousTermsFromStressHeatFluxMatrix(StressHeatFlux, uVal, vVal, 1);
			std::tie(ViscousTerm[0][1], ViscousTerm[1][1], ViscousTerm[2][1], ViscousTerm[3][1]) = math::viscousTerms::calcViscousTermsFromStressHeatFluxMatrix(StressHeatFlux, uVal, vVal, 2);
			return ViscousTerm;
		}

		void calcVolumeIntegralTerms(int element, std::vector<double> &VolIntTerm1, std::vector<double> &VolIntTerm2, std::vector<double> &VolIntTerm3, std::vector<double> &VolIntTerm4)
		{
			/*User's guide:
			All input array have form:
			- number of rows: orderElem*/

			std::vector<std::vector<double>> ViscousTerms(4, std::vector<double>(2, 0.0));
			std::vector<std::vector<double>> InviscidTerms(4, std::vector<double>(2, 0.0));
			//std::vector<double> VolInt(4, 0.0);

			std::vector<std::vector<double>>
				InvisGsVolX1(mathVar::nGauss + 1, std::vector<double>(mathVar::nGauss + 1, 0.0)),
				InvisGsVolY1(mathVar::nGauss + 1, std::vector<double>(mathVar::nGauss + 1, 0.0)),

				InvisGsVolX2(mathVar::nGauss + 1, std::vector<double>(mathVar::nGauss + 1, 0.0)),
				InvisGsVolY2(mathVar::nGauss + 1, std::vector<double>(mathVar::nGauss + 1, 0.0)),

				InvisGsVolX3(mathVar::nGauss + 1, std::vector<double>(mathVar::nGauss + 1, 0.0)),
				InvisGsVolY3(mathVar::nGauss + 1, std::vector<double>(mathVar::nGauss + 1, 0.0)),

				InvisGsVolX4(mathVar::nGauss + 1, std::vector<double>(mathVar::nGauss + 1, 0.0)),
				InvisGsVolY4(mathVar::nGauss + 1, std::vector<double>(mathVar::nGauss + 1, 0.0));

			std::vector<std::vector<double>>
				ViscGsVolX2(mathVar::nGauss + 1, std::vector<double>(mathVar::nGauss + 1, 0.0)),
				ViscGsVolY2(mathVar::nGauss + 1, std::vector<double>(mathVar::nGauss + 1, 0.0)),

				ViscGsVolX3(mathVar::nGauss + 1, std::vector<double>(mathVar::nGauss + 1, 0.0)),
				ViscGsVolY3(mathVar::nGauss + 1, std::vector<double>(mathVar::nGauss + 1, 0.0)),

				ViscGsVolX4(mathVar::nGauss + 1, std::vector<double>(mathVar::nGauss + 1, 0.0)),
				ViscGsVolY4(mathVar::nGauss + 1, std::vector<double>(mathVar::nGauss + 1, 0.0));

			double a(0.0), b(0.0);

			for (int na = 0; na <= mathVar::nGauss; na++)
			{
				for (int nb = 0; nb <= mathVar::nGauss; nb++)
				{
					std::tie(a, b) = auxUlti::getGaussCoor(na, nb);
					/*A INVISCID TERMS*/
					InviscidTerms = NSFEq::calcGaussInviscidTerm(element, a, b);
					/*A1. Inviscid term on Ox direction*/
					InvisGsVolX1[na][nb] = InviscidTerms[0][0];
					InvisGsVolX2[na][nb] = InviscidTerms[1][0];
					InvisGsVolX3[na][nb] = InviscidTerms[2][0];
					InvisGsVolX4[na][nb] = InviscidTerms[3][0];
					/*A2. Inviscid term on Oy direction*/
					InvisGsVolY1[na][nb] = InviscidTerms[0][1];
					InvisGsVolY2[na][nb] = InviscidTerms[1][1];
					InvisGsVolY3[na][nb] = InviscidTerms[2][1];
					InvisGsVolY4[na][nb] = InviscidTerms[3][1];

					/*B VISCOUS TERMS*/
					ViscousTerms = NSFEq::calcGaussViscousTerm(element, a, b);
					/*B1. Viscous term on Ox direction*/
					//ViscGsVolX1[na][nb] = ViscousTerms[0][0];
					ViscGsVolX2[na][nb] = ViscousTerms[1][0];
					ViscGsVolX3[na][nb] = ViscousTerms[2][0];
					ViscGsVolX4[na][nb] = ViscousTerms[3][0];
					/*B2. Viscous term on Oy direction*/
					//ViscGsVolY1[na][nb] = ViscousTerms[0][1];
					ViscGsVolY2[na][nb] = ViscousTerms[1][1];
					ViscGsVolY3[na][nb] = ViscousTerms[2][1];
					ViscGsVolY4[na][nb] = ViscousTerms[3][1];
				}
			}

			for (int order = 0; order <= mathVar::orderElem; order++)
			{
				/*CALCULATE INTEGRALS*/
				/*A INVISCID TERMS*/
				/*- integral term of (dB/dx)*(rho*u)*/
				VolIntTerm1[order] += process::volumeInte(element, InvisGsVolX1, order, 1);
				/*- integral term of (dB/dy)*(rho*v)*/
				VolIntTerm1[order] += process::volumeInte(element, InvisGsVolY1, order, 2);

				/*- integral term of (dB/dx)*(rho*u^2 + p)*/
				VolIntTerm2[order] += process::volumeInte(element, InvisGsVolX2, order, 1);
				/*- integral term of (dB/dy)*(rho*u*v)*/
				VolIntTerm2[order] += process::volumeInte(element, InvisGsVolY2, order, 2);

				/*- integral term of (dB/dx)*(rho*u*v)*/
				VolIntTerm3[order] += process::volumeInte(element, InvisGsVolX3, order, 1);
				/*- integral term of (dB/dy)*(rho*v^2 + p)*/
				VolIntTerm3[order] += process::volumeInte(element, InvisGsVolY3, order, 2);

				/*- integral term of (dB/dx)*(rho*totalE + p)*u*/
				VolIntTerm4[order] += process::volumeInte(element, InvisGsVolX4, order, 1);
				/*- integral term of (dB/dy)*(rho*totalE + p)*v*/
				VolIntTerm4[order] += process::volumeInte(element, InvisGsVolY4, order, 2);

				/*B VISCOUS TERMS*/
				/*- integral term of (dB/dx)*(0.0)*/
				VolIntTerm1[order] += 0.0;
				/*- integral term of (dB/dy)*(0.0)*/
				VolIntTerm1[order] += 0.0;

				/*- integral term of (dB/dx)*(tau_xx)*/
				VolIntTerm2[order] += process::volumeInte(element, ViscGsVolX2, order, 1);
				/*- integral term of (dB/dy)*(tau_xy)*/
				VolIntTerm2[order] += process::volumeInte(element, ViscGsVolY2, order, 2);

				/*- integral term of (dB/dx)*(tau_xy)*/
				VolIntTerm3[order] += process::volumeInte(element, ViscGsVolX3, order, 1);
				/*- integral term of (dB/dy)*(tau_yy)*/
				VolIntTerm3[order] += process::volumeInte(element, ViscGsVolY3, order, 2);

				/*- integral term of (dB/dx)*(u*tau_xx + v*tau_xy + Qx)*/
				VolIntTerm4[order] += process::volumeInte(element, ViscGsVolX4, order, 1);
				/*- integral term of (dB/dy)*(u*tau_xy + v*tau_yy + Qy)*/
				VolIntTerm4[order] += process::volumeInte(element, ViscGsVolY4, order, 2);
				/*End volume terms===========================================================================*/
			}
			//return VolInt;
		}

		void calcSurfaceIntegralTerms(int element, std::vector<double> &SurfIntTerm1, std::vector<double> &SurfIntTerm2, std::vector<double> &SurfIntTerm3, std::vector<double> &SurfIntTerm4)
		{
			/*User's guide:
			All input array have form:
			- number of rows: orderElem*/

			int elemType(auxUlti::checkType(element)), edgeName(0);
			int faceBcType(0);

			std::vector<double> inviscFlux1Temp(mathVar::nGauss + 1, 0.0),
				inviscFlux2Temp(mathVar::nGauss + 1, 0.0),
				inviscFlux3Temp(mathVar::nGauss + 1, 0.0),
				inviscFlux4Temp(mathVar::nGauss + 1, 0.0),
				
				ViscFlux1Temp(mathVar::nGauss + 1, 0.0),
				ViscFlux2Temp(mathVar::nGauss + 1, 0.0),
				ViscFlux3Temp(mathVar::nGauss + 1, 0.0),
				ViscFlux4Temp(mathVar::nGauss + 1, 0.0);

			std::vector<std::vector<double>> Fluxes(4, std::vector<double>(2, 0.0));

			std::vector<std::vector<double>>
				inviscFlux1(mathVar::nGauss + 1, std::vector<double>(4, 0.0)),
				inviscFlux2(mathVar::nGauss + 1, std::vector<double>(4, 0.0)),
				inviscFlux3(mathVar::nGauss + 1, std::vector<double>(4, 0.0)),
				inviscFlux4(mathVar::nGauss + 1, std::vector<double>(4, 0.0)),

				ViscFlux1(mathVar::nGauss + 1, std::vector<double>(4, 0.0)),
				ViscFlux2(mathVar::nGauss + 1, std::vector<double>(4, 0.0)),
				ViscFlux3(mathVar::nGauss + 1, std::vector<double>(4, 0.0)),
				ViscFlux4(mathVar::nGauss + 1, std::vector<double>(4, 0.0));

			//std::vector<double> SurInt(4, 0.0);

			for (int nface = 0; nface < elemType; nface++)
			{
				edgeName = meshVar::inedel[nface][element];
				faceBcType = auxUlti::getBCType(edgeName);

				if (faceBcType == 0)  //internal edge
				{
					for (int nGauss = 0; nGauss <= mathVar::nGauss; nGauss++)
					{
						Fluxes = process::NSFEq::getGaussVectorOfConserVarFluxesAtInternal(edgeName, element, nGauss);
						inviscFlux1[nGauss][nface] = Fluxes[0][0];
						inviscFlux2[nGauss][nface] = Fluxes[1][0];
						inviscFlux3[nGauss][nface] = Fluxes[2][0];
						inviscFlux4[nGauss][nface] = Fluxes[3][0];

						ViscFlux1[nGauss][nface] = Fluxes[0][1];
						ViscFlux2[nGauss][nface] = Fluxes[1][1];
						ViscFlux3[nGauss][nface] = Fluxes[2][1];
						ViscFlux4[nGauss][nface] = Fluxes[3][1];
					}
				}
				else  //boundary edge
				{
					for (int nGauss = 0; nGauss <= mathVar::nGauss; nGauss++)
					{
						Fluxes = NSFEqBCsImplement(element, edgeName, nGauss);
						inviscFlux1[nGauss][nface] = Fluxes[0][0];
						inviscFlux2[nGauss][nface] = Fluxes[1][0];
						inviscFlux3[nGauss][nface] = Fluxes[2][0];
						inviscFlux4[nGauss][nface] = Fluxes[3][0];

						ViscFlux1[nGauss][nface] = Fluxes[0][1];
						ViscFlux2[nGauss][nface] = Fluxes[1][1];
						ViscFlux3[nGauss][nface] = Fluxes[2][1];
						ViscFlux4[nGauss][nface] = Fluxes[3][1];
					}
				}
			}

			for (int order = 0; order <= mathVar::orderElem; order++)
			{
				for (int nface = 0; nface < elemType; nface++)
				{
					for (int nG = 0; nG <= mathVar::nGauss; nG++)
					{
						inviscFlux1Temp[nG] = inviscFlux1[nG][nface];
						inviscFlux2Temp[nG] = inviscFlux2[nG][nface];
						inviscFlux3Temp[nG] = inviscFlux3[nG][nface];
						inviscFlux4Temp[nG] = inviscFlux4[nG][nface];

						ViscFlux1Temp[nG] = ViscFlux1[nG][nface];
						ViscFlux2Temp[nG] = ViscFlux1[nG][nface];
						ViscFlux3Temp[nG] = ViscFlux1[nG][nface];
						ViscFlux4Temp[nG] = ViscFlux1[nG][nface];
					}
					SurfIntTerm1[order] += process::surfaceInte(element, edgeName, inviscFlux1Temp, order) + process::surfaceInte(element, edgeName, ViscFlux1Temp, order);
					SurfIntTerm2[order] += process::surfaceInte(element, edgeName, inviscFlux2Temp, order) + process::surfaceInte(element, edgeName, ViscFlux2Temp, order);
					SurfIntTerm3[order] += process::surfaceInte(element, edgeName, inviscFlux3Temp, order) + process::surfaceInte(element, edgeName, ViscFlux3Temp, order);
					SurfIntTerm4[order] += process::surfaceInte(element, edgeName, inviscFlux4Temp, order) + process::surfaceInte(element, edgeName, ViscFlux4Temp, order);
				}
			}
			//return SurInt;
		}

		std::vector<std::vector<double>> getGaussVectorOfConserVarFluxesAtInternal(int edgeName, int element, int nGauss)
		{
			std::vector<std::vector<double>> Fluxes(4, std::vector<double>(2, 0.0));
			double nx(auxUlti::getNormVectorComp(element, edgeName, 1)), ny(auxUlti::getNormVectorComp(element, edgeName, 2));
			double
				termX1P(0.0), termX1M(0.0),  //(rho*u)					or 0
				termX2P(0.0), termX2M(0.0),  //(rho*u^2 + p)			or tauxx
				termX3P(0.0), termX3M(0.0),  //(rho*u*v)				or tauxy
				termX4P(0.0), termX4M(0.0),  //(rho*totalE + p)*u		or tauxx*u + tauxy*v + Qx

				termY1P(0.0), termY1M(0.0),  //(rho*v)					or 0
				termY2P(0.0), termY2M(0.0),  //(rho*u*v)				or tauxy
				termY3P(0.0), termY3M(0.0),  //(rho*v^2 + p)			or tauyy
				termY4P(0.0), termY4M(0.0);  //(rho*totalE + p)*v		or tauxy*u + tauyy*v + Qy
			double rhoPlus(0.0), rhouPlus(0.0), rhovPlus(0.0), rhoEPlus(0.0), rhoMinus(0.0), rhouMinus(0.0), rhovMinus(0.0), rhoEMinus(0.0);

			/*INVISCID TERM*/
			//Get value
			std::tie(termX1P, termX1M) = process::NSFEq::getInternalValuesFromCalculatedArrays(edgeName, element, nGauss, 1, 1, 1);
			std::tie(termX2P, termX2M) = process::NSFEq::getInternalValuesFromCalculatedArrays(edgeName, element, nGauss, 1, 1, 2);
			std::tie(termX3P, termX3M) = process::NSFEq::getInternalValuesFromCalculatedArrays(edgeName, element, nGauss, 1, 1, 3);
			std::tie(termX4P, termX4M) = process::NSFEq::getInternalValuesFromCalculatedArrays(edgeName, element, nGauss, 1, 1, 4);

			std::tie(termY1P, termY1M) = process::NSFEq::getInternalValuesFromCalculatedArrays(edgeName, element, nGauss, 1, 2, 1);
			std::tie(termY2P, termY2M) = process::NSFEq::getInternalValuesFromCalculatedArrays(edgeName, element, nGauss, 1, 2, 2);
			std::tie(termY3P, termY3M) = process::NSFEq::getInternalValuesFromCalculatedArrays(edgeName, element, nGauss, 1, 2, 3);
			std::tie(termY4P, termY4M) = process::NSFEq::getInternalValuesFromCalculatedArrays(edgeName, element, nGauss, 1, 2, 4);

			std::tie(rhoPlus, rhoMinus) = process::getInternalValuesFromCalculatedArrays(edgeName, element, nGauss, 1);
			std::tie(rhouPlus, rhouMinus) = process::getInternalValuesFromCalculatedArrays(edgeName, element, nGauss, 2);
			std::tie(rhovPlus, rhovMinus) = process::getInternalValuesFromCalculatedArrays(edgeName, element, nGauss, 3);
			std::tie(rhoEPlus, rhoEMinus) = process::getInternalValuesFromCalculatedArrays(edgeName, element, nGauss, 4);

			/*Calculate fluxes*/
			Fluxes[0][0] = math::numericalFluxes::advectiveFlux(termX1P, termX1M, rhoPlus, rhoMinus, LxFConst[edgeName], nx) + math::numericalFluxes::advectiveFlux(termY1P, termY1M, rhoPlus, rhoMinus, LxFConst[edgeName], ny);
			Fluxes[1][0] = math::numericalFluxes::advectiveFlux(termX2P, termX2M, rhouPlus, rhouMinus, LxFConst[edgeName], nx) + math::numericalFluxes::advectiveFlux(termY2P, termY2M, rhouPlus, rhouMinus, LxFConst[edgeName], ny);
			Fluxes[2][0] = math::numericalFluxes::advectiveFlux(termX3P, termX3M, rhovPlus, rhovMinus, LxFConst[edgeName], nx) + math::numericalFluxes::advectiveFlux(termY3P, termY3M, rhovPlus, rhovMinus, LxFConst[edgeName], ny);
			Fluxes[3][0] = math::numericalFluxes::advectiveFlux(termX4P, termX4M, rhoEPlus, rhoEMinus, LxFConst[edgeName], nx) + math::numericalFluxes::advectiveFlux(termY4P, termY4M, rhoEPlus, rhoEMinus, LxFConst[edgeName], ny);
			
			/*VISCOUS TERM*/
			//Get value
			std::tie(termX1P, termX1M) = process::NSFEq::getInternalValuesFromCalculatedArrays(edgeName, element, nGauss, 2, 1, 1);
			std::tie(termX2P, termX2M) = process::NSFEq::getInternalValuesFromCalculatedArrays(edgeName, element, nGauss, 2, 1, 2);
			std::tie(termX3P, termX3M) = process::NSFEq::getInternalValuesFromCalculatedArrays(edgeName, element, nGauss, 2, 1, 3);
			std::tie(termX4P, termX4M) = process::NSFEq::getInternalValuesFromCalculatedArrays(edgeName, element, nGauss, 2, 1, 4);

			std::tie(termY1P, termY1M) = process::NSFEq::getInternalValuesFromCalculatedArrays(edgeName, element, nGauss, 2, 2, 1);
			std::tie(termY2P, termY2M) = process::NSFEq::getInternalValuesFromCalculatedArrays(edgeName, element, nGauss, 2, 2, 2);
			std::tie(termY3P, termY3M) = process::NSFEq::getInternalValuesFromCalculatedArrays(edgeName, element, nGauss, 2, 2, 3);
			std::tie(termY4P, termY4M) = process::NSFEq::getInternalValuesFromCalculatedArrays(edgeName, element, nGauss, 2, 2, 4);

			/*Calculate fluxes*/
			Fluxes[0][1] = math::numericalFluxes::diffusiveFlux(termX1M, termX1P, nx) + math::numericalFluxes::diffusiveFlux(termY1M, termY1P, ny);
			Fluxes[1][1] = math::numericalFluxes::diffusiveFlux(termX2M, termX2P, nx) + math::numericalFluxes::diffusiveFlux(termY2M, termY2P, ny);
			Fluxes[2][1] = math::numericalFluxes::diffusiveFlux(termX3M, termX3P, nx) + math::numericalFluxes::diffusiveFlux(termY3M, termY3P, ny);
			Fluxes[3][1] = math::numericalFluxes::diffusiveFlux(termX4M, termX4P, nx) + math::numericalFluxes::diffusiveFlux(termY4M, termY4P, ny);

			return Fluxes;
		}

		std::vector<double> solveTimeMarching(std::vector<double> &ddtArr, std::vector<double> &UnArr)
		{
			std::vector<double> OutArr(mathVar::orderElem + 1, 0.0);

			if (systemVar::ddtScheme==1) //Euler scheme
			{
				for (int order = 0; order <= mathVar::orderElem; order++)
				{
					OutArr[order] = dt * ddtArr[order] + UnArr[order];
				}
			}
			return OutArr;
		}
	}

	namespace Euler
	{
		void calcGlobalTimeStep()
		{
			std::vector<double>timeStepArr(meshVar::nelem2D, 1.0);
			for (int nelement = 0; nelement < meshVar::nelem2D; nelement++)
			{
				timeStepArr[nelement] = process::Euler::localTimeStep(nelement);
				
			}
			runTime += dt;
			dt = *std::min_element(timeStepArr.begin(), timeStepArr.end());  //find min value of vector
		}

		double localTimeStep(int element)
		{
			std::vector<double> vectorDeltaT;
			double deltaT(0.0), uVal(0.0), vVal(0.0), velocity(0.0), TVal(0.0), aSound(0.0), LocalMach(0.0), muVal(0.0),
				aG(0.0), bG(0.0), rhoVal(0.0), rhouVal(0.0), rhovVal(0.0), rhoEVal(0.0), size(meshVar::cellSize[element]);

			/*
			if (TVal<=0 || TVal != TVal) //
			{
				std::cout << "Unphysical T is detected at element " << element + meshVar::nelem1D + 1 << std::endl;
				//TVal = iniValues::TIni;
				//system("pause");
				double rhoVal(math::pointValue(element, aC, bC, 1, 2)),
					rhouVal(math::pointValue(element, aC, bC, 2, 2)),
					rhovVal(math::pointValue(element, aC, bC, 3, 2)),
					rhoEVal(math::pointValue(element, aC, bC, 4, 2));
				//std::cout << "value rho, rhou, rhov, rhoe = " << rhoVal << ", " << rhouVal << ", " << rhovVal << ", " << rhoEVal << std::endl;
			}
			*/
			for (int na = 0; na <= mathVar::nGauss; na++)
			{
				for (int nb = 0; nb <= mathVar::nGauss; nb++)
				{
					aG = mathVar::GaussPts[na][nb][0];
					bG = mathVar::GaussPts[na][nb][1];
					rhoVal = math::pointValue(element, aG, bG, 1, 2);
					rhouVal = math::pointValue(element, aG, bG, 2, 2);
					rhovVal = math::pointValue(element, aG, bG, 3, 2);
					rhoEVal = math::pointValue(element, aG, bG, 4, 2);

					uVal = rhouVal / rhoVal;
					vVal = rhovVal / rhoVal;
					velocity = sqrt(pow(uVal, 2) + pow(vVal, 2));
					TVal = math::CalcTFromConsvVar(rhoVal, rhouVal, rhovVal, rhoEVal);
					aSound = math::CalcSpeedOfSound(TVal);
					LocalMach = velocity / aSound;
					muVal = math::CalcVisCoef(TVal);
					vectorDeltaT.push_back((1.0 / pow((mathVar::orderElem + 1 + 1), 2))*(size*systemVar::CFL) / (fabs(velocity) + (aSound / LocalMach) + (muVal / size)));
				}
			}
			//deltaT = (1.0 / pow((mathVar::orderElem + 1 + 1), 2))*(size*systemVar::CFL) / (fabs(velocity) + (aSound / LocalMach) + (muVal / size));
			deltaT = *std::min_element(vectorDeltaT.begin(), vectorDeltaT.end());
			return deltaT;
		}

		double localErrorEstimate(int element, std::vector<double> &ddtArr)
		{
			double xC(-1.0 / 3.0), yC(1.0 / 3.0), errorVal(0.0);

			if (auxUlti::checkType(element) == 4)
			{
				xC = 0.0;
				yC = 0.0;
			}

			//Compute basis function
			math::basisFc(xC, yC);
			for (int order = 0; order <= mathVar::orderElem; order++)
			{
				errorVal += ddtArr[order] * mathVar::B[order];
			}

			return errorVal;
		}

		std::tuple <double, double, double, double> globalErrorEstimate(std::vector<double> &RhoError, std::vector<double> &RhouError, std::vector<double> &RhovError, std::vector<double> &RhoEError)
		{
			double rhoRes(*std::max_element(RhoError.begin(), RhoError.end())),
				rhouRes(*std::max_element(RhouError.begin(), RhouError.end())),
				rhovRes(*std::max_element(RhovError.begin(), RhovError.end())),
				rhoERes(*std::max_element(RhoEError.begin(), RhoEError.end()));
			return std::make_tuple(rhoRes, rhouRes, rhovRes, rhoERes);
		}
	}

	std::tuple<double, double> getInternalValuesFromCalculatedArrays(int edge, int element, int nG, int valType)
	{
		//valType = 1(rho*mu), 2(rhou*mu), 3(rhov*mu), 4(rhoE*mu)
		bool isMaster(auxUlti::checkMaster(element, edge));
		int locationPlus(-1), locationMinus(-1);
		double valPlus(0.0), valMinus(0.0);
		if (isMaster)
		{
			locationPlus = nG;
			locationMinus = nG + mathVar::nGauss + 1;
		}
		else
		{
			locationPlus = nG + mathVar::nGauss + 1;
			locationMinus = nG;
		}

		switch (valType)
		{
		case 1:
		{
			valPlus = interface_rho[edge][locationPlus];
			valMinus = interface_rho[edge][locationMinus];
		}
		break;
		case 2:
		{
			valPlus = interface_rhou[edge][locationPlus];
			valMinus = interface_rhou[edge][locationMinus];
		}
		break;
		case 3:
		{
			valPlus = interface_rhov[edge][locationPlus];
			valMinus = interface_rhov[edge][locationMinus];
		}
		break;
		case 4:
		{
			valPlus = interface_rhoE[edge][locationPlus];
			valMinus = interface_rhoE[edge][locationMinus];
		}
		break;
		default:
			break;
		}
		return std::make_tuple(valPlus, valMinus);
	}

	namespace limiter
	{
		void limiter()
		{
			if (systemVar::limiter==1)  //positivity preserving
			{
				for (int nelem = 0; nelem < meshVar::nelem2D; nelem++)
				{
					std::tie(theta1Arr[nelem], theta2Arr[nelem]) = limiter::Pp::quadratureCell::calcPpLimiterCoef(nelem);
				}
				if (limitVal::numOfLimitCell>0)
				{
					std::cout << "Limiter is applied at " << limitVal::numOfLimitCell << " cell(s)\n" ;
				}
			}
			else if (systemVar::limiter == 0)  //No limiter
			{
				for (int nelem = 0; nelem < meshVar::nelem2D; nelem++)
				{
					theta1Arr[nelem] = 1.0;
					theta2Arr[nelem] = 1.0;
				}
			}
		}

		namespace Pp
		{
			namespace triangleCell
			{
				//Function calculates coefficients of positivity preserving limiter
				std::tuple<double, double> calcPpLimiterCoef(int element)
				{
					double meanRho(0.0), minRho(0.0), theta1(0.0), theta2(0.0), omega(0.0);
					int elemType(auxUlti::checkType(element));

					double meanRhou(0.0), meanRhov(0.0), meanRhoE(0.0);

					//Find theta1
					minRho = math::limiter::triangleCell::calcMinRho(element);
					meanRho = rho[element][0];
					meanRhou = rhou[element][0];
					meanRhov = rhov[element][0];
					meanRhoE = rhoE[element][0];
					double meanT(math::CalcTFromConsvVar(meanRho, meanRhou, meanRhov, meanRhoE));
					double meanP(math::CalcP(meanT, meanRho));

					//Compute theta1
					std::tie(theta1, omega) = math::limiter::triangleCell::calcTheta1Coeff(meanRho, minRho, meanP);

					//Find theta2
					theta2 = math::limiter::triangleCell::calcTheta2Coeff(element, theta1, omega);

					//Reset limit flag
					limitVal::limitFlagLocal = false;
					return std::make_tuple(theta1, theta2);
				}
			}

			namespace quadratureCell
			{
				//Function calculates coefficients of positivity preserving limiter
				std::tuple<double, double> calcPpLimiterCoef(int element)
				{
					double meanRho(0.0), minRho(0.0), theta1(0.0), theta2(0.0), omega(0.0);
					int elemType(auxUlti::checkType(element));

					double meanRhou(0.0), meanRhov(0.0), meanRhoE(0.0), aG(0.0), bG(0.0);

					/*Note: according to Kontzialis et al, positivity preserving limiter for quadrilateral element, which is presented on Zhang's paper,
					shown a very good effect on results. Because of that, Zhang's limiter is used in this code for both triangular and quadrilateral elements*/

					//Find theta1
					minRho = math::limiter::quadratureCell::calcMinRhoQuad(element);
					meanRho = rho[element][0];

					meanRhou = rhou[element][0];
					meanRhov = rhov[element][0];
					meanRhoE = rhoE[element][0];
					double meanT(math::CalcTFromConsvVar(meanRho, meanRhou, meanRhov, meanRhoE));
					double meanP(math::CalcP(meanT, meanRho));

					//Compute theta1
					std::tie(theta1, omega) = math::limiter::quadratureCell::calcTheta1Coeff(meanRho, minRho, meanP);

					//Find theta2
					std::vector<double> vectort;
					/*
					int index(0);

					for (int na = 0; na <= mathVar::nGauss; na++)
					{
						for (int nb = 0; nb <= mathVar::nGauss; nb++)
						{
							vectort[index] = math::limiter::quadratureCell::calcTheta2Coeff(element, na, nb, theta1, omega, meanRho, meanRhou, meanRhov, meanRhoE, 1);
							index++;

							vectort[index] = math::limiter::quadratureCell::calcTheta2Coeff(element, na, nb, theta1, omega, meanRho, meanRhou, meanRhov, meanRhoE, 2);
							index++;
						}
					}
					*/

					//Compute t at all internal Gauss point
					for (int na = 0; na <= mathVar::nGauss; na++)
					{
						for (int nb = 0; nb <= mathVar::nGauss; nb++)
						{
							aG = mathVar::GaussPts[na][nb][0];
							bG = mathVar::GaussPts[na][nb][1];
							vectort.push_back(math::limiter::quadratureCell::calcTheta2Coeff(element, aG, bG, theta1, omega, meanRho, meanRhou, meanRhov, meanRhoE));
						}
					}
					//Compute t at edge DA
					aG = -1;
					for (int nG = 0; nG <= mathVar::nGauss; nG++)
					{
						bG = mathVar::xGauss[nG];
						vectort.push_back(math::limiter::quadratureCell::calcTheta2Coeff(element, aG, bG, theta1, omega, meanRho, meanRhou, meanRhov, meanRhoE));
					}
					//Compute t at edge BC
					aG = 1;
					for (int nG = 0; nG <= mathVar::nGauss; nG++)
					{
						bG = mathVar::xGauss[nG];
						vectort.push_back(math::limiter::quadratureCell::calcTheta2Coeff(element, aG, bG, theta1, omega, meanRho, meanRhou, meanRhov, meanRhoE));
					}
					//Compute t at edge AB
					bG = -1;
					for (int nG = 0; nG <= mathVar::nGauss; nG++)
					{
						aG = mathVar::xGauss[nG];
						vectort.push_back(math::limiter::quadratureCell::calcTheta2Coeff(element, aG, bG, theta1, omega, meanRho, meanRhou, meanRhov, meanRhoE));
					}
					//Compute t at edge CD
					bG = 1;
					for (int nG = 0; nG <= mathVar::nGauss; nG++)
					{
						aG = mathVar::xGauss[nG];
						vectort.push_back(math::limiter::quadratureCell::calcTheta2Coeff(element, aG, bG, theta1, omega, meanRho, meanRhou, meanRhov, meanRhoE));
					}

					theta2 = *std::min_element(vectort.begin(), vectort.end());  //find min value of vector
					if (limitVal::limitFlagLocal == true)
					{
						limitVal::numOfLimitCell++;
					}
					//Reset limit flag
					limitVal::limitFlagLocal = false;
					return std::make_tuple(theta1, theta2);
				}
			}
		}

	}

	double volumeInte(int elem, std::vector< std::vector<double> > &Ui, int order, int direction)
	{
		double Int(0.0);
		std::vector<std::vector<double>> A(mathVar::nGauss + 1, std::vector<double>(mathVar::nGauss + 1, 0.0));
		double dBi(0.0);

		for (int na = 0; na <= mathVar::nGauss; na++)
		{
			for (int nb = 0; nb <= mathVar::nGauss; nb++)
			{
				dBi = math::Calc_dBxdBy(elem, order, na, nb, direction);
				A[na][nb] = dBi * Ui[na][nb];
			}
		}
		Int = math::volumeInte(A, elem);

		return Int;
	}

	double surfaceInte(int elem, int edge, std::vector<double> &FluxVector, int order)
	{
		double Bi(0.0), a(0.0), b(0.0), inte(0.0);
		std::vector<double> F(mathVar::nGauss + 1, 0.0);

		//inte = 0.0;  //value of integral is reset for each order
		for (int nG = 0; nG <= mathVar::nGauss; nG++)
		{
			std::tie(a, b) = auxUlti::getGaussSurfCoor(edge, elem, nG);
			math::basisFc(a, b);  //mathVar::B is changed after this command line is excuted
			Bi = mathVar::B[order];
			F[nG] = Bi * FluxVector[nG];
		}
		inte = math::surfaceInte(F, edge, elem);

		return inte;
	}

	std::vector<std::vector<double>> calculateStiffMatrix(int element)
	{
		std::vector<std::vector<double>> StiffMatrix(mathVar::orderElem + 1, std::vector<double>(mathVar::orderElem + 1, 0.0));

		for (int order1 = 0; order1 <= mathVar::orderElem; order1++)
		{
			for (int order2 = 0; order2 <= mathVar::orderElem; order2++)
			{
				StiffMatrix[order1][order2] = process::calculateStiffMatrixElement(element, order1, order2);
			}
		}
		return StiffMatrix;
	}

	double calculateStiffMatrixElement(int element, int order1, int order2)
	{
		double B1(0.0), B2(0.0), Inte(0.0);
		std::vector<std::vector<double>> FMatrix(mathVar::nGauss + 1, std::vector<double>(mathVar::nGauss + 1, 0.0));
		for (int na = 0; na <= mathVar::nGauss; na++)
		{
			for (int nb = 0; nb <= mathVar::nGauss; nb++)
			{
				B1 = mathVar::BPts[order1][na][nb];
				B2 = mathVar::BPts[order2][na][nb];
				FMatrix[na][nb] = B1 * B2;
			}
		}
		Inte = math::volumeInte(FMatrix, element);
		return Inte;
	}

	bool checkRunningCond()
	{
		bool run(true);
		if (runTime >= systemVar::Ttime)
		{
			run = false;  //stop running if runtime bigger than Ttime (total time)
		}
		return run;
	}
}