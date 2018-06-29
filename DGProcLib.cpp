#include "DGProcLib.h"
#include "DGMath.h"
#include "varDeclaration.h"
//#include <math.h>
#include "DGAuxUltilitiesLib.h"
#include <tuple>
#include "DGBCsLib.h"

namespace meshParam
{
	void GaussParam()
	{
		math::Gauss(mathVar::nGauss);
		for (int na = 0; na <= mathVar::nGauss; na++)  //nGauss is started from 0
		{
			for (int nb = 0; nb <= mathVar::nGauss; nb++)
			{
				mathVar::GaussPts[na][nb][0] = mathVar::xGauss[na];
				mathVar::GaussPts[na][nb][1] = mathVar::xGauss[nb];

				mathVar::wGaussPts[na][nb][0] = mathVar::wGauss[na];
				mathVar::wGaussPts[na][nb][1] = mathVar::wGauss[nb];
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
				math::basisFc(a, b, mathVar::orderElem);
				math::dBasisFc(a, b, mathVar::orderElem);
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
}

namespace process
{
	void setIniValues()
	{
		iniValues::rhoIni = iniValues::pIni / (material::R*iniValues::TIni);
		material::Cp = material::R*material::gamma / (material::gamma - 1);
		material::Cv = material::Cp - material::R;
		iniValues::eIni = material::Cv*iniValues::TIni + 0.5*(pow(iniValues::uIni, 2) + pow(iniValues::vIni, 2) + pow(iniValues::wIni, 2));
		iniValues::muIni = math::CalcVisCoef(iniValues::TIni);

		/*rho*/
		std::vector<double> iniRho(mathVar::orderElem);
		iniRho = process::calcIniValues(iniValues::rhoIni);
		/*rhou*/
		std::vector<double> iniRhou(mathVar::orderElem);
		iniRhou = process::calcIniValues(iniValues::rhoIni*iniValues::uIni);
		/*rhov*/
		std::vector<double> iniRhov(mathVar::orderElem);
		iniRhov = process::calcIniValues(iniValues::rhoIni*iniValues::vIni);
		/*rhoE*/
		std::vector<double> iniRhoE(mathVar::orderElem);
		iniRhoE = process::calcIniValues(iniValues::rhoIni*iniValues::eIni);

		for (int nelem = 0; nelem < meshVar::nelem2D; nelem++)
		{
			for (int i = 0; i <= mathVar::orderElem; i++)
			{
				rho[nelem][i] = iniRho[i];
				rhou[nelem][i] = iniRhou[i];
				rhov[nelem][i] = iniRhov[i];
				rhoE[nelem][i] = iniRhoE[i];
			}
		}
	}

	std::vector<double> calcIniValues(double iniVal)
	{
		
		//std::vector<std::vector<double>> a(mathVar::orderElem + 1, std::vector<double>(mathVar::orderElem + 1, iniVal));
		//std::vector<double> b(mathVar::orderElem + 1, iniVal);
		std::vector<double> iniVector(mathVar::orderElem + 1, 0.0);
		double I(0.0);
		
		//int ptCount(0);
		/*
		for (int na = 0; na <= mathVar::orderElem; na++)
		{
			for (int nb = 0; nb <= mathVar::orderElem; nb++)
			{
				if (ptCount<=mathVar::orderElem)
				{
					for (int i = 0; i <= mathVar::orderElem; i++)
					{
						a[ptCount][i] = mathVar::BPts[i][na][nb];
					}
					ptCount++;
				}
				else
				{
					break;
				}
			}
		}

		for (int j = 0; j <= mathVar::orderElem; j++)
		{
			b[j] = iniVal;
		}

		iniVector = math::SolveSysEqs(a, b);
		*/
		for (int order = 0; order <= mathVar::orderElem; order++)
		{
			I = math::iniIntegral(order);
			iniVector[order] = I * iniVal;
		}
		return iniVector;
	}

	namespace auxEq
	{
		void solveAuxEquation()
		{
			std::vector<std::vector<double>> StiffMatrix(mathVar::orderElem + 1, std::vector<double>(mathVar::orderElem + 1, 0.0));
			std::vector<double> RHSOrderX(4, 0.0), RHSOrderY(4, 0.0);
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
				
				for (int order = 0; order <= mathVar::orderElem; order++)
				{
					//At each order, right hand side terms of all conservative variables are calculated

					//2) Calculate right hand side terms
					//2.1) Ox direction
					RHSOrderX = process::auxEq::CalcRHSTerm(nelement, order, 1);
					rhoRHSTermOxDir[order] = RHSOrderX[0];
					rhouRHSTermOxDir[order] = RHSOrderX[1];
					rhovRHSTermOxDir[order] = RHSOrderX[2];
					rhoERHSTermOxDir[order] = RHSOrderX[3];

					//2.2) Oy direction
					RHSOrderY = process::auxEq::CalcRHSTerm(nelement, order, 2);
					rhoRHSTermOyDir[order] = RHSOrderY[0];
					rhouRHSTermOyDir[order] = RHSOrderY[1];
					rhovRHSTermOyDir[order] = RHSOrderY[2];
					rhoERHSTermOyDir[order] = RHSOrderY[3];

					//4 vector RHS of 4 conservative variables for each direction is achieved after these steps have been finished
				}

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

		/*Function calculates right hand side terms of all conservative variables at ALL one order*/
		std::vector<std::vector<double>> CalcRHSTerm(int element, int dir)
		{
<<<<<<< HEAD
<<<<<<< HEAD
			/*NOTE: AUXILARY EQUATION:
			volumeInt(S_term) + volumeInt(U_term) - surfaceInt(U_term) = 0
			<=> volumeInt(S_term) = - volumeInt(U_term) + surfaceInt(U_term)
								   |				RHS term			   |*/
			
			std::vector<std::vector<double>> RHS(4, std::vector<double>(mathVar::orderElem + 1, 0.0));
			std::vector<std::vector<double>> SurfaceInt(4, std::vector<double>(mathVar::orderElem + 1, 0.0));
			//std::vector<double> VolInt(4, 0.0);
			//std::vector<double> SurInt(4, 0.0);
=======
			std::vector<double> VolInt(4, 0.0);
			std::vector<double> SurInt(4, 0.0);
>>>>>>> parent of 622c1e9... hoàn thành phần tính RHS, chờ công thức tính time step
=======
			std::vector<double> VolInt(4, 0.0);
			std::vector<double> SurInt(4, 0.0);
>>>>>>> parent of 622c1e9... hoàn thành phần tính RHS, chờ công thức tính time step
			std::vector<double> RHS(4, 0.0);
			std::vector<std::vector<double>> rhoGsVol(mathVar::nGauss + 1, std::vector<double>(mathVar::nGauss + 1, 0.0));
			std::vector<std::vector<double>> rhouGsVol(mathVar::nGauss + 1, std::vector<double>(mathVar::nGauss + 1, 0.0));
			std::vector<std::vector<double>> rhovGsVol(mathVar::nGauss + 1, std::vector<double>(mathVar::nGauss + 1, 0.0));
			std::vector<std::vector<double>> rhoEGsVol(mathVar::nGauss + 1, std::vector<double>(mathVar::nGauss + 1, 0.0));

			std::vector<double> Flux(4, 0.0);
			std::vector<double> rhoFlux(mathVar::nGauss + 1, 0.0);
			std::vector<double> rhouFlux(mathVar::nGauss + 1, 0.0);
			std::vector<double> rhovFlux(mathVar::nGauss + 1, 0.0);
			std::vector<double> rhoEFlux(mathVar::nGauss + 1, 0.0);
			int elemType(auxUlti::checkType(element)), edgeName(0);

			int faceBcType(0);
			double nVectorComp(0.0);

			/*Volume integral term===========================================================================*/
			//rho -------------------------------------------------------------------------------------------
			rhoGsVol = process::auxEq::getGaussMatrixOfConserVar(element, 1);
			VolInt[0] = process::volumeInte(element, rhoGsVol, order, dir);
			//rhou ------------------------------------------------------------------------------------------
			rhouGsVol = process::auxEq::getGaussMatrixOfConserVar(element, 2);
			VolInt[1] = process::volumeInte(element, rhouGsVol, order, dir);
			//rhov ------------------------------------------------------------------------------------------
			rhovGsVol = process::auxEq::getGaussMatrixOfConserVar(element, 3);
			VolInt[2] = process::volumeInte(element, rhovGsVol, order, dir);
			//rhou ------------------------------------------------------------------------------------------
			rhoEGsVol = process::auxEq::getGaussMatrixOfConserVar(element, 4);
			VolInt[3] = process::volumeInte(element, rhoEGsVol, order, dir);

			/*Surface integral term==========================================================================*/
			SurfaceInt = process::auxEq::calcSurfaceIntegralTerms(element, dir);

			for (int order = 0; order <= mathVar::orderElem; order++)
			{
				for (int i = 0; i < 4; i++)
				{
					RHS[0][order] += SurfaceInt[0][order];
				}
<<<<<<< HEAD
<<<<<<< HEAD
			}
=======
=======
>>>>>>> parent of 622c1e9... hoàn thành phần tính RHS, chờ công thức tính time step
				SurInt[0] += process::surfaceInte(element, edgeName, rhoFlux, order);
				SurInt[1] += process::surfaceInte(element, edgeName, rhouFlux, order);
				SurInt[2] += process::surfaceInte(element, edgeName, rhovFlux, order);
				SurInt[3] += process::surfaceInte(element, edgeName, rhoEFlux, order);
			}
			RHS = math::vectorSum(VolInt, SurInt);
<<<<<<< HEAD
>>>>>>> parent of 622c1e9... hoàn thành phần tính RHS, chờ công thức tính time step
=======
>>>>>>> parent of 622c1e9... hoàn thành phần tính RHS, chờ công thức tính time step
			return RHS;
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
					GaussMatrix[na][nb] = math::pointValue(element, a, b, 1, 2)*muGs;
				}
			}
			return GaussMatrix;
		}

		std::vector<double> getVectorOfConserVarFluxesAtInternal(int edge, int element, int nG, double nVectorComp)
		{
			std::vector<double> Flux(4, 0.0);
			double muGsP(0.0), muGsM(0.0), valP(0.0), valM(0.0);
			std::tie(muGsP, muGsM) = math::internalSurfaceValue(edge, element, nG, 7, 1);

			//rho flux
			std::tie(valP, valM) = math::internalSurfaceValue(edge, element, nG, 1, 2);
			Flux[0] = math::numericalFluxes::auxFlux(muGsM*valM, muGsP*valP, nVectorComp);

			//rhou flux
			std::tie(valP, valM) = math::internalSurfaceValue(edge, element, nG, 2, 2);
			Flux[1] = math::numericalFluxes::auxFlux(muGsM*valM, muGsP*valP, nVectorComp);

			//rhov flux
			std::tie(valP, valM) = math::internalSurfaceValue(edge, element, nG, 3, 2);
			Flux[2] = math::numericalFluxes::auxFlux(muGsM*valM, muGsP*valP, nVectorComp);

			//rhoE flux
			std::tie(valP, valM) = math::internalSurfaceValue(edge, element, nG, 4, 2);
			Flux[3] = math::numericalFluxes::auxFlux(muGsM*valM, muGsP*valP, nVectorComp);

			return Flux;
		}

		std::vector< std::vector<double> > getGaussVectorOfConserVarFlux(int element, int nGauss, int dir)
		{
			int elemType(auxUlti::checkType(element)), edgeName(0);

			/*Function calculates flux of ALL conservative variables at ONE Gauss point on ALL surfaces of element
			output array has form:
			- 4 rows: 4 fluxes of conservative variables at input nGauss
			- n column: number of faces of element*/
			std::vector<std::vector<double>> output(4, std::vector<double>(elemType, 0.0));

			int faceBcType(0);
			double nVectorComp(0.0);
			std::vector<double> Flux(4, 0.0);

			for (int nface = 0; nface < elemType; nface++)
			{
				edgeName = meshVar::inedel[nface][element];

				faceBcType = auxUlti::getBCType(edgeName);
				nVectorComp = (auxUlti::getNormVectorComp(element, edgeName, dir));

				if (faceBcType == 0)  //internal edge
				{
					Flux = process::auxEq::getVectorOfConserVarFluxesAtInternal(edgeName, element, nGauss, nVectorComp);
				}
				else  //boundary edge
				{
					Flux = auxEqBCsImplement(element, edgeName, nGauss, nVectorComp);
				}
				output[0][nface] = Flux[0];
				output[1][nface] = Flux[1];
				output[2][nface] = Flux[2];
				output[3][nface] = Flux[3];
			}
			return output;
		}

		std::vector<std::vector<double>> calcSurfaceIntegralTerms(int element, int dir)
		{
			int elemType(auxUlti::checkType(element)), edgeName(0);

			/*Surface integral array has form :
			- 4 rows : 4 fluxes of conservative variables at input nGauss
			- orderElem column : containt surface integrals at all order of accuracy*/
			std::vector<std::vector<double>> SurfInt(4, std::vector<double>(mathVar::orderElem + 1, 0.0));

			/*Fluxes array has form :
			- 4 rows : 4 fluxes of conservative variables at input nGauss
			- n column : number of faces of element*/
			std::vector<std::vector<double>> Fluxes(4, std::vector<double>(elemType, 0.0));

			/*rhoFlux, etc .. has the following form:
			- nGauss+1 rows: contain flux at all Gauss points
			- 4 columns: contain flux at all faces of element (default is 4 faces)*/
			std::vector<std::vector<double>> rhoFlux(mathVar::nGauss + 1, std::vector<double>(elemType, 0.0)),
				rhouFlux(mathVar::nGauss + 1, std::vector<double>(elemType, 0.0)),
				rhovFlux(mathVar::nGauss + 1, std::vector<double>(elemType, 0.0)),
				rhoEFlux(mathVar::nGauss + 1, std::vector<double>(elemType, 0.0));

			std::vector<double> rhoFluxFace(mathVar::nGauss + 1, 0.0),
				rhouFluxFace(mathVar::nGauss + 1, 0.0),
				rhovFluxFace(mathVar::nGauss + 1, 0.0),
				rhoEFluxFace(mathVar::nGauss + 1, 0.0);

			for (int nG = 0; nG <= mathVar::nGauss; nG++)
			{
				Fluxes = process::auxEq::getGaussVectorOfConserVarFlux(element, nG, dir);

				//distributes flux data to flux array
				for (int nface = 0; nface < elemType; nface++)
				{
					rhoFlux[nG][nface] = Fluxes[0][nface];
					rhouFlux[nG][nface] = Fluxes[1][nface];
					rhovFlux[nG][nface] = Fluxes[2][nface];
					rhoEFlux[nG][nface] = Fluxes[3][nface];
				}
			}

			for (int order = 0; order <= mathVar::orderElem; order++)
			{
				for (int nface = 0; nface < elemType; nface++)
				{
					for (int nG = 0; nG <= mathVar::nGauss; nG++)
					{
						rhoFluxFace[nG] = rhoFlux[nG][nface];
						rhouFluxFace[nG] = rhouFlux[nG][nface];
						rhovFluxFace[nG] = rhovFlux[nG][nface];
						rhoEFluxFace[nG] = rhoEFlux[nG][nface];
					}

					SurfInt[0][order] += process::surfaceInte(element, edgeName, rhoFluxFace, order);
					SurfInt[1][order] += process::surfaceInte(element, edgeName, rhouFluxFace, order);
					SurfInt[2][order] += process::surfaceInte(element, edgeName, rhovFluxFace, order);
					SurfInt[3][order] += process::surfaceInte(element, edgeName, rhoEFluxFace, order);
				}
			}
			return SurfInt;
		}
	}//end namespace auxEq

	namespace NSFEq
	{
<<<<<<< HEAD
<<<<<<< HEAD
		void solveNSFEquation()
		{
			std::vector<std::vector<double>> StiffMatrix(mathVar::orderElem + 1, std::vector<double>(mathVar::orderElem + 1, 0.0));
			std::vector<double> RHSOrder(4, 0.0);
			std::vector<double> rhoRHSTerm(mathVar::orderElem + 1, 0.0),
				rhouRHSTerm(mathVar::orderElem + 1, 0.0),
				rhovRHSTerm(mathVar::orderElem + 1, 0.0),
				rhoERHSTerm(mathVar::orderElem + 1, 0.0);

			std::vector<double> rhoVector(mathVar::orderElem + 1, 0.0),
				rhouVector(mathVar::orderElem + 1, 0.0),
				rhovVector(mathVar::orderElem + 1, 0.0),
				rhoEVector(mathVar::orderElem + 1, 0.0);

			for (int nelement = 0; nelement < meshVar::nelem2D; nelement++)
			{
				//1) Calculate Stiff matrix
				StiffMatrix = process::calculateStiffMatrix(nelement);

				for (int order = 0; order <= mathVar::orderElem; order++)
				{
					//At each order, right hand side terms of all conservative variables are calculated

					//2) Calculate right hand side terms
					RHSOrder = process::NSFEq::CalcRHSTerm(nelement, order);
					rhoRHSTerm[order] = RHSOrder[0];
					rhouRHSTerm[order] = RHSOrder[1];
					rhovRHSTerm[order] = RHSOrder[2];
					rhoERHSTerm[order] = RHSOrder[3];
					//4 vector RHS of 4 conservative variables for each direction is achieved after these steps have been finished
				}

				//3) Solve for conservative variables
				rhoVector = math::SolveSysEqs(StiffMatrix, rhoRHSTerm);
				rhouVector = math::SolveSysEqs(StiffMatrix, rhouRHSTerm);
				rhovVector = math::SolveSysEqs(StiffMatrix, rhovRHSTerm);
				rhoEVector = math::SolveSysEqs(StiffMatrix, rhoERHSTerm);

				//4) Save results to conservative variables array
				for (int order = 0; order <= mathVar::orderElem; order++)
				{
					rho[nelement][order] = rhoVector[order];
					rhou[nelement][order] = rhouVector[order];
					rhov[nelement][order] = rhovVector[order];
					rhoE[nelement][order] = rhoEVector[order];
				}
			}
		}

=======
>>>>>>> parent of 622c1e9... hoàn thành phần tính RHS, chờ công thức tính time step
=======
>>>>>>> parent of 622c1e9... hoàn thành phần tính RHS, chờ công thức tính time step
		/*Function calculates right hand side terms of all conservative variables at ONLY one order*/
		std::vector<double> CalcRHSTerm(int element, int order)
		{
			std::vector<double> RHS(4, 0.0);
			std::vector<double> VolInt(4, 0.0);
			/*Volume integral term===========================================================================*/
			VolInt = process::NSFEq::calcVolumeIntegralTerms(element, order);
			RHS[0] += VolInt[0];
			RHS[1] += VolInt[1];
			RHS[2] += VolInt[2];
			RHS[3] += VolInt[3];

			/*Surface integral term===========================================================================*/

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
			double uVal(math::pointValue(element, a, b, 2, 1)), vVal(math::pointValue(element, a, b, 3, 1)), muVal(math::pointValue(element, a, b, 7, 1));
			/*ViscousTerm 2D array has 4 row 2 column:
			- column 1: Ox direction
			- column 2: Oy direction*/
			std::vector<std::vector<double>> ViscousTerm(4, std::vector<double>(2, 0.0));
			std::vector<std::vector<double>> StressHeatFlux(2, std::vector<double>(3, 0.0));
			std::vector<double> vectorU(4, 0.0);
			std::vector<double> vectordUx(4, 0.0);
			std::vector<double> vectordUy(4, 0.0);

			/*calculate conservative and derivative variables*/

			vectorU[0] = math::pointValue(element, a, b, 1, 2);
			vectorU[1] = math::pointValue(element, a, b, 2, 2);
			vectorU[2] = math::pointValue(element, a, b, 3, 2);

			vectordUx[0] = math::pointAuxValue(element, a, b, 1, 1);
			vectordUy[0] = math::pointAuxValue(element, a, b, 1, 2);

			vectordUx[1] = math::pointAuxValue(element, a, b, 2, 1);
			vectordUy[1] = math::pointAuxValue(element, a, b, 2, 2);

			vectordUx[2] = math::pointAuxValue(element, a, b, 3, 1);
			vectordUy[2] = math::pointAuxValue(element, a, b, 3, 2);

			vectordUx[3] = math::pointAuxValue(element, a, b, 4, 1);
			vectordUy[3] = math::pointAuxValue(element, a, b, 4, 2);

			/*calculate stresses and heat fluxes*/
			StressHeatFlux = math::viscousTerms::calcStressTensorAndHeatFlux(muVal, vectorU, vectordUx, vectordUy);
			std::tie(ViscousTerm[0][0], ViscousTerm[1][0], ViscousTerm[2][0], ViscousTerm[3][0]) = math::viscousTerms::calcViscousTermsFromStressHeatFluxMatrix(StressHeatFlux, uVal, vVal, 1);
			std::tie(ViscousTerm[0][1], ViscousTerm[1][1], ViscousTerm[2][1], ViscousTerm[3][1]) = math::viscousTerms::calcViscousTermsFromStressHeatFluxMatrix(StressHeatFlux, uVal, vVal, 2);
			return ViscousTerm;
		}

		std::vector<double> calcVolumeIntegralTerms(int element, int order)
		{
			std::vector<std::vector<double>> ViscousTerms(4, std::vector<double>(2, 0.0));
			std::vector<std::vector<double>> InviscidTerms(4, std::vector<double>(2, 0.0));
			std::vector<double> VolInt(4, 0.0);

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

			/*CALCULATE INTEGRALS*/
			/*A INVISCID TERMS*/
			/*- integral term of (dB/dx)*(rho*u)*/
			VolInt[0] += process::volumeInte(element, InvisGsVolX1, order, 1);
			/*- integral term of (dB/dy)*(rho*v)*/
			VolInt[0] += process::volumeInte(element, InvisGsVolY1, order, 2);

			/*- integral term of (dB/dx)*(rho*u^2 + p)*/
			VolInt[1] += process::volumeInte(element, InvisGsVolX2, order, 1);
			/*- integral term of (dB/dy)*(rho*u*v)*/
			VolInt[1] += process::volumeInte(element, InvisGsVolY2, order, 2);

			/*- integral term of (dB/dx)*(rho*u*v)*/
			VolInt[2] += process::volumeInte(element, InvisGsVolX3, order, 1);
			/*- integral term of (dB/dy)*(rho*v^2 + p)*/
			VolInt[2] += process::volumeInte(element, InvisGsVolY3, order, 2);

			/*- integral term of (dB/dx)*(rho*totalE + p)*u*/
			VolInt[3] += process::volumeInte(element, InvisGsVolX4, order, 1);
			/*- integral term of (dB/dy)*(rho*totalE + p)*v*/
			VolInt[3] += process::volumeInte(element, InvisGsVolY4, order, 2);

			/*B VISCOUS TERMS*/
			/*- integral term of (dB/dx)*(0.0)*/
			VolInt[0] += 0.0;
			/*- integral term of (dB/dy)*(0.0)*/
			VolInt[0] += 0.0;

			/*- integral term of (dB/dx)*(tau_xx)*/
			VolInt[1] += process::volumeInte(element, ViscGsVolX2, order, 1);
			/*- integral term of (dB/dy)*(tau_xy)*/
			VolInt[1] += process::volumeInte(element, ViscGsVolY2, order, 2);

			/*- integral term of (dB/dx)*(tau_xy)*/
			VolInt[2] += process::volumeInte(element, ViscGsVolX3, order, 1);
			/*- integral term of (dB/dy)*(tau_yy)*/
			VolInt[2] += process::volumeInte(element, ViscGsVolY3, order, 2);

			/*- integral term of (dB/dx)*(u*tau_xx + v*tau_xy + Qx)*/
			VolInt[3] += process::volumeInte(element, ViscGsVolX4, order, 1);
			/*- integral term of (dB/dy)*(u*tau_xy + v*tau_yy + Qy)*/
			VolInt[3] += process::volumeInte(element, ViscGsVolY4, order, 2);
			/*End volume terms===========================================================================*/
			return VolInt;
		}

		std::vector<double> calcSurfaceIntegralTerms(int element, int order)
		{
			std::vector<std::vector<double>> Fluxes(4, std::vector<double>(2, 0.0));
			std::vector <double> inviscFlux1(mathVar::nGauss + 1, 0.0),
				inviscFlux2(mathVar::nGauss + 1, 0.0),
				inviscFlux3(mathVar::nGauss + 1, 0.0),
				inviscFlux4(mathVar::nGauss + 1, 0.0);

			std::vector <double> ViscFlux1(mathVar::nGauss + 1, 0.0),
				ViscFlux2(mathVar::nGauss + 1, 0.0),
				ViscFlux3(mathVar::nGauss + 1, 0.0),
				ViscFlux4(mathVar::nGauss + 1, 0.0);

			std::vector<double> SurInt(4, 0.0);

			int elemType(auxUlti::checkType(element)), edgeName(0);
			int faceBcType(0);

			for (int nface = 0; nface < elemType; nface++)
			{
				edgeName = meshVar::inedel[nface][element];
				faceBcType = auxUlti::getBCType(edgeName);

				if (faceBcType == 0)  //internal edge
				{
					for (int nGauss = 0; nGauss <= mathVar::nGauss; nGauss++)
					{
						Fluxes = process::NSFEq::getGaussVectorOfConserVarFluxesAtInternal(edgeName, element, nGauss);
						inviscFlux1[nGauss] = Fluxes[0][0];
						inviscFlux2[nGauss] = Fluxes[1][0];
						inviscFlux3[nGauss] = Fluxes[2][0];
						inviscFlux4[nGauss] = Fluxes[3][0];

						ViscFlux1[nGauss] = Fluxes[0][1];
						ViscFlux2[nGauss] = Fluxes[1][1];
						ViscFlux3[nGauss] = Fluxes[2][1];
						ViscFlux4[nGauss] = Fluxes[3][1];
					}
				}
				else  //boundary edge
				{
					for (int nGauss = 0; nGauss <= mathVar::nGauss; nGauss++)
					{
						Fluxes = NSFEqBCsImplement(element, edgeName, nGauss);
						inviscFlux1[nGauss] = Fluxes[0][0];
						inviscFlux2[nGauss] = Fluxes[1][0];
						inviscFlux3[nGauss] = Fluxes[2][0];
						inviscFlux4[nGauss] = Fluxes[3][0];

						ViscFlux1[nGauss] = Fluxes[0][1];
						ViscFlux2[nGauss] = Fluxes[1][1];
						ViscFlux3[nGauss] = Fluxes[2][1];
						ViscFlux4[nGauss] = Fluxes[3][1];
					}
				}
				SurInt[0] += process::surfaceInte(element, edgeName, inviscFlux1, order) + process::surfaceInte(element, edgeName, ViscFlux1, order);
				SurInt[1] += process::surfaceInte(element, edgeName, inviscFlux2, order) + process::surfaceInte(element, edgeName, ViscFlux2, order);
				SurInt[2] += process::surfaceInte(element, edgeName, inviscFlux3, order) + process::surfaceInte(element, edgeName, ViscFlux3, order);
				SurInt[3] += process::surfaceInte(element, edgeName, inviscFlux4, order) + process::surfaceInte(element, edgeName, ViscFlux4, order);
			}
			return SurInt;
		}

		std::vector<std::vector<double>> getGaussVectorOfConserVarFluxesAtInternal(int edgeName, int element, double nGauss)
		{
			std::vector<std::vector<double>> Fluxes(4, std::vector<double>(2, 0.0));
			std::vector<double> UMinus(4, 0.0),
				UPlus(4, 0.0),
				dUXMinus(4, 0.0),
				dUXPlus(4, 0.0),
				dUYMinus(4, 0.0),
				dUYPlus(4, 0.0),
				normalVector(2, 0.0);

			/*Normal vector*/
			normalVector[0] = auxUlti::getNormVectorComp(element, edgeName, 1);
			normalVector[1] = auxUlti::getNormVectorComp(element, edgeName, 2);

			/*INVISCID TERMS*/
			std::tie(UPlus[0], UMinus[0]) = math::internalSurfaceValue(edgeName, element, nGauss, 1, 2);  //rho
			std::tie(UPlus[1], UMinus[1]) = math::internalSurfaceValue(edgeName, element, nGauss, 2, 2);  //rhou
			std::tie(UPlus[2], UMinus[2]) = math::internalSurfaceValue(edgeName, element, nGauss, 3, 2);  //rhov
			std::tie(UPlus[3], UMinus[3]) = math::internalSurfaceValue(edgeName, element, nGauss, 4, 2);  //rhoE
			
			/*VISCOUS TERMS*/
			std::tie(dUXPlus[0], dUXMinus[0]) = math::internalSurfaceDerivativeValue(edgeName, element, nGauss, 1, 1);  //drhox
			std::tie(dUYPlus[0], dUYMinus[0]) = math::internalSurfaceDerivativeValue(edgeName, element, nGauss, 1, 2);  //drhoy

			std::tie(dUXPlus[1], dUXMinus[1]) = math::internalSurfaceDerivativeValue(edgeName, element, nGauss, 2, 1);  //drhoux
			std::tie(dUYPlus[1], dUYMinus[1]) = math::internalSurfaceDerivativeValue(edgeName, element, nGauss, 2, 2);  //drhoux

			std::tie(dUXPlus[2], dUXMinus[2]) = math::internalSurfaceDerivativeValue(edgeName, element, nGauss, 3, 1);  //drhovx
			std::tie(dUYPlus[2], dUYMinus[2]) = math::internalSurfaceDerivativeValue(edgeName, element, nGauss, 3, 2);  //drhovx

			std::tie(dUXPlus[3], dUXMinus[3]) = math::internalSurfaceDerivativeValue(edgeName, element, nGauss, 4, 1);  //drhoEx
			std::tie(dUYPlus[3], dUYMinus[3]) = math::internalSurfaceDerivativeValue(edgeName, element, nGauss, 4, 2);  //drhoEx
			
			/*Calculate fluxes*/
			Fluxes = math::numericalFluxes::NSFEqFluxFromConserVars(UPlus, UMinus, dUXPlus, dUXMinus, dUYPlus, dUYMinus, normalVector);
			return Fluxes;
		}
	}

	double volumeInte(int elem, std::vector< std::vector<double> > &Ui, int order, int direction)
	{
		double Int(0.0);
		std::vector<std::vector<double>> A(mathVar::nGauss + 1, std::vector<double>(mathVar::nGauss + 1, 0.0));
		double dBi(0.0), a(0.0), b(0.0);

		for (int na = 0; na <= mathVar::nGauss; na++)
		{
			for (int nb = 0; nb <= mathVar::nGauss; na++)
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
		std::vector<double> F(mathVar::nGauss, 0.0);
		double Int(0.0);

		inte = 0.0;  //value of integral is reset for each order
		for (int nG = 0; nG <= mathVar::nGauss; nG++)
		{
			std::tie(a, b) = auxUlti::getGaussSurfCoor(edge, elem, nG);
			math::basisFc(a, b, mathVar::orderElem);  //mathVar::B is changed after this command line is excuted
			Bi = mathVar::B[order];
			F[nG] = Bi * FluxVector[nG];
		}
		inte += math::surfaceInte(F, edge, elem);
		Int = inte;

		return Int;
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
}