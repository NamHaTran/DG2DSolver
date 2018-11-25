#include "DGController.h"
#include "DGIOLib.h"
#include "DGMessagesLib.h"
#include "DGMeshReaderLib.h"
#include "DGPostProcessLib.h"
#include "VarDeclaration.h"
#include "CommandCheck.h"
#include "DGProcLib.h"
#include "DGAuxUltilitiesLib.h"
#include <iostream>
#include <windows.h>

void Executer(std::string cmd)
{
	if (preProcessKey::checkUnvReader(cmd))
	{
		std::string UnvReaderLoc(systemVar::wD + "\\Ultilities\\MeshReader\\UnvReader.exe");
		auxUlti::openFileEXE(UnvReaderLoc);
	}
	else if (processKey::checkDGRun(cmd))
	{
		//PreProcessing();
		Processing();
	}
	else if (postProcessKey::checkExit(cmd))
	{
		systemVar::endKey = true;
	}
	else if (preProcessKey::checkUnvHelper(cmd))
	{
		message::UnvReaderHelp();
	}
	else if (preProcessKey::checkBCsHelper(cmd))
	{
		message::BCsHelp();
	}
	else if (preProcessKey::reSubmit(cmd))
	{
		IO::getCase();
		PreProcessing();
	}
	else if (preProcessKey::debug::checkElement(cmd))
	{
		int input(-1);
		std::cout << "Input element ID (ID supplied by SALOME): ";
		std::cin >> input;
		std::cout << " \n";
		debugTool::checkElemInfor(input);
	}
	else
	{
		std::cout << "	ERROR: unknow command <" << cmd << ">!!!\n";
	}
}

void Processing()
{
	/*SET INITIAL VALUES*/
	process::setIniValues();

	std::cout << " \n" << "Simulation is started\n";

	while (process::checkRunningCond())
	{
		systemVar::iterCount++;
		std::cout << "Iteration " << systemVar::iterCount << std::endl;

		//APPLY LIMITER
		limitVal::numOfLimitCell = 0;
		process::limiter::limiter();

		//CALCULATE TIME STEP
		process::Euler::calcGlobalTimeStep();

		//SOLVE AUXILARY EQUATION
		process::auxEq::solveAuxEquation();

		//SOLVE NSF EQUATION
		process::NSFEq::solveNSFEquation();

		//UPDATE VARIABLES
		process::NSFEq::updateVariables();
	
		systemVar::savingCout++;
		if (systemVar::savingCout == systemVar::wrtI) //
		{
			std::cout << "Saving case...\n" << std::endl;
			DG2Tecplot::exportCellCenteredData(systemVar::iterCount);
			systemVar::savingCout = 0;
		}
	}
}

void PreProcessing()
{
	/*LOAD MESH*/
	IO::loadMesh();

	/*LOAD CONSTANTS*/
	IO::loadConstants();

	/*LOAD p T U*/
	IO::loadpTU();
	//Check subsonic
	refValues::subsonic = auxUlti::checkSubSonic();

	/*PROCESS MESH*/
	MshReader::meshProcess();

	/*CALCULATE JACOBIAN, BASIS FUNCTION AND GAUSSIAN*/
	meshParam::GaussParam();
	meshParam::basisFcParam();
	meshParam::JacobianParam();

	/*CALCULATE CELL METRICS*/
	meshParam::calcCellMetrics();

	/*CALCULATE COORDINATES DERIVATIVES*/
	meshParam::derivCoordinates();

	/*RESIZE ARRAYS*/
	auxUlti::resizeDGArrays();

	auxUlti::mappingEdges();
}

void PostProcessing()
{
	
}