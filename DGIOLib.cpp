#include "DGIOLib.h"
#include "DGMessagesLib.h"
#include "ConstDeclaration.h"
#include "VarDeclaration.h"
#include "DGAuxUltilitiesLib.h"
#include <fstream>
#include <string>
#include <iostream>
#include <sstream>
#include <vector>

namespace IO
{
	void dispLogo()
	{
		std::string logoStr(" ");
		logoStr = message::headerFile();
		std::cout << logoStr << std::endl;
	}

	void getCase()
	{
		std::cout << "***Getting case's information***\n";
		systemVar::wD = auxUlti::workingdir();

		/*Get caseName from SubmitCase*/
		std::string submitLoc(" ");  //submitLoc contents location of submitingCase
		submitLoc = systemVar::wD + "\\CASES\\SubmitCase.txt";

		std::ifstream submitCaseFlux(submitLoc.c_str());
		if (submitCaseFlux)
		{
			std::string line, keyWord;
			int numLine(0), keyWFlag(0);
			while (std::getline(submitCaseFlux, line))  //Read submitCaseFlux line by line
			{
				std::istringstream line2str(line);
				std::vector<std::string> ptr;
				//Split <line2str> stringstream into array of words
				while ((line2str >> keyWord))
				{
					ptr.push_back(keyWord);
				}

				int numWrd = ptr.size();

				if (numWrd == 2)
				{
					std::string str1(ptr[0]), str2("caseName");
					if (str1.compare(str2) == 0)
					{
						systemVar::caseName = ptr[1];
						keyWFlag = 1;
						std::cout << "	Case " << systemVar::caseName << " has been submitted\n";
					}
				}
			}
			if (keyWFlag == 0)
			{
				std::cout << message::undfKeyW("caseName", submitLoc) << std::endl;
			}
		}
		else
		{
			message::writeLog((systemVar::wD + "\\CASES\\"), "", message::opFError("SubmitCase.txt", submitLoc));
		}
		systemVar::pwd = systemVar::wD + "\\CASES\\" + systemVar::caseName;
	}

	void loadMesh()
	{
		/*Declare loading locations*/
		std::string  Elem1DLoc = systemVar::pwd + "\\Constant\\Mesh\\Elements1D.txt";
		std::string  ptLoc = systemVar::pwd + "\\Constant\\Mesh\\Points.txt";
		std::string  Elem2DLoc = systemVar::pwd + "\\Constant\\Mesh\\Elements2D.txt";
		std::string  bcLoc = systemVar::pwd + "\\Constant\\boundaryPatch.txt";

		/*Load Points*/
		std::ifstream ptFlux(ptLoc.c_str());
		if (ptFlux)
		{
			std::string line(" ");
			double x, y, z;
			meshVar::npoin = 0;
			while (std::getline(ptFlux, line))
			{
				std::istringstream ptData(line);
				ptData >> meshVar::npoin >> x >> y >> z;
				meshVar::Points[meshVar::npoin - 1][0] = x;
				meshVar::Points[meshVar::npoin - 1][1] = y;
				meshVar::Points[meshVar::npoin - 1][2] = z;
			}
		}
		else
		{
			message::writeLog(systemVar::pwd, systemVar::caseName, message::opFError("Points.txt", ptLoc));
		}

		/*Load 1D elements*/
		std::ifstream Elem1DFlux(Elem1DLoc.c_str());
		if (Elem1DFlux)
		{
			int node1(0), node2(0), bcGrp(0);  //nelem1D is number of 1D elements
			std::string line(" ");
			while (std::getline(Elem1DFlux, line))
			{
				std::istringstream elData(line);
				elData >> meshVar::nelem1D >> node1 >> node2 >> bcGrp;
				meshVar::Elements1D[meshVar::nelem1D - 1][0] = node1 - 1;
				meshVar::Elements1D[meshVar::nelem1D - 1][1] = node2 - 1;
				meshVar::Elements1D[meshVar::nelem1D - 1][2] = bcGrp;
			}
		}
		else
		{
			message::writeLog(systemVar::pwd, systemVar::caseName, message::opFError("Elements1D.txt", Elem1DLoc));
		}

		/*Load 2D elements*/
		std::ifstream Elem2DFlux(Elem2DLoc.c_str());
		if (Elem2DFlux)
		{
			int temp(0), node1(0), node2(0), node3(0), node4(0);  //nelem2D is number of 2D elements
			std::string line(" ");
			meshVar::nelem2D = 0;
			while (std::getline(Elem2DFlux, line))
			{
				meshVar::nelem2D++;
				std::istringstream elData(line);
				elData >> temp >> node1 >> node2 >> node3 >> node4;
				meshVar::Elements2D[meshVar::nelem2D - 1][0] = node1 - 1;
				meshVar::Elements2D[meshVar::nelem2D - 1][1] = node2 - 1;
				meshVar::Elements2D[meshVar::nelem2D - 1][2] = node3 - 1;
				meshVar::Elements2D[meshVar::nelem2D - 1][3] = node4 - 1;
			}
		}
		else
		{
			message::writeLog(systemVar::pwd, systemVar::caseName, message::opFError("Elements2D.txt", Elem2DLoc));
		}

		/*Load boundary conditions*/
		std::ifstream bcFlux(bcLoc.c_str());
		if (bcFlux)
		{
			int bcType(0), boundIndex(0);
			std::string line(" "), keyWord(" ");
			while (std::getline(bcFlux, line))
			{
				std::istringstream line2str(line);
				std::vector<std::string> ptr;
				//Split <line2str> stringstream into array of words
				while ((line2str >> keyWord))
				{
					if ((keyWord.compare(" ") != 0))
					{
						ptr.push_back(keyWord);
					}
				}

				int numWrd = ptr.size();
				if (numWrd != 0)
				{
					std::string str1(ptr[0]);
					if (str1.compare("{") == 0)  //0 means str1 is the same as str2
					{
						//start group
						boundIndex++;
					}

					if ((numWrd >= 2)&(boundIndex != 0))
					{
						if ((boundIndex <= bcSize))
						{
							meshVar::BoundaryType[boundIndex - 1][0] = boundIndex;
							std::string str2 = ptr[1];
							if (str1.compare("Type") == 0)
							{
								if (str2.compare("wall") == 0)
								{
									meshVar::BoundaryType[boundIndex - 1][1] = 1;  //type 1 is wall
								}
								else if (str2.compare("patch") == 0)
								{
									meshVar::BoundaryType[boundIndex - 1][1] = 2;  //type 2 is patch
								}
								else if (str2.compare("symmetry") == 0)
								{
									meshVar::BoundaryType[boundIndex - 1][1] = 3;  //type 3 is symmetry
								}
								else
								{
									std::string errorStr = "Boundary type <" + ptr[1] + R"(> is unknown, available boundary types in file boundaryPatch are:
	wall
	patch
	symmetry)";
									message::writeLog(systemVar::pwd, systemVar::caseName, errorStr);
								}
							}
							else if (str1.compare("Method") == 0)
							{
								if (str2.compare("weakRiemann") == 0)
								{
									meshVar::BoundaryType[boundIndex - 1][2] = 1;  //type 1 weak Riemann
								}
								else if (str2.compare("weakPrescribed") == 0)
								{
									meshVar::BoundaryType[boundIndex - 1][2] = 2;  //type 2 is weak prescribed
								}
								else
								{
									std::string errorStr = "Boundary method <" + ptr[1] + R"(> is unknown, available boundary methods in file boundaryPatch are:
	weakRiemann
	weakPrescibed)";
									message::writeLog(systemVar::pwd, systemVar::caseName, errorStr);
								}
							}
						}
						else
						{
							std::string errorStr = "Maximum number of boundaries DG2D solver supports is only " + std::to_string(bcSize) + ". Make sure number of boundaries is less than " + std::to_string(bcSize);
							message::writeLog(systemVar::pwd, systemVar::caseName, errorStr);
						}
					}
				}
				meshVar::nBc = boundIndex;
			}
		}
		else
		{
			message::writeLog(systemVar::pwd, systemVar::caseName, message::opFError("boundaryPatch.txt", bcLoc));
		}
	}

	void SaveMeshInfor()
	{
		/*List of arrays need to be saved:
		- inedel: array contents information of edges belong to elements, column index is element index, each row in column contents index of edge belong to element, number of row is 4 because of default quad element.
		- ineled: array contents information of elements sharing edges, column index is edge index, each row in column contents index of element which edge is belong to, row 3 contents pointer.
		- inpoed: array contents information of points belong to edges, column 3 contents group which edge belongs to (group 0 is internal group), column 4 contents type of boundary (type 0 is internal edge)
		- normalVector: array contents information of normal vector of edges
		- MasterElemOfEdge: array content master element of iedge, use it with normalVector to get information of normal vector of edge*/

		std::string headerFile(message::headerFile());
		std::string  inedelLoc = systemVar::pwd + "\\Constant\\Mesh\\inedel.txt";
		std::string  ineledLoc = systemVar::pwd + "\\Constant\\Mesh\\ineled.txt";
		std::string  inpoedLoc = systemVar::pwd + "\\Constant\\Mesh\\inpoed.txt";
		std::string  normVectorLoc = systemVar::pwd + "\\Constant\\Mesh\\normalVector.txt";
		std::string  MasterElemOfEdgeLoc = systemVar::pwd + "\\Constant\\Mesh\\MasterElemOfEdge.txt";

		/*inedel*/
		std::ofstream Fluxinedel(inedelLoc.c_str());
		if (Fluxinedel)
		{
			Fluxinedel << headerFile << std::endl << "	inedel array\n" << "\n";
			for (int i = 0; i < 4; i++)
			{
				for (int j = 0; j < meshVar::nelem2D; j++)
				{
					Fluxinedel << meshVar::inedel[i][j] << " ";
				}
				Fluxinedel << "\n";
			}
		}
		else
		{
			message::writeLog(systemVar::pwd, systemVar::caseName, message::opFError("inedel.txt", inedelLoc));
		}

		/*ineled*/
		std::ofstream Fluxineled(ineledLoc.c_str());
		if (Fluxineled)
		{
			Fluxineled << headerFile << std::endl << "	ineled array\n" << "\n";
			for (int i = 0; i < 2; i++)
			{
				for (int j = 0; j < meshVar::inpoedCount; j++)
				{
					Fluxineled << meshVar::ineled[i][j] << " ";
				}
				Fluxineled << "\n";
			}
		}
		else
		{
			message::writeLog(systemVar::pwd, systemVar::caseName, message::opFError("ineled.txt", ineledLoc));
		}

		/*inpoed*/
		std::ofstream Fluxinpoed(inpoedLoc.c_str());
		if (Fluxinpoed)
		{
			Fluxinpoed << headerFile << std::endl << "	inpoed array\n" << "\n";
			for (int i = 0; i < 4; i++)
			{
				for (int j = 0; j < meshVar::inpoedCount; j++)
				{
					Fluxinpoed << meshVar::inpoed[i][j] << " ";
				}
				Fluxinpoed << "\n";
			}
		}
		else
		{
			message::writeLog(systemVar::pwd, systemVar::caseName, message::opFError("inpoed.txt", inpoedLoc));
		}

		/*normalVector*/
		std::ofstream FluxnormVector(normVectorLoc.c_str());
		if (FluxnormVector)
		{
			FluxnormVector << headerFile << std::endl << "	normalVector array\n" << "\n";
			for (int i = 0; i < 2; i++)
			{
				for (int j = 0; j < meshVar::inpoedCount; j++)
				{
					FluxnormVector << meshVar::normalVector[i][j] << " ";
				}
				FluxnormVector << "\n";
			}
		}
		else
		{
			message::writeLog(systemVar::pwd, systemVar::caseName, message::opFError("normalVector.txt", normVectorLoc));
		}

		/*inpoed*/
		std::ofstream FluxMasterElemOfEdge(MasterElemOfEdgeLoc.c_str());
		if (FluxMasterElemOfEdge)
		{
			FluxMasterElemOfEdge << headerFile << std::endl << "	MasterElemOfEdge array\n" << "\n";
			for (int i = 0; i < meshVar::inpoedCount; i++)
			{
				FluxMasterElemOfEdge << meshVar::MasterElemOfEdge[i] << " ";
				FluxMasterElemOfEdge << "\n";
			}
		}
		else
		{
			message::writeLog(systemVar::pwd, systemVar::caseName, message::opFError("MasterElemOfEdge.txt", MasterElemOfEdgeLoc));
		}
	}

	void loadConstants()
	{
		/*Read DGOptions*/
		std::string DGOptfileName("DGOptions.txt");
		std::string DGOptLoc(systemVar::wD + "\\CASES\\" + systemVar::caseName + "\\System");
		std::string DGOptkeyWordsDouble[2] = { "CourantNumber", "totalTime(s)" }, DGOptkeyWordsInt[3] = { "numberOfGaussPoints", "orderOfAccuracy", "writeInterval" }, DGOptkeyWordsBool[1] = { "writeLog" }, DGOptkeyWordsStr[1] = {};
		double DGOptoutDB[2] = {};
		int DGOptoutInt[3] = {};
		bool DGOptoutBool[1] = {};
		std::string DGOptoutStr[1] = {};

		readDataFile(DGOptfileName, DGOptLoc, DGOptkeyWordsDouble, DGOptkeyWordsInt, DGOptkeyWordsBool, DGOptkeyWordsStr, DGOptoutDB, DGOptoutInt, DGOptoutBool, DGOptoutStr, 2, 3, 1, 0);
		
		systemVar::CFL = DGOptoutDB[0];
		systemVar::Ttime = DGOptoutDB[1];
		mathVar::nGauss = DGOptoutInt[0];
		mathVar::orderElem = DGOptoutInt[1];
		systemVar::wrtI = DGOptoutBool[0];

		if (mathVar::orderElem>pow(mathVar::nGauss,2))
		{
			message::writeLog((systemVar::wD + "\\CASES\\" + systemVar::caseName), systemVar::caseName, message::nGaussOrderElemError());
		}
		
		/*Read Material*/
		std::string MatfileName("Material.txt");
		std::string MatLoc(systemVar::wD + "\\CASES\\" + systemVar::caseName + "\\Constant");
		std::string MatkeyWordsDouble[5] = { "gammaRatio", "gasConstant", "PrandtlNumber", "SutherlandAs", "SutherlandTs" }, MatkeyWordsInt[1] = {}, MatkeyWordsBool[1] = {}, MatkeyWordsStr[1] = {};
		double MatoutDB[5] = {};
		int MatoutInt[1] = {};
		bool MatoutBool[1] = {};
		std::string MatoutStr[1] = {};

		readDataFile(MatfileName, MatLoc, MatkeyWordsDouble, MatkeyWordsInt, MatkeyWordsBool, MatkeyWordsStr, MatoutDB, MatoutInt, MatoutBool, MatoutStr, 5, 0, 0, 0);
		material::gamma = MatoutDB[0];
		material::R = MatoutDB[1];
		material::Pr = MatoutDB[2];
		material::As = MatoutDB[3];
		material::Ts = MatoutDB[4];
	}

	void loadpTU()
	{
		/*Read U*/  //U must be read first of all
		readNonScalar();

		/*Read T*/
		readScalar("T");

		/*Read p*/
		readScalar("p");
	}

	void readDataFile(std::string fileName, std::string direction, std::string keyWordsDbl[], std::string keyWordsInt[], std::string keyWordsBool[], std::string keyWordsStr[], double *outDbl, int *outInt, bool *outBool, std::string *outStr, int numParamDbl, int numParamInt, int numParamBool, int numParamStr)  //Declaration of funciton which returns pointer of 1D array
	{
		/*User's guide:
		This function returns datas of type double and type int read from files.
		Input arguments:
		- fileName: name of file you want to read (file name and extension)
		- direction: working diectory contents data file (with no "\\" characters at the end)
		- keyWordsDbl: array contents keyWords of double values listed in file
		- keyWordsInt: array contents keyWords of int values listed in file
		- keyWordsBool: array contents keyWords of bool values listed in file
		- keyWordsStr: array contents keyWords of string values listed in file
		- outDbl: output array contents double values
		- outInt: output array contents int values
		- numParamDbl, numParamInt, numParamBool, numParamStr: number of double, int, bool, string parameters*/

		double dataDbl(0.0);
		int dataInt(0);
		std::string dataStr("abc");
		std::cout << "	Reading " << fileName <<"\n";
		std::string FileLoc(direction + "\\" + fileName);
		std::ifstream FileFlux(FileLoc.c_str());
		int indexDbl(0), indexInt(0), indexBool(0), indexStr(0);
		if (FileFlux)
		{
			std::string line, keyWord;
			while (std::getline(FileFlux, line))
			{
				std::istringstream line2str(line);
				std::vector<std::string> ptr;
				//Split <line2str> stringstream into array of words
				while ((line2str >> keyWord))
				{
					ptr.push_back(keyWord);
				}

				int numWrd = ptr.size();

				if (numWrd >= 2)
				{
					std::string str1(ptr[0]), str2(" "), str3(" "), str4(" "), str5(" ");
					std::istringstream strdata(ptr[1]);
					if (indexDbl<numParamDbl)
					{
						str2 = keyWordsDbl[indexDbl];
					}
					if (indexInt<numParamInt)
					{
						str3 = keyWordsInt[indexInt];
					}
					if (indexBool<numParamBool)
					{
						str4 = keyWordsBool[indexBool];
					}
					if (indexStr<numParamStr)
					{
						str5 = keyWordsStr[indexStr];
					}

					if (str1.compare(str2) == 0)  //double value
					{
						strdata >> dataDbl;
						outDbl[indexDbl]=dataDbl;
						indexDbl++;
					}
					else if ((str1.compare(str3) == 0))  //int value
					{
						strdata >> dataInt;
						outInt[indexInt] = dataInt;
						indexInt++;
					}
					else if ((str1.compare(str4) == 0))  //bool value
					{
						if ((ptr[1].compare("true") == 0) || (ptr[1].compare("yes") == 0))
						{
							outBool[indexBool] = true;
						}
						else if ((ptr[1].compare("false") == 0) || (ptr[1].compare("no") == 0))
						{
							outBool[indexBool] = false;
						}
						indexBool++;
					}
					else if ((str1.compare(str5) == 0))  //string value
					{
						outStr[indexStr] = ptr[1];
						indexStr++;
					}
				}
				if ((indexDbl>=numParamDbl) && (indexInt>=numParamInt) && (indexBool >= numParamBool) && (indexStr >= numParamStr))
				{
					break;
				}
			}
		}
		else
		{
			message::writeLog((systemVar::wD + "\\CASES\\" + systemVar::caseName), systemVar::caseName, message::opFError(fileName, FileLoc));
		}
	}

	void readNonScalar()
	{
		std::string fileName("U.txt"), tempStr("");
		std::string Loc(systemVar::wD + "\\CASES\\" + systemVar::caseName + "\\0");
		std::cout << "	Reading " << fileName << "\n";
		std::string FileLoc(Loc + "\\" + fileName);
		std::ifstream FileFlux(FileLoc.c_str());
		std::string str2("initialValue"), str3("noSlip"), str4("slip"), str5("fixedValue"), str6("Value"), str7("Group"), str8("Type"), str9("inletOutlet"), str10("inOutZeroGrad"), str11("symmetry");
		int bcGrp(0), type(0);

		if (FileFlux)
		{
			std::string line, keyWord;
			while (std::getline(FileFlux, line))
			{
				std::istringstream line2str(line);
				std::vector<std::string> ptr;
				//Split <line2str> stringstream into array of words
				while ((line2str >> keyWord))
				{
					ptr.push_back(keyWord);
				}

				int numWrd = ptr.size();
				if (numWrd != 0)
				{
					std::string str1(ptr[0]);
					if (str1.compare(str2) == 0)  //initial data
					{
						std::istringstream str_u(ptr[1]);
						std::istringstream str_v(ptr[2]);
						std::istringstream str_w(ptr[3]);
						str_u >> iniValues::uIni;
						str_v >> iniValues::vIni;
						str_w >> iniValues::wIni;
					}
					else if (str1.compare(str8) == 0)  //Bc Type
					{
						std::string str0(ptr[1]);

						if (meshVar::BoundaryType[bcGrp - 1][1] == 1)  //WALL
						{
							if ((str0.compare(str3) == 0))  //Type noSlip
							{
								bcValues::UBcType[bcGrp - 1] = 3;
								bcValues::uBC[bcGrp - 1] = 0.0;
								bcValues::vBC[bcGrp - 1] = 0.0;
								bcValues::wBC[bcGrp - 1] = 0.0;
							}
							else if ((str0.compare(str4) == 0))  //Type slip
							{
								//Use Maxwell-Smoluchovsky boundary condition
								bcValues::UBcType[bcGrp - 1] = 4;
								bcValues::uBC[bcGrp - 1] = 0.0;
								bcValues::vBC[bcGrp - 1] = 0.0;
								bcValues::vBC[bcGrp - 1] = 0.0;

								std::getline(FileFlux, line);
								std::istringstream fixedUStream(line);
								fixedUStream >> tempStr >> bcValues::uWall[bcGrp - 1] >> bcValues::vWall[bcGrp - 1] >> bcValues::wWall[bcGrp - 1];
							}
							else if ((str0.compare(str5) == 0))  //Type fixedValue
							{
								bcValues::UBcType[bcGrp - 1] = 2;
								std::getline(FileFlux, line);
								std::istringstream fixedUStream(line);
								fixedUStream >> tempStr >> bcValues::uBC[bcGrp - 1] >> bcValues::vBC[bcGrp - 1] >> bcValues::wBC[bcGrp - 1];
							}
							else
							{
								message::writeLog((systemVar::wD + "\\CASES\\" + systemVar::caseName), systemVar::caseName, message::undfBcType(str0, "U", "wall"));
							}
						}
						else if (meshVar::BoundaryType[bcGrp - 1][1] == 2)  //PATCH
						{
							if ((str0.compare(str10) == 0))  //Type inOutZeroGrad
							{
								bcValues::UBcType[bcGrp - 1] = 6;
								bcValues::uBC[bcGrp - 1] = 0.0;
								bcValues::vBC[bcGrp - 1] = 0.0;
								bcValues::wBC[bcGrp - 1] = 0.0;
							}
							else if ((str0.compare(str9) == 0))  //Type inletOutlet
							{
								bcValues::UBcType[bcGrp - 1] = 5;
								std::getline(FileFlux, line);
								std::istringstream fixedUStream(line);
								fixedUStream >> tempStr >> bcValues::uBC[bcGrp - 1] >> bcValues::vBC[bcGrp - 1] >> bcValues::wBC[bcGrp - 1];
							}
							else
							{
								message::writeLog((systemVar::wD + "\\CASES\\" + systemVar::caseName), systemVar::caseName, message::undfBcType(str0, "U", "patch"));
							}
						}
						else if (meshVar::BoundaryType[bcGrp - 1][1] == 3)  //SYMMETRY
						{
							if ((str0.compare(str11) == 0))  //Type symmetry
							{
								bcValues::UBcType[bcGrp - 1] = 7;
								bcValues::uBC[bcGrp - 1] = 0.0;
								bcValues::vBC[bcGrp - 1] = 0.0;
								bcValues::wBC[bcGrp - 1] = 0.0;
							}
							else
							{
								message::writeLog((systemVar::wD + "\\CASES\\" + systemVar::caseName), systemVar::caseName, message::undfBcType(str0, "U", "symmetry"));
							}
						}
					}
					else if ((str1.compare(str7) == 0))  //Group
					{
						std::istringstream str_bcGrp(ptr[1]);
						str_bcGrp >> bcGrp;
					}
				}
			}
		}
		else
		{
			message::writeLog((systemVar::wD + "\\CASES\\" + systemVar::caseName), systemVar::caseName, message::opFError(fileName, FileLoc));
		}
	}

	void readScalar(std::string fileName)
	{
		fileName = fileName + ".txt";
		std::string tempStr("");
		std::string Loc(systemVar::wD + "\\CASES\\" + systemVar::caseName + "\\0");
		std::cout << "	Reading " << fileName << "\n";
		std::string FileLoc(Loc + "\\" + fileName);
		std::ifstream FileFlux(FileLoc.c_str());
		std::string str2("initialValue"), str3("zeroNormalGrad"), str4("temperatureJump"), str5("fixedValue"), 
			str6("Value"), str7("Group"), str8("Type"), str9("inletOutlet"), str10("inOutZeroGrad"), str11("symmetry");
		int bcGrp(0), type(0);

		if (FileFlux)
		{
			std::string line, keyWord;
			if (fileName.compare("p.txt") == 0)
			{
				while (std::getline(FileFlux, line))
				{
					std::istringstream line2str(line);
					std::vector<std::string> ptr;
					//Split <line2str> stringstream into array of words
					while ((line2str >> keyWord))
					{
						ptr.push_back(keyWord);
					}

					int numWrd = ptr.size();
					if (numWrd != 0)
					{
						std::string str1(ptr[0]);
						if (str1.compare(str2) == 0)  //initial data
						{
							std::istringstream str_val(ptr[1]);
							str_val >> iniValues::pIni;;
						}
						else if (str1.compare(str8) == 0)  //Bc Type
						{
							std::string str0(ptr[1]);

							if (meshVar::BoundaryType[bcGrp - 1][1] == 1)  //WALL
							{
								if ((str0.compare(str3) == 0))  //Type zeroNormalGrad
								{
									bcValues::pBcType[bcGrp - 1] = 1;
									bcValues::pBC[bcGrp - 1] = iniValues::pIni;
								}
								else if ((str0.compare(str5) == 0))  //Type fixedValue
								{
									bcValues::pBcType[bcGrp - 1] = 2;
									std::getline(FileFlux, line);
									std::istringstream Stream(line);
									Stream >> tempStr >> bcValues::pBC[bcGrp - 1];
								}
								else
								{
									message::writeLog((systemVar::wD + "\\CASES\\" + systemVar::caseName), systemVar::caseName, message::undfBcType(str0, "p", "wall"));
								}
							}
							else if (meshVar::BoundaryType[bcGrp - 1][1] == 2)  //PATCH
							{
								if ((str0.compare(str10) == 0))  //Type inOutZeroGrad
								{
									bcValues::pBcType[bcGrp - 1] = 6;
									bcValues::pBC[bcGrp - 1] = iniValues::pIni;
								}
								else if ((str0.compare(str9) == 0))  //Type inletOutlet
								{
									bcValues::pBcType[bcGrp - 1] = 5;
									std::getline(FileFlux, line);
									std::istringstream Stream(line);
									Stream >> tempStr >> bcValues::pBC[bcGrp - 1];
								}
								else
								{
									message::writeLog((systemVar::wD + "\\CASES\\" + systemVar::caseName), systemVar::caseName, message::undfBcType(str0, "p", "patch"));
								}
							}
							else if (meshVar::BoundaryType[bcGrp - 1][1] == 3)  //SYMMETRY
							{
								if ((str0.compare(str11) == 0))  //Type symmetry
								{
									bcValues::pBcType[bcGrp - 1] = 7;
									bcValues::pBC[bcGrp - 1] = iniValues::pIni;
								}
								else
								{
									message::writeLog((systemVar::wD + "\\CASES\\" + systemVar::caseName), systemVar::caseName, message::undfBcType(str0, "p", "symmetry"));
								}
							}
						}
						else if ((str1.compare(str7) == 0))  //Group
						{
							std::istringstream str_bcGrp(ptr[1]);
							str_bcGrp >> bcGrp;
						}
					}
				}
			}
			else if (fileName.compare("T.txt") == 0)
			{
				while (std::getline(FileFlux, line))
				{
					std::istringstream line2str(line);
					std::vector<std::string> ptr;
					//Split <line2str> stringstream into array of words
					while ((line2str >> keyWord))
					{
						ptr.push_back(keyWord);
					}

					int numWrd = ptr.size();
					if (numWrd != 0)
					{
						std::string str1(ptr[0]);
						if (str1.compare(str2) == 0)  //initial data
						{
							std::istringstream str_val(ptr[1]);
							str_val >> iniValues::TIni;;
						}
						else if (str1.compare(str8) == 0)  //Bc Type
						{
							std::string str0(ptr[1]);

							if (meshVar::BoundaryType[bcGrp - 1][1] == 1)  //WALL
							{
								if ((str0.compare(str3) == 0))  //Type zeroNormalGrad
								{
									bcValues::TBcType[bcGrp - 1] = 1;
									bcValues::TBC[bcGrp - 1] = iniValues::TIni;
								}
								else if ((str0.compare(str5) == 0))  //Type fixedValue
								{
									bcValues::TBcType[bcGrp - 1] = 2;
									std::getline(FileFlux, line);
									std::istringstream Stream(line);
									Stream >> tempStr >> bcValues::TBC[bcGrp - 1];
								}
								else if ((str0.compare(str4) == 0))  //Type temperatureJump
								{
									bcValues::TBcType[bcGrp - 1] = 4;
									bcValues::TBC[bcGrp - 1] = iniValues::TIni;
									std::getline(FileFlux, line);
									std::istringstream Stream(line);
									Stream >> tempStr >> bcValues::TWall[bcGrp - 1];
								}
								else
								{
									message::writeLog((systemVar::wD + "\\CASES\\" + systemVar::caseName), systemVar::caseName, message::undfBcType(str0, "T", "wall"));
								}
							}
							else if (meshVar::BoundaryType[bcGrp - 1][1] == 2)  //PATCH
							{
								if ((str0.compare(str10) == 0))  //Type inOutZeroGrad
								{
									bcValues::TBcType[bcGrp - 1] = 6;
									bcValues::TBC[bcGrp - 1] = iniValues::TIni;
								}
								else if ((str0.compare(str9) == 0))  //Type inletOutlet
								{
									bcValues::TBcType[bcGrp - 1] = 5;
									std::getline(FileFlux, line);
									std::istringstream Stream(line);
									Stream >> tempStr >> bcValues::TBC[bcGrp - 1];
								}
								else
								{
									message::writeLog((systemVar::wD + "\\CASES\\" + systemVar::caseName), systemVar::caseName, message::undfBcType(str0, "T", "patch"));
								}
							}
							else if (meshVar::BoundaryType[bcGrp - 1][1] == 3)  //SYMMETRY
							{
								if ((str0.compare(str11) == 0))  //Type symmetry
								{
									bcValues::TBcType[bcGrp - 1] = 7;
									bcValues::TBC[bcGrp - 1] = iniValues::TIni;
								}
								else
								{
									message::writeLog((systemVar::wD + "\\CASES\\" + systemVar::caseName), systemVar::caseName, message::undfBcType(str0, "T", "symmetry"));
								}
							}
						}
						else if ((str1.compare(str7) == 0))  //Group
						{
							std::istringstream str_bcGrp(ptr[1]);
							str_bcGrp >> bcGrp;
						}
					}
				}
			}

		}
		else
		{
			message::writeLog((systemVar::wD + "\\CASES\\" + systemVar::caseName), systemVar::caseName, message::opFError(fileName, FileLoc));
		}

		for (int i = 0; i < meshVar::nBc; i++)
		{
			if ((bcValues::UBcType[i]!= bcValues::TBcType[i])&&(bcValues::UBcType[i]==4))
			{
				message::writeLog((systemVar::wD + "\\CASES\\" + systemVar::caseName), systemVar::caseName, message::SlipBcCompatibleError(i+1));
			}
		}
	}
}