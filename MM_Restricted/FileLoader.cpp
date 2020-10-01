/*************************************************************************************************************************/
// This source is provided for NON-COMMERCIAL RESEARCH PURPOSES only, and is provided “as is” WITHOUT ANY WARRANTY; 
// without even the implied warranty of fitness for a particular purpose. The redistribution of the code is not permitted.
//
// If you use the source or part of it in a publication, cite the following paper:
// 
// T. Bolkart, S. Wuhrer
// 3D Faces in Motion: Fully Automatic Registration and Statistical Analysis.
// Computer Vision and Image Understanding, 131:100-115, 2015
//
// Copyright (c) 2015 Timo Bolkart, Stefanie Wuhrer
/*************************************************************************************************************************/

#include "FileLoader.h"
#include "MathHelper.h"
#include "FileHandler.h"

#include <fstream>
#include <string>
#include <iostream>


bool FileLoader::loadFile(const std::string& sstrFileName, DataContainer& outData)
{
	if(!FileHandler::fileExist(sstrFileName))
	{
		return false;
	}

	const std::string sstrSuffix = FileHandler::getFileExtension(sstrFileName);
	if(sstrSuffix=="off")
	{
		if(!loadOFF(sstrFileName, outData))
		{
			std::cout << "Unable to load " << sstrFileName << std::endl;
			return false;
		}

		// Clean the mesh if there are any faces. If there are no faces, we have a point cloud which should not be cleaned.
		if(outData.getNumFaces() > 0)
		{
			MathHelper::cleanMesh(outData);
		}
	}
	else
	{
		//Extend this for use of other file formats.

		std::cout << "File format " << sstrSuffix << " not supported" << std::endl;
		return false;
	}

	return true;
}

bool FileLoader::loadLandmarks(const std::string& sstrFileName, std::vector<double>& landmarks, std::vector<bool>& loaded)
{
	if(!FileHandler::fileExist(sstrFileName))
	{
		return false;
	}

	std::string sstrSuffix = FileHandler::getFileExtension(sstrFileName);
	if(sstrSuffix=="txt")
	{
		return loadSimpleLandmarks(sstrFileName, landmarks, loaded);
	}

	return false;
}

bool FileLoader::loadRestrictedMultilinearModel(const std::string& sstrFileName, std::vector<size_t>& modeDims, std::vector<size_t>& truncModeDims, std::vector<double>& multModel
																, std::vector<double>& modeMean, std::vector<double>& neutralModelMean, std::vector<double>& mean)
{
	if(sstrFileName.empty())
	{
		return false;
	}

	const char* cstrFileName = sstrFileName.c_str();
	if(cstrFileName==NULL)
	{
		return false;
	}

	FILE* pFile(NULL);
	errno_t err = fopen_s(&pFile, cstrFileName, "r");
	if(pFile==NULL || err!=0)
	{
		return false;
	}

	char strOutput[1000];
	int numModes(0);
	if(!readNextNumber(pFile, INTEGER, strOutput, numModes))
	{
		return false;
	}

	for(int i = 0; i < numModes; ++i)
	{
		char strOutput[1000];
		size_t modeDim(0);
		if(!readNextNumber(pFile, INTEGER, strOutput, modeDim))
		{
			return false;
		}

		modeDims.push_back(modeDim);
	}

	for(int i = 0; i < numModes; ++i)
	{
		char strOutput[1000];
		size_t truncatedModeDimension(0);
		if(!readNextNumber(pFile, INTEGER, strOutput, truncatedModeDimension))
		{
			return false;
		}

		truncModeDims.push_back(truncatedModeDimension);
	}

	if(numModes<1)
	{
		return false;
	}

	size_t numTensorElements = modeDims[0];
	size_t modeMeanSize = 0;
	const size_t meanSize = modeDims[0];

	for(int i = 1; i < numModes; ++i)
	{
		const size_t d_i = modeDims[i];
		const size_t m_i = truncModeDims[i];

		numTensorElements *= d_i;
		modeMeanSize += m_i;
	}

	multModel.reserve(numTensorElements);
	modeMean.reserve(modeMeanSize);
	mean.reserve(meanSize);

	bool bEnd(false);
	while(!bEnd)
	{
		char strOutput[1000];
		
		while(!bEnd 
			&& !isEqual(strOutput, "MultilinearModel") 
			&& !isEqual(strOutput, "ModeMean")
			&& !isEqual(strOutput, "NeutralModeMean")
			&& !isEqual(strOutput, "Mean"))
		{
			bEnd = readNextNode(pFile, strOutput);
		}

		if(!bEnd && isEqual(strOutput, "MultilinearModel"))
		{
			processBlockStructure(pFile, strOutput, multModel);
		}

		if(!bEnd && isEqual(strOutput, "ModeMean"))
		{
			processBlockStructure(pFile, strOutput, modeMean);
		}

		if(!bEnd && isEqual(strOutput, "NeutralModeMean"))
		{
			processBlockStructure(pFile, strOutput, neutralModelMean);
		}

		if(!bEnd && isEqual(strOutput, "Mean"))
		{
			processBlockStructure(pFile, strOutput, mean);
		}
	}

	return fclose(pFile)==0;
}

bool FileLoader::loadDataFile(const std::string& sstrDataFileName, std::vector<double>& data)
{
	const char* cstrFileName = sstrDataFileName.c_str();
	if(cstrFileName==NULL)
	{
		return false;
	}

	FILE* pFile(NULL);
	errno_t err = fopen_s(&pFile, cstrFileName, "r");
	if(pFile==NULL || err!=0)
	{
		return false;
	}

	while(true)
	{
		double val(0.0);
		char strOutput[1000];
		if(!readNextNumber(pFile, FLOATING_POINT, strOutput, val))
		{
			break;
		}

		data.push_back(val);
	}

	return fclose(pFile)==0;
}

bool FileLoader::loadOFF(const std::string& sstrFileName, DataContainer& outData)
{
	const char* cstrFileName = sstrFileName.c_str();
	if(cstrFileName==NULL)
	{
		return false;
	}

	FILE* pFile(NULL);
	errno_t err = fopen_s(&pFile, cstrFileName, "r");
	if(pFile==NULL || err!=0)
	{
		return false;
	}

	char output[1000];
	readNextNode(pFile, output);

	bool bColorOff(false);
	if(strcmp(output, "OFF") == 0)
	{
		bColorOff = false;
	}
	else if(strcmp(output, "COFF") == 0)
	{
		bColorOff = true;
	}
	else
	{
		return false;
	}

	int numVertices(0);
	readNextNumber(pFile, INTEGER, output, numVertices);

	int numFaces(0);
	readNextNumber(pFile, INTEGER, output, numFaces);

	int numEdges(0);
	readNextNumber(pFile, INTEGER, output, numEdges);

	if(numVertices < 1 /*|| numFaces < 1*/)
	{
		return false;
	}

	std::vector<double> vertexList;
	std::vector<Vec3i*>& vertexIndexList = outData.getVertexIndexList();
	std::vector<Vec3d*>& vertexColors = outData.getVertexColorList();

	bool bEnd(false);
	for(int vertex = 0; vertex < numVertices; ++vertex)
	{
		char tmpOutput[1000];
	
		for(int i = 0; i < 3; ++i)
		{
			double number(0.0);		
			if(!readNextNumber(pFile, FLOATING_POINT, tmpOutput, number))
			{
				bEnd = true;
				break;
			}

			vertexList.push_back(number);
		}

		if(vertexList.size()%3 != 0)
		{
			std::cout << "Loaded vertex list of wrong dimension" << std::endl;
		}

		if(bColorOff)
		{
			Vec3d* pColor = new Vec3d(0.0,0.0,0.0);
			for(int i = 0; i < 3; ++i)
			{
				double colorValue(0);		
				if(!readNextNumber(pFile, FLOATING_POINT, tmpOutput, colorValue))
				{
					bEnd = true;
					break;
				}

				(*pColor)[i] = colorValue;
			}

			if(bEnd)
			{
				delete pColor;
			return false;
			}
			else
			{
				vertexColors.push_back(pColor);
			}

			double alphaValue(0);		
			if(!readNextNumber(pFile, FLOATING_POINT, tmpOutput, alphaValue))
			{
				return false;
			}
		}
	}

	outData.setVertexList(vertexList);

	for(int face = 0; face < numFaces; ++face)
	{
		char tmpOutput[1000];

		int numPolyPoints(0);		
		if(!readNextNumber(pFile, INTEGER, tmpOutput, numPolyPoints))
		{
			break;
		}

		if(numPolyPoints!=3)
		{
			std::cout << "loadOff() - only triangles supported" << std::endl; 
			return false;
		}

		Vec3i* pPoly = new Vec3i(0,0,0);
		for(int i = 0; i < 3; ++i)
		{
			int vertexIndex(0);		
			if(!readNextNumber(pFile, INTEGER, tmpOutput, vertexIndex))
			{
				bEnd = true;
				break;
			}

			(*pPoly)[i] = vertexIndex;
		}

		if(bEnd)
		{
			delete pPoly;
			break;
		}
		else
		{
			vertexIndexList.push_back(pPoly);
		}
	}

	return fclose(pFile)==0;
}

bool FileLoader::loadSimpleLandmarks(const std::string& sstrFileName, std::vector<double>& landmarks, std::vector<bool>& loaded)
{
	const char* cstrFileName = sstrFileName.c_str();
	if(cstrFileName==NULL)
	{
		return false;
	}

	FILE* pFile(NULL);
	errno_t err = fopen_s(&pFile, cstrFileName, "r");
	if(pFile==NULL || err!=0)
	{
		return false;
	}

	std::vector<double> tmpLandmarks;
	if(!loadDataFile(sstrFileName, tmpLandmarks))
	{
		std::cout << "Loading simple landmarks failed " << sstrFileName << std::endl;
		return false;
	}

	if(tmpLandmarks.size() % 3 != 0)
	{
		std::cout << "Wrong landmarks dimension" << std::endl;
		return false;
	}

	const size_t numLandmarks = tmpLandmarks.size()/3;

	loaded.clear();
	loaded.resize(numLandmarks);
	for(int i = 0; i < numLandmarks; ++i)
	{
		loaded[i] = true;
	}

	landmarks = tmpLandmarks;

	return fclose(pFile)==0;
}

bool FileLoader::readNextNode(FILE* pFile, char* cstrOutput)
{
	//Read string of characters until a whitespace (blank, newline, tab) is found.
	int ret = fscanf(pFile, "%1000s", cstrOutput);

	//Ignore comments marked by # at the beginning of the line
	while(ret!=EOF && cstrOutput[0]=='#')
	{
		//Read until the end of the line is reached
		do
		{
			ret=fgetc(pFile);
		}
		while(ret!=EOF && ret!='\n');

		ret=fscanf(pFile, "%1000s", cstrOutput);
	}

	// separate "{", "}", "[" and "]"
	size_t len = strlen(cstrOutput); // returns length of string s
	if (len > 1)
	{
		char currChar = cstrOutput[len-1];

		if(currChar == '{' || currChar == '}'
			|| currChar == '[' || currChar == ']')
		{
			ungetc(cstrOutput[len-1], pFile);
			cstrOutput[len-1] = 0;
		}
	}

	return ret==EOF;
}

void FileLoader::processBlockStructure(FILE* pFile, char* strOutput, std::vector<double>& values)
{
	bool bEnd = readNodeBlocksUntil(pFile, "{", strOutput);;
	while(!bEnd)
	{
		double number(0.0);
		if(!readNextNumber(pFile, FLOATING_POINT, strOutput, number))
		{
			break;
		}

		values.push_back(number);
	}
}