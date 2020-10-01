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

#include "FileWriter.h"
#include "FileHandler.h"
#include "DataContainer.h"

#include <vector>
#include <fstream>

bool FileWriter::saveFile(const std::string& sstrFileName, const DataContainer& data)
{
	if(sstrFileName.empty())
	{
		return false;
	}

	const std::string suffix = FileHandler::getFileExtension(sstrFileName);
	if(suffix==std::string("off"))
	{
		return FileWriter::writeOff(sstrFileName, data);
	}

	return false;
}

bool FileWriter::writeOff(const std::string& sstrFileName, const DataContainer& data)
{
	if(sstrFileName.empty())
	{
		return false;
	}

	const size_t numVertices = data.getNumVertices();
	const size_t numFaces = data.getVertexIndexList().size();
	const size_t numEdges = 0;

	const char* cstrFileName = sstrFileName.c_str();
	if(cstrFileName==NULL)
	{
		return false;
	}

	std::fstream output;
	output.open(cstrFileName, std::ios::out);
	output.precision(7);
	output.setf(std::ios::fixed,std::ios::floatfield);

	const std::vector<double>& vertices = data.getVertexList();
	const std::vector<Vec3d*>& vertexColors = data.getVertexColorList();

	bool bValidColors = vertices.size() == 3*vertexColors.size();

	if(bValidColors)
	{
		output << "COFF" << std::endl;
	}
	else
	{
		output << "OFF" << std::endl;
	}

	output << numVertices << " " << numFaces << " " << numEdges << std::endl;

	for(size_t i = 0; i < numVertices; ++i)
	{
		output << vertices[3*i] << " " << vertices[3*i+1] << " " << vertices[3*i+2] << " ";
		if(bValidColors)
		{
			const Vec3d& color = *(vertexColors[i]);
			output << color[0] << " " << color[1] << " " << color[2] << " " << "1.0";
		}

		output << std::endl;
	}

	const std::vector<Vec3i*>& vertexIndexList = data.getVertexIndexList();
	std::vector<Vec3i*>::const_iterator currFaceIter = vertexIndexList.begin();
	std::vector<Vec3i*>::const_iterator endFaceIter = vertexIndexList.end();
	for( ; currFaceIter != endFaceIter; ++currFaceIter)
	{
		const Vec3i* pFace = *currFaceIter;
		output << "3" << " " << (*pFace)[0] << " " << (*pFace)[1]  << " " << (*pFace)[2]  << std::endl;
	}

	output.close();
	return true;
}