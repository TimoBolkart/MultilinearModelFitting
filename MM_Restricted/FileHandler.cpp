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

#include "FileHandler.h"

#include <iostream>
#include <fstream>


std::string FileHandler::getFileExtension(const std::string& sstrFileName)
{
	if(sstrFileName.empty())
	{
		return "";
	}

	const size_t pos = sstrFileName.rfind(".");
	if(pos != std::string::npos)
	{
		return sstrFileName.substr(pos+1);
	}

	return "";
}

bool FileHandler::fileExist(const std::string& sstrFileName)
{
	if(sstrFileName.empty())
	{
		return false;
	}

	std::fstream inStream;
	inStream.open(sstrFileName, std::ios::in);

	if(inStream.is_open())
	{
		inStream.close();
		return true;
	}
	
	return false;
}

std::string FileHandler::getFileName(const std::string& sstrFullFileName)
{
	std::string sstrFileName = "";

	if(sstrFullFileName.empty())
	{
		return sstrFileName;
	}

	if(!FileHandler::fileExist(sstrFullFileName))
	{
		return sstrFileName;
	}

	size_t posEnd = sstrFullFileName.rfind(".");
	if(posEnd == std::string::npos)
	{
		posEnd = sstrFullFileName.size();
	}

	size_t posStart = sstrFullFileName.rfind("/");
	if(posStart != std::string::npos)
	{
		sstrFileName = sstrFullFileName.substr(posStart+1, posEnd-posStart-1);
		return sstrFileName;
	}

	posStart = sstrFullFileName.rfind("\\");
	if(posStart != std::string::npos)
	{
		sstrFileName = sstrFullFileName.substr(posStart+1, posEnd-posStart-1);
		return sstrFileName;
	}

	return sstrFileName;
}

std::string FileHandler::getFilePath(const std::string& sstrFileName)
{
	std::string sstrFilePath = "";

	if(sstrFileName.empty())
	{
		return sstrFilePath;
	}

	if(!FileHandler::fileExist(sstrFileName))
	{
		return sstrFilePath;
	}

	size_t pos = sstrFileName.rfind("/");
	if(pos != std::string::npos)
	{
		sstrFilePath = sstrFileName.substr(0, pos);
		return sstrFilePath;
	}

	pos = sstrFileName.rfind("\\");
	if(pos != std::string::npos)
	{
		sstrFilePath = sstrFileName.substr(0, pos);
		return sstrFilePath;
	}

	return sstrFilePath;
}


