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

#ifndef FILELOADER_H
#define FILELOADER_H

#include "Definitions.h"
#include "DataContainer.h"

class FileLoader
{
	enum Type
	{
		FLOATING_POINT,
		INTEGER
	};

public:
	//! Load geometry file. Extend this to load other file formats.
	//! \param sstrFileName			full file name
	//! \param outData				loaded geometry
	//! \return true if successful
	bool loadFile(const std::string& sstrFileName, DataContainer& outData);

	//! Load landmarks. Extend this to load other landmark file formats.
	//! \param sstrFileName			full landmark file name
	//! \param landmarks				vertices of loaded landmarks
	//! \param loaded					flags that indicate which vertices are successfully loaded (could be used in the presence of occlusions)
	//! \return true if successful
	bool loadLandmarks(const std::string& sstrFileName, std::vector<double>& landmarks, std::vector<bool>& loaded);

	//! Load the multilinear model.
	//! \param sstrFileName			full multilinear model file name
	//! \param modeDims				mode dimensions of the data d1 d2 d3
	//! \param truncModeDims		mode dimensions of the truncated multilinear model tensor d1 m2 m3
	//! \param multModel				multilinear model tensor
	//! \param modeMean				(m2+m3)-dimensional vector with mean coefficients of mode 2 and mode 3
	//! \param neutralModelMean	(m2+m3)-dimensional vector with mean coefficients of mode 2 and netural mean coefficients of mode 3
	//! \param mean					learned data mean
	//! \return true if successful
	bool loadRestrictedMultilinearModel(const std::string& sstrFileName, std::vector<size_t>& modeDims, std::vector<size_t>& truncModeDims, std::vector<double>& multModel
													, std::vector<double>& modeMean, std::vector<double>& neutralModelMean, std::vector<double>& mean);

	//! Load file of double  values.
	//! \param sstrDataFileName	full file name
	//! \param data					loaded data
	//! \return true if successful
	bool loadDataFile(const std::string& sstrDataFileName, std::vector<double>& data);

private:
	//! Specific loader for off files. Use similar structure to consider other file formats.
	//! \param sstrFileName		full file name
	//! \param outData			loaded geometry
	//! \return true if successful
	bool loadOFF(const std::string& sstrFileName, DataContainer& outData);

	//!
	//! \param sstrFileName
	//! \param landmarks
	//! \param loaded
	//! \return true if successful
	bool loadSimpleLandmarks(const std::string& sstrFileName, std::vector<double>& landmarks, std::vector<bool>& loaded);

	bool readNextNode(FILE* pFile, char* cstrOutput);

	//! return true if successful, false if not successful or eof is reached
	template<typename Unit>
	bool readNextNumber(FILE* pFile, const Type type, char* strOutput, Unit& number)
	{
		bool bEnd = readNextNode(pFile, strOutput);
		if(bEnd)
		{
			return false;
		}

		return convertToNumber(strOutput, type, number);
	}

	void processBlockStructure(FILE* pFile, char* strOutput, std::vector<double>& values);

	template<typename Unit>
	bool convertToNumber(const char* strIn, const Type& type, Unit& out)
	{
		char strBuffer[20];
		for(size_t i = 0; i < 20; ++i)
		{
			if(strIn[i] == ',')
			{
				strBuffer[i] = '\0';
				break;
			}
			else
			{
				strBuffer[i] = strIn[i];
			}
		}

		const char* cstrFormat = type==FLOATING_POINT ? "%lf%*c" : "%ld%*c";
		return sscanf(strBuffer, cstrFormat, &out)==1; 
	}

	template<typename Unit, size_t DIM>
	void processCoordinates(FILE* pFile, char* cstrOutput, std::vector<VecNX<Unit, DIM>*>& vertexList)
	{
		bool bEnd = readNodeBlocksUntil(pFile, "[", cstrOutput);
		assert(!bEnd && !isEqual(cstrOutput, "}"));

		do
		{
			VecNX<double, DIM>* pVec = new VecNX<double, DIM>();
			if(!readDataBlock(pFile, FLOATING_POINT, cstrOutput, *pVec))
			{
				delete pVec;
				break;
			}
			
			vertexList.push_back(pVec);
		}
		while(!isEqual(cstrOutput, "}"));
	}

	template<typename Unit, size_t DIM>
	void processIndices(FILE* pFile, char* /*cstrOutput*/, std::vector<VecNX<Unit, DIM>*>& indexList)
	{
		const char endChar[] = {'[', ']', '}', '\0'};

		char tmpOutput[1000];
		bool bEnd = readBlockUntil(pFile, endChar, tmpOutput);
		assert(!bEnd && !isEqual(tmpOutput, "}"));

		while(!bEnd)
		{
			std::vector<int> values;
			bEnd = readIndexBlock(pFile, values);
			if(bEnd || values.size() < DIM)
			{
				break;
			}

			VecNX<Unit, DIM>* pVec = new VecNX<Unit, DIM>();
			for(size_t i = 0; i < DIM; ++i)
			{
				(*pVec)[i] = values[i];
			}

			indexList.push_back(pVec);		

			if(values.size() == 4 && DIM == 3)
			{
				VecNX<Unit, DIM>* pVec2 = new VecNX<Unit, DIM>();
				(*pVec2)[0] = values[0];
				(*pVec2)[1] = values[2];
				(*pVec2)[2] = values[3];

				indexList.push_back(pVec2);		
			}
		}
	}

	template<typename Unit, size_t DIM>
	bool readDataBlock(FILE* pFile, const Type& type, char* cstrOutput, VecNX<Unit, DIM>& vec)
	{
		for(size_t i = 0; i < DIM; ++i)
		{
			Unit value = static_cast<Unit>(0.0);
			if(readNextNode(pFile, cstrOutput) || !convertToNumber(cstrOutput, type, value)
				|| isEqual(cstrOutput, "}") || isEqual(cstrOutput, "]"))
			{
				return false;
			}

			vec[i] = value;
		}

		return true;
	}

	bool readIndexBlock(FILE* pFile, std::vector<int>& values)
	{
		bool bEnd(false);
		while(!bEnd)
		{
			char cstrTmpOutput[20] = {0};

			bEnd = getIndexNumberString(pFile, cstrTmpOutput);

			size_t len = strlen(cstrTmpOutput);
			if(len>0)
			{
				int value(0);
				if(!convertToNumber(cstrTmpOutput, INTEGER, value))
				{
					return true;
				}

				if(value==-1)
				{
					return false;
				}

				values.push_back(value);
			}
			
			if(bEnd)
			{
				return true;
			}
		}

		return bEnd;
	}

	bool getIndexNumberString(FILE* pFile, char* cstrOutput)
	{
		char cstrTmpOutput[20] = {0};

		bool bEnd(false);
		while(!bEnd)
		{
			size_t len = strlen(cstrTmpOutput);

			char currChar = fgetc(pFile);
			if(currChar == ',' || currChar == ' ' || currChar == '\n')
			{
				if(len == 0)
				{
					continue;
				}
				else
				{
					break;
				}
			}
			else if(currChar == '}' || currChar == ']' || currChar == EOF)
			{
				bEnd = true;
				break;
			}

			cstrTmpOutput[len] = currChar;
		}

		size_t len = strlen(cstrTmpOutput);
		for(size_t i = 0; i < len; ++i)
		{
			cstrOutput[i] = cstrTmpOutput[i];
		}

		return bEnd;
	}

	//! Reads next block until it reaches the specified character.
	//! Stops reading if the end of file is reached.
	//! Stops reading if "}" or "]" is reached.
	bool readNodeBlocksUntil(FILE* pFile, const char* cstrEnd, char* cstrOutput)
	{
		bool bEnd(false);
		do
		{
			bEnd = readNextNode(pFile, cstrOutput);
		}
		while(!bEnd && !isEqual(cstrOutput, cstrEnd) 
				&& !isEqual(cstrOutput, "}") && !isEqual(cstrOutput, "]"));

		return bEnd;
	}

	bool readBlockUntil(FILE* pFile, const char* cstrEnd, char* cstrOutput)
	{
		const size_t numStopChar = strlen(cstrEnd);
		if(numStopChar == 0)
		{
			return false;
		}

		char cstrTmpOutput[1000] = {0};

		bool bEnd(false);

		bool bStop(false);
		while(!bStop)
		{
			const size_t len = strlen(cstrTmpOutput);

			char currChar = fgetc(pFile);
			for(size_t i = 0; i < numStopChar; ++i)
			{
				if(currChar == cstrEnd[i])
				{
					bStop = true;
					break;
				}
			}

			if(currChar!= EOF)
			{
				cstrTmpOutput[len] = currChar;
			}
			else
			{
				bEnd = true;
			}
		}

		size_t len = strlen(cstrTmpOutput);
		for(size_t i = 0; i < len; ++i)
		{
			cstrOutput[i] = cstrTmpOutput[i];
		}

		return bEnd;
	}

	bool isEqual(const char* cstrS1, const char* cstrS2)
	{
		return strcmp(cstrS1, cstrS2) == 0;
	}
};

#endif