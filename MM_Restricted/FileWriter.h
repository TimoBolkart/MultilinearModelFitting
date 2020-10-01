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

#ifndef FILEWRITER_H
#define FILEWRITER_H

#include "DataContainer.h"

#include <string>
#include <vector>

class FileWriter
{
public:
	//! Save geometry file. Extend this to support other file formats.
	//! \param sstrFileName		full file name of the exported file
	//! \param data				data that should be exported
	//! \return true if successful
	static bool saveFile(const std::string& sstrFileName, const DataContainer& data);

private:

	//! Specific writer for off files.
	//! \param sstrFileName		full file name of the exported file
	//! \param data				data that should be exported
	//! \return true if successful
	static bool writeOff(const std::string& sstrFileName, const DataContainer& data);
};

#endif