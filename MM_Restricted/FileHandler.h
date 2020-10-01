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

#ifndef FILEHANDLER_H
#define FILEHANDLER_H

#include <string>
class FileHandler
{
public:
	//! Get the file extension of a full file name.
	//! e.g. ...\...\Test.txt returns txt
	//! \param sstrFile			full file name
	//! \return file extension
	static std::string getFileExtension(const std::string& sstrFileName);

	//! Checks if a file name exists.
	//! \param sstrFile			full file name
	//! \return true if file exists
	static bool fileExist(const std::string& sstrFileName);

	//! Get the file name of a full file name.
	//! e.g. ...\...\Test.txt returns Test
	//! \param sstrFile			full file name
	//! \return file name if file exists
	static std::string getFileName(const std::string& sstrFullFileName);

	//! Get the file name of a full file name.
	//! e.g. ...\...\Test.txt returns ...\...
	//! \param sstrFile			full file name
	//! \return file path if file exists
	static std::string getFilePath(const std::string& sstrFileName);
};

#endif