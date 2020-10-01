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

#ifndef DATACONTAINER_H
#define DATACONTAINER_H

#include <vector>
#include "VectorNX.h"

//T					data type of vertex coordinates
//IndexDim			number of indices which define a face
//						e.g. 3 for triangle, 4 for quad,...
template<typename T, size_t IndexDim>
class ImportMesh
{
	//Indices always integer values
	typedef VecNX<int, IndexDim> IndexVec;

public:
	ImportMesh()
	{

	}

	ImportMesh(const std::vector<T>& vertexList, const std::vector<IndexVec*>& vertexIndexList)
	{
		copyPtrArray(vertexList, m_vertexList);
		copyPtrArray(vertexIndexList, m_vertexIndexList);
	}

	ImportMesh(const ImportMesh& mesh)
	{
		m_vertexList = mesh.getVertexList();
		copyPtrArray(mesh.getVertexIndexList(), m_vertexIndexList);
		copyPtrArray(mesh.getVertexColorList(), m_vertexColorList);
	}

	~ImportMesh()
	{
		clear();
	}

	ImportMesh& operator=(const ImportMesh& mesh)
	{
		clear();

		m_vertexList = mesh.getVertexList();
		copyPtrArray(mesh.getVertexIndexList(), m_vertexIndexList);
		copyPtrArray(mesh.getVertexColorList(), m_vertexColorList);

		return *this;
	}

	void clear()
	{
		m_vertexList.clear();
		clearPtrArray(m_vertexIndexList);
		clearPtrArray(m_vertexColorList);
	}

	size_t getNumVertices() const
	{
		return m_vertexList.size()/3;
	}

	size_t getNumFaces() const
	{
		return m_vertexIndexList.size();
	}

	void setVertexList(const std::vector<T>& vertices)
	{
		m_vertexList.clear();
		m_vertexList = vertices;
	}

	std::vector<T> getVertexList()
	{
		return m_vertexList;
	}

	const std::vector<T>& getVertexList() const
	{
		return m_vertexList;
	}

	std::vector<IndexVec*>& getVertexIndexList()
	{
		return m_vertexIndexList;
	}

	const std::vector<IndexVec*>& getVertexIndexList() const
	{
		return m_vertexIndexList;
	}

	void clearVertexColorList()
	{
		clearPtrArray(m_vertexColorList);
	}

	std::vector<Vec3d*>& getVertexColorList()
	{
		return m_vertexColorList;
	}

	const std::vector<Vec3d*>& getVertexColorList() const
	{
		return m_vertexColorList;
	}

private:
	template<typename ArrayType>
	void clearPtrArray(std::vector<ArrayType*>& vec)
	{
		std::vector<ArrayType*>::iterator currIter = vec.begin();
		std::vector<ArrayType*>::iterator endIter = vec.end();
		for(; currIter!=endIter; ++currIter)
		{
			if(*currIter != NULL)
			{
				delete *currIter;
			}
		}

		vec.clear();
	}

	template<typename ArrayType>
	void copyPtrArray(const std::vector<ArrayType*>& inVec, std::vector<ArrayType*>& outVec)
	{
		outVec.reserve(inVec.size());

		std::vector<ArrayType*>::const_iterator currIter = inVec.begin();
		std::vector<ArrayType*>::const_iterator endIter = inVec.end();
		for(; currIter!=endIter; ++currIter)
		{
			ArrayType* pNewType = new ArrayType(**currIter);
			outVec.push_back(pNewType);
		}
	}

	// Vertex data
	std::vector<T> m_vertexList;
	std::vector<IndexVec*> m_vertexIndexList;
	std::vector<Vec3d*> m_vertexColorList;
};

typedef ImportMesh<double,3> DataContainer;

#endif