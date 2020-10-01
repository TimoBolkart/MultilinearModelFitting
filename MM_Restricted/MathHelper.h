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

#ifndef MathHelper_H
#define MathHelper_H

#include "Definitions.h"
#include "DataContainer.h"
#include "VectorNX.h"

class MathHelper
{
public:
	//! Compute normals of a mesh, based on Max1999 - Weights for Computing Vertex Normals from Facet Normals.
	//! \param poly				mesh the normals should be computed for
	//! \param vertexNormals	computed vertex normals, empty if poly does not contain any faces
	static void calcVertexNormals(const DataContainer& poly, std::vector<double>& vertexNormals);

	//! Compute projection of point p1 into a plane defined by a point p2 and a normal n2.
	//! \param p1			point that should be projection
	//! \param p2			point wihin the tangential plane
	//! \param n2			normal vector of the trangential plane
	//! \param outPoint	projection of p1 within the tangential plane
	static void getPlaneProjection(const Vec3d& p1, const Vec3d& p2, const Vec3d& n2, Vec3d& outPoint);

	//! Initialize transformation by the identity transformation.
	//! \param s		scaling s=1
	//! \param R		3x3 identity matrix 
	//! \param t		zero length translation vector t=(0, 0, 0)
	static void initTransformation(double& s, std::vector<double>& R, std::vector<double>& t);

	//! Computes inverse transformation (x = invs*invR*y + invt) of given rigid transformation (y = s*M*x + v).
	//! \param s				scaling
	//! \param R				3x3 rotation matrix
	//! \param sstrRotOp		defines if M = R (sstrRotOp = "N") or M = R^T (sstrRotOp = "T")
	//! \param t				translation vector
	//! \param sstrTransOp	defines if v = +t (sstrTransOp = "+") or v = -t (sstrTransOp = "-")
	//! \param invs			scaling of the invert transformation
	//! \param invR			rotation matrix of the invert transformation
	//! \param invt			transformation of the invert transformation
	static void invertTransformation(const double s, const std::vector<double>& R, const std::string& sstrRotOp, const std::vector<double>& t, const std::string& sstrTransOp, double& invs, std::vector<double>& invR, std::vector<double>& invt);

	//! Transform mesh vertices (y = s*M*x + v). 
	//! \param s				scaling
	//! \param R				3x3 rotation matrix
	//! \param sstrRotOp		defines if M = R (sstrRotOp = "N") or M = R^T (sstrRotOp = "T")
	//! \param t				translation vector
	//! \param sstrTransOp	defines if v = +t (sstrTransOp = "+") or v = -t (sstrTransOp = "-")
	//! \param mesh			mesh to be tranformed
	static void transformMesh(const double s, const std::vector<double>& R, const std::string& sstrRotOp, const std::vector<double>& t, const std::string& sstrTransOp, DataContainer& mesh);

	//! Transform data vector (y = s*M*x + v). 
	//! \param s				scaling
	//! \param R				3x3 rotation matrix
	//! \param sstrRotOp		defines if M = R (sstrRotOp = "N") or M = R^T (sstrRotOp = "T")
	//! \param t				translation vector
	//! \param sstrTransOp	defines if v = +t (sstrTransOp = "+") or v = -t (sstrTransOp = "-")
	//! \param data			data vector to be tranformed
	static void transformData(const double s, const std::vector<double>& R, const std::string& sstrRotOp, const std::vector<double>& t, const std::string& sstrTransOp, std::vector<double>& data);

	//! Transform data into local coordinate system of model using landmarks (+ some rigid ICP refinement steps, if enabled).
	//! \param modelData						vertices in the model coordinate system
	//! \param modelNormals					normals of the model
	//! \param modelLandmarks				landmarks at the model coordinate system
	//! \param modelLandmarksLoaded		flags which model landmarks are loaded
	//! \param data							data that should be transformed
	//! \param dataNormals					data normals
	//! \param dataLandmarks				landmarks of the data
	//! \param dataLandmarksLoaded		flags which data landmarks are loaded
	//! \param s								computed scaling
	//! \param R								computed 3x3 rotation matrix
	//! \param t								computed translaten vector
	//! \param bICPRefinement				flag if ICP should be used to refine the computed transformation
	//! \param bScaling						flag if the computed transformation should contain scaling
	//! \return true if alignment computation was successful
	static bool alignData(const std::vector<double>& modelData, const std::vector<double>& modelNormals, const std::vector<double>& modelLandmarks, const std::vector<bool>& modelLandmarksLoaded, std::vector<double>& data, std::vector<double>& dataNormals
								, const std::vector<double>& dataLandmarks, const std::vector<bool>& dataLandmarksLoaded, double& s, std::vector<double>& R, std::vector<double>& t, bool bICPRefinement = true, bool bScaling = true);

	//! Computes a rigid transformation between two sets of corresponding landmarks, that transforms lmks2 into the coordinate system of lmks1.
	//! \param lmks1							first set of landmarks
	//! \param lmks2Loaded					flags which landmarks are loaded 
	//! \param lmks2							second set of landmarks
	//! \param lmks2Loaded					flags which landmarks are loaded
	//! \param s								computed scaling
	//! \param R								computed 3x3 rotation matrix
	//! \param t								computed translaten vector
	//! \return true if alignment computation was successful
	static bool calcRigidLandmarkAlignment(const std::vector<double>& lmks1, const std::vector<bool>& lmks1Loaded, const std::vector<double>& lmks2, const std::vector<bool>& lmks2Loaded
														, double& s, std::vector<double>& R, std::vector<double>& t, bool bScaling = true);

	//! Computes scaling, rotation and translation for modelData = s*R*sequenceData+t
	static bool calcAlignmentTrafo(const std::vector<double>& modelData, const std::vector<double>& sequenceData, double& s, std::vector<double>& R, std::vector<double>& t, bool bScaling = true);
	
	//! Compute the ICP refinement for a given rigid transformation (|s*R*data+t - modelData|^2 -> min).
	//! \param modelData			vertices of the model
	//! \param modelNormals		normals of the model
	//! \param data				vertices of the data that should be transformed
	//! \param dataNormals		normals of the data that should be transformed
	//! \param s					input: initial scaling factor			output: refined scaling
	//! \param R					input: initial rotation matrix		output: refined rotation matrix
	//! \param t					input: initial translation vector	output: refined translation vector
	//! \return	true if refinement computation was successful
	static bool calcICPRefinenment(const std::vector<double>& modelData, const std::vector<double>& modelNormals, const std::vector<double>& data, const std::vector<double>& dataNormals, const double maxIPCPointDist
											, double& s, std::vector<double>& R, std::vector<double>& t, bool bScaling = true);

	//! Compute scaling that scales the source vertices to the target vertices.
	//! \param source			source data
	//! \param target			target data
	//! \param s				computed scaling factor
	//! \return true if successful
	static bool computeScaling(const std::vector<double>& source, const std::vector<double>& target, double& s);

	//! Calculating least squares rotation matrix R, that R*source = target
	//! \param source			source data
	//! \param target			target data
	//! \param R				computed rotation matrix
	//! \return true if successful
	static bool computeRotationAlignmentMatrix(const std::vector<double>& source, const std::vector<double>& target, std::vector<double>& R);

	//! Computes the 3d mean vector of the data and subtracts the mean of each 3d data vertex.
	//! \param data			data that are centered
	//! \param mean			computed mean
	static void centerData(std::vector<double>& data, std::vector<double>& mean);

	//! Computes the dataDim-dimensional mean vector of the data and subtracts the mean of each dataDim-dimensional data point.
	//! \param data				data that should be centered
	//! \param dataDim			dimension of one data point		
	//! \param centeredData		computed centered data
	//! \param mean				computed mean
	static void centerData(const std::vector<double>& data, const size_t dataDim, std::vector<double>& centeredData, std::vector<double>& mean);

	//! Rotates all vertices of a mesh.
	//! \param R			3x3 rotation matrix
	//! \param sstrOp		defines if multiplied by R (sstrOp = "N") or R^T (sstrOp = "T")
	//! \param mesh		rotated mesh
	static void rotateMesh(const std::vector<double>& R, const std::string& sstrOp, DataContainer& mesh);

	//! Rotates all vertices of a data vector.
	//! \param R			3x3 rotation matrix
	//! \param sstrOp		defines if multiplied by R (sstrOp = "N") or R^T (sstrOp = "T")
	//! \param mesh		rotated data
	static void rotateData(const std::vector<double>& R, const std::string& sstrOp, std::vector<double>& data);

	//! Translates mesh.
	//! \param t			3d translation vector
	//! \param sstrOp		defines if multiplied translated by t (sstrOp = "+")  or -t (sstrOp = "-")
	//! \param mesh		translated mesh
	static void translateMesh(const std::vector<double>& t, const std::string& sstrOp, DataContainer& mesh);

	//! Translates all vertices of a data vector.
	//! \param t			3d translation vector
	//! \param sstrOp		defines if multiplied translated by t (sstrOp = "+")  or -t (sstrOp = "-")
	//! \param mesh		translated data
	static void translateData(const std::vector<double>& t, const std::string& sstrOp, std::vector<double>& data);

	//! Scales mesh.
	//! \param factor		scaling factor
	//! \param mesh		scaled mesh
	static void scaleMesh(const double factor, DataContainer& mesh);

	//! Scales all vertices of a data vector.
	//! \param factor		scaling factor
	//! \param mesh		scaled data
	static void scaleData(const double factor, std::vector<double>& data);

	//! Invert squared matrix with dimension dim x dim
	//! \param matrix			squared matrix that should be inverted
	//! \param dim				dimension of the matrix
	//! \param invMatrix		inverted output matrix
	//! \return true if successful
	static bool invertMatrix(const std::vector<double>& matrix, const size_t dim, std::vector<double>& invMatrix);

	static void matrixMult(const std::vector<double>& MA, const size_t numARows, const size_t numACols, const std::string sstrTransA, const std::vector<double>& MB, const size_t numBCols, std::vector<double>& AB);

	//! Remove all vertices that are not part of any triangle.
	//! \param mesh		mesh that should be cleaned
	static void cleanMesh(DataContainer& mesh);
};

#endif