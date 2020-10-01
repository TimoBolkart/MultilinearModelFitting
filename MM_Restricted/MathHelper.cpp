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

#include "MathHelper.h"
#include "KDTree3.h"

#include <vector>
#include <map>
#include <set>

#include <iostream>
#include <sstream>

namespace clapack
{
	extern "C"
	{
		#include "blaswrap.h"
		#include "f2c.h"
		extern int dgemm_(char *transa, char *transb, integer *m, integer *n
								, integer *k, doublereal *alpha, doublereal *a, integer *lda
								, doublereal *b, integer *ldb, doublereal *beta, doublereal *c
								, integer *ldc);
		extern int dgels_(char *trans, integer *m, integer *n, integer *nrhs, doublereal *a
								, integer *lda, doublereal *b, integer *ldb, doublereal *work
								, integer *lwork, integer *info);
		extern int dgesdd_(char *jobz, integer *m, integer *n, doublereal *a
								, integer *lda, doublereal *s, doublereal *u, integer *ldu
								, doublereal *vt, integer *ldvt, doublereal *work, integer *lwork
								, integer *iwork, integer *info);
		extern int dgetrf_(integer *m, integer *n, doublereal *a, integer *lda
								, integer *ipiv, integer *info);
		extern int dgetri_(integer *n, doublereal *a, integer *lda, integer *ipiv
								, doublereal *work, integer *lwork, integer *info);	
	}
}

//Compute normals, based on Max1999 - Weights for Computing Vertex Normals from Facet Normals
void MathHelper::calcVertexNormals(const DataContainer& poly, std::vector<double>& vertexNormals)
{
	vertexNormals.clear();

	const std::vector<double>& vertexList = poly.getVertexList();
	const std::vector<Vec3i*>& vertexIndexList = poly.getVertexIndexList();
	if(vertexIndexList.empty())
	{
		return;
	}

	std::map<int, std::vector<Vec3i*>> pointPolygonMap;
	{
		std::vector<Vec3i*>::const_iterator currIndexIter = vertexIndexList.begin();
		std::vector<Vec3i*>::const_iterator endIndexIter = vertexIndexList.end();
		for( ; currIndexIter!=endIndexIter; ++currIndexIter)
		{
			Vec3i* indexVec = *currIndexIter;

			for(size_t i = 0; i < 3; ++i)
			{
				int currIndex = (*indexVec)[i];

				std::map<int, std::vector<Vec3i*>>::iterator mapIter = pointPolygonMap.find(currIndex);
				if(mapIter!=pointPolygonMap.end())
				{
					std::vector<Vec3i*>& vec = mapIter->second;
					vec.push_back(indexVec);
				}
				else
				{
					std::vector<Vec3i*> vec;
					vec.push_back(indexVec);
					pointPolygonMap.insert(std::make_pair(currIndex, vec));
				}
			}
		}
	}

	const size_t vertexNormalDim = vertexList.size();
	vertexNormals.clear();
	vertexNormals.reserve(vertexNormalDim);

	for(size_t i = 0; i < vertexNormalDim; ++i)
	{
		vertexNormals.push_back(0.0);
	}
	
	std::map<int, std::vector<Vec3i*>>::const_iterator currPointPolyIter = pointPolygonMap.begin();
	std::map<int, std::vector<Vec3i*>>::const_iterator endPointPolyIter = pointPolygonMap.end();
	for(; currPointPolyIter != endPointPolyIter; ++currPointPolyIter)
	{
		const int currIndex = currPointPolyIter->first;
		const std::vector<Vec3i*>& neighborPolygons = currPointPolyIter->second;

		const Vec3d currVertex(vertexList[3*currIndex], vertexList[3*currIndex+1], vertexList[3*currIndex+2]);
		Vec3d vertexNormal(0.0, 0.0, 0.0);

		for(size_t i = 0; i < neighborPolygons.size(); ++i)
		{
			const Vec3i& currPoly = (*neighborPolygons[i]);

			int tmpPos(-1);
			for(size_t j = 0; j < 3; ++j)
			{
				tmpPos = currPoly[j] == currIndex ? static_cast<int>(j) : tmpPos;
			}

			if(tmpPos == -1 || currPoly[tmpPos] != currIndex)
			{
				std::cout << "Wrong neighboring index while computing normal vector" << std::endl;
				continue;
			}

			const int prevIndex = currPoly[(tmpPos+2)%3];
			const int nextIndex = currPoly[(tmpPos+1)%3];

			const Vec3d prevVertex(vertexList[3*prevIndex], vertexList[3*prevIndex+1], vertexList[3*prevIndex+2]);
			const Vec3d nextVertex(vertexList[3*nextIndex], vertexList[3*nextIndex+1], vertexList[3*nextIndex+2]);

			const Vec3d v1 = nextVertex-currVertex;
			const Vec3d v2 = prevVertex-currVertex;

			const double v1SqrLength = v1.sqrLength();
			const double v2SqrtLength = v2.sqrLength();

			if(v1SqrLength < DBL_EPSILON || v2SqrtLength < DBL_EPSILON)
			{
				std::cout << "Zero length edge while computing normal vector" << std::endl;
				continue;
			}

			const double factor = 1.0/(v1SqrLength*v2SqrtLength);

			Vec3d v1xv2;
			v1.crossProduct(v2, v1xv2);
			vertexNormal += v1xv2*factor;
		}
		
		if(!vertexNormal.normalize())
		{
			std::cout << "Zero length normal vector while computing normal vector" << std::endl;
			continue;
		}
		
		vertexNormals[3*currIndex] = vertexNormal[0];
		vertexNormals[3*currIndex+1] = vertexNormal[1];
		vertexNormals[3*currIndex+2] = vertexNormal[2];
	}
}

void MathHelper::getPlaneProjection(const Vec3d& p1, const Vec3d& p2, const Vec3d& n2, Vec3d& outPoint)
{
	Vec3d diff = p2-p1;
	double d = diff.dotProduct(n2);

	outPoint = n2;
	outPoint.scalarMult(d);

	outPoint = p1 + outPoint;
}

void MathHelper::initTransformation(double& s, std::vector<double>& R, std::vector<double>& t)
{
	s = 1.0;

	R.clear();
	R.reserve(9);

	for(size_t i = 0; i < 9; ++i)
	{
		R.push_back(i%4 == 0 ? 1.0 : 0.0);
	}

	t.clear();
	t.reserve(3);

	for(size_t i = 0; i < 3; ++i)
	{
		t.push_back(0.0);
	}
}

void MathHelper::invertTransformation(const double s, const std::vector<double>& R, const std::string& sstrRotOp, const std::vector<double>& t, const std::string& sstrTransOp
									, double& invs, std::vector<double>& invR, std::vector<double>& invt)
{
	const std::string sstrTmpTransOp = sstrTransOp == "+" ? "-" : "+";
	const std::string sstrTmpRotOp = sstrRotOp == "N" ? "T" : "N";

	MathHelper::initTransformation(invs, invR, invt);

	// s* = 1/s
	invs *= 1/s;
	
	// R* = R^T
	MathHelper::rotateData(R, sstrTmpRotOp, invR);

	// t* = - 1/s * R^T * t
	MathHelper::translateData(t, sstrTmpTransOp, invt);
	MathHelper::rotateData(R, sstrTmpRotOp, invt);
	MathHelper::scaleData(1/s, invt);
}

void MathHelper::transformMesh(const double s, const std::vector<double>& R, const std::string& sstrRotOp, const std::vector<double>& t, const std::string& sstrTransOp, DataContainer& mesh)
{
	std::vector<double> vertices =  mesh.getVertexList();

	MathHelper::transformData(s, R, sstrRotOp, t, sstrTransOp, vertices);

	mesh.setVertexList(vertices);
}

void MathHelper::transformData(const double s, const std::vector<double>& R, const std::string& sstrRotOp, const std::vector<double>& t, const std::string& sstrTransOp, std::vector<double>& data)
{
	MathHelper::rotateData(R, sstrRotOp, data);
	MathHelper::scaleData(s, data);
	MathHelper::translateData(t, sstrTransOp, data);
}

bool MathHelper::alignData(const std::vector<double>& modelData, const std::vector<double>& modelNormals, const std::vector<double>& modelLandmarks, const std::vector<bool>& modelLandmarksLoaded, std::vector<double>& data
									, std::vector<double>& dataNormals, const std::vector<double>& dataLandmarks, const std::vector<bool>& dataLandmarksLoaded, double& s, std::vector<double>& R, std::vector<double>& t, bool bICPRefinement, bool bScaling)
{
	double sLmk(0.0);
	std::vector<double> RLmk;
	std::vector<double> tLmk;
	if(!MathHelper::calcRigidLandmarkAlignment(modelLandmarks, modelLandmarksLoaded, dataLandmarks, dataLandmarksLoaded, sLmk, RLmk, tLmk, bScaling))
	{
		return false;
	}

	s = sLmk;

	t.clear();
	t = tLmk;

	R.clear();
	R = RLmk;

	if(bICPRefinement)
	{
		if(!MathHelper::calcICPRefinenment(modelData, modelNormals, data, dataNormals, ICP_MAX_POINT_DISTANCE, s, R, t, bScaling))
		{
			return false;
		}
	}

	MathHelper::transformData(s, R, "N", t, "+", data);

	if(!dataNormals.empty())
	{
		MathHelper::rotateData(R, "N", dataNormals);
	}

	return true;
}

bool MathHelper::calcRigidLandmarkAlignment(const std::vector<double>& lmks1, const std::vector<bool>& lmks1Loaded, const std::vector<double>& lmks2, const std::vector<bool>& lmks2Loaded
															, double& s, std::vector<double>& R, std::vector<double>& t, bool bScaling)
{
	if(lmks1.size() != lmks2.size())
	{
		std::cout << "MathHelper::calcRigidLandmarkAlignment(...) - Number of landmarks does not match " << lmks1.size() << " != " << lmks2.size() << std::endl;
		return false;
	}

	std::vector<double> modelLandmarks;
	std::vector<double> targetLandmarks;

	//Use the first 8 landmarks for rigid alignment
	for(size_t i = 0; i < 8; ++i)
	{
		if(lmks1Loaded[i] && lmks2Loaded[i])
		{
			for(size_t j = 0; j < 3; ++j)
			{
				modelLandmarks.push_back(lmks1[3*i+j]);
				targetLandmarks.push_back(lmks2[3*i+j]);
			}
		}
	}

	return MathHelper::calcAlignmentTrafo(modelLandmarks, targetLandmarks, s, R, t, bScaling);
}

bool MathHelper::calcAlignmentTrafo(const std::vector<double>& modelData, const std::vector<double>& sequenceData, double& s, std::vector<double>& R, std::vector<double>& t, bool bScaling)
{
	std::vector<double> tmpModelData = modelData;
	std::vector<double> tmpSequenceData = sequenceData;

	std::vector<double> modelMean;
	MathHelper::centerData(tmpModelData, modelMean);

	std::vector<double> sequenceMean;
	MathHelper::centerData(tmpSequenceData, sequenceMean);
	
	s = 1.0;

	if(bScaling)
	{
		if(!computeScaling(tmpSequenceData, tmpModelData, s))
		{
			return false;
		}

		MathHelper::scaleData(s, tmpSequenceData);
	}

	if(!MathHelper::computeRotationAlignmentMatrix(tmpSequenceData, tmpModelData, R))
	{
		return false;
	}

	// t = -s*R*modelMean + sequenceMean
	t.clear();
	t = sequenceMean;

	MathHelper::transformData(-s, R, "N", modelMean, "+", t);

	return true;
}

bool MathHelper::calcICPRefinenment(const std::vector<double>& modelData, const std::vector<double>& modelNormals, const std::vector<double>& data, const std::vector<double>& dataNormals, const double maxIPCPointDist
													, double& s, std::vector<double>& R, std::vector<double>& t, bool bScaling)
{
	if(fabs(s) < DBL_EPSILON || R.empty() || t.empty())
	{
		std::cout << "No refinement step possible due to wrong input transformation" << std::endl;
		return false;
	}

	//If normals are not valid, they are ignored
	bool bHasValidNormals = (modelData.size() == modelNormals.size()) && (data.size() == dataNormals.size());

	const size_t modelDataSize = modelData.size();
	const size_t numModelPoints = modelDataSize/3;

	double invs(1.0);
	std::vector<double> invR;
	std::vector<double> invt;
	MathHelper::invertTransformation(s, R, "N", t, "+", invs, invR, invt);

	KDTree3 kdTree(data);

	//No need to scale the max point distance while iterating, since we compute nearest neighbors in target data system.
	//Since the maximum point distance is specified in the model system, it need to be scaled
	const double sqrMaxDist(pow(invs*maxIPCPointDist, 2));

	for(int iter = 0; iter < ICP_NUMBER_ALIGNMENT_ITERATIONS; ++iter)
	{
		std::vector<double> copyModelData = modelData;
		std::vector<double> copyModelNormals = modelNormals;

		//Transform model points and normals into the local coordinate system of the target data.
		MathHelper::transformData(invs, invR, "N", invt, "+", copyModelData);
		MathHelper::rotateData(invR, "N", copyModelNormals);

		std::vector<double> modelPointSelection;
		std::vector<double> nnDataPointSelection;

		for(size_t j = 0; j < numModelPoints; ++j)
		{
			const Vec3d sourcePoint(copyModelData[3*j], copyModelData[3*j+1], copyModelData[3*j+2]);

			std::vector<double> queryPoint;
			queryPoint.push_back(sourcePoint[0]);
			queryPoint.push_back(sourcePoint[1]);
			queryPoint.push_back(sourcePoint[2]);

			int pointIndex(0);
			double sqrPointDist(0);
			kdTree.getNearestPoint(queryPoint, pointIndex, sqrPointDist);

			const int targetIndex = 3*pointIndex;

			bool bPointValid = (sqrPointDist <= sqrMaxDist);
			if(bHasValidNormals)
			{
				const Vec3d sourceNormal(copyModelNormals[3*j], copyModelNormals[3*j+1], copyModelNormals[3*j+2]);
				const Vec3d nnTargetNormal(dataNormals[targetIndex], dataNormals[targetIndex+1], dataNormals[targetIndex+2]);
				const double angle = sourceNormal.angle(nnTargetNormal);
				bPointValid &= (angle <= ICP_MAX_NORMAL_ANGLE);
			}

			if(bPointValid)
			{
				modelPointSelection.push_back(sourcePoint[0]);
				modelPointSelection.push_back(sourcePoint[1]);
				modelPointSelection.push_back(sourcePoint[2]);

				const Vec3d nnTargetPoint(data[targetIndex], data[targetIndex+1], data[targetIndex+2]);

				if(bHasValidNormals)
				{
					const Vec3d nnTargetNormal(dataNormals[targetIndex], dataNormals[targetIndex+1], dataNormals[targetIndex+2]);

					Vec3d planeProjectionPoint;					
					MathHelper::getPlaneProjection(sourcePoint, nnTargetPoint, nnTargetNormal, planeProjectionPoint);

					nnDataPointSelection.push_back(planeProjectionPoint[0]);
					nnDataPointSelection.push_back(planeProjectionPoint[1]);
					nnDataPointSelection.push_back(planeProjectionPoint[2]);
				}
				else
				{
					nnDataPointSelection.push_back(nnTargetPoint[0]);
					nnDataPointSelection.push_back(nnTargetPoint[1]);
					nnDataPointSelection.push_back(nnTargetPoint[2]);
				}
			}
		}

		if(modelPointSelection.size() < 12)
		{
			std::cout << "MathHelper::alignData() - No refinement step possible" << std::endl;
			return false;
		}

		double sModel(1.0);
		std::vector<double> RModel;
		std::vector<double> tModel;
		if(!MathHelper::calcAlignmentTrafo(nnDataPointSelection, modelPointSelection, sModel, RModel, tModel, bScaling))
		{
			return false;
		}

		invs*=sModel;
		MathHelper::rotateData(RModel, "N", invR);
		MathHelper::transformData(sModel, RModel, "N", tModel, "+", invt);
	}

	MathHelper::invertTransformation(invs, invR, "N", invt, "+", s, R, t);

	return true;
}

bool MathHelper::computeScaling(const std::vector<double>& source, const std::vector<double>& target, double& s)
{
	if(source.size() != target.size() || source.size()%3 != 0)
	{
		s = 1.0;
		return false;
	}
		
	const size_t dataSize = source.size();
	const size_t numPoints = source.size()/3;

	double sourceDist(0.0);
	double targetDist(0.0);

/**/
	for(size_t i = 0; i < dataSize; ++i)
	{
		sourceDist += pow(source[i], 2);
		targetDist += pow(target[i], 2);
	}

	sourceDist = sqrt(sourceDist/static_cast<double>(numPoints));
	targetDist = sqrt(targetDist/static_cast<double>(numPoints));

	s = sourceDist > 0.0 ? targetDist / sourceDist : 1.0;
/**/

/*
	for(size_t i = 0; i < numPoints; ++i)
	{
		double sqrModelPointDist(0.0);
		double sqrtSequencePointDist(0.0);
		for(size_t j = 0; j < 3; ++j)
		{
			sqrModelPointDist += pow(tmpModelData[i], 2);
			sqrtSequencePointDist += pow(tmpSequenceData[i], 2);
		}

		modelDist += sqrt(sqrModelPointDist);
		sequenceDist += sqrt(sqrtSequencePointDist);
	}

	//modelDist /= static_cast<double>(numPoints);
	//sequenceDist /= static_cast<double>(numPoints);

	s = sequenceDist > 0.0 ? modelDist / sequenceDist : 1.0;
	MathHelper::scaleData(s, tmpSequenceData);
*/
	return true;
}

bool MathHelper::invertMatrix(const std::vector<double>& matrix, const size_t dim, std::vector<double>& invMatrix)
{
	const size_t matrixSize = matrix.size();
	
	double* A = new double[matrixSize];
	
	for(size_t i = 0; i < matrixSize; ++i)
	{
		A[i] = matrix[i];
	}

	long int n = static_cast<long int>(dim);
	long int* ipiv = new long int [n];
	double* work = new double[n];
	long int info;
	clapack::dgetrf_(&n, &n, A, &n, ipiv, &info);
	if(info != 0)
	{
		delete [] A;
		delete [] ipiv;
		delete [] work;
		return false;
	}

	clapack::dgetri_(&n, A, &n, ipiv, work, &n, &info);
	if(info != 0)
	{
		delete [] A;
		delete [] ipiv;
		delete [] work;
		return false;
	}

	invMatrix.clear();
	invMatrix.resize(matrixSize);
	for(size_t i = 0; i < matrixSize; ++i)
	{
		invMatrix[i] = A[i];
	}

	delete [] A;
	delete [] ipiv;
	delete [] work;
	return true;
}

void MathHelper::matrixMult(const std::vector<double>& MA, const size_t numARows, const size_t numACols, const std::string sstrTransA, const std::vector<double>& MB, const size_t numBCols, std::vector<double>& AB)
{
	char transA = sstrTransA == "N" ? 'N' : 'T';
	char transB = 'N';

	//long int m = sstrTransA == "N" ? static_cast<long int>(numARows) : static_cast<long int>(numACols);
	long int m = static_cast<long int>(numARows);
	long int n = static_cast<long int>(numBCols);
	//long int k = sstrTransA == "N" ? static_cast<long int>(numACols) : static_cast<long int>(numARows);
	long int k = static_cast<long int>(numACols);
	double alpha = 1.0;

	long int lda = sstrTransA == "N" ? m : k;
	long int ldb = k;
	long int ldc = m;

	double* A = new double[m*k];
	for(int i = 0; i < m*k; ++i)
	{
		A[i] = MA[i];
	}


	double* B = new double[k*n];
	for(int i = 0; i < k*n; ++i)
	{
		B[i] = MB[i];
	} 


	double* C = new double[m*n];

	double beta = 0.0;

	clapack::dgemm_(&transA, &transB, &m, &n, &k, &alpha, A, &lda, B, &ldb, &beta, C, &ldc);

	AB.reserve(m*n);
	for(int i = 0; i < m*n; ++i)
	{
		AB.push_back(C[i]);
	} 

	delete [] A;
	delete [] B;
	delete [] C;
}

bool MathHelper::computeRotationAlignmentMatrix(const std::vector<double>& source, const std::vector<double>& target, std::vector<double>& R)
{
	const int dim = 3;
	const size_t minSampleSize = min(source.size(), target.size());
	const int minNumPoints = static_cast<int>(minSampleSize)/dim;

	char trans = 'N';
	long int m = static_cast<long int>(minNumPoints);
	long int n = static_cast<long int>(dim);

	double* X = new double[n*n];

	{
		long int nrhs = n;

		double* A = new double[m*n];
		double* B = new double[m*nrhs];

		for(int i = 0; i < m; i++)
		{
			A[i] = source[dim*i];
			A[m+i] = source[dim*i+1];
			A[2*m+i] = source[dim*i+2];

			B[i] = target[dim*i];
			B[m+i] = target[dim*i+1];
			B[2*m+i] = target[dim*i+2];
		}

		long int lda = m;
		long int lbd = m;

		long int lwork = max(1, 2*m*n);
		double* work = new double[lwork];
	
		long int info = 0;

		//Estimate rotation matrix by A*X=B (points are rows of A and B)
		clapack::dgels_(&trans, &m, &n, &nrhs, A, &lda, B, &lbd, work, &lwork, &info);
		if(info != 0)
		{
			std::cout<<"Calculation of rotation matrix not successful "<< info << std::endl;
		
			delete [] A;
			delete [] B;
			delete [] work;
			return false;
		}

		for(int i = 0; i < n; ++i)
		{
			for(int j = 0; j < n; ++j)
			{
				X[i*n+j] = B[i*m+j];
			}
		}

		delete [] A;
		delete [] B;
		delete [] work;
	}

	double* U = new double[n*n];
	double* VT = new double[n*n];
	
	{
		char jobz = 'S';

		double* S = new double[n];

		long int ldu = n;
		long int ldvt = n;

		long int lwork = 7*n*n+4*n;
		double* work = new double[lwork];

		long int * iwork = new long int[8*n];
		long int info = 0;

		//Singular value decomposition of estimated rotation matrix (X=U*S*VT)
		clapack::dgesdd_(&jobz, &n, &n, X, &n, S, U, &ldu, VT, &ldvt, work, &lwork, iwork, &info);
		if(info != 0)
		{
			std::cout<<"SVD decomposition of R not successful "<< info << std::endl;

			delete [] U;
			delete [] VT;
			delete [] S;
			delete [] work;
			delete [] iwork;
			delete [] X;
			return false;
		}

		delete [] S;
		delete [] work;
		delete [] iwork;
	}

	{
		char transU = 'N';
		char transVT = 'N';

		long int ldu = n;
		long int ldvt = n;
		long int ldr = n;

		double alpha = 1.0;
		double beta = 0.0;
		
		//X=U*VT
		clapack::dgemm_(&transU, &transVT, &n, &n, &n, &alpha, U, &ldu, VT, &ldvt, &beta, X, &ldr);
	}

	//disallow reflections
   const double det = X[0]*X[4]*X[8]+X[3]*X[7]*X[2]+X[6]*X[1]*X[5]-X[2]*X[4]*X[6]-X[5]*X[7]*X[0]-X[8]*X[1]*X[3];
	if(det < 0)
	{
		std::cout<<"Determinant is "<< det << std::endl;
		delete [] X;
		delete [] U;
		delete [] VT;
		return false;
	}

	R.clear();
	R.resize(n*n);

	//Transpose tmpR matrix because we return rotation matrix which is multiplied from left to the data 
	for(int i = 0; i < n; ++i)
	{
		for(int j = 0; j < n; ++j)
		{
			R[i*n+j] = X[j*n+i];
		}
	}

	delete [] X;
	delete [] U;
	delete [] VT;

	return true;
}

void MathHelper::centerData(std::vector<double>& data, std::vector<double>& mean)
{
	mean.clear();
	mean.push_back(0.0);
	mean.push_back(0.0);
	mean.push_back(0.0);

	const size_t numPoints = data.size()/3;

	for(size_t i = 0; i < numPoints; ++i)
	{
		for(size_t j = 0; j < 3; ++j)
		{
			mean[j] += data[i*3+j];
		}
	}

	const double invNumPoints(1.0/static_cast<double>(numPoints));
	for(size_t i = 0; i < 3; ++i)
	{
		mean[i] *= invNumPoints;
	}

	MathHelper::translateData(mean, "-", data);
}

void MathHelper::centerData(const std::vector<double>& data, const size_t dataDim, std::vector<double>& centeredData, std::vector<double>& mean)
{
	if(data.size() % dataDim != 0)
	{
		std::cout << "Data dimension for centering not correct" << std::endl;

		return;
	}

	mean.clear();
	mean.reserve(dataDim);
	for(int i = 0; i < dataDim; ++i)
	{
		mean.push_back(0.0);
	}

	const size_t numSamples = static_cast<size_t>(static_cast<double>(data.size())/static_cast<double>(dataDim));
	if(numSamples == 0)
	{
		return;
	}

	for(size_t i = 0; i < numSamples; ++i)
	{
		const size_t startIndex = i*dataDim;

		for(size_t j = 0; j < dataDim; ++j)
		{
			const size_t currIndex = startIndex+j;
			mean[j] += data[currIndex];
		}
	}

	const double factor = 1.0 / static_cast<double>(numSamples);
	for(size_t i = 0; i < dataDim; ++i)
	{
		mean[i] *= factor;
	}

	centeredData.reserve(numSamples*dataDim);
	for(size_t i = 0; i < numSamples; ++i)
	{
		const size_t startIndex = i*dataDim;

		for(size_t j = 0; j < dataDim; ++j)
		{
			const size_t currIndex = startIndex+j;
			centeredData.push_back(data[currIndex]-mean[j]);
		}
	}
}

void MathHelper::rotateMesh(const std::vector<double>& R, const std::string& sstrOp, DataContainer& mesh)
{
	std::vector<double> vertices = mesh.getVertexList();

	MathHelper::rotateData(R, sstrOp, vertices);

	mesh.setVertexList(vertices);
}

void MathHelper::rotateData(const std::vector<double>& R, const std::string& sstrOp, std::vector<double>& data)
{
	if(R.size() != 9)
	{
		return;
	}

	const int numPoints = static_cast<int>(data.size()/3);

	char transA = sstrOp == "N" ? 'N' : 'T';
	char transB = 'N';

	long int m = static_cast<long int>(3);
	long int n = static_cast<long int>(numPoints);
	long int k = static_cast<long int>(3);	
	double alpha = 1.0;

	long int lda = m;
	long int ldb = k;
	long int ldc = m;

	double* A = new double[m*m];
	for(int i = 0; i < m*m; ++i)
	{
		A[i] = R[i];
	}

	double* B = new double[m*n];
	for(int i = 0; i < m*n; ++i)
	{
		B[i] = data[i];
	} 

	double* C = new double[m*n];

	double beta = 0.0;

	clapack::dgemm_(&transA, &transB, &m, &n, &k, &alpha, A, &lda, B, &ldb, &beta, C, &ldc);

	data.clear();
	data.reserve(m*n);
	for(int i = 0; i < m*n; ++i)
	{
		data.push_back(C[i]);
	} 
	
	delete [] A;
	delete [] B;
	delete [] C;
}

void MathHelper::translateMesh(const std::vector<double>& t, const std::string& sstrOp, DataContainer& mesh)
{
	std::vector<double> vertices = mesh.getVertexList();

	MathHelper::translateData(t, sstrOp, vertices);

	mesh.setVertexList(vertices);
}

void MathHelper::translateData(const std::vector<double>& t, const std::string& sstrOp, std::vector<double>& data)
{
	if(t.size() != 3)
	{
		return;
	}

	const size_t numPoints = data.size() / 3;
	for(size_t i = 0; i < numPoints; ++i)
	{
		const size_t startIndex = 3*i;

		for(size_t j = 0; j < 3; ++j)
		{
			const double value = sstrOp == "+" ? data[startIndex+j]+t[j] : data[startIndex+j]-t[j];
			data[startIndex+j] = value;
		}
	}
}

void MathHelper::scaleMesh(const double factor, DataContainer& mesh)
{
	std::vector<double> vertices = mesh.getVertexList();

	MathHelper::scaleData(factor, vertices);

	mesh.setVertexList(vertices);
}

void MathHelper::scaleData(const double factor, std::vector<double>& data)
{
	for(size_t i = 0; i < data.size(); ++i)
	{
		data[i] *= factor;
	}
}

void MathHelper::cleanMesh(DataContainer& mesh)
{
	std::set<size_t> validPointIndices;

	const std::vector<Vec3i*> triangles = mesh.getVertexIndexList();
	for(size_t i = 0; i < triangles.size(); ++i)
	{
		const Vec3i& currTriangle = *(triangles[i]);
		if((currTriangle[0] == currTriangle[1])
			|| (currTriangle[1] == currTriangle[2])
			|| (currTriangle[2] == currTriangle[0]))
		{
			continue;
		}

		validPointIndices.insert(currTriangle[0]);
		validPointIndices.insert(currTriangle[1]);
		validPointIndices.insert(currTriangle[2]);
	}

	const std::vector<double>& meshVertices = mesh.getVertexList();

	std::vector<std::pair<int, int>> oldNewMap;

	size_t newId = 0;
	for(size_t id = 0; id < mesh.getNumVertices(); ++id)
	{
		if(validPointIndices.find(id) == validPointIndices.end())
		{
			oldNewMap.push_back(std::make_pair(static_cast<int>(id), -1));
		}
		else
		{
			oldNewMap.push_back(std::make_pair(static_cast<int>(id), static_cast<int>(newId)));
			++newId;
		}
	}

	DataContainer cleanMesh;
	std::vector<double> cleanMeshVertices;
	for(size_t oldId = 0; oldId < mesh.getNumVertices(); ++oldId)
	{
		const int newPointId = oldNewMap[oldId].second;
		if(newPointId != -1)
		{
			cleanMeshVertices.push_back(meshVertices[3*oldId]);
			cleanMeshVertices.push_back(meshVertices[3*oldId+1]);
			cleanMeshVertices.push_back(meshVertices[3*oldId+2]);
		}
	}
	
	cleanMesh.setVertexList(cleanMeshVertices);

	std::vector<Vec3i*>& cleanMeshTriangles = cleanMesh.getVertexIndexList();
	for(size_t i = 0; i < triangles.size(); ++i)
	{
		const Vec3i& oldTriangle = *(triangles[i]);
		if((oldTriangle[0] == oldTriangle[1])
			|| (oldTriangle[1] == oldTriangle[2])
			|| (oldTriangle[2] == oldTriangle[0]))
		{
			continue;
		}

		const int newAId = oldNewMap[oldTriangle[0]].second;
		const int newBId = oldNewMap[oldTriangle[1]].second;
		const int newCId = oldNewMap[oldTriangle[2]].second;
		if(newAId == -1 || newBId == -1 || newCId == -1)
		{
			std::cout << "Error" << std::endl;
		}

		Vec3i* newTriangle = new Vec3i(newAId, newBId, newCId);
		cleanMeshTriangles.push_back(newTriangle);
	}

	if(cleanMesh.getNumVertices() != mesh.getNumVertices())
	{
		std::cout << "Removed " << mesh.getNumVertices() - cleanMesh.getNumVertices() << " vertices" << std::endl;
		std::cout << "Removed " << mesh.getNumFaces() - cleanMesh.getNumFaces() << " faces" << std::endl;
	}

	mesh = cleanMesh;
}
