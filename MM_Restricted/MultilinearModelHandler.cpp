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

#include "MultilinearModelHandler.h"	
#include "FileLoader.h"
#include "FileWriter.h"
#include "MathHelper.h"
#include "MMProjectionCostFunction.h"
#include "KDTree3.h"

#include <vnl/vnl_cost_function.h>
#include <vnl/algo/vnl_lbfgsb.h>


MultilinearModelHandler::MultilinearModelHandler()
{

}

MultilinearModelHandler::~MultilinearModelHandler()
{

}

void MultilinearModelHandler::clear()
{
	m_multilinearModel.clear();
	m_modeDimensions.clear();
	m_truncModeDimensions.clear();
	m_mean.clear();
	m_meanWeights.clear();
	m_neutralMeanWeights.clear();
}

void MultilinearModelHandler::reconstructForVariations(const std::vector<double>& variations, std::vector<double>& outPoints)
{
	if(m_truncModeDimensions.size()!=3)
	{
		std::cout << "Wrong number of mode diemensions " << m_truncModeDimensions.size() << " != " << 3 << std::endl;
		return;
	}

	const size_t m2 = m_truncModeDimensions[1];
	const size_t m3 = m_truncModeDimensions[2];
	if(variations.size() != m2+m3)
	{
		std::cout << "Wrong input variation vector dimension " << variations.size() << " != " << m2+m3 << std::endl;
		return;
	}

	std::vector<double> weightVector;
	for(size_t i = 0; i < m2; ++i)
	{
		double tmpWeight = getWeightForVariation(variations[i], m_meanWeights[i], 2);
		weightVector.push_back(tmpWeight);
	}

	for(size_t i = 0; i < m3; ++i)
	{
		double tmpWeight = getWeightForVariation(variations[m2+i], m_meanWeights[m2+i], 3);
		weightVector.push_back(tmpWeight);
	}

	reconstructForWeights(weightVector, outPoints);
}

void MultilinearModelHandler::reconstructForWeights(const std::vector<double>& weightVector, std::vector<double>& outPoints)
{
	if(m_truncModeDimensions.size()!=3)
	{
		std::cout << "Wrong number of mode diemensions " << m_truncModeDimensions.size() << " != " << 3 << std::endl;
		return;
	}

	const size_t m2 = m_truncModeDimensions[1];
	const size_t m3 = m_truncModeDimensions[2];
	if(weightVector.size() != m2+m3)
	{
		std::cout << "Wrong input weight vector dimension " << weightVector.size() << " != " << m2+m3 << std::endl;
		return;
	}

	std::vector<double> w2;
	w2.reserve(m_truncModeDimensions[1]);

	for(size_t i = 0; i < m_truncModeDimensions[1]; ++i)
	{
		w2.push_back(weightVector[i]);
	}

	Tensor t2;
	m_multilinearModel.modeMultiply(w2, "T", m_truncModeDimensions[1], 1, 2, t2);

	std::vector<double> w3;
	w3.reserve(m_truncModeDimensions[2]);

	for(size_t i = 0; i < m_truncModeDimensions[2]; ++i)
	{
		const size_t index = m_truncModeDimensions[1]+i;
		w3.push_back(weightVector[index]);
	}

	Tensor result;
	t2.modeMultiply(w3, "T", m_truncModeDimensions[2], 1, 3, result);

	outPoints.clear();
	outPoints.reserve(m_truncModeDimensions[0]);

	for(size_t i = 0; i < m_truncModeDimensions[0]; ++i)
	{
		const double value = result.getElement(i, 0, 0, 0);
		outPoints.push_back(value);
	}

	for(size_t i = 0; i < m_mean.size(); ++i)
	{
		outPoints[i] += m_mean[i];
	}
}

void MultilinearModelHandler::getDataMeanVertices(std::vector<double>& dataMean)
{
	reconstructForWeights(m_meanWeights, dataMean);
}

void MultilinearModelHandler::getNeutralDataMeanVertices(std::vector<double>& neutralDataMean)
{
	reconstructForWeights(m_neutralMeanWeights, neutralDataMean);
}

void MultilinearModelHandler::projectRegisteredToMM(const std::vector<double>& correspondingTargetPoints, const std::vector<bool>& validPoints, const std::vector<size_t>& modelLmkIndices, const std::vector<double>& targetLmks, const double landmarkWeight
																	, const std::vector<double>& initWeights, std::vector<double>& outWeights)
{
	//Center data
	const size_t d1(m_modeDimensions[0]);
	if(correspondingTargetPoints.size() != d1)
	{
		std::cout << "Wrong number of target points " << correspondingTargetPoints.size() << " != " << d1 << std::endl;
		return;
	}

	if(3*validPoints.size() != d1)
	{
		std::cout << "Number of valid point flags wrong " << validPoints.size() << " != " << d1/3 << std::endl;
		return;
	}

	std::vector<double> centeredTargetPoints;
	centeredTargetPoints.resize(d1);

	for(size_t i = 0; i < d1; ++i)
	{
		centeredTargetPoints[i] = correspondingTargetPoints[i] - m_mean[i];
	}

	//Center landmarks
	const size_t numLmks = std::min<size_t>(modelLmkIndices.size(), targetLmks.size()/3);
	
	std::vector<double> centeredTargetLmks;
	centeredTargetLmks.resize(3*numLmks);

	for(size_t i = 0; i < numLmks; ++i)
	{
		const size_t currLmkIndex = modelLmkIndices[i];
		for(size_t j = 0; j < 3; ++j)
		{
			centeredTargetLmks[3*i+j] = targetLmks[3*i+j]-m_mean[3*currLmkIndex+j];
		}
	}

	projectRegisteredToMultilinearModelNonLinear(centeredTargetPoints, validPoints, modelLmkIndices, centeredTargetLmks, landmarkWeight, initWeights, PROJECTION_MODE2_PRIOR_BOX_SIZE, PROJECTION_MODE3_PRIOR_BOX_SIZE, outWeights);
}

void MultilinearModelHandler::projectToMM(const DataContainer& modelMesh, const std::vector<size_t>& modelLmkIndices, const std::vector<double>& targetPoints, const std::vector<double>& targetNormals, const std::vector<double>& targetLmks
														, std::vector<double>& outWeights)
{
	KDTree3 targetKDTree(targetPoints);
	DataContainer currModel = modelMesh;

	outWeights = m_meanWeights;

	double landmarkWeight = PROJECTION_LMK_WEIGHT;

	//Iterative computation of a correspondence using nearest neighbors and projection to the multilinear model space
	for(size_t iter = 0; iter < PROJECTION_NUMBER_ITERATIONS; ++iter)
	{
		std::cout << "Iteration " << iter+1 << std::endl;

#ifdef USE_LANDMARKS_JUST_FOR_INITIALIZATION
		if(iter > 0)
		{
			landmarkWeight = 0.0;
		}
#endif

		//Reconstruct the model of the current iteration
		std::vector<double> currModelVertices;
		reconstructForWeights(outWeights, currModelVertices);

		currModel.setVertexList(currModelVertices);

		const std::vector<double>& modelMeshVertices = currModel.getVertexList();

		std::vector<double> modelMeshNormals;
		MathHelper::calcVertexNormals(currModel, modelMeshNormals);

		//Compute the nearest neighbors of the current model and compute its corresponding points at the target data
		std::vector<double> currNearestNeighbors;
		std::vector<bool> currValidValues;
		getNearestNeighbors(modelMeshVertices, modelMeshNormals, targetKDTree, targetPoints, targetNormals, PROJECTION_MAX_POINT_DISTANCE, PROJECTION_MAX_NORMAL_ANGLE, currNearestNeighbors, currValidValues);

		//Project the corresponding vertices to the multilinear model space
		std::vector<double> currOutWeights;
		projectRegisteredToMM(currNearestNeighbors, currValidValues, modelLmkIndices, targetLmks, landmarkWeight, outWeights, currOutWeights);

		outWeights = currOutWeights;
	}
}

bool MultilinearModelHandler::importMultilinearModel(const std::string& sstrFileName)
{
	clear();

	std::vector<double> multModel;

	FileLoader loader;
	if(!loader.loadRestrictedMultilinearModel(sstrFileName, m_modeDimensions, m_truncModeDimensions, multModel, m_meanWeights, m_neutralMeanWeights, m_mean))
	{
		return false;
	}

	const size_t d1 = m_modeDimensions[0];
	const size_t d2 = m_modeDimensions[1];
	const size_t d3 = m_modeDimensions[2];

	const size_t m1 = m_truncModeDimensions[0];
	const size_t m2 = m_truncModeDimensions[1];
	const size_t m3 = m_truncModeDimensions[2];

	m_multilinearModel.init(multModel, d1, m2, m3);

	std::cout << std::endl;
	std::cout << "*************************************************" << std::endl;
	std::cout << "Data dimension" << std::endl;
	std::cout << "Mode1: " << d1 << std::endl;
	std::cout << "Mode2: " << d2 << std::endl;
	std::cout << "Mode3: " << d3 << std::endl;
	std::cout << std::endl;
	std::cout << "Multilinear Model dimension (truncated)" << std::endl;
	std::cout << "Mode1: " << m1 << std::endl;
	std::cout << "Mode2: " << m2 << std::endl;
	std::cout << "Mode3: " << m3 << std::endl;
	std::cout << "*************************************************" << std::endl;
	std::cout << std::endl;

	return true;
}

void MultilinearModelHandler::projectRegisteredToMultilinearModelNonLinear(const std::vector<double>& centeredTargetPoints, const std::vector<bool>& validPoints, const std::vector<size_t>& modelLmkIndices, const std::vector<double>& centeredTargetLmks
																									, const double landmarkWeight, const std::vector<double>& initWeight, const double maxMode2Variation, const double maMode3Variaton, std::vector<double>& outWeight)
{
	const size_t d1(m_modeDimensions[0]);
	const size_t m2(m_truncModeDimensions[1]);
	const size_t m3(m_truncModeDimensions[2]);

	if(d1 != centeredTargetPoints.size())
	{
		std::cout << "MultilinearModelHandler::projectSequenceToMultilinearModelNonLinear(...) - Dimension of target points wrong " <<  d1 << " != " << centeredTargetPoints.size() << std::endl;
		return;
	}

	if(3*validPoints.size() != centeredTargetPoints.size())
	{
		std::cout << "MultilinearModelHandler::projectSequenceToMultilinearModelNonLinear(...) - Number of valid point flags wrong " <<  validPoints.size() << " != " << centeredTargetPoints.size()/3 << std::endl;
		return;
	}

	if(3*modelLmkIndices.size() != centeredTargetLmks.size())
	{
		std::cout << "MultilinearModelHandler::projectSequenceToMultilinearModelNonLinear(...) - Number of landmarks for fitting not equal " <<  modelLmkIndices.size() << " != " << centeredTargetLmks.size()/3 << std::endl;
		return;
	}

	const size_t numWeights = m2+m3;
	if(initWeight.size() != numWeights)
	{
		std::cout << "MultilinearModelHandler::projectSequenceToMultilinearModelNonLinear(...) - Initial weights of wrong size " <<  initWeight.size() << " != " << m2+m3 << std::endl;
		return;
	}

	vnl_vector<long> boundSelection(numWeights, 2);

	vnl_vector<double> lowerBounds(numWeights, 0.0);
	vnl_vector<double> upperBounds(numWeights, 0.0);

	for(size_t i = 0; i < m2; ++i)
	{
		const double lowerBound = getWeightForVariation(-maxMode2Variation, m_meanWeights[i], 2);
		const double upperBound = getWeightForVariation(maxMode2Variation, m_meanWeights[i], 2);

		lowerBounds[i] = lowerBound;
		upperBounds[i] = upperBound;
	}

	for(size_t i = 0; i < m3; ++i)
	{
		const size_t sourceIndex = m2+i;

		const double lowerBound = getWeightForVariation(-maMode3Variaton, m_meanWeights[sourceIndex], 3);
		const double upperBound = getWeightForVariation(maMode3Variaton, m_meanWeights[sourceIndex], 3);

		lowerBounds[m2+i] = lowerBound;
		upperBounds[m2+i] = upperBound;
	}

	vnl_vector<double> x(numWeights, 0.0);
	for(size_t i = 0; i < numWeights; ++i)
	{
		x[i] = initWeight[i];
	}

	MMProjectionCostFunction fkt(&m_multilinearModel, centeredTargetPoints, validPoints, modelLmkIndices, centeredTargetLmks, landmarkWeight);
	vnl_lbfgsb minimizer(fkt);
	minimizer.set_cost_function_convergence_factor(1000); //for higher precision: 10000000
	minimizer.set_projected_gradient_tolerance(0.001);		//for higher precision: 0.00001
	minimizer.set_max_function_evals(50);						//for higher presicion: 1000

	minimizer.set_bound_selection(boundSelection);
	minimizer.set_lower_bound(lowerBounds);
	minimizer.set_upper_bound(upperBounds);
#ifdef TRACE_NONLINEAR_PROJECTION_METHOD
	minimizer.set_trace(true);
#endif
	minimizer.minimize(x);

	outWeight.clear();
	outWeight.resize(numWeights);

	for(size_t i = 0; i < numWeights; ++i)
	{
		outWeight[i] = x[i];
	}
}

void MultilinearModelHandler:: getNearestNeighbors(const std::vector<double>& sourcePoints, const std::vector<double>& sourceNormals, KDTree3& targetKDTree, const std::vector<double>& targetPoints, const std::vector<double>& targetNormals
									, const double maxProjectionDist, const double maxAngle, std::vector<double>& nearestNeighbors, std::vector<bool>& validValues)
{
	const size_t numSourceVertices = sourcePoints.size()/3;

	//If normals are not valid, they are ignored
	bool bHasValidNormals = (sourcePoints.size() == sourceNormals.size()) && (targetPoints.size() == targetNormals.size());

	int validPoints(0);
	int invalidPoints(0);

	for(size_t i = 0; i < numSourceVertices; ++i)
	{
		const Vec3d sourcePoint(sourcePoints[3*i], sourcePoints[3*i+1], sourcePoints[3*i+2]);

		std::vector<double> querySourcePoint;
		querySourcePoint.push_back(sourcePoint[0]);
		querySourcePoint.push_back(sourcePoint[1]);
		querySourcePoint.push_back(sourcePoint[2]);

		int nnPointIndex(0);
		double nnSqrPointDist(0);
		targetKDTree.getNearestPoint(querySourcePoint, nnPointIndex, nnSqrPointDist);

		bool bPointValid = (sqrt(nnSqrPointDist) <= maxProjectionDist);
		if(bHasValidNormals)
		{
			const Vec3d sourceNormal(sourceNormals[3*i], sourceNormals[3*i+1], sourceNormals[3*i+2]);
			const Vec3d targetNormal(targetNormals[3*nnPointIndex], targetNormals[3*nnPointIndex+1], targetNormals[3*nnPointIndex+2]);
			const double angle = sourceNormal.angle(targetNormal);

			bPointValid &= (angle <= maxAngle);
		}

		validValues.push_back(bPointValid);

		if(bPointValid)
		{
			Vec3d nnPoint(targetPoints[3*nnPointIndex], targetPoints[3*nnPointIndex+1], targetPoints[3*nnPointIndex+2]);

			if(bHasValidNormals)
			{
				const Vec3d targetNormal(targetNormals[3*nnPointIndex], targetNormals[3*nnPointIndex+1], targetNormals[3*nnPointIndex+2]);

				Vec3d planeProjectionPoint;					
				MathHelper::getPlaneProjection(sourcePoint, nnPoint, targetNormal, planeProjectionPoint);

				nearestNeighbors.push_back(planeProjectionPoint[0]);
				nearestNeighbors.push_back(planeProjectionPoint[1]);
				nearestNeighbors.push_back(planeProjectionPoint[2]);
			}
			else
			{
				nearestNeighbors.push_back(nnPoint[0]);
				nearestNeighbors.push_back(nnPoint[1]);
				nearestNeighbors.push_back(nnPoint[2]);
			}

			++validPoints;
		}
		else
		{
			nearestNeighbors.push_back(0.0);
			nearestNeighbors.push_back(0.0);
			nearestNeighbors.push_back(0.0);

			++invalidPoints;
		}
	}
	
	std::cout << "getNearestNeighbors(...) - Valid points: " << validPoints << " Invalid points: " << invalidPoints << std::endl;
}