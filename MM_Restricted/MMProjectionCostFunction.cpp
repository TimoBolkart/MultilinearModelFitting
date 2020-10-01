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

#include "MMProjectionCostFunction.h"
#include "MultilinearModel.h"
#include "Definitions.h"

#include <iostream>

MMProjectionCostFunction::MMProjectionCostFunction(const Tensor* pMM, const std::vector<double>& centeredTargetPoints, const std::vector<bool>& validPoints, const std::vector<size_t>& modelLmkIndices
																	, const std::vector<double>& centeredTargetLmks, const double landmarkWeight)
: vnl_cost_function(pMM->getModeDimension(2)+pMM->getModeDimension(3))
{
	init(pMM, centeredTargetPoints, validPoints, modelLmkIndices, centeredTargetLmks, landmarkWeight);
}

MMProjectionCostFunction::~MMProjectionCostFunction()
{

}

void MMProjectionCostFunction::compute(const vnl_vector<double>& x, double *f, vnl_vector<double>* g)
{
	const size_t m2(m_pMultilinearModel->getModeDimension(2));
	const size_t m3(m_pMultilinearModel->getModeDimension(3));

	std::vector<double> w2;
	w2.resize(m2);

	for(size_t i = 0; i < m2; ++i)
	{
		w2[i] = x[i];
	}

	Tensor M3;
	m_pMultilinearModel->modeMultiply(w2, "T", m2, 1, 2, M3);

	std::vector<double> w3;
	w3.resize(m3);

	for(size_t i = 0; i < m3; ++i)
	{
		w3[i] = x[m2+i];
	}

	Tensor M2;
	m_pMultilinearModel->modeMultiply(w3, "T", m3, 1, 3, M2);

	Tensor currData;
	M3.modeMultiply(w3, "T", m3, 1, 3, currData);

	double E_data(0.0);
	std::vector<double> gradE_data;
	calcDataEnergy(M3, M2, currData, E_data, gradE_data);

	double E_lmk(0.0);
	std::vector<double> gradE_lmk;
	calcLmkEnergy(M3, M2, currData, E_lmk, gradE_lmk);

	*f = E_data+m_landmarkWeight*E_lmk;

	for(size_t i = 0; i < m2+m3; ++i)
	{
		(*g)[i] = gradE_data[i]+m_landmarkWeight*gradE_lmk[i];
	}
}

void MMProjectionCostFunction::init(const Tensor* pMM, const std::vector<double>& centeredTargetPoints, const std::vector<bool>& validPoints, const std::vector<size_t>& modelLmkIndices, const std::vector<double>& centeredTargetLmks
												, const double landmarkWeight)
{
	if(!checkInputDataDimensions(pMM, centeredTargetPoints, validPoints, modelLmkIndices, centeredTargetLmks))
	{
		std::cout << "MMProjectionCostFunction::init(...) - Unable to initialize" << std::endl;
		return;
	}

	m_pMultilinearModel = pMM;
	m_centeredTargetPoints = centeredTargetPoints;
	m_lmkIndices = modelLmkIndices;
	m_centeredTargetLmks = centeredTargetLmks;	
	m_landmarkWeight = landmarkWeight;

	const size_t numPoints = validPoints.size();
	m_pointWeights.resize(numPoints);

	for(size_t i = 0; i < numPoints; ++i)
	{
		m_pointWeights[i] = validPoints[i] ? 1.0 : 0.0;
	}
}

void MMProjectionCostFunction::calcDataEnergy(const Tensor& M3, const Tensor& M2, const Tensor& currData, double& E_data, std::vector<double>& gradE_data)
{
	const size_t d1(m_pMultilinearModel->getModeDimension(1));
	const size_t m2(m_pMultilinearModel->getModeDimension(2));
	const size_t m3(m_pMultilinearModel->getModeDimension(3));

	E_data = 0.0;
	initVector(gradE_data, m2+m3, 0.0);

	std::vector<double> tmpDiffVec;
	tmpDiffVec.resize(d1);

	const size_t numPoints = d1/3;

	double factor(0.0);

	for(size_t i = 0; i < numPoints; ++i)
	{
		const double currFactor = m_pointWeights[i];
		factor += currFactor;

		const size_t startIndex = 3*i;
		for(size_t j = 0; j < 3; ++j)
		{
			const double tmpDiff = currFactor > DBL_EPSILON ? currFactor*(currData.getElement(startIndex+j, 0, 0, 0) - m_centeredTargetPoints[startIndex+j]) : 0.0;
			E_data += pow(tmpDiff, 2);

			tmpDiffVec[startIndex+j] = tmpDiff;
		}
	}

	const double multFactor = factor > DBL_EPSILON ? 1.0/factor : 1.0;
	E_data *= multFactor;

	// calc mode 2 gradient
	for(size_t i = 0; i < m2; ++i)
	{
		double out = 0.0;
		for(size_t j = 0; j < d1; ++j)
		{
			out += M2.getElement(j, i, 0, 0)*tmpDiffVec[j]; 
		}

		gradE_data[i] += 2.0*multFactor*out;
	}

	// calc mode 3 gradient
	for(size_t i = 0; i < m3; ++i)
	{
		double out = 0.0;
		for(size_t j = 0; j < d1; ++j)
		{
			out += M3.getElement(j, 0, i, 0)*tmpDiffVec[j]; 
		}

		gradE_data[m2+i] += 2.0*multFactor*out;
	}
}

void MMProjectionCostFunction::calcLmkEnergy(const Tensor& M3, const Tensor& M2, const Tensor& currData, double& E_lmk, std::vector<double>& gradE_lmk)
{
	const size_t m2(m_pMultilinearModel->getModeDimension(2));
	const size_t m3(m_pMultilinearModel->getModeDimension(3));

	E_lmk = 0.0;
	initVector(gradE_lmk, m2+m3, 0.0);

	if(m_lmkIndices.empty() || m_centeredTargetLmks.empty())
	{
		return;
	}

	const size_t numLmks = m_lmkIndices.size();

	std::vector<double> tmpDiffVec;
	initVector(tmpDiffVec, 3*numLmks, 0.0);

	for(size_t i = 0; i < numLmks; ++i)
	{
		const size_t currLmkIndex = m_lmkIndices[i];
		
		for(size_t j = 0; j < 3; ++j)
		{
			const double tmpDiff = currData.getElement(3*currLmkIndex+j, 0, 0, 0) - m_centeredTargetLmks[3*i+j];
			E_lmk += pow(tmpDiff, 2);

			tmpDiffVec[3*i+j] = tmpDiff;
		}
	}

	const double multFactor = 1.0 / static_cast<double>(numLmks);
	E_lmk *= multFactor;

	// calc mode 2 gradient
	for(size_t i = 0; i < m2; ++i)
	{
		double out = 0.0;
		for(size_t j = 0; j < numLmks; ++j)
		{
			const size_t currLmkIndex = m_lmkIndices[j];

			for(size_t k = 0; k < 3; ++k)
			{
				out += M2.getElement(3*currLmkIndex+k, i, 0, 0)*tmpDiffVec[3*j+k]; 
			}
		}

		gradE_lmk[i] += 2.0*multFactor*out;
	}

	// calc mode 3 gradient
	for(size_t i = 0; i < m3; ++i)
	{
		double out = 0.0;
		for(size_t j = 0; j < numLmks; ++j)
		{
			const size_t currLmkIndex = m_lmkIndices[j];

			for(size_t k = 0; k < 3; ++k)
			{
				out += M3.getElement(3*currLmkIndex+k, 0, i, 0)*tmpDiffVec[3*j+k]; 
			}
		}

		gradE_lmk[m2+i] += 2.0*multFactor*out;
	}
}

bool MMProjectionCostFunction::checkInputDataDimensions(const Tensor* pMM, const std::vector<double>& centeredTargetPoints, const std::vector<bool>& validPoints
																			, const std::vector<size_t>& modelLmkIndices, const std::vector<double>& centeredTargetLmks)
{
	if(pMM == NULL)
	{
		std::cout << "Multilinear model instance not valid" << std::endl;
		return false;
	}

	const size_t d1(pMM->getModeDimension(1));
	if(centeredTargetPoints.size() != d1)
	{
		std::cout << "Target point dimension not valid " << centeredTargetPoints.size() << " != " << d1 << std::endl;
		return false;
	}

	if(3*validPoints.size() != d1)
	{
		std::cout << "Number of valid target point flags not valid " << validPoints.size() << " != " << d1/3;
		return false;
	}

	if(3*modelLmkIndices.size() != centeredTargetLmks.size())
	{
		std::cout << "Number of landmarks for fitting not valid " << modelLmkIndices.size() << " != " << centeredTargetLmks.size()/3;
		return false;
	}

	return true;
}