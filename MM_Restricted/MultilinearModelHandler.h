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

#ifndef MULTILINEARMODELHANLDER
#define MULTILINEARMODELHANLDER

#include "Definitions.h"
#include "DataContainer.h"
#include "MultilinearModel.h"

#include <vector>
#include <map>

class KDTree3;

class MultilinearModelHandler
{
public:
	//! Handler for all multilinear model related stuff computation stuff. 
	MultilinearModelHandler();

	~MultilinearModelHandler();

	//! Clear all data structures.
	void clear();

	//! Compute the vertices for variations. The variation parameters are centered at 0 and and scaled. 
	//! E.g. To get the variations where some principal components are r*stdev far from the mean, call this function with a vector, where corresponding entries are r.
	//! \param variations		(m2+m3)-dimensional input vector where the first m2 entries represent the variations of the second mode, and the following m3 entries the variation of third mode
	//! \param outPoints			vertices of the reconstruction
	void reconstructForVariations(const std::vector<double>& variations, std::vector<double>& outPoints);

	//! Compute the vertices for mode coefficients.
	//! \param weightVector		(m2+m3)-dimensional input vector where the first m2 entries represent the coefficients of the second mode, and the following m3 entries the cofficients of third mode
	//! \param outPoints			vertices of the reconstruction
	void reconstructForWeights(const std::vector<double>& weightVector, std::vector<double>& outPoints);

	//! Get the vertices of the learned data mean.
	//! \param dataMean			learned data mean
	void getDataMeanVertices(std::vector<double>& dataMean);

	//! Get the vertices of the learned neutral data mean.
	//! \param neutralDataMean	learned neutral data mean
	void getNeutralDataMeanVertices(std::vector<double>& neutralDataMean);

	//! Project data to the multilinear model space where all data vertices are in full correspondence to the template data.
	//! \param correspondingTargetPoints	vertices that are in full correspondence to the template vertices
	//! \param validPoints						flags that indicate which vertices are valid
	//! \param modelLmkIndices					landmarks of the template model used for fitting
	//! \param targetLmks						landmarks of the target model used for fitting
	//! \param landmarkWeight					weight of landmark energy (0: not used, <1: lower influence than vertex data, 1: same influence than vertex data, >1: higher influence than vertex data)
	//! \param initWeights						(m2+m3)-dimensional vector with initial coefficients in multilinear model space
	//! \param outWeights						coefficients of the projection 
	void projectRegisteredToMM(const std::vector<double>& correspondingTargetPoints, const std::vector<bool>& validPoints, const std::vector<size_t>& modelLmkIndices, const std::vector<double>& targetLmks, const double landmarkWeight
										, const std::vector<double>& initWeights, std::vector<double>& outWeights);

	//! Project data to the multilinear model space.
	//! \param modelMesh					template mesh of the model
	//! \param modelLmkIndices			indices of the template data landmarks used for fitting
	//! \param targetPoints				vertices of the data that projected to the multilinear model space
	//! \param targetNormals			vertex normals of the projected data (if empty, no normals are used)
	//! \param targetLmks				corresponding landmarks of the projected data
	//! \param outWeights				coefficients fo the projection
	void projectToMM(const DataContainer& modelMesh, const std::vector<size_t>& modelLmkIndices, const std::vector<double>& targetPoints, const std::vector<double>& targetNormals, const std::vector<double>& targetLmks, std::vector<double>& outWeights);

	//! Import a multilinear model.
	//! \param sstrFileName			full file name of the model
	//! \return true if import was successful
	bool importMultilinearModel(const std::string& sstrFileName);

private:

	//! Convert from centered and scaled variation parameter to model coefficient.
	//! \param variation			parameter in the centered and scaled variation space
	//! \param meanWeight		mean weight of the model of corresponding coefficient
	//! \param mode				mode of the model
	//! \return	the converted model coefficient
	double getWeightForVariation(const double variation, const double meanWeight, const size_t mode)
	{
		if((mode<1) || (m_modeDimensions.size() < mode-1))
		{
			return 0.0;
		}

		const double tmp = sqrt(static_cast<double>(m_modeDimensions[mode-1]));
		return variation/tmp+meanWeight;
	}

	//! Convert from the model coefficient to the centered and scaled variation parameter.
	//! \param weight				model coefficient
	//! \param meanWeight		mean weight of the model of corresponding coefficient
	//! \param mode				mode of the model
	//! \return	the converted parameter in the centered and scaled variation space
	double getVariationForWeight(const double weight, const double meanWeight, const size_t mode)
	{
		if((mode<1) || (m_modeDimensions.size() < mode-1))
		{
			return 0.0;
		}

		const double tmp = sqrt(static_cast<double>(m_modeDimensions[mode-1]));
		return (weight-meanWeight)*tmp;
	}

	//! Helper function to project data to the multilinear model space.
	//! \param centeredTargetPoints		vertices of the projection target that are in full correspondence with the template data and with subtracted learned model mean
	//! \param validPoints					flags that indicate which vertices are valid
	//! \param modelLmkIndices				landmarks of the template model used for fitting
	//! \param centeredTargetLmks			landmarks of the target model used for fitting with subtracted learned model mean
	//! \param landmarkWeight				weight of landmark energy (0: not used, <1: lower influence than vertex data, 1: same influence than vertex data, >1: higher influence than vertex data)
	//! \param initWeight					(m2+m3)-dimensional vector with initial coefficients in multilinear model space
	//! \param maxMode2Variation			size of the prior-box of the second mode (maximum distance from the model mode mean of the second mode)
	//! \param maxMode3Variaton			size of the prior-box of the third mode (maximum distance from the model mode mean of the third mode)
	//! \param outWeight						coefficients fo the projection
	void projectRegisteredToMultilinearModelNonLinear(const std::vector<double>& centeredTargetPoints, const std::vector<bool>& validPoints, const std::vector<size_t>& modelLmkIndices, const std::vector<double>& centeredTargetLmks
																	, const double landmarkWeight, const std::vector<double>& initWeight, const double maxMode2Variation, const double maxMode3Variaton, std::vector<double>& outWeight);
	
	//! Computes nearest neighbors for a set of vertices.
	//! \param sourcePoints					set of vertices the nearest neighbors are computed for
	//! \param sourceNormals				vertex normals of input vertices
	//! \param targetKDTree					KDTree of the target vertices
	//! \param targetPoints					target vertices
	//! \param targetNormals				target vertex normals (not used if empty)
	//! \param maxProjectionDist			maximum distance between a source point and its computed nearest neighbor
	//! \param maxAngle						maximum angle between the normals of a source point and its computed nearest neighbor (just used if normals are not empty)
	//! \param nearestNeighbors			output set of vertices that are corresponding to the source points w.r.t. the nearest neighbor measure
	//! \param validValues					flags that indicate which vertices are valid
	void getNearestNeighbors(const std::vector<double>& sourcePoints, const std::vector<double>& sourceNormals, KDTree3& targetKDTree, const std::vector<double>& targetPoints, const std::vector<double>& targetNormals
									, const double maxProjectionDist, const double maxAngle, std::vector<double>& nearestNeighbors, std::vector<bool>& validValues);

	//! Multilinear model tensor of dimension d1 x m2 x m3
	Tensor m_multilinearModel;
	//! Dimension of the original data d1 d2 d3
	std::vector<size_t> m_modeDimensions;
	//! Dimensions of the truncated model d1 m2 m3
	std::vector<size_t> m_truncModeDimensions;
	//! Learned data mean
	std::vector<double> m_mean;
	//! (m2+m3)-dimensional vector with mean coefficients of mode 2 and mode 3
	std::vector<double> m_meanWeights;
	//! (m2+m3)-dimensional vector with mean coefficients of mode 2 and netural mean coefficients of mode 3
	std::vector<double> m_neutralMeanWeights;
};
#endif