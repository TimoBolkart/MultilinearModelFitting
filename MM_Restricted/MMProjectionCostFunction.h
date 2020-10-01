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

#ifndef MMPROJECTIONCOSTFUNCTION_H
#define MMPROJECTIONCOSTFUNCTION_H

#include <vnl/vnl_vector.h>
#include <vnl/vnl_cost_function.h>

#include <vector>

class Tensor;

class MMProjectionCostFunction : public vnl_cost_function
{
public:
	//! Specific projection cost function.
	//! \param pMM								pointer to the multilinear model used for error minimization
	//! \param centeredTargetPoints		vertices of the centered target object
	//! \param validPoints					flags if vertices are valid or not
	//! \param modelLmkIndices				indices of the landmarks on the template
	//! \param centeredTargetLmks			landmark vertices for the target object
	//! \param landmarkWeight				weight of the landmark energy term
	MMProjectionCostFunction(const Tensor* pMM, const std::vector<double>& centeredTargetPoints, const std::vector<bool>& validPoints, const std::vector<size_t>& modelLmkIndices, const std::vector<double>& centeredTargetLmks, const double landmarkWeight);

	~MMProjectionCostFunction();

	//! Derived function to compute function value and gradient.
	//! \param x		initial values for current iteration
	//! \param f		function value evaluated on x
	//! \param g		gradient on x
	virtual void compute(const vnl_vector<double>& x, double *f, vnl_vector<double>* g);

private:
	void init(const Tensor* pMM, const std::vector<double>& centeredTargetPoints, const std::vector<bool>& validPoints, const std::vector<size_t>& modelLmkIndices, const std::vector<double>& centeredTargetLmks, const double landmarkWeight);

	//! Compute the value of the data energy term and its gradient.
	//! \param M3				matrix obtained by mode multiplication of the multilinear model tensor with the coefficients of the second mode
	//! \param M2				matrix obtained by mode multiplication of the multilinear model tensor with the coefficients of the third mode
	//! \param currData		coordinates of the data for values of the last iteration
	//! \param E_data			output data energy value
	//! \param gradE_data	output data energy gradient
	void calcDataEnergy(const Tensor& M3, const Tensor& M2, const Tensor& currData, double& E_data, std::vector<double>& gradE_data);

	//! Compute the value of the landmark energy term and its gradient.
	//! \param M3				matrix obtained by mode multiplication of the multilinear model tensor with the coefficients of the second mode
	//! \param M2				matrix obtained by mode multiplication of the multilinear model tensor with the coefficients of the third mode
	//! \param currData		coordinates of the data for values of the last iteration	
	//! \param E_data			output landmark energy value
	//! \param gradE_data	output data energy gradient
	void calcLmkEnergy(const Tensor& M3, const Tensor& M2, const Tensor& currData, double& E_lmk, std::vector<double>& gradE_lmk);

	//! Checks if the input data are of the right dimensions.
	//! \param pMM								pointer to the multilinear model used for error minimization
	//! \param centeredTargetPoints		vertices of the centered target object
	//! \param validPoints					flags if vertices are valid or not
	//! \param modelLmkIndices				indices of the landmarks on the template
	//! \param centeredTargetLmks			landmark vertices for the target object
	//! \return									true if the data are of right dimensions
	bool checkInputDataDimensions(const Tensor* pMM, const std::vector<double>& centeredTargetPoints, const std::vector<bool>& validPoints, const std::vector<size_t>& modelLmkIndices, const std::vector<double>& centeredTargetLmks);

	template<typename T>
	void initVector(std::vector<T>& vec, const size_t dim, const T& value)
	{
		vec.clear();
		vec.resize(dim);

		for(size_t i = 0; i < dim; ++i)
		{
			vec[i] = value;
		}
	}

	const Tensor* m_pMultilinearModel;
	std::vector<double> m_centeredTargetPoints;					
	std::vector<double> m_pointWeights;
	std::vector<size_t> m_lmkIndices;					
	std::vector<double> m_centeredTargetLmks;		
	double m_landmarkWeight;
};

#endif