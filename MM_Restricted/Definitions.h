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

#ifndef DEFINITIONS_H
#define DEFINITIONS_H

/*
Parameter for the refinement of the initial landmark-based alignment using a rigid ICP
*/

//Number of ICP alignment iterations.
const int ICP_NUMBER_ALIGNMENT_ITERATIONS = 1;
//Maximum point distance for ICP (unit of the learned model).
const double ICP_MAX_POINT_DISTANCE = 10.0;
//Maximum angle between normals of points used for ICP.
const double ICP_MAX_NORMAL_ANGLE = 45.0;


/*
Parameter for the fitting
*/

//Weighting of the landmarks for fitting. (0: not used, <1: lower influence than vertex data, 1: same influence than vertex data, >1: higher influence than vertex data).
const double PROJECTION_LMK_WEIGHT = 2.0;
//If enabled, landmarks are just used for initialization while the first fitting iteration. Especially recommended if landmarks are not so reliable.
#define USE_LANDMARKS_JUST_FOR_INITIALIZATION

//Size of the prior-box of the second mode (distance from the mean within the second mode that are allowed while fitting).
const double PROJECTION_MODE2_PRIOR_BOX_SIZE = 1.0;
//Size of the prior-box of the third mode (distance from the mean within the third mode that are allowed while fitting).
const double PROJECTION_MODE3_PRIOR_BOX_SIZE = 1.0;

//Number of iteration steps while fitting.
const size_t PROJECTION_NUMBER_ITERATIONS = 6;

//Maximum point distance for fitting (unit of the learned model).
const double PROJECTION_MAX_POINT_DISTANCE = 10.0;
//Maximum angle between normals of points used for fitting.
const double PROJECTION_MAX_NORMAL_ANGLE = 45.0;

//If enabled, non-linear optimization outputs function and gradient in each step.
//#define TRACE_NONLINEAR_PROJECTION_METHOD

#endif