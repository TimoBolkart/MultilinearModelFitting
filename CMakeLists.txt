cmake_minimum_required(VERSION 2.7)

INCLUDE(ConfigureANN.cmake)
INCLUDE(ConfigureCLAPACK.cmake)
INCLUDE(ConfigureITK.cmake)

IF(ANN_FOUND AND CLAPACK_FOUND AND ITK_FOUND)
  ADD_SUBDIRECTORY(MM_Restricted)
ELSE(ANN_FOUND AND CLAPACK_FOUND AND ITK_FOUND)
  Message("PROBLEM: One of the required libraries not found. MM_Restricted will not be compiled.")  
ENDIF(ANN_FOUND AND CLAPACK_FOUND AND ITK_FOUND) 