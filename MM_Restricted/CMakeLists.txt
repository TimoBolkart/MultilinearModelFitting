SET( MM_Files
	FileHandler.cpp
	FileLoader.cpp
	FileWriter.cpp
	KDTree3.cpp
	MathHelper.cpp
	MMProjectionCostFunction.cpp
	MultilinearModel.cpp
	MultilinearModelHandler.cpp
	Main.cpp
)

INCLUDE_DIRECTORIES(${ANN_INCLUDE_DIR})
INCLUDE_DIRECTORIES(${CLAPACK_INCLUDE_DIR}) 
INCLUDE_DIRECTORIES(${ITK_INCLUDES}) 
   
ADD_EXECUTABLE(MM_Restricted ${MM_Files})
TARGET_LINK_LIBRARIES(MM_Restricted ${ITK_LIBRARIES} ${ANN_LIBRARIES} ${CLAPACK_LIBRARIES})   