target_sources(${PROJECT_NAME} PRIVATE 
	${CMAKE_CURRENT_SOURCE_DIR}/MicrostructureSpecificationBase.cpp
	${CMAKE_CURRENT_SOURCE_DIR}/MicrostructureGenerator.cpp
	${CMAKE_CURRENT_SOURCE_DIR}/ShearLoopGenerator.cpp
	${CMAKE_CURRENT_SOURCE_DIR}/ShearLoopDensitySpecification.cpp
	${CMAKE_CURRENT_SOURCE_DIR}/ShearLoopIndividualSpecification.cpp
	${CMAKE_CURRENT_SOURCE_DIR}/PeriodicDipoleDensitySpecification.cpp
	${CMAKE_CURRENT_SOURCE_DIR}/PeriodicDipoleIndividualSpecification.cpp
	${CMAKE_CURRENT_SOURCE_DIR}/PeriodicDipoleGenerator.cpp
	${CMAKE_CURRENT_SOURCE_DIR}/PrismaticLoopDensitySpecification.cpp
	${CMAKE_CURRENT_SOURCE_DIR}/PrismaticLoopIndividualSpecification.cpp
	${CMAKE_CURRENT_SOURCE_DIR}/PrismaticLoopGenerator.cpp
	${CMAKE_CURRENT_SOURCE_DIR}/FrankLoopsDensitySpecification.cpp
	${CMAKE_CURRENT_SOURCE_DIR}/FrankLoopsIndividualSpecification.cpp
	${CMAKE_CURRENT_SOURCE_DIR}/FrankLoopsGenerator.cpp
	${CMAKE_CURRENT_SOURCE_DIR}/StackingFaultTetrahedraDensitySpecification.cpp
	${CMAKE_CURRENT_SOURCE_DIR}/StackingFaultTetrahedraIndividualSpecification.cpp
	${CMAKE_CURRENT_SOURCE_DIR}/StackingFaultTetrahedraGenerator.cpp
	${CMAKE_CURRENT_SOURCE_DIR}/PlanarLoopIndividualSpecification.cpp
	${CMAKE_CURRENT_SOURCE_DIR}/PlanarLoopGenerator.cpp
	${CMAKE_CURRENT_SOURCE_DIR}/SphericalInclusionDensitySpecification.cpp
	${CMAKE_CURRENT_SOURCE_DIR}/SphericalInclusionIndividualSpecification.cpp
	${CMAKE_CURRENT_SOURCE_DIR}/SphericalInclusionsGenerator.cpp
	${CMAKE_CURRENT_SOURCE_DIR}/PolyhedronInclusionIndividualSpecification.cpp
	${CMAKE_CURRENT_SOURCE_DIR}/PolyhedronInclusionsGenerator.cpp
#	${CMAKE_CURRENT_SOURCE_DIR}/VTKGenerator.cpp
			)


