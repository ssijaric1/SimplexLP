set(LPBENCH_NAME lpBench)				#console benchmark target

file(GLOB LPBENCH_SOURCES   ${CMAKE_CURRENT_LIST_DIR}/cli/*.cpp)
file(GLOB LPBENCH_INC_CORE  ${CMAKE_CURRENT_LIST_DIR}/src/core/*.h)
file(GLOB LPBENCH_INC_SPARSE ${NATID_SDK_INC}/sparse/*.h)

# add executable
add_executable(${LPBENCH_NAME} ${LPBENCH_SOURCES} ${LPBENCH_INC_CORE}
				${LPBENCH_INC_SPARSE})

source_group("src"          FILES ${LPBENCH_SOURCES})
source_group("src\\core"    FILES ${LPBENCH_INC_CORE})
source_group("inc\\sparse"  FILES ${LPBENCH_INC_SPARSE})

target_link_libraries(${LPBENCH_NAME}
		debug ${MU_LIB_DEBUG}     optimized ${MU_LIB_RELEASE}
		debug ${MATRIX_LIB_DEBUG} optimized ${MATRIX_LIB_RELEASE})

setIDEPropertiesForExecutable(${LPBENCH_NAME})
