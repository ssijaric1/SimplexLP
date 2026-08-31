set(SIMPLEXLP_NAME simplexLP)			#GUI application target

#collect source files (src + subfolders)
file(GLOB SIMPLEXLP_SOURCES    ${CMAKE_CURRENT_LIST_DIR}/src/*.cpp)
file(GLOB SIMPLEXLP_INCS       ${CMAKE_CURRENT_LIST_DIR}/src/*.h)
file(GLOB SIMPLEXLP_INC_CORE   ${CMAKE_CURRENT_LIST_DIR}/src/core/*.h)
file(GLOB SIMPLEXLP_INC_VIZ    ${CMAKE_CURRENT_LIST_DIR}/src/viz/*.h)
file(GLOB SIMPLEXLP_INC_BENCH  ${CMAKE_CURRENT_LIST_DIR}/src/bench/*.h)
set(SIMPLEXLP_PLIST            ${CMAKE_CURRENT_LIST_DIR}/src/Info.plist)

file(GLOB SIMPLEXLP_INC_SPARSE ${NATID_SDK_INC}/sparse/*.h)
file(GLOB SIMPLEXLP_INC_DENSE  ${NATID_SDK_INC}/dense/*.h)

# add executable
add_executable(${SIMPLEXLP_NAME} ${SIMPLEXLP_INCS} ${SIMPLEXLP_SOURCES}
				${SIMPLEXLP_INC_CORE} ${SIMPLEXLP_INC_VIZ} ${SIMPLEXLP_INC_BENCH}
				${SIMPLEXLP_INC_SPARSE} ${SIMPLEXLP_INC_DENSE})

source_group("src"            FILES ${SIMPLEXLP_SOURCES} ${SIMPLEXLP_INCS})
source_group("src\\core"      FILES ${SIMPLEXLP_INC_CORE})
source_group("src\\viz"       FILES ${SIMPLEXLP_INC_VIZ})
source_group("src\\bench"     FILES ${SIMPLEXLP_INC_BENCH})
source_group("inc\\sparse"    FILES ${SIMPLEXLP_INC_SPARSE})
source_group("inc\\dense"     FILES ${SIMPLEXLP_INC_DENSE})

target_link_libraries(${SIMPLEXLP_NAME}
		debug ${MU_LIB_DEBUG}      optimized ${MU_LIB_RELEASE}
		debug ${NATGUI_LIB_DEBUG}  optimized ${NATGUI_LIB_RELEASE}
		debug ${MATRIX_LIB_DEBUG}  optimized ${MATRIX_LIB_RELEASE})

setTargetPropertiesForGUIApp(${SIMPLEXLP_NAME} ${SIMPLEXLP_PLIST})

setIDEPropertiesForGUIExecutable(${SIMPLEXLP_NAME} ${CMAKE_CURRENT_LIST_DIR})

setPlatformDLLPath(${SIMPLEXLP_NAME})

setAppIcon(${SIMPLEXLP_NAME} ${CMAKE_CURRENT_LIST_DIR})
