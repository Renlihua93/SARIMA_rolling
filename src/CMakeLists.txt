# SPDX-License-Identifier: BSD-3-Clause
# SARIMA library (see README.md).
include_directories(${CMAKE_CURRENT_SOURCE_DIR})
include_directories(${CMAKE_SOURCE_DIR}/include)

set(SOURCE_FILES
	sarima.c
	fit.c
	ts.c
	optim.c
	linalg.c
	cstats.c
	load_matrix.c
	linear_fit.c
	mean_fit.c
	mean_fit.h
)

add_library(sarimalib STATIC ${SOURCE_FILES}
		mean_fit.c
		mean_fit.h)
target_include_directories(sarimalib PUBLIC
	${CMAKE_SOURCE_DIR}/include
	${CMAKE_CURRENT_SOURCE_DIR}
)

set_target_properties(sarimalib PROPERTIES
	ARCHIVE_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/lib"
)
