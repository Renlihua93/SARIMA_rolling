# SPDX-License-Identifier: BSD-3-Clause
add_executable(sarimatest sarimatest.c
		../src/mean_fit.c
		../src/mean_fit.h)

target_link_libraries(sarimatest PRIVATE sarimalib)
add_dependencies(sarimatest sarimalib)

if(NOT MSVC)
	target_link_libraries(sarimatest PRIVATE m)
endif()

target_include_directories(sarimatest PRIVATE ${CMAKE_SOURCE_DIR}/include)


target_compile_definitions(sarimatest PRIVATE
	CTSA_BUILT_WITH_CMAKE=1
	CTSA_DATA_FILE="${CTSA_DATA_FILE}"
)

set_target_properties(sarimatest PROPERTIES RUNTIME_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}")

file(COPY ${CMAKE_SOURCE_DIR}/data/
		DESTINATION ${CMAKE_BINARY_DIR}/data)
