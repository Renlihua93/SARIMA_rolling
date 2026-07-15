# SPDX-License-Identifier: BSD-3-Clause
cmake_minimum_required(VERSION 3.16)

project(ctsa_sarima_simple C)

set(CMAKE_C_STANDARD 99)
set(CMAKE_C_STANDARD_REQUIRED ON)

if(MSVC)
	add_compile_options(/W3)
else()
	add_compile_options(-Wall -Wextra)
endif()

add_subdirectory(src)
add_subdirectory(test)
