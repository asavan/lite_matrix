cmake_minimum_required(VERSION 3.16)
get_filename_component(ProjectId ${CMAKE_CURRENT_SOURCE_DIR} NAME)
string(REPLACE " " "_" ProjectId ${ProjectId})
project(${ProjectId})
set(CMAKE_CXX_STANDARD 20)
set(CMAKE_DISABLE_PRECOMPILE_HEADERS ON)

if (MSVC)
    add_compile_options(/W4 /WX)
else()
    add_compile_options(-Wall -Wextra -pedantic -Werror)
endif()

add_library(MatrixLib reverse_matrix.cpp)
target_sources(MatrixLib
    PRIVATE
        matrix.hpp
    PUBLIC
        matrix.h
		reverse.h
)

add_executable(${ProjectId} main.cpp)
target_link_libraries(${ProjectId} PRIVATE MatrixLib)
# target_link_libraries(tictactoe user32 gdi32 kernel32)
set_property(DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR} PROPERTY VS_STARTUP_PROJECT ${ProjectId})

file(COPY matrix.txt DESTINATION "${CMAKE_CURRENT_BINARY_DIR}")
