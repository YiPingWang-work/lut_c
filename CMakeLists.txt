cmake_minimum_required(VERSION 3.30)
project(lut_c C)

set(CMAKE_C_STANDARD 17)

add_executable(lut_c main.c
        lut/mul_mat_with_lut.c
)
