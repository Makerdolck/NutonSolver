file(GLOB headers
        "includes/*.h"
        "includes/*.cpp")

add_executable(Solver Solver.cpp ${headers})
# target_include_directories(Solver PUBLIC ./includes)