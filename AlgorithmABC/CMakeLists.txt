cmake_minimum_required(VERSION 2.8)


project(frontiers23)

set(CMAKE_BUILD_TYPE Release)
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -lpthread -Ofast")


set(CMAKE_CXX_STANDARD 17)


add_executable(main
    main.cpp
    ABC.cpp
    model.cpp)

include_directories(main PUBLIC
    ./
    model/
    vendor/
    )
