# frontiers23: Socioeconomic determinants of Stay-at-Home Policies during the First COVID-19 Wave

## Build instructions for the C++ code

The code works on Windows, MacOS and Linux.

### Using CMake

CMake should automatically detect the machine's compiler. Minimum version of CMake required: 2.8

On MacOS and Linux:

1. Create folder inside `/AlgorithmABC` named `build`
2. Enter folder `build`
3. Create makefiles: run `cmake ..`
4. Build the program: run `make`
5. Run the program: `./main`

On Windows:

1. Create folder inside `/AlgorithmABC` named `build`
2. Enter folder `build`
3. Configure CMake files: run `cmake ..`
4. Build the program: `cmake --build . --config release`
5. Run the program: `.\Release\main.exe`

### Using g++

On MacOS and Linux:

1. Create folder inside `/AlgorithmABC` named `build`
2. Enter folder `build`
3. Command: `g++ ABC.cpp main.cpp model.cpp -I . -o build/main -I vendor -lpthread -std=c++17 -Ofast`
4. Run the program: `./main`

On Windows:

1. Create folder inside `/AlgorithmABC` named `build`
2. Enter folder `build`
3. Command: `g++ ABC.cpp main.cpp model.cpp -I . -o build/main.exe -I vendor -lpthread -std=c++17 -Ofast`
4. Run the program: `.\main.exe`
