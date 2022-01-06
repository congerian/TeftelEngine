mkdir .\build\
cd ./build && ^
cmake ../ -G Ninja -DCMAKE_CUDA_COMPILER="C:/Program Files/NVIDIA GPU Computing Toolkit/CUDA/v11.5/bin/nvcc.exe" && ^
ninja -j 16 && ^
cd "..\bin\win32\release" && "Test TEFEN.exe"
PAUSE