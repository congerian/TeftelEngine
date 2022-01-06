mkdir .\build\
cd ./build && ^
cmake ../ -G Ninja -DCMAKE_CUDA_COMPILER="C:/Program Files/NVIDIA GPU Computing Toolkit/CUDA/v11.5/bin/nvcc.exe" && ^
ninja -j 16
PAUSE