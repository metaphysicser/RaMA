# RaMA

## Requirements
- GCC 9+ or any C++17 standard-compliant compiler (**Required**)
- CMake 3.1+ (**Required**)

## Get Started
Git clone the program. Please check Alignment/WFA2-lib has been downloaded.
~~~sh
git clone --recursive https://github.com/metaphysicser/RaMA.git
~~~
Complile the program 
~~~sh
cd RaMA 
mkdir build && cd build
cmake ..
make -j$(nproc)
~~~
There are some flags that can be used. For instance:
~~~sh
cmake .. -DUSE_M64=ON
cmake .. -DEXTRA_FLAGS="-mavx2"
~~~

## How to use RaMA
simple case:
~~~sh
./RaMA -i /path/to/input.fasta -o /path/to/output_dir
~~~
