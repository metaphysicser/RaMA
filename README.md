# RaMA

## Requirements
Before installing RaMA, ensure your **Linux** system meets the following requirements:
- GCC version 9 or higher. (**Required**)
- CMake version 3.1 or higher. (**Required**)

## Get Started
### Clone the Repository
Start by cloning the RaMA repository along with its submodules. Ensure that `Alignment/WFA2-lib` is downloaded properly:
~~~sh
git clone --recursive https://github.com/metaphysicser/RaMA.git
~~~
### Compilation
Navigate to the RaMA directory, create a build directory, and compile the program:
~~~sh
cd RaMA 
mkdir build && cd build
cmake .. -DCMAKE_BUILD_TYPE=Release
make -j$(nproc)
~~~
Compilation Flags
Customize your build with additional flags if necessary:
~~~sh
# Enable 64-bit mode
cmake .. -DCMAKE_BUILD_TYPE=Release -DUSE_M64=ON

# Use extra compiler flags, e.g., for AVX2 support
cmake .. -DCMAKE_BUILD_TYPE=Release -DEXTRA_FLAGS="-mavx2"
~~~

## How to use RaMA
### General Usage
Execute RaMA by specifying the input FASTA file and the output directory:
~~~sh
./RaMA -i /path/to/input.fasta -o /path/to/output_dir
~~~

### Detailed Usage
Control RaMA's behavior with the following arguments:
~~~plaintext
Arguments:
    -h, --help               Show this help message and exit 
    -i, --input              Input FASTA file path containing the sequences for alignment.
    -o, --output             Output directory path for saving alignment results and additional files.
    -t, --threads            Number of threads for the alignment process. Defaults to the number of available cores if unspecified.
    -s, --save               Saves anchor binary files to the output directory for future use, including SA, LCP, and Linear Sparse Table.
    -l, --load               Loads existing anchor binary files from the output directory to skip SA, LCP, and Linear Sparse Table construction.
    -c, --max_match_count    Maximum number of rare matches to use for anchor finding. Altering this value is generally not recommended.
    -m, --match              Match score for sequence alignment. Lower values favor matching characters. Default is 0.
    -x, --mismatch           Mismatch penalty. Higher values penalize mismatches more. Default is 3.
    -g, --gap_open1          Penalty for initiating a short gap. Key for handling different gap lengths. Default is 4.
    -e, --gap_extension1     Penalty for extending a short gap. Less severe than gap opening penalty. Default is 2.
    -G, --gap_open2          Penalty for initiating a long gap. Aims to manage long gaps strategically. Default is 12.
    -E, --gap_extension2     Penalty for extending a long gap. Provides a lenient approach to long gap management. Default is 1.
~~~
**Note**: The -i (input) and -o (output) parameters are required. For most users, the default settings are recommended, but feel free to adjust them as needed.

## License
[Apache 2.0](https://github.com/metaphysicser/RaMA/blob/master/LICENSE) © [[MALABZ_UESTC](https://github.com/malabz) [Pinglu Zhang](https://github.com/metaphysicser)]

## Contact

RaMA is actively developed and maintained by [ZOU's Lab](https://github.com/malabz). For any questions, feedback, or support, we encourage you to connect with us:

- **Issue Tracker**: For reporting bugs, suggesting enhancements, or requesting features, please use our [GitHub Issues](https://github.com/metaphysicser/RaMA/issues).
- **Email**: For direct communication, you can reach out to Pinglu Zhang at [pingluzhang@outlook.com](mailto:pingluzhang@outlook.com).

