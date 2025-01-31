# RaMA
![Platforms](https://anaconda.org/malab/rama/badges/platforms.svg) ![last_release](https://anaconda.org/malab/rama/badges/latest_release_date.svg) ![license](https://anaconda.org/malab/rama/badges/license.svg)

**RaMA** (Rare Match Aligner) is a high-performance tool specialized in **pairwise centromere alignment**, optimized for ultra-long tandem repeats in complex genomic regions. Leveraging rare matches as anchors and a 2-piece affine gap cost model, RaMA is specifically designed for aligned centromere sequences, delivering highly accurate and efficient pairwise alignments for **assembled centromere regions**.

## Get Started
RaMA can be installed using two methods: via conda or by compiling the source code. 
- If you want to quickly get started with RaMA, we highly recommend using the conda installation method. 
- if you prefer to customize RaMA, such as enabling M64 mode to handle ultra-large-scale data (not generally recommended) or using different instruction sets to accelerate alignment, we suggest installing from the source code.
### Install via conda
We highly recommand you install RaMA via conda. Simply run the following command:
~~~
conda create -n RaMA
conda activate RaMA
conda install -c malab -c conda-forge rama
RaMA -h
~~~
### Install via source code

#### Requirements
Before installing RaMA, ensure your **Linux** system meets the following requirements:
- GCC version 9 or higher. (**Required**)
- CMake version 3.1 or higher. (**Required**)
#### Clone the Repository
Start by cloning the RaMA repository along with its submodules. Ensure that `Alignment/WFA2-lib` is downloaded properly:
~~~sh
git clone --recursive https://github.com/metaphysicser/RaMA.git
~~~
#### Compilation
Navigate to the RaMA directory, create a build directory, and compile the program:
~~~sh
cd RaMA 
mkdir build && cd build
cmake .. -DCMAKE_BUILD_TYPE=Release
make -j$(nproc)
./RaMA -h
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
RaMA -r /path/to/ref.fasta -q /path/to/query.fasta -o /path/to/output_dir
~~~
Note! Both ref.fasta and query.fasta should contain only one sequence.
### Detailed Usage
Control RaMA's behavior with the following arguments:
~~~plaintext
Arguments:
    -h, --help               Show this help message and exit 
    -r, --reference          Reference FASTA file path containing the reference sequences for alignment.
    -q, --query              Query FASTA file path containing the query sequences for alignment.
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

    -a, --sam_output        Whether to output in SAM format. If specified, results will be saved in SAM format.
    -p, --paf_output        Whether to output in PAF format. If specified, results will be saved in PAF format.
~~~
**Note**: The -r -q (input) and -o (output) parameters are required. For most users, the default settings are recommended, but feel free to adjust them as needed.

## Output Files
After running RaMA, the following output files will be generated in the specified output directory:

1. **output.fasta**: Alignment results in FASTA format.
2. **cigar.txt**: Alignment results in CIGAR format.
3. reliable_region.csv: Contains information about the alignment reliability and rare matches. We define a rare match as an alignment with high reliability. Additionally, if the alignment result between rare matches consists of a simple insertion, simple deletion, or a perfect match, this interval is also considered part of the Rare Match-based Reliable Region. Example content:
   ```plaintext
   cigar    reliablity    rare match
   1178I    1             0
   326=     1             1
   1X       1             0
   1073=    1             1
   1X       1             0
- confidence: Indicates whether the alignment is considered reliable.
- rare match: Indicates if the corresponding CIGAR is a rare match.
4. first_anchor.csv: Rare match anchors obtained during the first iteration.
5. final_anchor.csv: All rare match anchors after the final iteration.
6. intervals_need_align.csv: Regions that require wavefront alignment.
7. RaMA.log: Contains information about the alignment process.
8. output.sam: If the -a option is selected, the result will be saved in SAM format. 
9. output.paf: If the -p option is selected, the result will be saved in PAF format.

## License
[Apache 2.0](https://github.com/metaphysicser/RaMA/blob/master/LICENSE) © [[MALABZ_UESTC](https://github.com/malabz) [Pinglu Zhang](https://metaphysicser.github.io/)]

## Contact

RaMA is actively developed and maintained by [ZOU's Lab](https://github.com/malabz). For any questions, feedback, or support, we encourage you to connect with us:

- **Issue Tracker**: For reporting bugs, suggesting enhancements, or requesting features, please use our [GitHub Issues](https://github.com/metaphysicser/RaMA/issues).
- **Email**: For direct communication, you can reach out to Pinglu Zhang at [pingluzhang@outlook.com](mailto:pingluzhang@outlook.com).

