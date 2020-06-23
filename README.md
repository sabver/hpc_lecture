# my student id
19M38171 

# about final report

## my final version code
1. ***final_report*/cuda_cavity_one_block.cu**

2. other code files in final_report almost are test codes or failure version codes.

## how to verify my code's correctness 
1. final_report/u_json.txt,final_report/v_json.txt,final_report/p_json.txt is the result of my code with parameter THREAD_NUM = 1024

2. or you could run my code to get the file u_json.txt,v_json.txt and p_json.txt,but you should modify the path of output file in code's function write_string_to_file

3. run final_report/cuda_cavity.py to see the picture with u_json.txt and v_json.txt and p_json

# hpc_lecture

|          | Topic                                | Sample code               |
| -------- | ------------------------------------ | ------------------------- |
| Class 1  | Introduction to parallel programming |                           |
| Class 2  | Shared memory parallelization        | 02_openmp                 |
| Class 3  | Distributed memory parallelization   | 03_mpi                    |
| Class 4  | SIMD parallelization                 | 04_simd                   |
| Class 5  | GPU programming                      | 05_cuda,05_openacc        |
| Class 6  | Parallel programing models           | 06_starpu                 |
| Class 7  | Cache blocking                       | 07_cache_cpu,07_cache_gpu |
| Class 8  | High Performance Python              | 08_cython                 |
| Class 9  | I/O libraries                        | 09_io                     |
| Class 10 | Parallel debugger                    | 10_debugger               |
| Class 11 | Parallel profiler                    | 11_profiler               |
| Class 12 | Containers                           |                           |
| Class 13 | Scientific computing                 | 13_fdm,13_solver          |
| Class 14 | Deep Learning                        | 14_dl                     |
