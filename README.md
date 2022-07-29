# Algorithms for the Escherization problem
This code is a C++ implementation of the algorithms developed in the following papers. 

[1] Yuichi Nagata and Shinji Imahori, An Efficient Exhaustive Search Algorithm for the Escherization Problem, Algorithmica, Vol.82, No.9, 2502-2534, 2020. [[link]](https://link.springer.com/article/10.1007/s00453-020-00695-6)

[2] Yuichi Nagata and Shinji Imahori, Escherization with Large Deformations Based on As-Rigid-As-Possible Shape Modeling, ACM Transactions on Graphics, Vol.41, No.2:11, 1-16, 2022. [[link]](https://dl.acm.org/doi/full/10.1145/3487017)

The following algorithms are implemented.
- EST with the Euclidean (or AD) distance 
- EST with the E_I (or E_IR) distance 
- Heuristic search with the E_I (or E_IR) distance 

# Directory structure
.
|-- README.md
|-- Program (source codes of the algorithms)
|-- Data (example date of goal figures)

# Environment
The author ran the programs on Ubuntu 16.04 (or 20.04), and the following preparations may be required to compile the source code. 

- Make the Eigen library (https://eigen.tuxfamily.org) available on your PC. 
- Install the X11 library on your PC. 
  sudo apt install libx11-dev 
