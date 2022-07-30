# Escherization with Large Deformations Based on As-Rigid-As-Possible Shape Modeling
This code is a C++ implementation of the algorithms developed in the following papers. 

[1] Yuichi Nagata and Shinji Imahori, An Efficient Exhaustive Search Algorithm for the Escherization Problem, Algorithmica, Vol.82, No.9, 2502-2534, 2020. [[link]](https://link.springer.com/article/10.1007/s00453-020-00695-6)

[2] Yuichi Nagata and Shinji Imahori, Escherization with Large Deformations Based on As-Rigid-As-Possible Shape Modeling, ACM Transactions on Graphics, Vol.41, No.2:11, 1-16, 2022. [[link]](https://dl.acm.org/doi/full/10.1145/3487017)

The following algorithms are implemented.
- EST with the Euclidean or AD distance 
- EST with the E_I or E_IR distance 
- Heuristic search with the E_I or E_IR distance 

# Environment
The author ran the programs on Ubuntu 16.04 (or 20.04), and the following preparations may be required to compile the source code. 
- Make the Eigen library (https://eigen.tuxfamily.org) available on your PC. 
- Install the X11 library on your PC. `$ sudo apt install libx11-dev`

# Goal figures
Directory Data contains some goal figure files, which are divided into two types. 
- Polygonal type: e.g. pegasus_60.dat 
- Mesh type: e.g. pegasus_60_27.dat 

See document.txt in directory Data for the data format.  

# How to use 
It is recommended to run EST with the Euclidean distance first because this is the most basic one. 

## EST with the Euclidean or AD distance 
Before compiling the program, rewrite parts of the program directly. 
- Set EST in line 4 of env_E.h.
- Set E or AD in line 4 of search_E_EST.h, depending on whether the Euclidean or AD distance is used.  

Compile:
```
$ ./build_E.exe
$ mv jikken jikken_E
```
Then, a executable file jikken_E is created.

Execution:  
```
$ ./jikken_E <string1> <string2> <integer1> <integer2>
```
  
&nbsp;　\<string1\> : file name of a goal figure (polygonal type)  
&nbsp;　\<string2\> : file name to which results are written  
&nbsp;　\<integer1\> : results are displayed?  0:no, 1:yes  
&nbsp;　\<integer2\> : the number of top solutions stored  

(Example)
```
$ ./jikken_E pegasus_60.dat ABC 1 20
```
   
Result:  
The specified number of top solutions are saved in a file (ABC_pegasus_60.dat.tile in this example) and are displayed in the terminal, as shown in the example below 
```
order =  0 eval=18994.244366: IH= 6 (18 4 18 7 4 3 ) 22  
order =  1 eval=21028.815657: IH= 6 (5 14 5 12 14 4 ) 48  
order =  2 eval=21868.573515: IH= 4 (13 4 7 13 14 3 ) 20  
....  
```

Each line represent the following values:  
&nbsp; \- order: ranking of the solution  
&nbsp; \- eval : the distance value  
&nbsp; \- IH   : Isohedral type  
&nbsp; \- (...) : the numbers of points assigned to the tiling edges of the template  
&nbsp; \- bright-most number : the value of j (see [1] or [2])  

If \<integer1\> = 1, the top solutions (tile figures) are displayed. Press the return key to display the next solution. 

The top solutions stored in a file can also be displayed (see "Display Result"). 


## Display Results

The solutions saved in a file (e.g. ABC_pegasus_60.dat.tile) can be displayed using this program. 

Compile:  
```
./build_display.exe
```
Then, the executable file "display.exe" is created.

Execution:  
```
$ ./display.exe <string1> 
```
&nbsp; \<string1\> : file name of a result file 

Press the return key to display the next solution. 


## EST with the E_I or E_IR distance 
Before compiling the program, rewrite parts of the program directly. 
- Set EST in line 4 of env_I.h.
- Set I or IR in line 4 of search_I_EST.h, depending on whether the E_I or E_IR distance is used.

Compile:
```
$ ./build_I.exe
$ mv jikken jikken_I
```
Then, a executable file jikken_I is created.

Execution:  
```
$ ./jikken_I <string1> <string2> <integer1> <integer2> <double1>
```
  
&nbsp;　\<string1\> : file name of a goal figure (mesh type)  
&nbsp;　\<string2\> : file name to which results are written  
&nbsp;　\<integer1\> : results are displayed?  0:no, 1:yes  
&nbsp;　\<integer2\> : the number of top solutions stored  
&nbsp;　\<double1\> : the value of $\alpha$ (see [2]) (0.5 is recommended)

(Example)
```
$ ./jikken_I pegasus_60_27.dat ABC 1 20 0.5
```
   
Result:  
Same as EST with the Euclidean distance.

## Heuristic search with the E_I or E_IR distance 

**Step 1: Obtaining promising configurations**

Before compiling the program, rewrite a part of the program directly. 
- Set CONF in line 4 of env_E.h.

Compile:
```
$ ./build_E.exe
$ mv jikken jikken_get_conf
```
Then, a executable file jikken_get_conf is created.

Execution:  
```
$ ./jikken_get_conf <string1> <string2> <integer1> <integer2>
```
  
&nbsp;　\<string1\> : file name of a goal figure (polygonal type)  
&nbsp;　\<string2\> : file name to which results are written  
&nbsp;　\<integer1\> : set 0 or 1 (a parameter not used)  
&nbsp;　\<integer2\> : the number of top configurations stored (10000-50000 is recommended)

(Example)
```
$ ./jikken_get_conf pegasus_60.dat ABC 0 50000
```
Then, the specified number of promising configurations are saved in a file (ABC_pegasus_60.dat.conf in this example).

**Step 2: Perform optimization for the specified configurations**

Before compiling the program, rewrite parts of the program directly. 
- Set CONF in line 4 of env_I.h.
- Set I or IR in line 4 of search_I_EST.h, depending on whether the E_I or E_IR distance is used.  

Compile:
```
$ ./build_I.exe
$ mv jikken jikken_I_conf
```
Then, a executable file jikken_I_conf is created.

Execution:  
```
$ ./jikken_I_conf <string1> <string2> <integer1> <integer2> <double1> <string3>
```
  
&nbsp;　\<string1\> : file name of a goal figure (mesh type)  
&nbsp;　\<string2\> : file name to which results are written  
&nbsp;　\<integer1\> : results are displayed?  0:no, 1:yes  
&nbsp;　\<integer2\> : the number of top solutions stored  
&nbsp;　\<double1\> : the value of $\alpha$ (see [2]) (0.5 is recommended)  
&nbsp;　\<string3\> : file name of a promising configuration data  

(Example)
```
$ ./jikken_I_conf pegasus_60_27.dat ABC 1 20 0.5 ABC_pegasus_60.dat.conf
```
   
Result:  
Same as EST with the Euclidean distance.


