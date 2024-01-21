# uBlade:Efficient Batch Processing for Uncertain Graph Queries
Organization
--------

This repository contains code for our paper "uBlade:Efficient Batch Processing for Uncertain Graph Queries". It is implemented based on Julian Shun's **Ligra** https://github.com/jshun/ligra


Compilation
--------

Compilation is done from within the apps/ directory.

Compilers

* g++ &gt;= 5.3.0 with OpenMP


to compile, run

```
$ make -j  #compiles with all threads
$ g++ -std=c++14 -O3 -fopenmp -DBYTERLE -o weighted weighted.C
```

The following commands cleans the directory:
```
$ make clean #removes all executables
$ make cleansrc #removes all executables and linked files from the ligra/ directory
```
Running code
-------

```
$ ./weighted -w -src 2745 -tar 5292 -delta 1000 -nb 1024 -batch 20 -sample 100 -thread 10 ../input/gnu_in.txt
```
The input graph is always treated as an weighted graph for edge probabiliies, hence flag "-w" is necessary.
Flags "-src" and "-tar" stand for the source node and target node number. The code we present are weighted and unweighted examples of graph traversal, target node only depends the output.
Flags "-delta" and "-nb" are delta-stepping parameters delta and number of buckets.
Flags "-batch", "-sample" and "-thread" stand for batch_size, total sample number and thread numbers.


Input Formats
-----------
We support the adjacency graph format used by the [Problem Based Benchmark suite](http://www.cs.cmu.edu/~pbbs/benchmarks/graphIO.html) and [Ligra](https://github.com/jshun/ligra).

The adjacency graph format starts with a sequence of offsets one for each vertex, followed by a sequence of directed edges ordered by their source vertex. The m edge weights should be stored after all of the edge targets in the .adj file. The offset for a vertex i refers to the location of the start of a contiguous block of out edges for vertex i in the sequence of edges. The block continues until the offset of the next vertex, or the end if i is the last vertex. All vertices and offsets are 0 based and represented in decimal. The specific format is as follows:

```
WeightedAdjacencyGraph
<n>
<m>
<o0>
<o1>
...
<o(n-1)>
<e0>
<e1>
...
<e(m-1)>
<w0>
<w1>
...
<w(m-1)>
```

The edge weight &lt;w> is a compressed integer, we use three decimal places to express the probability of edges, hence w = (length + probability) * 1000. For example, an edge &lt;e> with length 26 and probability 0.250 has &lt;w> = 26250.

This file is represented as plain text.

User Define Applications
-----------------------
We provide basic usage for both weighted and unweighted graphs, with and without edge weights. The output of the applications is the number of possible worlds that reach the inputted target nodes with the most probable shortest distance (if weighted). All sampling results are stored in the array "distances," which is an n * sample_number-sized array where all nodes' reliability and distance are stored consecutively. All applications listed in our experiments can be performed by collecting information from this array. For example, for SSSD, the most probable shortest distance from source to any node t can be obtained by collecting information from distances[i], where i varies from t×sample_number to (t×sample_number−1).

