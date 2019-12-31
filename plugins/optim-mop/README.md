# Ibex-MOP

*ibexMop* is an interval branch & bound solver for **Nonlinear BiObjective Optimization** problems.

*ibexMop* constructs an **upper envelope** for the non-dominated set
by following a branch & bound strategy starting with an initial *box* (containing the variable domains)
and building a search tree. In each iteration of the algorithm,
a node is selected and treated by *filtering*, *upper-bounding*
and *splitting* techniques.

The solver returns a set of vectors Y' guaranteeing a maximal distance *eps* between
any *non-dominated* feasible vector and the returned set Y'. It includes 
some methods to take into account the *upper envelope* 
for filtering dominated solutions (e.g., dominance peeler).

*ibexMop* includes three upper bounding methods:

  * (inner_polytope) The upper envelope Y' is represented by using a dominace-free set of *vectors*. An
  [inner polytope algorithm](http://citeseerx.ist.psu.edu/viewdoc/download?doi=10.1.1.653.5777&rep=rep1&type=pdf)
  is used for finding feasible solutions. The algorithm constructs a feasible and convex polytope and then it finds
two feasible vectors inside this polytope by minimizing a linearization of each one of the objective functions.
Then it finds a set of N-2 (by default N=50) equidistant feasible solutions.

  * (ub1) The upper envelope Y' is represented by using a set of *upper line segments*. 
  It is warrantied that every vector in any segment is epsilon-dominated 
  by at least one feasible solution. For finding upper line segments, the 
  algorithm first finds a feasible line segment in the domain space using an inner polytope algorithm 
  and then it generates a line segment passing over the image of this feasible line in the objective space. 
  
  * (ub2) The upper envelope Y' is also represented by using a set of upper line segments. 
  The algorithm performs the same procedure than ub1 for finding upper line segments, however
  it goes further and it is able to find the whole upper envelope related to the feasible line segment.

![ub1](https://i.imgur.com/H6zAwpO.png)
Example of an upper line segment.
  
![Cy Comparison](https://i.imgur.com/Wphf10d.png)
Example of using ub2 for finding an upper envelope for the blue feasible curve in the objective space.

## Download
````
git clone https://github.com/INFPUCV/ibex-lib.git
git checkout -t origin/mop-server
````

## Installation

````
    sudo apt install python3 flex bison gcc g++ make pkg-config libz-dev zlib1g-dev python3-tk
    pip3 install -r requirements.txt
    
    ./waf configure --with-optim --with-optim-mop --with-affine --prefix=. --gaol-dir= --lp-lib=soplex
    ./waf install
````      
       
## Main options
```
./__build__/plugins/optim-mop/ibexmop {OPTIONS} [filename]
```
  OPTIONS:

      -h, --help                        Display this help menu
      -f[string], --filt=[string]       the filtering method (default: acidhc4)
      --linear-relax=[string]           the linear relaxation method (default:
                                        compo)
      -b[string], --bis=[string]        the bisection method (default:
                                        largestfirst)
      -s[string], --search=[string]     the search strategy (default: NDSdist)
      --eps=[float]                     the desired precision (default: 0.01)
      -t[float], --timelimit=[float]    timelimit (default: 100)
      --cy-contract                     Contract using the box y+cy, w_ub=+inf.
      --cy-contract-full                Contract using the additional constraint
                                        cy.
      --eps-contract                    Contract using eps.
      --no-bisecty                      Do not bisect y variables.
      --ub=[string]                     Upper bounding strategy (default: ub2).
      --rh=[float]                      Termination criteria for the ub2
                                        algorithm (dist < rh*ini_dist)
      --server_mode                     Server Mode (some options are discativated).
      --server_out=[string]             Server Output File
      --server_in=[string]              Server Instructions File
      -v, --verbose                     Verbose output. Shows the dominance-free
                                        set of solutions obtained by the solver.
      --trace                           Activate trace. Updates of loup/uplo are
                                        printed while minimizing.
      --plot                            Save a file to be plotted by plot.py.
      filename                          The name of the MINIBEX file.
      "--" can be used to terminate flag options and force all following
      arguments to be treated as positional options

### Filtering Method (-f):
 + hc4
 + acidhc4
### Linear Relaxation Method (--linear-relax):
 + no
 + compo (a method combining two linearization techniques: AF2 and XNewton)
### Search Strategy (-s):
 + weighted_sum (or the [OC search strategy](http://www.sciencedirect.com/science/article/pii/S0377221716303824))
 + NDSdist
### Bisection Method (-b):
 + largestfirst
 + roundrobin
### Upper bounding method (-ub):
 + inner_polytope
 + ub1
 + ub2

## Run an example:

     ./__build__/plugins/optim-mop/ibexmop plugins/optim-mop/benchs/osy.txt  --cy-contract-full --eps-contract --ub=ub2 --eps=0.0001

## Run an example (server mode):

    ./__build__/plugins/optim-mop/ibexmop plugins/optim-mop/benchs/tan.txt --cy-contract-full --eps_r=0.01 --ub=ub2 --server_mode --server_in=intructions.txt --server_out=output2.txt
    
    ./__build__/plugins/optim-mop/ibexmop plugins/optim-mop/benchs/osy.txt --cy-contract-full --eps_r=0.001 --ub=ub2 --server_mode --server_in=intructions.txt --server_out=output2.txt

The output file (server_out) contains two lists of points (segments) representing an upper and lower envelope for the optimal solutions.

The input file (server_in) allows us to give instructions to the solver. For the moment two commands:
* zoom_in y1_lb y1_ub y2_lb y2_ub
* zoom_out y1_lb y1_ub y2_lb y2_ub
* get_solution output_file y1 y2  (the point y=(y1,y2) belonging to a segment, the solver create a file output_file with a solution vector x and its image f(x) which dominates y)

For plotting the non-dominated vectors returned by the solver

     python3 plugins/optim-mop/main/plot3.py



## Format of the instances (Minibex):

Instances can be written in the [Minibex language](http://www.ibex-lib.org/doc/minibex.html),
considering that the objectives *must correspond* to the first two constraints with the following syntax:
```
Constraints
<expression of the first objective function> = z1;
<expression of the second objective function> = z2;
// other constraints...
```
You can see the set of instances in [plugins/optim-mop/benchs](https://github.com/INFPUCV/ibex-lib/tree/master-mop/plugins/optim-mop/benchs).


## Authors:
 - Ignacio Araya - <ignacio.araya@pucv.cl>
