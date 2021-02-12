# Marzone 4.0
Marzone 4.0 has been largely refactored from its original state. It now is mostly written in C++ and utilizes the object oriented aspects of the language for better encapsulation of different functionalities. 

## Latest Release
### **Marzone 4 is here for you to try**
 
Marzone has been rebuilt in C++ and we are ready for you to test the new version 4.  
 
**Upgrades include:**  
- support for multithreading and parallelization (full utilisation of CPU)  
- improved computational speed and efficiency  
- additional error reporting  
 
Download for Windows / MacOS / Linux: [TBC]()  

You can access the code for the new version in a separate branch of this repo: /marzone4  
We will officially announce the release of the new version and merge changes to main soon, and are looking for your feedback in the meantime.  
Please feel free to use our [Google Group](https://groups.google.com/g/marxan). You can also create GitHub issues on this repo or e-mail marxancloud@gmail.com  


## Overview of Code Structure
If you want to start contributing to the codebase, we encourage you to familiarise yourself with the type of problems that Marzone is intended to solve, and be familiar with the input files and their roles within the program. 
- The Marzone manual is located at the [marxan website](https://marxansolutions.org/)

The original codebase has been rewritten and major parts of it separated into classes. Here is a brief overview of the classes, their dependencies and their intended functionality:

| File/Class      | Dependencies | Description |
|-----------------|--------------|-------------|
| common.hpp      | None      | Contains all common data structures for marzone, used by all components of the program |
| util.hpp        | None      | Contains stand-alone utility functions used by most files in the program |
| logger.hpp      | common    | Logging functions for debug printing, error/warning printing and writing to log files. |
| costs.hpp/Costs | common, util     | Parses cost files and stores information relating to costs, cost names and values. |
| species.hpp/Species         | common, util     | Parses species file and stores mapping information relating to species lookup, targets and penalties |
| pu.hpp/Pu | species, costs, logger| Parses pu, pulock, puzone, pu connection/bound and puvspr files and stores pu mapping information, pu costs and species amounts in each pu. Also responsible for computations relating to these values |
| zones.hpp/Zones | pu | Parses zone information, zone costs, zone targets, zone boundary, zone contributions and other input files relating to zones. Stores the zone mapping and values relating to zones. Performs computation on these values when needed.|
| reserve.hpp/Reserve | zones | Holds information specific to a single reserve configuration, including its species amounts, shortfall, objective value, connection and costs. Also responsible for computing changes in reserve units.|



\
The solvers have been moved to **/solver** folder. Currently Marzone supports the following:
- Simulated Annealing 
- Iterative Improvement
- General Heuristic methods specific to Marzone.

All solvers depend on **reserve, pu, zones and species** components. To add a new solver, study and follow the existing pattern of the solvers under /solver. 


## Build command (windows and linux)
```
g++ -O3 -std=c++17 -static -fopenmp marzone.cpp -o bin/marzone
```

## Build command (Mac)
```
TBC
```

## How to run tests (WSL/Linux only)
Prerequisites - install [cpputest](https://cpputest.github.io/manual.html). 
Please ensure any changes to computation does not break the existing unit tests. 

- cd test
- make
- ./test
