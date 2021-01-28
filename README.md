# marzone
Marzone 4.0 has been largely refactored from its original state. It now is mostly written in C++ and utilizes the object oriented aspects of the language for better encapsulation of different functionalities. 

Diagram overview TBC. 

Solvers have been moved to /solver.

## Build command (windows)
```
g++ -O3 -std=c++17 -static -fopenmp marzone.cpp -o bin/marzone
```

## How to run tests (WSL/Linux only)
Prerequisites - install [cpputest](https://cpputest.github.io/manual.html)

- cd test
- make
- run ./test
