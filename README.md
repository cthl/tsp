# tsp
Traveling salesman code based on Gurobi using branch and cut

Simply run `make` to compile the program, after adjusting the path to Gurobi in Makefile.
The program reads in a TSP from stdin, e.g., `cat tests/hk48.txt | ./tsp`.
