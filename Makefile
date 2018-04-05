GUROBIDIR = /usr/site/gurobi-5.5.0/
INC       = $(GUROBIDIR)/include/
CLIB      = -Wl,-rpath,$(GUROBIDIR)/lib/ -lgurobi55
CCLIB     = -L$(GUROBIDIR)/lib/ -lgurobi_c++ 

tsp: tsp.cc
	g++ -m64 -std=c++11 -g -o tsp tsp.cc -I$(INC) $(CLIB) $(CCLIB) -lpthread -lm

