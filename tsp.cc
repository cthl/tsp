// tsp.cc - traveling salesman code based on Gurobi using branch and cut

#include <string>
#include <iostream>
#include <iomanip>
#include <vector>
#include <cmath>
#include <chrono>
#include <vector>
#include <algorithm>
// Gurobi
#include "gurobi_c++.h"

// Data structure to represent an edge of the input graph
struct Edge
{
  int end1;
  int end2;
  double weight;
};

// Function declarations
void solve_TSP(int                     num_nodes,
               int                     num_edges,
               const std::vector<Edge> &edges,
               std::vector<double>     &x);

int main(int argc, char **argv)
{
  // Read problem from stdin.
  std::string line;
  // Read problem size.
  std::getline(std::cin, line);
  // Remove leading spaces.
  while (line[0] == ' ') line = line.substr(1);
  const int num_nodes = std::stoi(line.substr(0, line.find(" ")));
  const int num_edges = std::stoi(line.substr(line.find(" ") + 1));
  // Read graph.
  std::vector<Edge> edges;
  edges.resize(num_edges);
  for (int e = 0; e < num_edges; e++) {
    std::getline(std::cin, line);
    // Remove leading spaces.
    while (line[0] == ' ') line = line.substr(1);
    edges[e].end1 = std::stoi(line.substr(0, line.find(" ")));
    line = line.substr(line.find(" ") + 1);
    edges[e].end2 = std::stoi(line.substr(0, line.find(" ")));
    line = line.substr(line.find(" ") + 1);
    edges[e].weight = std::stod(line);
  }
  std::cout << "Loaded TSP with " << num_nodes << " nodes and "
                                  << num_edges << " edges.\n";

  // Allocate memory for the solution.
  std::vector<double> x;
  x.resize(num_edges);

  std::cout << "Computation begins.\n";
  // Start timer.
  const auto t_start = std::chrono::high_resolution_clock::now();
  // Solve TSP using Gurobi (for the LPs).
  try {
    solve_TSP(num_nodes, num_edges, edges, x);
  }
  catch (const GRBException &e) {
    std::cerr << "Gurobi exception: " << e.getMessage() << std::endl;
    std::exit(1);
  }
  // Stop timer.
  const auto t_end = std::chrono::high_resolution_clock::now();
  const std::chrono::duration<double> dtime = t_end - t_start;

  std::cout << "Computation finished (" 
            << std::fixed << std::setprecision(3)
            << dtime.count() << "s).\n";

  // Print optimal solution.
  std::cout << "The best tour is:\n";
  double c_optimal = 0.0;
  for (int e = 0; e < num_edges; e++) {
    // See if the edge is used.
    if (x[e] > 0.0) {
      std::cout << edges[e].end1 << " "
                << edges[e].end2 << " "
                << std::setprecision(1)
                << edges[e].weight << std::endl;
      c_optimal += x[e]*edges[e].weight;
    }
  }
  std::cout << "The cost of the best tour is: " << c_optimal << std::endl;

  return 0;
}

void solve_TSP(int                     num_nodes,
               int                     num_edges,
               const std::vector<Edge> &edges,
               std::vector<double>     &x)
{
  // Set up environment.
  GRBEnv env;
  // Create initial model.
  GRBModel model(env);
  // Add variables.
  model.addVars(num_edges, GRB_CONTINUOUS);
  model.update();
  GRBVar *vars = model.getVars();
  // Set up objective function.
  GRBLinExpr obj;
  for (int e = 0; e < num_edges; e++) {
    obj += GRBLinExpr(vars[e], edges[e].weight);
  }
  model.setObjective(obj, GRB_MINIMIZE);
  // Add initial constraints.
  for (int n = 0; n < num_nodes; n++) {
    GRBLinExpr lhs;
    // Add all edges that are adjacent to the current node.
    for (int e = 0; e < num_edges; e++) {
      if (edges[e].end1 == n || edges[e].end2 == n) {
        lhs += GRBLinExpr(vars[e], 1.0);
      }
    }
    model.addConstr(lhs, GRB_EQUAL, GRBLinExpr(2.0));
  }
  model.update();

  // Replace this with branch and cut.
  model.optimize();

  // Extract solution.
  for (int e = 0; e < num_edges; e++) {
    x[e] = vars[e].get(GRB_DoubleAttr_X);
  }
}

