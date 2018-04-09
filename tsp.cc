// tsp.cc - traveling salesman code based on Gurobi using branch and cut

#include <string>
#include <iostream>
#include <iomanip>
#include <vector>
#include <cmath>
#include <chrono>
#include <vector>
#include <queue>
#include <limits>
// Gurobi
#include "gurobi_c++.h"

// Auxiliary function to decide whether or not a number is integral.
bool is_integral(double d)
{
  return std::fabs(d - std::round(d)) < 1.0e-9;
}

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
               std::vector<double>     &x,
               int                     &lp_solves);

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
  int lp_solves;
  try {
    solve_TSP(num_nodes, num_edges, edges, x, lp_solves);
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

  // Print number of LP solves.
  std::cout << "Solved a total of " << lp_solves << " LPs." << std::endl;

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
               std::vector<double>     &x,
               int                     &lp_solves)
{
  // Set up environment.
  GRBEnv env;
  // Create initial model.
  GRBModel initial_model(env);
  // Add variables.
  initial_model.addVars(num_edges, GRB_CONTINUOUS);
  initial_model.update();
  GRBVar *vars = initial_model.getVars();
  // Set up objective function.
  GRBLinExpr obj;
  for (int e = 0; e < num_edges; e++) {
    obj += GRBLinExpr(vars[e], edges[e].weight);
  }
  initial_model.setObjective(obj, GRB_MINIMIZE);
  // Add initial constraints.
  for (int n = 0; n < num_nodes; n++) {
    GRBLinExpr lhs;
    // Add all edges that are adjacent to the current node.
    for (int e = 0; e < num_edges; e++) {
      if (edges[e].end1 == n || edges[e].end2 == n) {
        lhs += GRBLinExpr(vars[e], 1.0);
      }
    }
    initial_model.addConstr(lhs, GRB_EQUAL, GRBLinExpr(2.0));
  }
  for (int e = 0; e < num_edges; e++) {
    initial_model.addConstr(GRBLinExpr(vars[e], 1.0),
                            GRB_LESS_EQUAL,
                            GRBLinExpr(1.0));
  }
  initial_model.update();
  // Make sure that we don't use the wrong set of variables later.
  delete[] vars;
  vars = NULL;

  lp_solves = 0;

  // Branch and cut.
  double cost;
  double opt_cost = std::numeric_limits<double>::infinity();
  std::queue<GRBModel> problems;
  problems.push(initial_model);
  while (problems.size() > 0) {
    // Get next problem in queue.
    GRBModel &model = problems.front();

    // Solve current model.
    model.optimize();
    lp_solves++;

    if (model.get(GRB_IntAttr_Status) != GRB_OPTIMAL) {
      // Do not continue branch if problem is infeasible.
      problems.pop();
      continue;
    }

    // Check cost.
    cost = model.get(GRB_DoubleAttr_ObjVal);
    if (cost >= opt_cost) {
      // Cost is too high; stop following this branch.
      problems.pop();
      continue;
    }

    // Get current solution.
    vars = model.getVars();
    for (int e = 0; e < num_edges; e++) {
      x[e] = vars[e].get(GRB_DoubleAttr_X);
    }

    // Add subtour elimination here and go back to the LP solve...

    // Branch using a fractional variable.
    bool integral_sol = true;
    for (int e = 0; e < num_edges; e++) {
      if (!is_integral(x[e])) {
        integral_sol = false;
        // Add model with <= constraint for fractional variable.
        GRBModel model_le(model);
        model_le.addConstr(GRBLinExpr(vars[e], 1.0),
                           GRB_LESS_EQUAL,
                           GRBLinExpr(std::floor(x[e])));
        model_le.update();
        problems.push(model_le);
        // Add model with >= constraint for fractional variable.
        GRBModel model_ge(model);
        model_ge.addConstr(GRBLinExpr(vars[e], 1.0),
                           GRB_GREATER_EQUAL,
                           GRBLinExpr(std::ceil(x[e])));
        model_ge.update();
        problems.push(model_ge);
        // Print information about the queue.
        std::cout << "Branching; there are now "
                  << problems.size()
                  << " models in the queue" << std::endl;
        // Stop after creating one branch!
        break;
      }
    }

    // Update optimal cost if integral solution was found.
    if (integral_sol) {
      opt_cost = cost;
    }

    // Clean up and remove current problem from queue.
    delete[] vars;
    vars = NULL;
    problems.pop();
  }
}

