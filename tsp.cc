// tsp.cc - traveling salesman code based on Gurobi using branch and cut

#include <string>
#include <iostream>
#include <iomanip>
#include <vector>
#include <cmath>
#include <chrono>
#include <vector>
#include <deque>
#include <limits>
// Gurobi
#include "gurobi_c++.h"

// Data structure to represent an edge of the input graph
struct Edge
{
  int end1;
  int end2;
  double weight;
};

// If memory is limited, a soft limit for the maximum number of LPs in the
// queue can be set for the branch and cut algorithm.
// This is effectively disabled by default.
const int lp_soft_limit = 1000000;

// Tolerance to determine whether or not a number is integral
const double tol = 1.0e-9;

// Auxiliary function to decide whether or not a number is integral
inline bool is_integral(double d)
{
  return std::fabs(d - std::round(d)) < tol;
}

void find_components(int                       num_nodes,
                     int                       num_edges,
                     std::vector<Edge>         edges,
                     const std::vector<double> &x,
                     int                       &num_components,
                     std::vector<int>          &components)
{
  // Mark all unassigned nodes with -1.
  components.resize(num_nodes);
  for (int n = 0; n < num_nodes; n++) {
    components[n] = -1;
  }

  // Component index
  int c = 0;
  // Node indices
  int n1, n2;
  // Find all connected components.
  while (true) {
    // Find an unassigned node.
    n1 = -1;
    for (int n = 0; n < num_nodes; n++) {
      if (components[n] == -1) {
        n1 = n;
        break;
      }
    }
    if (n1 == -1) {
      // All nodes have been assigned.
      break;
    }

    // Assign node to current component.
    components[n1] = c;
    // Mark the entire connected component.
    while (true) {
      // Find an unassigned connected node.
      n2 = -1;
      for (int e = 0; e < num_edges; e++) {
        // Skip edges that are not used in the current solution.
        if (x[e] < tol) continue;
        if (edges[e].end1 == n1 && components[edges[e].end2] == -1) {
          n2 = edges[e].end2;
          break;
        }
        if (edges[e].end2 == n1 && components[edges[e].end1] == -1) {
          n2 = edges[e].end1;
          break;
        }
      }
      // No connected node found. Continue with next component.
      if (n2 == -1) break;

      // Assign the connected node to the current component.
      components[n2] = c;
      // Merge the two nodes.
      for (int e = 0; e < num_edges; e++) {
        if (edges[e].end1 == n2) {
          edges[e].end1 = n1;
        }
        if (edges[e].end2 == n2) {
          edges[e].end2 = n1;
        }
      }
    }
    c++;
  }
  num_components = c;
}

void solve_TSP(int                     num_nodes,
               int                     num_edges,
               const std::vector<Edge> &edges,
               bool                    integral_weights,
               std::vector<double>     &x_opt,
               int                     &lp_solves,
               int                     &subtour_constraints)
{
  // Allocate memory for the solution.
  std::vector<double> x;
  x.resize(num_edges);

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

  lp_solves = 0;
  subtour_constraints = 0;

  // Branch and cut.
  double cost;
  double cost_opt = std::numeric_limits<double>::infinity();
  int num_components;
  std::vector<int> components;
  std::deque<GRBModel> problems;
  problems.push_back(initial_model);
  while (problems.size() > 0) {
    // Get next problem in queue.
    GRBModel model = problems.front();
    problems.pop_front();

    // In this loop, the LP is solved repeatedly until a solution without
    // subtours is found.
    bool skipped = false;
    while (true) {
      // Solve current model.
      model.optimize();
      lp_solves++;

      if (model.get(GRB_IntAttr_Status) != GRB_OPTIMAL) {
        // Do not continue branch if problem is infeasible.
        skipped = true;
        break;
      }

      // Check cost. If it is too high, stop following this branch.
      cost = model.get(GRB_DoubleAttr_ObjVal);
      // We can cut off branches more aggressively when the weights,
      // and thus the optimal cost, are integral.
      if (integral_weights && cost > cost_opt - 1.0 + tol) {
          skipped = true;
          break;
      }
      if (cost > cost_opt) {
        skipped = true;
        break;
      }

      // Get current solution.
      delete[] vars;
      vars = model.getVars();
      for (int e = 0; e < num_edges; e++) {
        x[e] = vars[e].get(GRB_DoubleAttr_X);
      }

      // Find connected components of the solution and eliminate subtours.
      find_components(num_nodes,
                      num_edges,
                      edges,
                      x,
                      num_components,
                      components);

      if (num_components == 1) {
        // There are no more subtours that could be eliminated.
        break;
      }

      // We will add one constraint per connected component.
      std::vector<GRBLinExpr> lhs;
      lhs.resize(num_components);
      for (int e = 0; e < num_edges; e++) {
        // Identify the component this edge belongs to.
        const int c1 = components[edges[e].end1];
        const int c2 = components[edges[e].end2];
        // Skip edges that connect two different components.
        if (c1 != c2) continue;
        // Add edge to the subtour elimination constraint.
        lhs[c1] += GRBLinExpr(vars[e], 1.0);
      }
      // Compute the size of each component.
      // This is required for the right-hand side of the constraints.
      std::vector<int> component_sizes;
      component_sizes.resize(num_components);
      for (int c = 0; c < num_components; c++) {
        component_sizes[c] = 0;
      }
      for (int n = 0; n < num_nodes; n++) {
        component_sizes[components[n]]++;
      }
      // Add constraints to model.
      for (int c = 0; c < num_components; c++) {
        model.addConstr(lhs[c],
                        GRB_LESS_EQUAL,
                        GRBLinExpr(component_sizes[c] - 1.0));
        subtour_constraints++;
      }
      std::cout << "Added " << num_components
                << " subtour elimination constraints." << std::endl;
      model.update();
    }

    if (skipped) continue;

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
        // Add model with >= constraint for fractional variable.
        GRBModel model_ge(model);
        model_ge.addConstr(GRBLinExpr(vars[e], 1.0),
                           GRB_GREATER_EQUAL,
                           GRBLinExpr(std::ceil(x[e])));
        model_ge.update();
        // Check if the soft limit for the number of LPs is hit.
        if (problems.size() < lp_soft_limit) {
          // Add new problems at the end of the queue.
          // This corresponds to breadth-first search.
          problems.push_back(model_le);
          problems.push_back(model_ge);
        }
        else {
          // Add new problems at the beginning of the queue.
          // This corresponds to depth-first search.
          problems.push_front(model_le);
          problems.push_front(model_ge);
        }
        // Print information about the queue.
        std::cout << "Branching; there are now "
                  << problems.size()
                  << " models in the queue."
                  << std::endl << std::endl;
        // Stop after creating one branch!
        break;
      }
    }

    // Update optimal cost and optimal solution if integral solution was found.
    if (integral_sol) {
      cost_opt = cost;
      x_opt = x;
    }
  }
}

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

  // Check if the edge weights are integral.
  // If so, we can optimize a bit more aggressively in some places.
  bool integral_weights = true;
  for (int e = 0; e < num_edges; e++) {
    if (!is_integral(edges[e].weight)) {
      integral_weights = false;
      break;
    }
  }

  std::cout << "Computation begins.\n";
  // Start timer.
  const auto t_start = std::chrono::high_resolution_clock::now();
  // Solve TSP using Gurobi (for the LPs).
  std::vector<double> x_opt;
  int lp_solves;
  int subtour_constraints;
  try {
    solve_TSP(num_nodes,
              num_edges,
              edges,
              integral_weights,
              x_opt,
              lp_solves,
              subtour_constraints);
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

  // Print additional information.
  std::cout << "Solved a total of " << lp_solves << " LPs." << std::endl;
  std::cout << "Added a total of " << subtour_constraints
            << " subtour elimination constraints." << std::endl;

  // Print optimal solution.
  std::cout << "The best tour is:\n";
  double c_optimal = 0.0;
  // Set output format.
  if (integral_weights) {
    std::cout << std::setprecision(0);
  }
  else {
    std::cout << std::setprecision(1);
  }
  for (int e = 0; e < num_edges; e++) {
    // See if the edge is used.
    if (x_opt[e] > 0.0) {
      std::cout << edges[e].end1 << " "
                << edges[e].end2 << " "
                << edges[e].weight << std::endl;
      c_optimal += x_opt[e]*edges[e].weight;
    }
  }
  std::cout << "The cost of the best tour is: " << c_optimal << std::endl;

  return 0;
}

