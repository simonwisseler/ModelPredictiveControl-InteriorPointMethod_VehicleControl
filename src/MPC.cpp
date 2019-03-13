#include "MPC.h"
#include <cppad/cppad.hpp>
#include <cppad/ipopt/solve.hpp>
#include <iostream>
#include <string>
#include <vector>
#include "Eigen-3.3/Eigen/Core"

using CppAD::AD;
using Eigen::VectorXd;

// predicting one second into the future
size_t N = 10;
double dt = 0.1;

// This value assumes the model presented in the classroom is used.
//
// It was obtained by measuring the radius formed by running the vehicle in the
//   simulator around in a circle with a constant steering angle and velocity on
//   a flat terrain.
//
// Lf was tuned until the radius formed by simulating the model
//   presented in the classroom matched the previous radius.
//
// This is the length from car's front to center of gravity that has a similar radius.
const double Lf = 2.67;

// starting indices for block of x-variables, y-variables, ... within vars vector
// for each state variable, there are N variables (one for each timestep)
// for each control variable, there are N-1 variables (one for each timestep up to N-1, but none for timestep N)
unsigned int x_start = 0;
unsigned int y_start = x_start + N;
unsigned int psi_start = y_start + N;
unsigned int v_start = psi_start + N;
unsigned int cte_start = v_start + N;
unsigned int psi_error_start = cte_start + N;
unsigned int delta_start = psi_error_start + N;
unsigned int a_start = delta_start + N - 1;

// reference velocity to be used in cost function (deviations from reference velocity are penalized)
double v_ref = 70;

class FG_eval {
public:
    // fitted polynomial coefficients
    VectorXd coeffs;
    FG_eval(VectorXd coeffs) { this->coeffs = coeffs; }
    
    typedef CPPAD_TESTVECTOR(AD<double>) ADvector;
    void operator()(ADvector& fg, const ADvector& vars) {
        /**
         * `fg` is a vector of the cost constraints, `vars` is a vector of variable
         *   values (state & actuators)
         */
        // cost function
        fg[0] = 0;
        
        for (unsigned int i = 0; i < N; i++) {
            // cost from cross-track errors
            fg[0] += 3000 * CppAD::pow(vars[cte_start + i], 2);
            // cost from heading errors
            fg[0] += 3000 * CppAD::pow(vars[psi_error_start + i], 2);
            // Cost from deviations from reference velocity
            fg[0] += 1 * CppAD::pow(vars[v_start + i] - v_ref, 2);
        }
        
        for (unsigned int i = 0; i < N - 1; i++) {
            // Cost from acceleration
            fg[0] += 5 * CppAD::pow(vars[delta_start + i], 2);
            // Cost from steering
            fg[0] += 5 * CppAD::pow(vars[a_start + i], 2);
            // Cost from high combined acceleration and steering
            fg[0] += 700 * CppAD::pow(vars[delta_start + i] * vars[v_start+i], 2);
        }
        
        for (unsigned int i = 0; i < N - 2; i++) {
            // Cost from changes in steering value (make the ride less bumpy!)
            fg[0] += 200 * CppAD::pow(vars[delta_start + i + 1] - vars[delta_start + i], 2);
            // Cost from changes in acceleration (make the ride less bumpy!)
            fg[0] += 10 * CppAD::pow(vars[a_start + i + 1] - vars[a_start + i], 2);
        }
        
        // initial state constraints
        // (add 1 to each of the starting indices as total cost at fg[0])
        fg[1 + x_start] = vars[x_start];
        fg[1 + y_start] = vars[y_start];
        fg[1 + psi_start] = vars[psi_start];
        fg[1 + v_start] = vars[v_start];
        fg[1 + cte_start] = vars[cte_start];
        fg[1 + psi_error_start] = vars[psi_error_start];
        
        // kinematic bicycle model constraints
        // (loop from t=1 to t=N-1 as values for timestep 0 already defined under initial state constraints)
        // [variable]1: [variable] at timestep t, [variable]0: [variable] at timestep t-1
        for (unsigned int t = 1; t < N; t++) {
            AD<double> x1 = vars[x_start + t];
            AD<double> x0 = vars[x_start + t - 1];
            AD<double> y1 = vars[y_start + t];
            AD<double> y0 = vars[y_start + t - 1];
            AD<double> psi1 = vars[psi_start + t];
            AD<double> psi0 = vars[psi_start + t - 1];
            AD<double> v1 = vars[v_start + t];
            AD<double> v0 = vars[v_start + t - 1];
            AD<double> cte1 = vars[cte_start + t];
            AD<double> cte0 = vars[cte_start + t - 1];
            AD<double> psi_error1 = vars[psi_error_start + t];
            AD<double> psi_error0 = vars[psi_error_start + t - 1];
            
            AD<double> delta = vars[delta_start + t - 1];
            AD<double> a = vars[a_start + t - 1];
            // to account for latency, use previous control values
            if (t > 1) {
                a = vars[a_start + t - 2];
                delta = vars[delta_start + t - 2];
            }
            
            AD<double> f0 = coeffs[0] + coeffs[1] * x0 + coeffs[2] * CppAD::pow(x0, 2) + coeffs[3] * CppAD::pow(x0, 3);
            // psi_des0 = atan(rise/run)
            AD<double> psi_des0 = CppAD::atan(coeffs[1] + 2 * coeffs[2] * x0 + 3 * coeffs[3] * CppAD::pow(x0, 2)); // returns value in radians
            
            // as above, add 1 to each of starting indices as total cost at fg[0] = 0
            fg[1 + x_start + t] = x1 - (x0 + v0 * CppAD::cos(psi0) * dt);
            fg[1 + y_start + t] = y1 - (y0 + v0 * CppAD::sin(psi0) * dt);
            fg[1 + psi_start + t] = psi1 - (psi0 - v0/Lf * CppAD::tan(delta) * dt);
            fg[1 + v_start + t] = v1 - (v0 + a * dt);
            fg[1 + cte_start + t] = cte1 - ((f0 - y0) + (v0 * CppAD::sin(psi_error0) * dt));
            fg[1 + psi_error_start + t] = psi_error1 - ((psi0 - psi_des0) - v0/Lf * CppAD::tan(delta) * dt);
        }
    }
};

//
// MPC class definition implementation.
//
MPC::MPC() {}
MPC::~MPC() {}

std::vector<double> MPC::Solve(const VectorXd &state, const VectorXd &coeffs) {
    bool ok = true;
    typedef CPPAD_TESTVECTOR(double) Dvector;
    
    // N variables for each state variable and (N-1) variables for each control variable
    unsigned int n_vars = N * 6 + (N-1) * 2;
    
    // one constraint for each state variable and timestep
    unsigned int n_constraints = N * 6;
    
    
    // initial value of the independent variables
    // (all 0 besides initial state)
    Dvector vars(n_vars);
    
    for (int i = 0; i < n_vars; i++) {
        vars[i] = 0;
    }
    
    vars[x_start] = state[0];
    vars[y_start] = state[1];
    vars[psi_start] = state[2];
    vars[v_start] = state[3];
    vars[cte_start] = state[4];
    vars[psi_error_start] = state[5];
    
    
    // lower and upper limits for variables
    // (all [max negative, max positive] besides actuators)
    Dvector vars_lowerbound(n_vars);
    Dvector vars_upperbound(n_vars);
    
    for (unsigned int i = 0; i < delta_start; i++) {
        vars_lowerbound[i] = -1.0e19;
        vars_upperbound[i] = 1.0e19;
    }
    
    for (unsigned int i = delta_start; i < a_start; i++) {
        vars_lowerbound[i] = -25 * M_PI / 180;
        vars_upperbound[i] = 25 * M_PI / 180;
    }
    
    for (unsigned int i = a_start; i < n_vars; i++) {
        vars_lowerbound[i] = -1.0;
        vars_upperbound[i] = 1.0;
    }
    
    
    // Lower and upper limits for the constraints
    // (all 0 besides initial state)
    Dvector constraints_lowerbound(n_constraints);
    Dvector constraints_upperbound(n_constraints);
    
    for (unsigned int i = 0; i < n_constraints; i++) {
        constraints_lowerbound[i] = 0;
        constraints_upperbound[i] = 0;
    }
    
    constraints_lowerbound[x_start] = state[0];
    constraints_lowerbound[y_start] = state[1];
    constraints_lowerbound[psi_start] = state[2];
    constraints_lowerbound[v_start] = state[3];
    constraints_lowerbound[cte_start] = state[4];
    constraints_lowerbound[psi_error_start] = state[5];
    
    constraints_upperbound[x_start] = state[0];
    constraints_upperbound[y_start] = state[1];
    constraints_upperbound[psi_start] = state[2];
    constraints_upperbound[v_start] = state[3];
    constraints_upperbound[cte_start] = state[4];
    constraints_upperbound[psi_error_start] = state[5];
    
    
    // object that computes objective and constraints
    FG_eval fg_eval(coeffs);
    
    // NOTE: You don't have to worry about these options
    // options for IPOPT solver
    std::string options;
    // Uncomment this if you'd like more print information
    options += "Integer print_level  0\n";
    // NOTE: Setting sparse to true allows the solver to take advantage
    //   of sparse routines, this makes the computation MUCH FASTER. If you can
    //   uncomment 1 of these and see if it makes a difference or not but if you
    //   uncomment both the computation time should go up in orders of magnitude.
    options += "Sparse  true        forward\n";
    options += "Sparse  true        reverse\n";
    // NOTE: Currently the solver has a maximum time limit of 0.5 seconds.
    // Change this as you see fit.
    options += "Numeric max_cpu_time          0.5\n";
    
    // place to return solution
    CppAD::ipopt::solve_result<Dvector> solution;
    
    // solve the problem
    // (see https://www.coin-or.org/CppAD/Doc/doxydoc/html/namespaceCppAD_1_1ipopt_a2b3f613ffa4b230d7c407d3c04c64dd4.html )
    CppAD::ipopt::solve<Dvector, FG_eval>(
                                          options, vars, vars_lowerbound, vars_upperbound, constraints_lowerbound,
                                          constraints_upperbound, fg_eval, solution);
    
    // Check some of the solution values
    ok &= solution.status == CppAD::ipopt::solve_result<Dvector>::success;
    
    // Cost
    auto cost = solution.obj_value;
    std::cout << "Cost " << cost << std::endl;
    
    
    // save control values for first timestep in vector to be returned
    std::vector<double> result;
    result.push_back(solution.x[delta_start]);
    result.push_back(solution.x[a_start]);
    
    // save planned x, y coordinates for all planned timesteps in vector to be returned
    for (unsigned int t = 1; t < N; t++) {
        result.push_back(solution.x[x_start + t]);
        result.push_back(solution.x[y_start + t]);
    }
    
    return result;
}

