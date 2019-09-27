#include "MPC.h"
#include <cppad/cppad.hpp>
#include <cppad/ipopt/solve.hpp>
#include <iostream>
#include <string>
#include <vector>
#include "Eigen-3.3/Eigen/Core"

using CppAD::AD;
using Eigen::VectorXd;


const size_t N = 10;
// Index where the first element was stored in the vector
const int Index_px_start = 0;
const int Index_py_start = Index_px_start + N;
const int Index_psi_start = Index_py_start + N; //orientation
const int Index_v_start = Index_psi_start + N;
// cross-track error. The difference between the trajectory defined by 
// the waypoints and the current vehicle position y in the coordinate space of the vehicle.
const int Index_cte_start = Index_v_start + N; 
const int Index_epsi_start = Index_cte_start + N; // orientation error
const int Index_delta_start = Index_epsi_start + N; // steering angle
const int Index_a_start = Index_delta_start + N - 1;


// weights for cost function
const double W_cte = 3.0;
const double W_epsi = 1000.0;
const double W_v = 1.0;
const double W_delta = 5000.0;
const double W_a = 1.0;
const double W_diff_delta = 15000.0; // weight cost for high difference between consecutive steering actuations
const double W_diff_a = 15.0; // weight cost for high difference between consecutive acceleration actuations

const double ref_cte = 0;
const double ref_epsi = 0;
const double ref_v = 60;

class FG_eval {
    public:
    double Lf = 2.67;
    double dt = 0.1; 
    int time_step = 10;
    // Fitted polynomial coefficients
    VectorXd coeffs;
    FG_eval(VectorXd coeffs) { this->coeffs = coeffs; }

    typedef CPPAD_TESTVECTOR(AD<double>) ADvector;

    void operator()(ADvector& fg, const ADvector& vars) {
        /**
        * TODO: implement MPC
        * `fg` is a vector of the cost constraints, `vars` is a vector of variable 
        *   values (state & actuators)
        * NOTE: You'll probably go back and forth between this function and
        *   the Solver function below.
        */
        // fg[0] cost function value(fx), fg[>0] variables(gx)

        fg[0] = 0.0;
        for (int i = 0; i < time_step; i++) {
            fg[0] += W_cte * CppAD::pow(vars[Index_cte_start + i], 2);
            fg[0] += W_epsi * CppAD::pow(vars[Index_epsi_start + i], 2);
            fg[0] += W_v * CppAD::pow(vars[Index_v_start + i] - ref_v, 2);
        }

        for (int i = 0; i < time_step - 1; i++) {
            fg[0] += W_delta * CppAD::pow(vars[Index_delta_start + i], 2);
            fg[0] += W_a * CppAD::pow(vars[Index_a_start + i], 2);
        }

        for (int i = 0; i < time_step - 2; i++) {
            fg[0] += W_diff_delta * CppAD::pow(vars[Index_delta_start + i + 1] - vars[Index_delta_start + i], 2);
            fg[0] += W_diff_a * CppAD::pow(vars[Index_a_start + i + 1] - vars[Index_a_start + i], 2);
        }

        fg[1 + Index_px_start] = vars[Index_px_start];
        fg[1 + Index_py_start] = vars[Index_py_start];
        fg[1 + Index_psi_start] = vars[Index_psi_start];
        fg[1 + Index_v_start] = vars[Index_v_start];
        fg[1 + Index_cte_start] = vars[Index_cte_start];
        fg[1 + Index_epsi_start] = vars[Index_epsi_start];

        for (int i = 1; i < time_step; i++) {
            // state at time t+1
            AD<double> x1 = vars[Index_px_start + i];
            AD<double> y1 = vars[Index_py_start + i];
            AD<double> psi1 = vars[Index_psi_start + i];
            AD<double> v1 = vars[Index_v_start + i];
            AD<double> cte1 = vars[Index_cte_start + i];
            AD<double> epsi1 = vars[Index_epsi_start + i];
            // state at time t
            AD<double> x0 = vars[Index_px_start + i - 1];
            AD<double> y0 = vars[Index_py_start + i - 1];
            AD<double> psi0 = vars[Index_psi_start + i - 1];
            AD<double> v0 = vars[Index_v_start + i - 1];
            AD<double> cte0 = vars[Index_cte_start + i - 1];
            AD<double> epsi0 = vars[Index_epsi_start + i - 1];
            AD<double> delta0 = vars[Index_delta_start + i - 1];
            AD<double> a0 = vars[Index_a_start + i - 1];

            AD<double> f0 = coeffs[3]*x0*x0*x0 + coeffs[2]*x0*x0 + coeffs[1]*x0 + coeffs[0];
            AD<double> psides0 = CppAD::atan(coeffs[1] + (2*coeffs[2]*x0) + (3*coeffs[3]*x0*x0));
            // NOTE: Handle the latency here: since the dt is same as the latency 100ms, so
            fg[1 + Index_px_start + i] = x1 - (x0 + v0 * CppAD::cos(psi0) * dt);
            fg[1 + Index_py_start + i] = y1 - (y0 + v0 * CppAD::sin(psi0) * dt);
            fg[1 + Index_psi_start + i] = psi1 - (psi0 - v0 * delta0 / Lf * dt);
            fg[1 + Index_v_start + i] = v1 - (v0 + a0 * dt);
            fg[1 + Index_cte_start + i] = cte1 - ((f0 - y0) + (v0 * CppAD::sin(epsi0) * dt));
            fg[1 + Index_epsi_start + i] = epsi1 - ((psi0 - psides0) - v0 * delta0 / Lf * dt);
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

    /**
    * TODO: Set the number of model variables (includes both states and inputs).
    * For example: If the state is a 4 element vector, the actuators is a 2
    *    element vector and there are 10 timesteps. The number of variables is:
    *    4 * 10 + 2 * 9
    */
    // number of independent variables (domain dimension for f and g)
    size_t n_vars = num_of_states * N + num_of_actutions * (N - 1);
    /**
    * TODO: Set the number of constraints
    */
    // number of constraints (range dimension for g)
    size_t n_constraints = num_of_states * N;

    // Initial states
    double x = state[0];
    double y = state[1];
    double psi = state[2];
    double v = state[3];
    double cte = state[4];
    double epsi = state[5];

    // Initial value of the independent variables.
    // SHOULD BE 0 besides initial state.
    Dvector vars(n_vars);
    for (size_t i = 0; i < n_vars; ++i) {
        vars[i] = 0;
    }

    vars[Index_px_start] = x;
    vars[Index_py_start] = y;
    vars[Index_psi_start] = psi;
    vars[Index_v_start] = v;
    vars[Index_cte_start] = cte;
    vars[Index_epsi_start] = epsi;

    Dvector vars_lowerbound(n_vars);
    Dvector vars_upperbound(n_vars);
    /**
    * TODO: Set lower and upper limits for variables.
    */

    // Set all non-actuators upper and lowerlimits
    // to the max negative and positive values.
    for (size_t i = 0; i < Index_delta_start; i++) {
        vars_lowerbound[i] = -1.0e19;
        vars_upperbound[i] = 1.0e19;
    }

    // The upper and lower limits of delta are set to -25 and 25 degrees (values in radians).
    for (size_t i = Index_delta_start; i < Index_a_start; i++) {
        vars_lowerbound[i] = -0.436332;
        vars_upperbound[i] = 0.436332;
    }

    // Acceleration/decceleration upper and lower limits.
    for (size_t i = Index_a_start; i < n_vars; i++) {
        vars_lowerbound[i] = -1.0;
        vars_upperbound[i] =  1.0;
    }

    // Lower and upper limits for the constraints
    // Should be 0 besides initial state.
    Dvector constraints_lowerbound(n_constraints);
    Dvector constraints_upperbound(n_constraints);
    for (size_t i = 0; i < n_constraints; ++i) {
        constraints_lowerbound[i] = 0;
        constraints_upperbound[i] = 0;
    }  

    constraints_lowerbound[Index_px_start] = x;
    constraints_lowerbound[Index_py_start] = y;
    constraints_lowerbound[Index_psi_start] = psi;
    constraints_lowerbound[Index_v_start] = v;
    constraints_lowerbound[Index_cte_start] = cte;
    constraints_lowerbound[Index_epsi_start] = epsi;

    constraints_upperbound[Index_px_start] = x;
    constraints_upperbound[Index_py_start] = y;
    constraints_upperbound[Index_psi_start] = psi;
    constraints_upperbound[Index_v_start] = v;
    constraints_upperbound[Index_cte_start] = cte;
    constraints_upperbound[Index_epsi_start] = epsi;

    // object that computes objective and constraints
    FG_eval fg_eval(coeffs);

    // NOTE: You don't have to worry about these options
    // options for IPOPT solver
    std::string options;
    // Uncomment this if you'd like more print information
    options += "Integer print_level 0\n";
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
    CppAD::ipopt::solve<Dvector, FG_eval>(
        options, vars, vars_lowerbound, vars_upperbound, constraints_lowerbound,
        constraints_upperbound, fg_eval, solution);

    // Check some of the solution values
    ok &= solution.status == CppAD::ipopt::solve_result<Dvector>::success;

    // Cost
    auto cost = solution.obj_value;
    std::cout << "Cost " << cost << std::endl;

    /**
    * TODO: Return the first actuator values. The variables can be accessed with
    *   `solution.x[i]`.
    *
    * {...} is shorthand for creating a vector, so auto x1 = {1.0,2.0}
    *   creates a 2 element double vector.
    */
    double steering = solution.x[Index_delta_start];
    double acceleration = solution.x[Index_a_start];
    std::vector<double> outputs = {steering, acceleration};

    // attach the predicted route to display
    for (size_t i = 0; i < N; i++) {
        outputs.push_back(solution.x[Index_px_start + i]);
        outputs.push_back(solution.x[Index_py_start + i]);
    }
    return outputs;
}