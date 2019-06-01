// Various time-stepping methods for solving coupled ODE systems

#ifndef TIMESTEPPING_HPP_
#define TIMESTEPPING_HPP_

#include <cmath>
#include <iostream>
#include <vector>
#include <utility>
#include <limits>
#include <algorithm>

namespace timestepping {

// ODE solver using some time-stepping scheme
// Returns (time values, state values)
// DerivFn(double t, vector<double> y) -> vector<double> dy
// Timestepper(double t, vector<double> y, DerivFn deriv) -> (new_t, new_y)
template<typename DerivFn, typename Timestepper>
std::vector< std::pair< double, std::vector<double> > > odesolve(
    DerivFn deriv, std::vector<double> y0, double t_final, Timestepper step) {

    std::vector< std::pair< double, std::vector<double> > >
        odesolution{{0, y0}};
    while(odesolution.back().first < t_final) {
        try {
            odesolution.push_back(step(
                odesolution.back().first, odesolution.back().second, deriv));
        } catch(const std::exception &e) {
            std::cout << e.what() << std::endl;
            break;
        }
    }
    return odesolution;
}

// Supporting function
// c1*v1 + c2*v2
std::vector<double> vlincombo(std::vector<double> v1, std::vector<double> v2,
    double c1=1, double c2=1) {
    std::vector<double> vsum;
    vsum.reserve(v1.size());
    for(unsigned i = 0; i < v1.size(); ++i) {
        vsum.push_back(c1*v1[i] + c2*v2[i]);
    }
    return vsum;
}

// 2nd-order Runge-Kutta scheme
class RK2 {
    private:
        double dt;
    public:
        RK2(double dt):dt(dt) {}

        template<typename DerivFn>
        std::pair< double, std::vector<double> > operator()(
            double t, std::vector<double> y, DerivFn deriv) {
            std::vector<double> dy = deriv(t, y);
            // Make a full step using half-step peek values
            dy = deriv(t + dt/2, vlincombo(dy, y, dt/2));
            return std::make_pair(t + dt, vlincombo(dy, y, dt));
        }
};

// 4th-order Runge-Kutta scheme
class RK4 {
    private:
        double dt;
    public:
        RK4(double dt):dt(dt) {}

        void set_dt(double dt) {
            this->dt = dt;
        }

        template<typename DerivFn>
        std::pair< double, std::vector<double> > operator()(
            double t, std::vector<double> y, DerivFn deriv) {
            // Estimate derivative as quadrature of four points:
            // dy1, dy2, dy3, dy4
            std::vector<double> dy1 = deriv(t, y);
            std::vector<double> dy2 = deriv(t + dt/2, vlincombo(dy1, y, dt/2));
            std::vector<double> dy3 = deriv(t + dt/2, vlincombo(dy2, y, dt/2));
            std::vector<double> dy4 = deriv(t + dt, vlincombo(dy3, y, dt));
            // Make a full step using the quadrature derivative estimate
            for(unsigned i = 0; i < y.size(); ++i) {
                y[i] += dt/6 * (dy1[i] + 2*(dy2[i] + dy3[i]) + dy4[i]);
            }
            return std::make_pair(t + dt, y);
        }
};

// Adaptive 4/5-th order Runge-Kutta scheme
class AdaptiveRK {
    private:
        double tol;
        double dt;
        double dt_shrink;   // Shrink factor on time step adjustment
        double dt_adjust_lim;   // Max factor of adjustment in a single iteration
        unsigned max_dt_adjusts;

        RK4 rk4stepper;
    public:
        AdaptiveRK(double tol, double dt=1e-3,
            double dt_shrink=0.9, double dt_adjust_lim=4,
            unsigned max_dt_adjusts=100)
            :tol(tol), dt(dt),
            dt_shrink(dt_shrink), dt_adjust_lim(dt_adjust_lim),
            max_dt_adjusts(max_dt_adjusts), rk4stepper(dt) {}

        template<typename DerivFn>
        std::pair< double, std::vector<double> > operator()(
            double t, std::vector<double> y, DerivFn deriv) {
            
            for(unsigned i = 0; i < max_dt_adjusts; ++i) {
                // Time after the step
                double t_new = t + dt;

                // Two half-steps
                rk4stepper.set_dt(dt/2);
                std::vector<double> y_small = rk4stepper(t + dt/2,
                    rk4stepper(t, y, deriv).second, deriv).second;

                // One full step
                rk4stepper.set_dt(dt);
                std::vector<double> y_big = rk4stepper(t, y, deriv).second;

                // Estimate maximum relative truncation error of any component
                double error_ratio = 0;
                for(unsigned cmp = 0; cmp < y_small.size(); ++cmp) {
                    // Buffer by the numeric double precision limit
                    double desired_err = tol *
                        (fabs(y_small[cmp]) + fabs(y_big[cmp]))/2
                        + std::numeric_limits<double>::epsilon();
                    error_ratio = std::max(error_ratio,
                        fabs(y_small[cmp] - y_big[cmp]) / desired_err);
                }

                // Estimate better time step
                // This persists for the next time step if no more adjustments
                // are needed
                dt = std::max(dt / dt_adjust_lim,
                    std::min(dt * dt_adjust_lim,
                        dt * dt_shrink * pow(error_ratio, -0.2)
                    ));

                // Check if the error is within the desired tolerance
                if(error_ratio < 1) {
                    // Return the calculation with the smaller time step
                    return std::make_pair(t_new, y_small);
                }
            }

            // Give up after too many adjustments
            throw std::runtime_error(
                "Maximum number of time step adjustments exceeded");
        }
};

}

#endif