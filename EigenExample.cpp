#define _USE_MATH_DEFINES

#include <iostream>
#include <eigen3/Eigen/Dense>
#include <cmath>

// Function to normalize y1 and y2 to ensure the constraint y1^2 + y2^2 = 1 is satisfied
void normalize(Eigen::Vector2d& y) {
    double r = std::sqrt(y(0) * y(0) + y(1) * y(1));
    y(0) /= r;
    y(1) /= r;
}

// Function to compute the derivatives (dy/dt) based on the equation dy1/dt = y2
Eigen::Vector2d derivatives(const Eigen::Vector2d& y) {
    Eigen::Vector2d dydt;
    dydt(0) = y(1); // dy1/dt = y2
    dydt(1) = 0.0;  // No differential equation for y2
    return dydt;
}

// Runge-Kutta 4th order method for one time step
Eigen::Vector2d RK4_step(const Eigen::Vector2d& y, double dt) {
    Eigen::Vector2d k1, k2, k3, k4;

    // Calculate k1
    k1 = derivatives(y);

    // Calculate k2
    k2 = derivatives(y + 0.5 * dt * k1);

    // Calculate k3
    k3 = derivatives(y + 0.5 * dt * k2);

    // Calculate k4
    k4 = derivatives(y + dt * k3);

    // Update the solution
    Eigen::Vector2d y_new = y + (dt / 6.0) * (k1 + 2 * k2 + 2 * k3 + k4);

    // Normalize the result to ensure the constraint is satisfied
    normalize(y_new);

    return y_new;
}

int main() {
    // Initial conditions
    Eigen::Vector2d y;
    y(0) = 1.0 / std::sqrt(2); // y1(0) = 1/sqrt(2)
    y(1) = 1.0 / std::sqrt(2); // y2(0) = 1/sqrt(2)

    double y1_true; 
    double y2_true; 

    // Time step and number of steps
    double dt = 0.01;
    int num_steps = 100; // For example, 100 time steps

    // Time loop for integrating the equations
    for (int step = 0; step < num_steps; ++step) {
        y1_true = std::cos(step * dt + M_PI/4);
        y2_true = std::sin(step * dt + M_PI/4);

        // Print the current state (y1 and y2)
        std::cout << "Step " << step << ": y1 = " << y(0) << ", y2 = " << y(1) << ": y1_true = " << y1_true << ", y2_true = " << y2_true << std::endl;

        // Take a Runge-Kutta 4th order step
        y = RK4_step(y, dt);
    }

    return 0;
}
