#pragma once

#include <cstddef>
#include <vector>

namespace goy {

/// Fixed-step classical Runge-Kutta 4.
///
/// Templated on the system rather than on an abstract base class: the
/// rhs() call is inlined, so swapping integrators costs nothing at runtime.
/// Scratch buffers live in the object and are sized once, so a step performs
/// no allocation -- this matters when stepping 10^8 times.

template <typename System>
class RungeKutta4 {
public:
    explicit RungeKutta4(std::size_t n)
        : k1_(n), k2_(n), k3_(n), k4_(n), tmp_(n) {}

    static constexpr const char* name() { return "rk4"; }
    static constexpr int order() { return 4; }

    void step(const System& sys, std::vector<double>& u, double dt) {
        const std::size_t n = u.size();

        sys.rhs(u, k1_);
        for (std::size_t i = 0; i < n; ++i) tmp_[i] = u[i] + 0.5 * dt * k1_[i];

        sys.rhs(tmp_, k2_);
        for (std::size_t i = 0; i < n; ++i) tmp_[i] = u[i] + 0.5 * dt * k2_[i];

        sys.rhs(tmp_, k3_);
        for (std::size_t i = 0; i < n; ++i) tmp_[i] = u[i] + dt * k3_[i];

        sys.rhs(tmp_, k4_);
        const double sixth = dt / 6.0;
        for (std::size_t i = 0; i < n; ++i) {
            u[i] += sixth * (k1_[i] + 2.0 * (k2_[i] + k3_[i]) + k4_[i]);
        }
    }

private:
    std::vector<double> k1_, k2_, k3_, k4_, tmp_;
};

/// Second-order Heun (explicit trapezoidal). Present so that a result can be
/// cross-checked at a different order: if a feature of the solution survives
/// a change of integrator, it is far more likely to be physical.
template <typename System>
class Heun {
public:
    explicit Heun(std::size_t n) : k1_(n), k2_(n), tmp_(n) {}

    static constexpr const char* name() { return "heun"; }
    static constexpr int order() { return 2; }

    void step(const System& sys, std::vector<double>& u, double dt) {
        const std::size_t n = u.size();

        sys.rhs(u, k1_);
        for (std::size_t i = 0; i < n; ++i) tmp_[i] = u[i] + dt * k1_[i];

        sys.rhs(tmp_, k2_);
        const double half = 0.5 * dt;
        for (std::size_t i = 0; i < n; ++i) u[i] += half * (k1_[i] + k2_[i]);
    }

private:
    std::vector<double> k1_, k2_, tmp_;
};

}  // namespace goy
