#include "goy/model.hpp"

#include <cmath>
#include <algorithm>
#include <stdexcept>

namespace goy {

Model::Model(const Parameters& p)
    : n_(static_cast<std::size_t>(p.n_shells)),
      nu_(p.nu),
      k_(static_cast<std::size_t>(p.n_shells)),
      k2_(static_cast<std::size_t>(p.n_shells)),
      inv_k_(static_cast<std::size_t>(p.n_shells)),
      a_(static_cast<std::size_t>(p.n_shells), 0.0),
      b_(static_cast<std::size_t>(p.n_shells), 0.0),
      c_(static_cast<std::size_t>(p.n_shells), 0.0),
      ekman_(static_cast<std::size_t>(p.n_shells), 0.0),
      forcing_(static_cast<std::size_t>(p.n_shells), 0.0) {
    if (n_ < 5) {
        throw std::invalid_argument("GOY model needs at least 5 shells");
    }

    // Wavenumbers grow geometrically: k_n = k0 * lambda^n.
    k_[0] = p.k0;
    for (std::size_t i = 1; i < n_; ++i) {
        k_[i] = k_[i - 1] * p.lambda;
    }
    for (std::size_t i = 0; i < n_; ++i) {
        k2_[i]    = k_[i] * k_[i];
        inv_k_[i] = 1.0 / k_[i];
    }

    // a_n = c0 * lambda^n.
    a_[0] = p.c0;
    for (std::size_t i = 1; i < n_; ++i) {
        a_[i] = a_[i - 1] * p.lambda;
    }

    // b and c follow from requiring exact conservation of energy and
    // enstrophy by the nonlinear terms.
    for (std::size_t i = 1; i < n_ - 1; ++i) {
        b_[i] = -a_[i - 1] * (k2_[i - 1] - k2_[i + 1]) / (k2_[i] - k2_[i + 1]);
    }
    for (std::size_t i = 2; i < n_; ++i) {
        c_[i] = a_[i - 2] * (k2_[i - 2] - k2_[i - 1]) / (k2_[i - 1] - k2_[i]);
    }

    // Boundary shells have no neighbours to couple to.
    a_[n_ - 2] = 0.0;
    a_[n_ - 1] = 0.0;
    b_[0]      = 0.0;
    b_[n_ - 1] = 0.0;
    c_[0]      = 0.0;
    c_[1]      = 0.0;

    // Ekman friction on the largest scales only.
    const std::size_t ekman_count =
    std::min(n_, static_cast<std::size_t>(std::max(0, p.ekman_shells)));
    for (std::size_t i = 0; i < ekman_count; ++i) {
        ekman_[i] = p.nu_ekman * inv_k_[i];
    }

    // Forcing at two shells.
    auto apply_force = [&](int shell) {
        if (shell >= 0 && static_cast<std::size_t>(shell) < n_) {
            forcing_[static_cast<std::size_t>(shell)] += p.force_amplitude;
        }
    };
    apply_force(p.force_shell_a);
    apply_force(p.force_shell_b);
}

void Model::rhs(const std::vector<double>& u, std::vector<double>& out) const {
    // Interior shells first: the branch-free bulk of the work. The original
    // code re-tested i==0, i==1, i==N-2, i==N-1 on every shell of every RK
    // stage; here the boundary coefficients are already zero, so the same
    // expression is valid everywhere and the branches disappear entirely.
    for (std::size_t i = 0; i < n_; ++i) {
        const double u_p1 = (i + 1 < n_) ? u[i + 1] : 0.0;
        const double u_p2 = (i + 2 < n_) ? u[i + 2] : 0.0;
        const double u_m1 = (i >= 1) ? u[i - 1] : 0.0;
        const double u_m2 = (i >= 2) ? u[i - 2] : 0.0;

        out[i] = a_[i] * u_p1 * u_p2
               + b_[i] * u_m1 * u_p1
               + c_[i] * u_m1 * u_m2
               - nu_ * k2_[i] * u[i]
               - ekman_[i] * u[i]
               + forcing_[i];
    }
}

double Model::energy(const std::vector<double>& u) const {
    // Kahan compensated summation to preserve 
    // small contributions across many decades
    double sum = 0.0, comp = 0.0;
    for (std::size_t i = 0; i < n_; ++i) {
        const double term = u[i] * u[i] - comp;
        const double next = sum + term;
        comp = (next - sum) - term;
        sum  = next;
    }
    return 0.5 * sum;
}

double Model::enstrophy(const std::vector<double>& u) const {
    double sum = 0.0, comp = 0.0;
    for (std::size_t i = 0; i < n_; ++i) {
        const double term = k2_[i] * u[i] * u[i] - comp;
        const double next = sum + term;
        comp = (next - sum) - term;
        sum  = next;
    }
    return 0.5 * sum;
}

double Model::enstrophy_dissipation(const std::vector<double>& spectrum) const {
    double sum = 0.0, comp = 0.0;
    for (std::size_t i = 0; i < n_; ++i) {
        const double term = 2.0 * nu_ * spectrum[i] * k2_[i] * k2_[i] - comp;
        const double next = sum + term;
        comp = (next - sum) - term;
        sum  = next;
    }
    return sum;
}

double Model::dissipation_wavenumber(double eta) const {
    if (eta <= 0.0) return 0.0;
    return std::pow(eta, 1.0 / 6.0) * std::pow(nu_, -0.5);
}

double Model::normalise_spectrum(double E, double eta) const {
    if (eta <= 0.0) return 0.0;
    return E * std::pow(nu_, -1.5) * std::pow(eta, -1.0 / 6.0);
}

}  // namespace goy
