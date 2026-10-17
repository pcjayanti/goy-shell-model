#pragma once

#include <cstddef>
#include <vector>

#include "goy/parameters.hpp"

namespace goy {

/// The GOY shell model for 2D turbulence.
///
///   du_n/dt =  a_n u_{n+1} u_{n+2}
///            + b_n u_{n-1} u_{n+1}
///            + c_n u_{n-1} u_{n-2}
///            - nu k_n^2 u_n
///            - nu_ekman u_n / k_n   (large scales only)
///            + f_n
///
/// The coefficients a, b, c are constructed so that the nonlinear terms
/// conserve both energy sum(u_n^2) and enstrophy sum(k_n^2 u_n^2) exactly.
/// That property is what the conservation test exploits.
class Model {
public:
    explicit Model(const Parameters& p);

    /// Evaluate du/dt into `out`. Hot path: no allocation, no pow().
    void rhs(const std::vector<double>& u, std::vector<double>& out) const;

    /// sum(u_n^2) / 2 -- the first inviscid invariant in 2D.
    [[nodiscard]] double energy(const std::vector<double>& u) const;

    /// sum(k_n^2 u_n^2) / 2 -- the second inviscid invariant in 2D.
    [[nodiscard]] double enstrophy(const std::vector<double>& u) const;

    /// Enstrophy dissipation rate, 2 nu sum(k_n^4 u_n^2 / (2 k_n)).
    [[nodiscard]] double enstrophy_dissipation(
        const std::vector<double>& mean_spectrum) const;

    /// Kolmogorov (dissipation) wavenumber, eta^(1/6) nu^(-1/2).
    [[nodiscard]] double dissipation_wavenumber(double eta) const;

    /// Normalised spectrum E * nu^(-3/2) * eta^(-1/6), per Eq. (2.11).
    [[nodiscard]] double normalise_spectrum(double E, double eta) const;

    [[nodiscard]] const std::vector<double>& k() const { return k_; }
    [[nodiscard]] std::size_t size() const { return n_; }

private:
    std::size_t n_;
    double nu_;

    std::vector<double> k_;        ///< wavenumbers
    std::vector<double> k2_;       ///< k^2, precomputed
    std::vector<double> inv_k_;    ///< 1/k, precomputed
    std::vector<double> a_, b_, c_;///< coupling coefficients
    std::vector<double> ekman_;    ///< nu_ekman / k_n, or 0
    std::vector<double> forcing_;  ///< f_n
};

}  // namespace goy
