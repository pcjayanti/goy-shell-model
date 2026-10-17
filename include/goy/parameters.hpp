#pragma once

#include <cstdint>
#include <string>
#include <vector>

namespace goy {

/// Every quantity that defines a run.
struct Parameters {
    // --- discretisation -----------------------------------------------
    int    n_shells   = 29;
    double lambda     = 2.0;              ///< geometric ratio k_{n+1}/k_n
    double k0         = 0.0009765625;     ///< 2^-10, smallest wavenumber
    double c0         = 0.0009765625;     ///< coupling prefactor

    // --- dissipation ---------------------------------------------------
    double nu         = 1.0e-10;          ///< kinematic viscosity
    double nu_ekman   = 0.009;            ///< Ekman friction coefficient
    int    ekman_shells = 6;              ///< friction acts on shells [0, this)

    // --- forcing --------------------------------------------------------
    double force_amplitude = 0.002;
    int    force_shell_a   = 19;          ///< 0-based shell index
    int    force_shell_b   = 20;

    // --- time stepping --------------------------------------------------
    double   dt = 1.0e-4;
    long long total_steps = 480000000;    ///< 48000 simulated seconds

    /// Simulated times (seconds) at which to emit a spectrum.
    std::vector<double> checkpoint_times = {
        500, 1000, 1500, 2000, 2500, 3000, 4000, 5000,
        6000, 8000, 10000, 12000, 24000, 48000};

    // --- reproducibility --------------------------------------------------
    std::uint64_t seed = 20160513;        ///< thesis submission date
    double initial_scale = 1.0e-3;

    // --- bookkeeping -------------------------------------------------------
    std::string run_name = "default";
    std::string output_dir = "runs";
    std::string integrator = "rk4";

    /// Total step count implied by the last checkpoint.
    [[nodiscard]] long long steps_for(double seconds) const {
        return static_cast<long long>(seconds / dt + 0.5);
    }
};

/// Thesis section 3.3.1 / Table 3.1: forcing at shells 10 and 11 (1-based),
/// Ekman friction on the first 10 shells. Reproduces Figs. 3.1-3.8.
Parameters preset_ohkitani_yamada();

/// Thesis section 3.3.2: forcing moved to shells 19 and 20 (1-based), Ekman
/// friction confined to the first 5 shells and strengthened by 1e4. This is
/// the configuration actually printed in Appendix A. Reproduces Figs. 3.16-3.17.
Parameters preset_sharp_inverse_cascade();

/// Inviscid, unforced. Energy and enstrophy are exact invariants of the
/// coupling coefficients, so any drift is purely numerical. Used by the
/// conservation test to decide whether the blow-up reported in the thesis
/// is physical or an artefact of fixed-step RK4 over 1e8 steps.
Parameters preset_conservation_test();

}  // namespace goy
