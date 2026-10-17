#pragma once

#include <algorithm>
#include <random>
#include <vector>

#include "goy/model.hpp"
#include "goy/parameters.hpp"
#include "goy/recorder.hpp"

namespace goy {

/// Running mean of u_n^2 over the whole integration, accumulated with Kahan
/// compensation.

class MeanSquare {
public:
    explicit MeanSquare(std::size_t n) : sum_(n, 0.0), comp_(n, 0.0) {}

    void accumulate(const std::vector<double>& u) {
        for (std::size_t i = 0; i < u.size(); ++i) {
            const double term = u[i] * u[i] - comp_[i];
            const double next = sum_[i] + term;
            comp_[i] = (next - sum_[i]) - term;
            sum_[i]  = next;
        }
        ++count_;
    }

    /// Time-averaged energy spectrum E(k_n) = <u_n^2> / (2 k_n).
    [[nodiscard]] std::vector<double> spectrum(
        const std::vector<double>& k) const {
        std::vector<double> E(sum_.size(), 0.0);
        if (count_ == 0) return E;
        for (std::size_t i = 0; i < sum_.size(); ++i) {
            E[i] = sum_[i] / (2.0 * static_cast<double>(count_) * k[i]);
        }
        return E;
    }

private:
    std::vector<double> sum_, comp_;
    long long count_ = 0;
};

/// Draw the initial condition.
///
/// The 2016 code called srand(time(NULL)) inside the per-shell random
/// function, so all 29 calls within the same wall-clock second reseeded to
/// the same value and returned the same magnitude. The initial state was not
/// random across shells: it was a constant vector with alternating signs.
/// Here a single generator is seeded once from Parameters::seed, which makes
/// runs both genuinely random across shells and exactly reproducible.
inline std::vector<double> make_initial_condition(const Parameters& p) {
    std::mt19937_64 rng(p.seed);
    std::uniform_real_distribution<double> dist(-1.0, 1.0);

    std::vector<double> u(static_cast<std::size_t>(p.n_shells));
    for (auto& value : u) {
        value = p.initial_scale * dist(rng);
    }
    return u;
}

/// Outcome of a run, returned so callers can assert on it in tests.
struct RunResult {
    std::vector<double> u_final;
    double energy_initial    = 0.0;
    double energy_final      = 0.0;
    double enstrophy_initial = 0.0;
    double enstrophy_final   = 0.0;
    long long steps          = 0;
    bool diverged            = false;

    [[nodiscard]] double energy_drift() const {
        return energy_initial > 0.0
                   ? std::abs(energy_final - energy_initial) / energy_initial
                   : 0.0;
    }
    [[nodiscard]] double enstrophy_drift() const {
        return enstrophy_initial > 0.0
                   ? std::abs(enstrophy_final - enstrophy_initial) /
                         enstrophy_initial
                   : 0.0;
    }
};

/// Integrate the model, emitting a spectrum at each requested time.
/// `recorder` may be null, in which case nothing is written (used by tests).
template <template <typename> class Stepper>
RunResult run(const Parameters& p, Model& model, Recorder* recorder) {
    const std::size_t n = static_cast<std::size_t>(p.n_shells);

    std::vector<double> u = make_initial_condition(p);
    Stepper<Model> stepper(n);
    MeanSquare mean(n);

    RunResult result;
    result.energy_initial    = model.energy(u);
    result.enstrophy_initial = model.enstrophy(u);

    // Checkpoints as step counts, sorted and de-duplicated.
    std::vector<long long> checkpoints;
    checkpoints.reserve(p.checkpoint_times.size());
    for (double t : p.checkpoint_times) checkpoints.push_back(p.steps_for(t));
    std::sort(checkpoints.begin(), checkpoints.end());
    checkpoints.erase(std::unique(checkpoints.begin(), checkpoints.end()),
                      checkpoints.end());

    const long long last =
        checkpoints.empty() ? p.total_steps : checkpoints.back();
    std::size_t next = 0;

    for (long long step = 1; step <= last; ++step) {
        stepper.step(model, u, p.dt);
        mean.accumulate(u);

        if (next < checkpoints.size() && step == checkpoints[next]) {
            const double time = static_cast<double>(step) * p.dt;
            const std::vector<double> E = mean.spectrum(model.k());
            const double eta = model.enstrophy_dissipation(E);
            const double k_d = model.dissipation_wavenumber(eta);

            if (recorder != nullptr) {
                std::vector<double> E_norm(n);
                for (std::size_t i = 0; i < n; ++i) {
                    E_norm[i] = model.normalise_spectrum(E[i], eta);
                }
                recorder->write_spectrum(time, model.k(), E, E_norm, k_d, u);
                recorder->write_scalars(time, model.energy(u),
                                        model.enstrophy(u), eta, k_d);
            }
            ++next;
        }

        // Fixed-step explicit schemes on this model can blow up outright.
        // Detect it rather than filling the output with inf and nan.
        if (!std::isfinite(u[0])) {
            result.diverged = true;
            result.steps = step;
            break;
        }
    }

    if (!result.diverged) result.steps = last;
    result.u_final         = u;
    result.energy_final    = model.energy(u);
    result.enstrophy_final = model.enstrophy(u);
    return result;
}

}  // namespace goy
