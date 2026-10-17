// Self-contained test suite

#include <cmath>
#include <cstdio>
#include <iostream>
#include <string>
#include <vector>

#include "goy/integrator.hpp"
#include "goy/model.hpp"
#include "goy/parameters.hpp"
#include "goy/simulation.hpp"

namespace {

int g_failures = 0;
int g_checks = 0;

void check(bool condition, const std::string& what) {
    ++g_checks;
    if (condition) {
        std::cout << "  PASS  " << what << "\n";
    } else {
        std::cout << "  FAIL  " << what << "\n";
        ++g_failures;
    }
}

void check_near(double value, double expected, double tolerance,
                const std::string& what) {
    const bool ok = std::abs(value - expected) <= tolerance;
    ++g_checks;
    if (ok) {
        std::cout << "  PASS  " << what << "  (" << value << ")\n";
    } else {
        std::cout << "  FAIL  " << what << "  got " << value << ", expected "
                  << expected << " +/- " << tolerance << "\n";
        ++g_failures;
    }
}

// ---------------------------------------------------------------------------

/// The coupling coefficients exist to make the nonlinear terms conserve both
/// energy and enstrophy exactly. We can verify that directly: with dissipation
/// and forcing switched off, d/dt of each invariant must vanish to round-off.
void test_nonlinear_terms_conserve_invariants() {
    std::cout << "nonlinear terms conserve energy and enstrophy\n";

    goy::Parameters p = goy::preset_conservation_test();
    goy::Model model(p);

    std::vector<double> u = goy::make_initial_condition(p);
    std::vector<double> du(u.size());
    model.rhs(u, du);

    // dE/dt = sum(u_n du_n),  dOmega/dt = sum(k_n^2 u_n du_n)
    double dE = 0.0, dOmega = 0.0, scale = 0.0;
    const std::vector<double>& k = model.k();
    for (std::size_t i = 0; i < u.size(); ++i) {
        dE     += u[i] * du[i];
        dOmega += k[i] * k[i] * u[i] * du[i];
        scale  += std::abs(u[i] * du[i]);
    }

    check(scale > 0.0, "the right-hand side is not trivially zero");
    check_near(dE / scale, 0.0, 1e-12, "relative dE/dt vanishes");

    double omega_scale = 0.0;
    for (std::size_t i = 0; i < u.size(); ++i) {
        omega_scale += std::abs(k[i] * k[i] * u[i] * du[i]);
    }
    check_near(dOmega / omega_scale, 0.0, 1e-12,
               "relative d(enstrophy)/dt vanishes");
}

/// The question the thesis left open: is the unbounded energy growth physical,
/// or is it numerical drift in fixed-step RK4 over a very long integration?
/// We integrate the inviscid, unforced system and measure how much the invariants
/// move. If they hold steady, drift is not the explanation.
void test_invariants_survive_long_integration() {
    std::cout << "invariants survive long inviscid integration\n";

    goy::Parameters p = goy::preset_conservation_test();
    p.checkpoint_times = {100.0};
    goy::Model model(p);

    goy::RunResult r = goy::run<goy::RungeKutta4>(p, model, nullptr);

    check(!r.diverged, "inviscid run does not diverge");
    check(r.energy_drift() < 1e-8,
          "energy drift below 1e-8 over " + std::to_string(r.steps) + " steps");
    check(r.enstrophy_drift() < 1e-6,
          "enstrophy drift below 1e-6 over " + std::to_string(r.steps) +
              " steps");
    std::printf("        energy drift    %.3e\n", r.energy_drift());
    std::printf("        enstrophy drift %.3e\n", r.enstrophy_drift());
}

/// Same seed must give the same answer. This is what the 2016 srand() bug
/// destroyed, and it is the property every other test depends on.
void test_runs_are_reproducible() {
    std::cout << "runs are reproducible from the seed\n";

    goy::Parameters p = goy::preset_conservation_test();
    p.checkpoint_times = {1.0};

    goy::Model m1(p), m2(p);
    goy::RunResult a = goy::run<goy::RungeKutta4>(p, m1, nullptr);
    goy::RunResult b = goy::run<goy::RungeKutta4>(p, m2, nullptr);

    bool identical = true;
    for (std::size_t i = 0; i < a.u_final.size(); ++i) {
        if (a.u_final[i] != b.u_final[i]) identical = false;
    }
    check(identical, "two runs with the same seed agree bit for bit");

    goy::Parameters q = p;
    q.seed = p.seed + 1;
    goy::Model m3(q);
    goy::RunResult c = goy::run<goy::RungeKutta4>(q, m3, nullptr);
    check(c.u_final[0] != a.u_final[0], "a different seed gives a different run");
}

/// The initial condition must actually vary across shells. Under the 2016 bug,
/// every shell received the same magnitude; this test would have caught it.
void test_initial_condition_varies_across_shells() {
    std::cout << "initial condition is random across shells\n";

    goy::Parameters p;
    const std::vector<double> u = goy::make_initial_condition(p);

    int distinct = 0;
    for (std::size_t i = 1; i < u.size(); ++i) {
        if (std::abs(std::abs(u[i]) - std::abs(u[0])) > 1e-15) ++distinct;
    }
    check(distinct == static_cast<int>(u.size()) - 1,
      "every shell magnitude differs from shell 0");

    bool in_range = true;
    for (double value : u) {
        if (std::abs(value) > p.initial_scale) in_range = false;
    }
    check(in_range, "all magnitudes are within the requested scale");
}

/// Halving the step size should reduce RK4's error by about 2^4.
/// Confirms the integrator really is fourth order and not silently degraded.
void test_rk4_converges_at_fourth_order() {
    std::cout << "RK4 shows fourth-order convergence\n";

    goy::Parameters p = goy::preset_conservation_test();
    p.checkpoint_times = {0.5};

    auto final_state = [&](double dt) {
        goy::Parameters q = p;
        q.dt = dt;
        goy::Model model(q);
        return goy::run<goy::RungeKutta4>(q, model, nullptr).u_final;
    };

    const std::vector<double> coarse = final_state(4.0e-4);
    const std::vector<double> medium = final_state(2.0e-4);
    const std::vector<double> fine   = final_state(1.0e-4); // treat the "fine" result as the reference

    double e_coarse = 0.0, e_medium = 0.0;
    for (std::size_t i = 0; i < fine.size(); ++i) {
        e_coarse = std::max(e_coarse, std::abs(coarse[i] - fine[i]));
        e_medium = std::max(e_medium, std::abs(medium[i] - fine[i]));
    }

    const double ratio = e_medium > 0.0 ? e_coarse / e_medium : 0.0;
    std::printf("        error ratio %.2f (expect ~16 for 4th order)\n", ratio);
    check(ratio > 8.0, "error shrinks by more than 8x when dt is halved");
    // we expect the ratio to be 16 for 4th order, and 4 for 2nd order.
    // But to account for the fact that we chose "fine" as the reference, 
    // we do not choose a threshold that is too close to 16.
}

/// Heun is second order, so it must be measurably worse than RK4 at the same
/// step size. We can compare them in terms of the energy drift in the inviscid case. 
/// If they agreed, one of them would be wrong.
void test_integrators_are_interchangeable_but_differ() {
    std::cout << "integrators are interchangeable via the template\n";

    goy::Parameters p = goy::preset_conservation_test();
    p.checkpoint_times = {1.0};

    goy::Model m1(p), m2(p);
    goy::RunResult rk4  = goy::run<goy::RungeKutta4>(p, m1, nullptr);
    goy::RunResult heun = goy::run<goy::Heun>(p, m2, nullptr);

    check(!rk4.diverged && !heun.diverged, "both integrators complete");
    check(rk4.energy_drift() < heun.energy_drift(),
          "RK4 conserves energy better than Heun");
    std::printf("        rk4 %.3e   heun %.3e\n", rk4.energy_drift(),
                heun.energy_drift());
}


}  // namespace

int main() {
    std::cout << "GOY shell model -- test suite\n\n";

    test_nonlinear_terms_conserve_invariants();
    test_invariants_survive_long_integration();
    test_runs_are_reproducible();
    test_initial_condition_varies_across_shells();
    test_rk4_converges_at_fourth_order();
    test_integrators_are_interchangeable_but_differ();

    std::cout << "\n" << (g_checks - g_failures) << "/" << g_checks
              << " checks passed\n";
    return g_failures == 0 ? 0 : 1;
}
