#include <algorithm>
#include <atomic>
#include <chrono>
#include <cmath>
#include <cstdlib>
#include <iostream>
#include <mutex>
#include <sstream>
#include <string>
#include <thread>
#include <vector>

#include "goy/integrator.hpp"
#include "goy/model.hpp"
#include "goy/parameters.hpp"
#include "goy/recorder.hpp"
#include "goy/simulation.hpp"

namespace {

std::mutex g_console;

void usage() {
    std::cout << R"(goy -- GOY shell model for 2D turbulence

Usage: goy [options]

  --preset NAME        oy | sharp | conservation      (default: sharp)
  --integrator NAME    rk4 | heun                     (default: rk4)
  --time SECONDS       final simulated time           (default: preset)
  --nu VALUE           kinematic viscosity
  --nu-ekman VALUE     Ekman friction coefficient
  --ekman-shells N     number of Ekman-active shells
  --force-a N          first forced shell  (0-based)
  --force-b N          second forced shell (0-based)
  --dt VALUE           time step
  --seed N             RNG seed for the initial condition
  --shells N           number of shells                (default: 29)
  --name NAME          run name (output subdirectory)
  --out DIR            output directory                (default: runs)
  --sweep-nu A,B,C     run these viscosities in parallel
  --threads N          worker threads for a sweep      (default: hardware)
  --help

Output lands in <out>/<name>/ as spectra.csv, scalars.csv, parameters.txt.
)";
}

std::vector<double> parse_list(const std::string& text) {
    std::vector<double> values;
    std::stringstream ss(text);
    std::string item;
    while (std::getline(ss, item, ',')) {
        if (!item.empty()) values.push_back(std::stod(item));
    }
    return values;
}

/// Run one configuration and report timing. Returns false if it diverged.
bool execute(const goy::Parameters& p, const std::string& integrator) {
    goy::Model model(p);
    goy::Recorder recorder(p);

    const auto t0 = std::chrono::steady_clock::now();
    goy::RunResult result =
        (integrator == "heun")
            ? goy::run<goy::Heun>(p, model, &recorder)
            : goy::run<goy::RungeKutta4>(p, model, &recorder);
    const auto t1 = std::chrono::steady_clock::now();

    const double seconds =
        std::chrono::duration<double>(t1 - t0).count();
    const double rate =
        seconds > 0.0 ? static_cast<double>(result.steps) / seconds : 0.0;

    std::lock_guard<std::mutex> lock(g_console);
    std::cout << "[" << p.run_name << "] " << result.steps << " steps in "
              << seconds << " s  (" << static_cast<long long>(rate)
              << " steps/s)\n"
              << "           energy " << result.energy_initial << " -> "
              << result.energy_final << "\n"
              << "           output " << recorder.directory() << "\n";
    if (result.diverged) {
        std::cout << "           DIVERGED at step " << result.steps << "\n";
    }
    return !result.diverged;
}

}  // namespace

int main(int argc, char** argv) {
    std::string preset = "sharp";
    std::string integrator = "rk4";
    std::string sweep_nu;
    unsigned threads = 0;
    double final_time = -1.0;

    // Deferred overrides: the preset is chosen first, then these applied.
    std::vector<std::pair<std::string, std::string>> overrides;

    for (int i = 1; i < argc; ++i) {
        const std::string arg = argv[i];
        auto next = [&]() -> std::string {
            if (i + 1 >= argc) {
                std::cerr << "missing value for " << arg << "\n";
                std::exit(2);
            }
            return argv[++i];
        };

        if (arg == "--help" || arg == "-h") { usage(); return 0; }
        else if (arg == "--preset")       preset = next();
        else if (arg == "--integrator")   integrator = next();
        else if (arg == "--sweep-nu")     sweep_nu = next();
        else if (arg == "--threads")      threads = static_cast<unsigned>(std::stoul(next()));
        else if (arg == "--time")         final_time = std::stod(next());
        else if (arg.rfind("--", 0) == 0) overrides.emplace_back(arg, next());
        else { std::cerr << "unknown argument: " << arg << "\n"; return 2; }
    }

    goy::Parameters base;
    if (preset == "oy")                base = goy::preset_ohkitani_yamada();
    else if (preset == "sharp")        base = goy::preset_sharp_inverse_cascade();
    else if (preset == "conservation") base = goy::preset_conservation_test();
    else { std::cerr << "unknown preset: " << preset << "\n"; return 2; }

    for (const auto& [key, value] : overrides) {
        if      (key == "--nu")            base.nu = std::stod(value);
        else if (key == "--nu-ekman")      base.nu_ekman = std::stod(value);
        else if (key == "--ekman-shells")  base.ekman_shells = std::stoi(value);
        else if (key == "--force-a")       base.force_shell_a = std::stoi(value);
        else if (key == "--force-b")       base.force_shell_b = std::stoi(value);
        else if (key == "--dt")            base.dt = std::stod(value);
        else if (key == "--seed")          base.seed = std::stoull(value);
        else if (key == "--shells")        base.n_shells = std::stoi(value);
        else if (key == "--name")          base.run_name = value;
        else if (key == "--out")           base.output_dir = value;
        else { std::cerr << "unknown option: " << key << "\n"; return 2; }
    }
    base.integrator = integrator;

    if (final_time > 0.0) {
        std::vector<double> kept;
        for (double t : base.checkpoint_times) {
            if (t <= final_time) kept.push_back(t);
        }
        if (kept.empty() || kept.back() < final_time) kept.push_back(final_time);
        base.checkpoint_times = kept;
        base.total_steps = base.steps_for(final_time);
    }

    try {
        if (sweep_nu.empty()) {
            return execute(base, integrator) ? 0 : 1;
        }

        const std::vector<double> values = parse_list(sweep_nu);
        std::vector<goy::Parameters> jobs;
        jobs.reserve(values.size());
        for (double nu : values) {
            goy::Parameters p = base;
            p.nu = nu;
            std::ostringstream name;
            name << base.run_name << "_nu" << nu;
            p.run_name = name.str();
            jobs.push_back(p);
        }

        if (threads == 0) {
            threads = std::max(1u, std::thread::hardware_concurrency());
        }
        threads = std::min<unsigned>(threads, static_cast<unsigned>(jobs.size()));
        std::cout << "sweeping " << jobs.size() << " configurations on "
                  << threads << " threads\n";

        std::atomic<std::size_t> cursor{0};
        std::atomic<bool> all_ok{true};
        std::vector<std::thread> pool;
        for (unsigned t = 0; t < threads; ++t) {
            pool.emplace_back([&] {
                for (;;) {
                    const std::size_t index = cursor.fetch_add(1);
                    if (index >= jobs.size()) return;
                    if (!execute(jobs[index], integrator)) all_ok = false;
                }
            });
        }
        for (auto& worker : pool) worker.join();
        return all_ok ? 0 : 1;

    } catch (const std::exception& error) {
        std::cerr << "error: " << error.what() << "\n";
        return 1;
    }
}
