#include "goy/recorder.hpp"

#include <filesystem>
#include <iomanip>
#include <stdexcept>

namespace goy {

namespace fs = std::filesystem;

Recorder::Recorder(const Parameters& p) {
    const fs::path dir = fs::path(p.output_dir) / p.run_name;
    fs::create_directories(dir);
    dir_ = dir.string();

    spectra_.open(dir / "spectra.csv");
    scalars_.open(dir / "scalars.csv");
    if (!spectra_ || !scalars_) {
        throw std::runtime_error("cannot open output files in " + dir_);
    }

    spectra_ << "time,shell,k,E,E_norm,k_norm,u\n";
    scalars_ << "time,energy,enstrophy,eta,k_d\n";

    spectra_ << std::setprecision(12);
    scalars_ << std::setprecision(12);

    std::ofstream manifest(dir / "parameters.txt");
    manifest << std::setprecision(17)
             << "run_name        " << p.run_name        << '\n'
             << "n_shells        " << p.n_shells        << '\n'
             << "lambda          " << p.lambda          << '\n'
             << "k0              " << p.k0              << '\n'
             << "c0              " << p.c0              << '\n'
             << "nu              " << p.nu              << '\n'
             << "nu_ekman        " << p.nu_ekman        << '\n'
             << "ekman_shells    " << p.ekman_shells    << '\n'
             << "force_amplitude " << p.force_amplitude << '\n'
             << "force_shell_a   " << p.force_shell_a   << '\n'
             << "force_shell_b   " << p.force_shell_b   << '\n'
             << "dt              " << p.dt              << '\n'
             << "integrator      " << p.integrator      << '\n'
             << "total_steps     " << p.total_steps     << '\n'
             << "seed            " << p.seed            << '\n'
             << "initial_scale   " << p.initial_scale   << '\n';
}

void Recorder::write_spectrum(double time,
                              const std::vector<double>& k,
                              const std::vector<double>& E,
                              const std::vector<double>& E_norm,
                              double k_d,
                              const std::vector<double>& u) {
    for (std::size_t i = 0; i < k.size(); ++i) {
        spectra_ << time << ',' << i << ',' << k[i] << ',' << E[i] << ','
                 << E_norm[i] << ','
                 << (k_d > 0.0 ? k[i] / k_d : 0.0) << ','
                 << u[i] << '\n';
    }
    spectra_.flush();
}

void Recorder::write_scalars(double time, double energy, double enstrophy,
                             double eta, double k_d) {
    scalars_ << time << ',' << energy << ',' << enstrophy << ',' << eta << ','
             << k_d << '\n';
    scalars_.flush();
}

}  // namespace goy
