#pragma once

#include <fstream>
#include <string>
#include <vector>

#include "goy/parameters.hpp"

namespace goy {

/// Writes a run's output into its own directory:
///
///   <output_dir>/<run_name>/spectra.csv    time,shell,k,E,E_norm,k_norm,u
///   <output_dir>/<run_name>/scalars.csv    time,energy,enstrophy,eta,k_d
///   <output_dir>/<run_name>/parameters.txt every value that defined the run
///
///  1. Raw values are written, not log10 of them. Logs are a presentation
///     choice and belong in the analysis layer.
///  2. One directory per run, with the parameters recorded alongside the
///     data.
///
class Recorder {
public:
    explicit Recorder(const Parameters& p);

    void write_spectrum(double time,
                        const std::vector<double>& k,
                        const std::vector<double>& E,
                        const std::vector<double>& E_norm,
                        double k_d,
                        const std::vector<double>& u);

    void write_scalars(double time, double energy, double enstrophy,
                       double eta, double k_d);

    [[nodiscard]] const std::string& directory() const { return dir_; }

private:
    std::string   dir_;
    std::ofstream spectra_;
    std::ofstream scalars_;
};

}  // namespace goy
