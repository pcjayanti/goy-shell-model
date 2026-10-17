#include "goy/parameters.hpp"

namespace goy {

Parameters preset_ohkitani_yamada() {
    Parameters p;
    p.run_name      = "ohkitani_yamada";
    p.nu            = 1.0e-10;
    p.nu_ekman      = 9.0e-7;   // the pre-1e4 value used in section 3.3.1
    p.ekman_shells  = 10;
    p.force_shell_a = 9;        // shells 10 and 11, 1-based
    p.force_shell_b = 10;
    return p;
}

Parameters preset_sharp_inverse_cascade() {
    Parameters p;
    p.run_name      = "sharp_inverse_cascade";
    p.nu            = 1.0e-10;
    p.nu_ekman      = 0.009;    // raised by 1e4, per section 3.3.2
    p.ekman_shells  = 6;
    p.force_shell_a = 19;       // shells 20 and 21, 1-based
    p.force_shell_b = 20;
    return p;
}

Parameters preset_conservation_test() {
    Parameters p;
    p.run_name        = "conservation_test";
    p.nu              = 0.0;
    p.nu_ekman        = 0.0;
    p.ekman_shells    = 0;
    p.force_amplitude = 0.0;
    p.dt              = 1.0e-4;
    p.total_steps     = 1000000;
    p.checkpoint_times = {10.0, 50.0, 100.0};
    return p;
}

}  // namespace goy
