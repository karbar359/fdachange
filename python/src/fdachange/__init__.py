"""fdachange: epidemic change-point tests for functional data.

Python twin of the fdachange R package (the reference implementation);
the two are held equal by a cross-language parity test suite.
"""

from ._weights import holder_weight, holder_weights, critical_value
from ._kernels import kernel_gaussian, kernel_min, kernel_2max, gram
from .wilcoxon import wilcoxon_statistic, wilcoxon_scan
from .lp import lp_statistic, lp_scan
from .rkhs import rkhs_statistic, rkhs_scan
from .far1 import far1_statistic, far1_scan, far1_betas
from .epidemic_test import epidemic_test, EpidemicTestResult
from .scenarios import (
    generate_s1, generate_s2, generate_s3, generate_s4, generate_s5,
    generate_f1, generate_f2, generate_f3, generate_f4, generate_f5,
    scenario_cdf,
)
from .projections import (
    project_l2, project_ito, project_wiener, project_c_max_var,
    project_pc1, project_pc_d, projection_fn, projection_multiplicity,
)
from .reconstruction import reconstruct
from .splits import scan_over_splits, make_run
from .nulls import (
    null_cache_dir, save_null_draws,
    wilcoxon_null_draws, wilcoxon_critical_value,
    lp_null_draws, lp_critical_value,
    rkhs_null_draws, rkhs_critical_value,
    far1_null_draws, far1_critical_value,
)
from .simulate import (
    simulate_wilcoxon_null, simulate_lp_null,
    simulate_rkhs_null, simulate_far1_null,
)
from .datasets import TelcoCluster, load_telco_far1

__version__ = "1.0.0"

__all__ = [
    "holder_weight", "holder_weights", "critical_value",
    "kernel_gaussian", "kernel_min", "kernel_2max", "gram",
    "wilcoxon_statistic", "wilcoxon_scan",
    "lp_statistic", "lp_scan",
    "rkhs_statistic", "rkhs_scan",
    "far1_statistic", "far1_scan", "far1_betas",
    "epidemic_test", "EpidemicTestResult",
    "generate_s1", "generate_s2", "generate_s3", "generate_s4",
    "generate_s5", "generate_f1", "generate_f2", "generate_f3",
    "generate_f4", "generate_f5", "scenario_cdf",
    "project_l2", "project_ito", "project_wiener", "project_c_max_var",
    "project_pc1", "project_pc_d", "projection_fn",
    "projection_multiplicity",
    "reconstruct", "scan_over_splits", "make_run",
    "null_cache_dir", "save_null_draws",
    "wilcoxon_null_draws", "wilcoxon_critical_value",
    "lp_null_draws", "lp_critical_value",
    "rkhs_null_draws", "rkhs_critical_value",
    "far1_null_draws", "far1_critical_value",
    "simulate_wilcoxon_null", "simulate_lp_null",
    "simulate_rkhs_null", "simulate_far1_null",
    "TelcoCluster", "load_telco_far1",
]
