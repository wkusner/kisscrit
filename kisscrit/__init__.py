"""kisscrit: critical configurations of the kissing / Tammes min-distance function on spheres.

Status: untrusted (exploration or tooling; no part of the proof).

Points x_1..x_N on S^{d-1}; objective r(X) (equivalently the max cosine). Tools: Clarke-criticality
via balanced stresses, co-index via the Lagrangian Hessian, the square KKT system of a contact graph,
monotone ascent, contact-graph surgery, high-precision refinement and minimal polynomials, and (d=3)
enumeration of hull triangulations.
"""
import warnings as _w
_w.filterwarnings('ignore', message='The hashes produced for graphs')
from .geometry import normalize, r_of_c, c_of_r, theta_of_c, contacts, random_config
from .stress import analyze, balanced_stress
from .kkt import solve_kkt, check_solution
from .ascent import ascend
from .search import Census, polish, random_ascent, surgery
from .symmetry import symmetries, group_orders
__version__ = '0.1.0'
