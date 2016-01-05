"""Provides tools to build and analyze overlap graphs on sequencing reads."""

from .graph import OverlapGraph

from .discovery import SeedExtensionParams, extend_segments, most_signitifcant_shift
from .plots import plot_shift_signifiance_discrimination, plot_all_seeds, plot_seed_extension_rws, plot_num_seeds_discrimination

from .builder import OverlapBuilder

__all__ = ['OverlapGraph', 'OverlapBuilder', 'SeedExtensionParams',
    'extend_segments', 'most_signitifcant_shift',
    'plot_shift_signifiance_discrimination', 'plot_all_seeds',
    'plot_seed_extension_rws'
]
