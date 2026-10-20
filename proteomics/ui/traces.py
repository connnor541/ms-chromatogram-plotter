"""Lazy, cached construction of exact / binned / smoothed traces per file."""
import streamlit as st

from proteomics.analysis.traces import apply_smoothing_pipeline, build_binned_trace, build_exact_trace


class TraceCache:
    # Traces are only built for files that need them (chromatogram / stacked / waterfall)
    def __init__(self, datasets, settings):
        self._datasets = datasets
        self._settings = settings
        self._cache = {}

    def get(self, label):
        if label not in self._cache:
            ds = self._datasets[label]
            s = self._settings
            all_exact, all_binned = {}, {}
            with st.expander(f"Trace & smoothing log - {label}", expanded=False):
                for fraction in ds['fractions']:
                    exact = build_exact_trace(ds['df_clean'], fraction, s.combine_mode)
                    binned_raw = build_binned_trace(ds['df_clean'], fraction, s.bin_width,
                                                    ds['x_min'], ds['x_max'], s.combine_mode)
                    all_exact[fraction] = exact
                    all_binned[fraction] = apply_smoothing_pipeline(binned_raw, s.smooth_mode, s.window, s.bin_width)
            self._cache[label] = (all_exact, all_binned)
        return self._cache[label]
