"""Proteomics Visualizer - Streamlit entry point.  Run with:  streamlit run app.py"""
import streamlit as st

from proteomics.ui.downloads import FigureDisplay, render_download_center
from proteomics.ui.loading import load_datasets, show_loaded_overview, stop_if_no_files
from proteomics.ui.orthogonality_section import render_orthogonality
from proteomics.ui.sections import render_boxplots, render_chromatograms, render_combined_visualizations
from proteomics.ui.sidebar import render_sidebar
from proteomics.ui.traces import TraceCache


st.set_page_config(page_title="Proteomics Visualizer", layout="wide")

st.title("Proteomics Data Visualization Dashboard")

settings = render_sidebar()

stop_if_no_files(settings)
datasets = load_datasets(settings)
if not datasets:
    st.stop()
show_loaded_overview(datasets)

traces = TraceCache(datasets, settings)
display = FigureDisplay()

render_chromatograms(datasets, settings, traces, display)
render_combined_visualizations(datasets, settings, traces, display)
render_boxplots(datasets, display)
render_orthogonality(datasets, display)
render_download_center(display.download_pngs)
