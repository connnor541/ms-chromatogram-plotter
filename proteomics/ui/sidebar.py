#Sidebar: file upload, fractionation method detection, processing and plot settings.
from dataclasses import dataclass, field

import streamlit as st

import proteomics.analysis.method_detection as md


@dataclass
class SidebarSettings:
    is_plgs: bool
    active: list            # list of (uploaded_file, method)
    labels: list            # unique label per active file
    smooth_mode: str
    bin_width: float
    window: float
    bar_width: float
    show_filter: bool
    combine_mode: str
    filter_quant: bool
    separate_modified: bool
    rt_tolerance: float
    show_chromo: list = field(default_factory=list)
    show_stacked: list = field(default_factory=list)
    show_waterfall: list = field(default_factory=list)
    show_overlap: list = field(default_factory=list)
    elev: int = 23
    azim: int = -83


# ------------------------------------------------------------------
# 1. Sidebar: files + fractionation method detection
# ------------------------------------------------------------------
def _render_file_inputs():
    st.sidebar.header("Data & Filter Settings")

    #prompt the user to specify data input type Progenesis or PLGS
    input_type = st.sidebar.radio("Input type", ['Progenesis', 'PLGS'], index=0, horizontal=True)
    is_plgs = input_type == 'PLGS'

    #makes the field for dropping input files
    if is_plgs:
        uploaded_files = st.sidebar.file_uploader(
            "Upload PLGS Excel file(s)", type=['xlsx'], accept_multiple_files=True,
            key="uploader_plgs"
        )
        st.sidebar.caption("One Excel file per fractionation method. Each sheet is one fraction, "
                           "numbered by sheet order (first sheet = fraction 1).")
    else:
        uploaded_files = st.sidebar.file_uploader(
            "Upload Proteomics CSV file(s)", type=['csv'], accept_multiple_files=True,
            key="uploader_progenesis"
        )

    active = []  # list of (uploaded_file, method)
    if uploaded_files:
        auto_methods, notes = md.auto_assign_methods([f.name for f in uploaded_files])
        options = md.METHOD_OPTIONS + [md.SKIP]

        with st.sidebar.expander("Fractionation methods", expanded=True):
            st.caption("Detected from the file name. Change it if it is wrong, "
                       "or choose 'Skip file' to leave a file out.")
            for note in notes:
                st.warning(note)
            for i, (f, detected) in enumerate(zip(uploaded_files, auto_methods)):
                default = detected if detected in md.METHOD_OPTIONS else md.SKIP
                # the detected default is part of the key so the dropdown resets when detection changes
                choice = st.selectbox(f.name, options, index=options.index(default),
                                      key=f"method_{i}_{f.name}_{default}")
                if choice != md.SKIP:
                    active.append((f, choice))

    labels = md.make_unique_labels([m for _, m in active])
    return is_plgs, active, labels


# ------------------------------------------------------------------
# 2. Sidebar: processing settings
# ------------------------------------------------------------------
def _render_processing_settings(is_plgs):
    smooth_mode = st.sidebar.selectbox("Smoothing Mode", ['none', 'savgol', 'gaussian', 'decimate'])
    bin_width = st.sidebar.slider("Bin Width (min)", 0.01, 0.5, 0.1, 0.01)
    window = st.sidebar.slider("Smoothing Window (min)", 0.1, 5.0, 1.0, 0.1)
    bar_width = st.sidebar.slider("Bar Width", 0.01, 0.2, 0.05, 0.01)
    show_filter = st.sidebar.checkbox("Show Filter Overlay", value=False)
    combine_mode = st.sidebar.radio("Aggregation Method", ['max', 'sum'])
    if is_plgs:
        filter_quant = False   # PLGS files have no 'Use in quantitation' column
    else:
        filter_quant = st.sidebar.checkbox("Enforce 'Use in quantitation' Filter", value=False)
    separate_modified = st.sidebar.checkbox(
        "Keep modified forms separate", value=False,
        help="Off: a modified peptide is merged with the other forms of the same sequence. "
             "On: every modification pattern of a sequence stays its own entry.")
    #Progenesis keeps its previous behaviour (identical RT only) by default, PLGS RTs of the same peak differ slightly
    #RT tolerance makes it that same sequence with slightly different RT can merge
    rt_tolerance = st.sidebar.number_input(
        "RT tolerance for merging (min)", min_value=0.0, max_value=5.0,
        value=0.05 if is_plgs else 0.0, step=0.01, format="%.3f",
        key="rt_tolerance_plgs" if is_plgs else "rt_tolerance_progenesis",
        help="Rows of the same sequence whose retention times differ by at most this value are merged. "
             "0 merges identical retention times only.")
    return (smooth_mode, bin_width, window, bar_width, show_filter, combine_mode, filter_quant,
            separate_modified, rt_tolerance)

# ------------------------------------------------------------------
# 3. Sidebar: which files to show per plot type (all off by default)
# ------------------------------------------------------------------
def _render_plot_selection(labels):
    st.sidebar.header("Per-file Plots")
    st.sidebar.caption("Pick the files (fractionation methods) to display. Nothing is shown by default.")
    show_chromo = st.sidebar.multiselect("Chromatograms", labels, default=[])
    show_stacked = st.sidebar.multiselect("Stacked Fractions", labels, default=[])
    show_waterfall = st.sidebar.multiselect("3D Waterfall", labels, default=[])
    show_overlap = st.sidebar.multiselect("Peptide Overlap (UpSet)", labels, default=[])
    return show_chromo, show_stacked, show_waterfall, show_overlap


def _render_waterfall_view():
    #change elevation and azimuth values to change perspective of waterfall plot. Elevation to move vertically and azimuth to move horizontally
    st.sidebar.header("3D Waterfall View")
    elev = st.sidebar.slider("Elevation", 0, 90, 23)
    azim = st.sidebar.slider("Azimuth", -180, 180, -83)
    return elev, azim


def render_sidebar() -> SidebarSettings:
    """Draw the whole sidebar (in the original widget order) and return everything the user chose."""
    is_plgs, active, labels = _render_file_inputs()
    (smooth_mode, bin_width, window, bar_width,
     show_filter, combine_mode, filter_quant,
     separate_modified, rt_tolerance) = _render_processing_settings(is_plgs)
    show_chromo, show_stacked, show_waterfall, show_overlap = _render_plot_selection(labels)
    elev, azim = _render_waterfall_view()

    return SidebarSettings(
        is_plgs=is_plgs, active=active, labels=labels,
        smooth_mode=smooth_mode, bin_width=bin_width, window=window, bar_width=bar_width,
        show_filter=show_filter, combine_mode=combine_mode, filter_quant=filter_quant,
        separate_modified=separate_modified, rt_tolerance=rt_tolerance,
        show_chromo=show_chromo, show_stacked=show_stacked,
        show_waterfall=show_waterfall, show_overlap=show_overlap,
        elev=elev, azim=azim,
    )
