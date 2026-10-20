https://ms-chromatogram-plotter.streamlit.app/

# Proteomics Visualizer

    pip install -r requirements.txt
    streamlit run app.py

## Layout

    app.py                          entry point (wires the sections together)
    proteomics/
      data/                         DATA READ-IN
        loaders.py                  raw Progenesis CSV / PLGS Excel -> DataFrame
        cleaning.py                 DataFrame -> tidy table (Fraction, Retention_time, Sequence, Intensity)
      analysis/                     computation, no plotting
        method_detection.py         fractionation method from file name
        aa_properties.py            pI / GRAVY / MW per peptide
        peptide_stats.py            per-fraction stats, fraction overlap sets
        traces.py                   exact / binned / smoothed traces
        orthogonality.py            entropy + mutual information
      plotting/                     DATA VISUALIZATION (matplotlib only)
        base.py                     thread lock, figure factory, PNG rendering, labels
        chromatograms.py            single + stacked chromatograms
        waterfall.py                3D waterfall
        overlap.py                  UpSet plot
        boxplots.py                 pI / GRAVY / intensity boxplot grids
        orthogonality_plots.py      entropy / heat-map figures
      ui/                           Streamlit layer
        sidebar.py                  all sidebar widgets -> SidebarSettings
        loading.py                  load + clean uploaded files -> datasets
        traces.py                   lazy trace cache
        sections.py                 chromatogram / combined / boxplot sections
        orthogonality_section.py    orthogonality section
        downloads.py                figure display + Download Center
