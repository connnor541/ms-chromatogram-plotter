"""Figure display plus the Download Center (PNG per plot and one ZIP)."""
import io
import re
import zipfile

import streamlit as st

from proteomics.plotting.base import fig_to_png


def slug(label):
    #function to replace faulty filename to safe one
    return re.sub(r'[^A-Za-z0-9]+', '_', label).strip('_')


class FigureDisplay:
    def __init__(self):
        self.download_pngs = {}   # maps file name (no extension) -> PNG bytes

    def show_fig(self, fig, name):
        #Render the figure once (thread-safe) and use the same PNG for display and download.
        png = fig_to_png(fig)
        st.image(png, width='stretch')
        self.download_pngs[name] = png


# ------------------------------------------------------------------
# 8. Download Center (only the plots that are currently displayed)
# ------------------------------------------------------------------
def render_download_center(download_pngs):
    st.divider()
    st.subheader("Download Center")

    if download_pngs:
        cols = st.columns(4)
        for i, (name, png) in enumerate(download_pngs.items()):
            cols[i % 4].download_button(
                label=f"Download {name}",
                data=png,
                file_name=f"{name.lower()}.png",
                mime="image/png"
            )

        zip_buffer = io.BytesIO()
        with zipfile.ZipFile(zip_buffer, "w") as zf:
            for name, png in download_pngs.items():
                zf.writestr(f"{name.lower()}.png", png)

        st.download_button(
            label="📦 Download ALL Displayed Plots (ZIP)",
            data=zip_buffer.getvalue(),
            file_name="all_plots.zip",
            mime="application/zip",
            width='stretch'
        )
    else:
        st.caption("No plots to download yet.")
