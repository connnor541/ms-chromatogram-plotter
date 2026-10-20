"""3D waterfall plot of all fractions of one file."""
import matplotlib
import numpy as np
import pandas as pd
import streamlit as st
from mpl_toolkits.mplot3d import proj3d

from matplotlib.colors import to_hex
from proteomics.plotting.base import fmt_fraction, locked, new_figure, with_method

def default_waterfall_colors(fractions):
    """Default hex color per fraction: the 'brg' colormap spread evenly over the fractions."""
    return [to_hex(c) for c in matplotlib.colormaps['brg'](np.linspace(0, 1, len(fractions)))]

def compute_zlabel_rotation(ax, flip=False):
    """
    Calculate the on-screen rotation for the z-axis label. We calculate the angle of rotation
    so that the z-label is perpendicular with the axis. Normally this angle is calculated internally by Matplotlib. 
    However, Matplotlib will always put the text in top to bottom orientation, bottom to top orientation looks better graphically. 
    So we need to calculate the angle rotation ourselves and then set the rotation of the z-label to that angle and then 
    flip the text upside down. 
    """
    renderer = ax.figure.canvas.get_renderer()
    ax.draw(renderer) 

    zaxis = ax.zaxis
    mins, maxs, tc, highs = zaxis._get_coord_info()
    minmax = np.where(highs, maxs, mins)
    maxmin = np.where(~highs, maxs, mins)

    edgep1s, edgep2s, _ = zaxis._get_all_axis_line_edge_points(
        minmax, maxmin, zaxis._label_position)
    edgep1, edgep2 = edgep1s[0], edgep2s[0]

    # ax.M will exist because we called ax.draw()
    pep = np.asarray(proj3d._proj_trans_points([edgep1, edgep2], ax.M))
    
    dx, dy = (ax.transAxes.transform([pep[0:2, 1]]) -
              ax.transAxes.transform([pep[0:2, 0]]))[0]

    angle = np.degrees(np.arctan2(dy, dx))
    if flip:
        angle += 180
    return angle

def place_zaxis_multiplier_text(ax, text, fontsize=10, pad_points=10):
    """
    Places a small multiplier label (e.g. r'$\\times 10^6$') right next to the
    TOP of the z-axis, following it correctly for ANY elev/azim.
    """
    
    fig = ax.get_figure()
    fig.canvas.draw()  # finalize the current projection matrix (ax.M) and layout

    zaxis = ax.zaxis
    mins, maxs, tc, highs = zaxis._get_coord_info()
    minmax = np.where(highs, maxs, mins)
    maxmin = np.where(~highs, maxs, mins)

    edgep1s, edgep2s, _ = zaxis._get_all_axis_line_edge_points(
        minmax, maxmin, zaxis._label_position)
    edgep_top = np.asarray(edgep2s[0], dtype=float)  # 3D point: top of the z-axis line

    # Outward direction in DATA space: away from the box's center, in the
    # (x, y) plane only (the corner's x,y is fixed along the whole vertical
    # z-axis edge, so this is well defined regardless of azimuth/elevation).
    x0, x1 = ax.get_xlim3d()
    y0, y1 = ax.get_ylim3d()
    cx, cy = (x0 + x1) / 2.0, (y0 + y1) / 2.0

    out_x, out_y = edgep_top[0] - cx, edgep_top[1] - cy
    norm_xy = np.hypot(out_x, out_y)
    if norm_xy < 1e-9:
        out_x, out_y = 1.0, 0.0  # degenerate fallback
    else:
        out_x, out_y = out_x / norm_xy, out_y / norm_xy

    scale = 0.08 * max(x1 - x0, y1 - y0)
    outward_point_3d = np.array([edgep_top[0] + out_x * scale,
                                  edgep_top[1] + out_y * scale,
                                  edgep_top[2]])

    # Project the axis-top point and the outward point into 2D display space
    proj = np.asarray(proj3d._proj_trans_points([edgep_top, outward_point_3d], ax.M))
    p_top_disp = ax.transData.transform(proj[0:2, 0])
    p_out_disp = ax.transData.transform(proj[0:2, 1])

    direction = p_out_disp - p_top_disp
    norm = np.hypot(*direction)
    direction = direction / norm if norm > 1e-6 else np.array([0.0, 1.0])

    pad_px = pad_points * fig.dpi / 72.0  # points -> pixels, AT THE CURRENT dpi
    text_disp = p_top_disp + direction * pad_px

    # Convert to DPI-independent figure-fraction coordinates before placing,
    # so a later savefig(dpi=...) at a different resolution can't shift it.
    text_fig_frac = fig.transFigure.inverted().transform(text_disp)

    fig.text(text_fig_frac[0], text_fig_frac[1], text, fontsize=fontsize,
              ha='center', va='center', transform=fig.transFigure)


@locked
def plot_waterfall_3d(all_binned, smooth_mode, elev=23, azim=-83, method=None, custom_colors=None):
    fractions = sorted(list(all_binned.keys()))
    if len(fractions) == 0:
        st.warning("No data available for the waterfall plot.")
        return None

    # Find max intensity across all fractions for scaling
    all_intensities = pd.concat([df['Processed'] for df in all_binned.values()])
    max_val = all_intensities.max()
    
    # Calculate scale factor: nearest power of 1000 (e.g., 10^3, 10^6, 10^9)
    # If max_val is 0, default to 1 to avoid division errors
    exponent = int(np.log10(max_val) // 3 * 3) if max_val > 0 else 0
    scale_factor = 10**exponent

    defaults = default_waterfall_colors(fractions)
    custom_colors = custom_colors or {}

    fig = new_figure(figsize=(10, 7))
    ax = fig.add_subplot(111, projection='3d')

    for i, fraction in enumerate(fractions):
        binned = all_binned[fraction]
        x = binned['bin_center'].values
        y_val = i * 2.5
        depth = np.full_like(x, y_val, dtype=float)
        
        # 2. Apply Dynamic Scale
        z = binned['Processed'].values / scale_factor 
        
        ax.plot(x, depth, z, color=custom_colors.get(fraction, defaults[i]), linewidth=1.8, alpha=0.8, label=fmt_fraction(fraction))

    ax.zaxis._axinfo['juggled'] = (1, 2, 0)
    ax.set_xlabel('Retention time (min)', fontsize=12, labelpad=10)
    ax.set_zlabel('Intensity (a.u.)', fontsize=12, labelpad=6)
    ax.zaxis.set_rotate_label(False)
    ax.set_yticks([])
    ax.set_ylabel('')

    for axis in (ax.xaxis, ax.yaxis, ax.zaxis):
        axis.pane.set_facecolor((1, 1, 1, 1))
        axis.pane.set_edgecolor('lightgray')
        axis._axinfo['grid']['color'] = (0.85, 0.85, 0.85, 1)

    ax.view_init(elev=elev, azim=azim)
    if method:
        ax.set_title(with_method('3D Waterfall', method), fontsize=14, fontweight='bold')
    ax.set_box_aspect((1.5, 2.0, 0.8))

    # Apply the rotation function
    zlabel_angle = compute_zlabel_rotation(ax, flip=False)
    ax.zaxis.get_label().set_rotation(zlabel_angle)    
    
    ax.legend(loc='upper right', frameon=True, fontsize=10)
    fig.tight_layout()

    # 3. Dynamic Scaling Label (e.g., x 10^6)
    place_zaxis_multiplier_text(ax, rf'$\times 10^{{{exponent}}}$', fontsize=10, pad_points=48)   

    return fig
