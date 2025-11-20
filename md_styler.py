# md_styler.py
from __future__ import annotations

import colorsys
from dataclasses import dataclass, field
from typing import Dict, List, Tuple, Optional, Iterable

import matplotlib as mpl
import matplotlib.pyplot as plt
from matplotlib.ticker import ScalarFormatter, AutoMinorLocator
from mpl_toolkits.axes_grid1 import Divider, Size

# ---------- utilities

def _rgb01_to_hex(rgb: Tuple[float, float, float]) -> str:
    r, g, b = [max(0, min(1, x)) for x in rgb]
    return "#{:02x}{:02x}{:02x}".format(int(round(r * 255)),
                                        int(round(g * 255)),
                                        int(round(b * 255)))

def _matte_rgb(rgb: Tuple[float, float, float],
               sat_scale: float = 0.75,
               light_shift: float = -0.02) -> Tuple[float, float, float]:
    """
    Convert RGB->HLS, reduce saturation for AO 'chalky' look, optionally tweak lightness.
    Values are conservative to preserve print reliability.
    """
    r, g, b = rgb
    h, l, s = colorsys.rgb_to_hls(r, g, b)
    s = max(0.0, min(1.0, s * sat_scale))
    l = max(0.0, min(1.0, l + light_shift))
    r2, g2, b2 = colorsys.hls_to_rgb(h, l, s)
    return (r2, g2, b2)

# ---------- main class

@dataclass
class MDStyler:
    # global matte tuning
    sat_scale: float = 0.75   # lower = more matte
    light_shift: float = -0.02  # slight darken helps AO feel & print contrast

    # geometry (axes/data rectangle) in inches
    axes_square: Tuple[float, float] = (2.2, 2.2)      # 1:1
    axes_horizontal: Tuple[float, float] = (3.0, 2.0)  # 3:2

    # dpi & fonts
    dpi: int = 300
    base_fontsize: float = 9.0
    font_family: str = "Helvetica"  # will fall back to Arial/DejaVu Sans if unavailable

    # AA/CG line styles
    ls_aa: str = "-"
    ls_cg: Tuple[int, Tuple[int, ...]] = (0, (6, 2))

    # expose after init
    colors: Dict[str, str] = field(init=False)
    palette6: List[str] = field(init=False)

    def __post_init__(self):
        # ---- VMD base RGB (0-1) from your console (AOChalky/Glass cues)
        base = {
            "cyan":   (0.0, 0.88, 1.0),            # cyan2
            "orange": (1.0, 0.5, 0.0),             # orange
            "black":  (0.0, 0.0, 0.0),             # AO 'gray' was used visually; keep true black for text/axes
            "gray":   (0.35, 0.35, 0.35),          # AO gray
            "green":  (0.0, 1.0, 0.0),             # green
            "red":    (0.81, 0.0, 0.0),            # red3
            "box":    (0.02, 0.38, 0.67),          # PBC box (blue2)
            "blue3":  (0.01, 0.04, 0.93),          # optional accent
        }

        # apply matte transform (except 'black' which we keep literal)
        matte = {}
        for k, rgb in base.items():
            if k == "black":
                matte[k] = rgb
            else:
                matte[k] = _matte_rgb(rgb, self.sat_scale, self.light_shift)

        self.colors = {k: _rgb01_to_hex(v) for k, v in matte.items()}

        # ---- build a 6-color, CVD-safe list harmonized to your VMD hues (also matte)
        #   order: high-contrast first four, then two accents
        #   cyan, orange, deep blue (from box), vermillion/red3, olive-green (from green), purple (derived)
        # create a subdued purple close to Okabe–Ito's but nudged toward your blue range
        oi_purple = _matte_rgb((0.59, 0.44, 0.84), self.sat_scale, self.light_shift)
        # olive-ish green toned from your green
        olive = _matte_rgb((0.33, 0.56, 0.0), self.sat_scale, self.light_shift)

        self.palette6 = [
            self.colors["cyan"],             # 1
            self.colors["orange"],           # 2
            self.colors["box"],              # 3 (deep blue)
            self.colors["red"],              # 4 (vermillion-ish)
            _rgb01_to_hex(olive),            # 5
            _rgb01_to_hex(oi_purple),        # 6
        ]

    # ---------- public API

    def apply(self) -> "MDStyler":
        """Apply rcParams for a JCTC-friendly, AO-matte style. Call at the top of analysis scripts."""
        mpl.rcdefaults()

        # Fonts & math
        mpl.rcParams.update({
            "font.family": "sans-serif",
            "font.sans-serif": [self.font_family, "Arial", "Helvetica", "DejaVu Sans"],
            "mathtext.fontset": "dejavusans",
            "text.kerning_factor": 0,
            # Base sizes
            "font.size": self.base_fontsize,
            "axes.titlesize": self.base_fontsize,
            "axes.labelsize": self.base_fontsize,
            "legend.fontsize": self.base_fontsize * 0.95,
            "xtick.labelsize": self.base_fontsize * 0.9,
            "ytick.labelsize": self.base_fontsize * 0.9,
            # Lines/markers
            "lines.linewidth": 2.0,
            "lines.markersize": 5.5,
            # Spines/ticks/grid
            "axes.spines.top": False,
            "axes.spines.right": False,
            "xtick.direction": "in",
            "ytick.direction": "in",
            "xtick.major.width": 0.9,
            "ytick.major.width": 0.9,
            "xtick.minor.width": 0.8,
            "ytick.minor.width": 0.8,
            "xtick.major.size": 3.5,
            "ytick.major.size": 3.5,
            "xtick.minor.size": 2.5,
            "ytick.minor.size": 2.5,
            "grid.linewidth": 0.6,
            "grid.color": "#d0d0d0",
            "axes.grid": False,
            # Figure
            "figure.dpi": self.dpi,
            "savefig.dpi": self.dpi,
            "savefig.transparent": False,
            "savefig.facecolor": "white",
            "figure.facecolor": "white",
            "axes.facecolor": "white",
            # Color cycle
            "axes.prop_cycle": mpl.cycler(color=self.palette6),
            # Errorbars
            "errorbar.capsize": 3,
        })

        return self

    # --- fixed data-rectangle creators (identical axes size regardless of labels)

    def _fixed_axes(self, axes_w: float, axes_h: float,
                    margins: Tuple[float, float, float, float] = (0.35, 0.15, 0.1, 0.25),
                    facecolor: str = "white"):
        """
        Create a figure with a *fixed-size axes area* in inches using Divider.
        margins = (left, right, bottom, top) in inches reserved for labels/colorbars/etc.
        Tweak margins per-plot if needed; the axes itself remains constant.
        """
        left, right, bottom, top = margins
        fig_w = axes_w + left + right
        fig_h = axes_h + bottom + top

        fig = plt.figure(figsize=(fig_w, fig_h), facecolor=facecolor)
        # Define sizes: fixed margins + fixed axes
        h = [Size.Fixed(left),  Size.Fixed(axes_w), Size.Fixed(right)]
        v = [Size.Fixed(bottom), Size.Fixed(axes_h), Size.Fixed(top)]
        divider = Divider(fig, (0, 0, 1, 1), h, v, aspect=False)
        ax = plt.axes(divider.get_position())
        ax.set_axes_locator(divider.new_locator(nx=1, ny=1))

        # sensible defaults: scientific notation & minor ticks
        self._apply_axes_defaults(ax)
        return fig, ax

    def fig_square(self, margins: Tuple[float, float, float, float] = (0.35, 0.15, 0.1, 0.25)):
        return self._fixed_axes(*self.axes_square, margins=margins)

    def fig_horizontal(self, margins: Tuple[float, float, float, float] = (0.35, 0.15, 0.1, 0.25)):
        return self._fixed_axes(*self.axes_horizontal, margins=margins)

    # --- accessors

    def get_color(self, name: str) -> str:
        return self.colors[name]

    def get_palette(self, n: int = 6) -> List[str]:
        if n <= 6:
            return self.palette6[:n]
        # Repeat in order if more needed; user can override
        times = (n + 5) // 6
        return (self.palette6 * times)[:n]

    # --- helpers

    def as_aa(self, **plot_kwargs):
        """Return default kwargs for AA (solid)."""
        return dict(linestyle=self.ls_aa, **plot_kwargs)

    def as_cg(self, **plot_kwargs):
        """Return default kwargs for CG (dashed)."""
        return dict(linestyle=self.ls_cg, **plot_kwargs)

    def hist_kde_compare(self,
                         ax: Optional[mpl.axes.Axes],
                         aa: Iterable[float],
                         cg: Iterable[float],
                         aa_color: Optional[str] = None,
                         cg_color: Optional[str] = None,
                         bins: Optional[int] = None,
                         alpha_hist: float = 0.35,
                         kde_bw: Optional[str] = "scott",
                         label_aa: str = "AA",
                         label_cg: str = "CG"):
        """
        AA: filled histogram; CG: KDE curve (dashed).
        Uses seaborn if available; otherwise falls back to matplotlib (and scipy if present).
        """
        if ax is None:
            _, ax = self.fig_horizontal()

        aa_color = aa_color or self.palette6[0]  # cyan-ish
        cg_color = cg_color or self.palette6[2]  # deep blue by default

        try:
            import seaborn as sns  # type: ignore
            # style-neutral call so our rcParams dominate
            sns.histplot(aa, ax=ax, bins=bins, stat="density",
                         color=aa_color, alpha=alpha_hist, edgecolor="none", label=label_aa)
            sns.kdeplot(cg, ax=ax, bw_method=kde_bw, color=cg_color,
                        linestyle=self.ls_cg, linewidth=2.0, label=label_cg)
        except Exception:
            # Matplotlib fallback
            import numpy as np
            aa = np.asarray(aa)
            cg = np.asarray(cg)
            ax.hist(aa, bins=bins or "auto", density=True, color=aa_color,
                    alpha=alpha_hist, edgecolor="none", label=label_aa)
            try:
                from scipy.stats import gaussian_kde  # type: ignore
                kde = gaussian_kde(cg, bw_method=kde_bw)
                xs = np.linspace(np.nanmin(cg), np.nanmax(cg), 512)
                ax.plot(xs, kde(xs), color=cg_color, linestyle=self.ls_cg,
                        linewidth=2.0, label=label_cg)
            except Exception:
                # crude line via histogram if scipy missing
                hist, edges = np.histogram(cg, bins=bins or "auto", density=True)
                xs = 0.5 * (edges[:-1] + edges[1:])
                ax.plot(xs, hist, color=cg_color, linestyle=self.ls_cg,
                        linewidth=2.0, label=label_cg)

        ax.legend(frameon=False)
        return ax

    def save(self, fig: mpl.figure.Figure, path: str,
             transparent: bool = False, tight: bool = True, bbox_inches: Optional[str] = None):
        if tight and bbox_inches is None:
            bbox_inches = "tight"
        fig.savefig(path, dpi=self.dpi, transparent=transparent, bbox_inches=bbox_inches)

    # ---------- internal

    def _apply_axes_defaults(self, ax: mpl.axes.Axes):
        ax.minorticks_on()
        ax.xaxis.set_minor_locator(AutoMinorLocator())
        ax.yaxis.set_minor_locator(AutoMinorLocator())

        # scientific/SI-ish notation
        fmt = ScalarFormatter(useMathText=True)
        fmt.set_powerlimits((-3, 3))  # switch to sci outside 1e-3..1e3
        ax.xaxis.set_major_formatter(fmt)
        ax.yaxis.set_major_formatter(fmt)

        # default label colors & spine colors to matte gray/black
        ax.tick_params(colors=self.colors["black"])
        for spine in ax.spines.values():
            spine.set_color("#000000")
            spine.set_linewidth(0.9)
