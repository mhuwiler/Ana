"""Flexible PDF report builder for analysis plots.

Usage:
    from utils.pdf import PDFReport

    pdf = PDFReport("output/plots/coi_study")
    for var in pdf.variables("b_coi*"):
        pdf.page_grid(var, "stacked", cols=2)
        pdf.page_grid(var, "sig", cols=2)
    pdf.save("report.pdf")
"""

import fnmatch
import glob
import math
import os

from fpdf import FPDF


class PDFReport:
    """Build a PDF report from analysis plot PNGs."""

    def __init__(self, plot_dir, orientation="L"):
        self._dir = plot_dir
        self._pdf = FPDF(orientation=orientation, unit="mm", format="A4")
        self._pdf.set_auto_page_break(auto=False)
        if orientation == "L":
            self._pw, self._ph = 297, 210
        else:
            self._pw, self._ph = 210, 297
        self._margin = 5

        # Auto-discover
        self._plots = {}    # (trigger, var, plot_type) → filepath
        self._triggers = set()
        self._variables = set()
        self._plot_types = set()
        self._discover()

    def _discover(self):
        """Scan mc/{type}/{trigger}/{var}.png structure."""
        for png in glob.glob(os.path.join(self._dir, "mc", "*", "*", "*.png")):
            rel = os.path.relpath(png, self._dir)
            parts = rel.split(os.sep)
            if len(parts) != 4:
                continue
            _, ptype, trig, fname = parts
            var = fname.replace(".png", "")
            self._triggers.add(trig)
            self._variables.add(var)
            self._plot_types.add(ptype)
            self._plots[(trig, var, ptype)] = png

        # Also discover cm/ plots
        for png in glob.glob(os.path.join(self._dir, "mc", "cm", "*.png")):
            var = os.path.basename(png).replace(".png", "")
            self._variables.add(f"cm_{var}")
            self._plot_types.add("cm")
            self._plots[("_all", f"cm_{var}", "cm")] = png

    def triggers(self):
        """Return sorted list of discovered trigger names."""
        return sorted(self._triggers)

    def plot_types(self):
        """Return sorted list of discovered plot type names."""
        return sorted(self._plot_types)

    def variables(self, pattern="*"):
        """Return sorted list of variable names matching pattern(s).

        pattern can be a string or list of strings (fnmatch patterns).
        """
        if isinstance(pattern, str):
            patterns = [pattern]
        else:
            patterns = pattern
        return sorted(v for v in self._variables
                      if any(fnmatch.fnmatch(v, p) for p in patterns))

    def _get_image(self, trigger, var, plot_type):
        """Get image path for a specific (trigger, var, plot_type)."""
        return self._plots.get((trigger, var, plot_type))

    def page_grid(self, var, plot_type, triggers=None, cols=2, title=None):
        """Add a page with a grid of trigger plots for one variable + one plot type.

        Parameters
        ----------
        var : str
            Variable name.
        plot_type : str
            Plot type (stacked, shape, sig, eff, etc.).
        triggers : list[str] or None
            Trigger names. None = all discovered triggers.
        cols : int
            Number of columns in the grid.
        title : str or None
            Page title. None = auto-generate from var + plot_type.
        """
        use_triggers = triggers or self.triggers()
        images = []
        labels = []
        for trig in use_triggers:
            path = self._get_image(trig, var, plot_type)
            images.append(path)
            labels.append(trig)

        if not any(images):
            return

        if title is None:
            title = f"{var} - {plot_type}"

        self._add_grid_page(images, labels, cols, title)

    def page_single(self, image_path, title=None):
        """Add a page with a single full-page image."""
        if not image_path or not os.path.exists(image_path):
            return
        self._pdf.add_page()
        m = self._margin
        if title:
            self._pdf.set_font("Helvetica", "B", 12)
            self._pdf.set_xy(m, m)
            self._pdf.cell(self._pw - 2 * m, 6, title, align="C")
            y0 = m + 8
        else:
            y0 = m
        self._pdf.image(image_path, x=m, y=y0,
                        w=self._pw - 2 * m, h=self._ph - y0 - m)

    def page_custom(self, images, labels=None, cols=2, title=None):
        """Add a page with arbitrary images in a grid.

        Parameters
        ----------
        images : list[str]
            Image file paths (None entries = empty cell).
        labels : list[str] or None
            Labels for each image (shown above each cell).
        cols : int
            Number of columns.
        title : str or None
            Page title.
        """
        if labels is None:
            labels = [""] * len(images)
        self._add_grid_page(images, labels, cols, title)

    def save(self, output_path):
        """Write PDF to disk."""
        os.makedirs(os.path.dirname(os.path.abspath(output_path)), exist_ok=True)
        self._pdf.output(output_path)
        sz = os.path.getsize(output_path) / 1e6
        n = self._pdf.pages_count
        print(f"Saved: {output_path} ({sz:.1f} MB, {n} pages)")

    def _add_grid_page(self, images, labels, cols, title):
        """Internal: add a page with images in a grid."""
        self._pdf.add_page()
        m = self._margin

        # Title
        y0 = m
        if title:
            self._pdf.set_font("Helvetica", "B", 12)
            self._pdf.set_xy(m, m)
            self._pdf.cell(self._pw - 2 * m, 6, title, align="C")
            y0 = m + 8

        n = len(images)
        rows = math.ceil(n / cols)
        usable_w = self._pw - 2 * m
        usable_h = self._ph - y0 - m
        cell_w = usable_w / cols
        cell_h = usable_h / rows

        has_labels = any(labels)
        label_h = 5 if has_labels else 0

        for i, img_path in enumerate(images):
            ri = i // cols
            ci = i % cols
            x = m + ci * cell_w
            # Each cell has: label (5mm) + image
            cell_total = (usable_h) / rows
            img_h = cell_total - label_h
            y = y0 + ri * cell_total

            # Label above each image
            if has_labels and i < len(labels) and labels[i]:
                self._pdf.set_font("Helvetica", "B", 9)
                self._pdf.set_xy(x, y)
                self._pdf.cell(cell_w, label_h, labels[i], align="C")

            if img_path and os.path.exists(img_path):
                pad = 1
                self._pdf.image(img_path, x=x + pad, y=y + label_h + pad,
                                w=cell_w - 2 * pad, h=img_h - 2 * pad)
