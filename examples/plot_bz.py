#!/usr/bin/env python3
"""Overlay of chi0 on the Brillouin zone, from the files written by chi_bz.x, as a PDF.

usage:  python3 plot_bz.py [<prefix>] [re|im] [experiment.dat] [options]

The plot is the map of chi0 on a plane of the Brillouin zone (chi_bz.x), a colour map in the shape
of the zone. 1D plots are cross sections of it: the last question of the map asks for one (constant x
or y in r.l.u.), and its 1D plot follows. (The line of a 1D system from chi_bz.x is drawn as a 1D
plot directly.) Without a prefix, or with the name of the chi file (e.g. chi_6.4521_300), the script
asks which map (of those in the folder). The questions:
  2D  part, planes, region, view, colour bar, colour map, colour scale (percentiles, gamma,
      Gamma left out), contours, zone scale, font, labels, experiment, output, and a 1D cross section
      (cut) at constant x or y of the map in r.l.u.; it is drawn as a solid line on the map, its
      values are written to <prefix>_cut_y0.25.dat, and its 1D plot follows with the 1D questions
      (and the range of its x axis)
  1D  part, lines (+c / -c: together = their average / panels / one), x axis range (r.l.u.), x axis
      unit (rlu or invA = 1/Angstrom), x tick spacing, y tick spacing, y axis range, mark the maximum
      (red dashed line), font, labels, experiment, output

  <prefix>        as printed by chi_bz.x, e.g. chi_6.4521_300_ab: reads <prefix>_map.dat,
                  <prefix>_bz.dat and <prefix>_points.dat. For a plane at k = c /= 0 chi_bz.x writes
                  the planes +c and -c into one file (chi_6.4521_300_ab_L+-0.1667); give it, or
                  just chi_6.4521_300_ab. By default the two planes are drawn together as one map,
                  their average at each H, K (--planes together); --planes panels puts them side by
                  side, each labelled with its plane, with one colour scale; --planes + or - one.
                  The same for a cross section: together (one curve, the average), panels or one
  re | im         Re chi0 (default) or Im chi0/w
  experiment.dat  optional: measured points in 1/Angstrom in the axes of the map (plane: columns
                  qx qy; line: column q). q-vectors in H K L are best given as points in bz.in.

The script asks for the settings in the terminal (Enter keeps the value in brackets), writes the
PDF, and then offers a numbered list of the settings: change one or several and the PDF is written
again (a PDF viewer that reloads the file shows each version). It prints the command that gives
the same figure without questions. No display is needed.

By default (view 'zone') the figure is just the first Brillouin zone filled with chi0: the colour
map fills the zone up to its border (values outside the zone in the chi_bz.x map are filled in
from the periodicity of the lattice) and stops there, with the zone border, the labelled points
and the colour bar; no axes, no grid lines. After each PDF the terminal gives the largest value
(away from Gamma) and its place: the fraction of the way from Gamma to the border, the direction
(toward a corner or an edge centre), q in 1/Angstrom and H K L.
On request, the zone can carry a scale: dashed copies of the border at fractions of the way from
Gamma to the border (along Gamma-M they are fractions of Gamma-M, along Gamma-K of Gamma-K),
guide lines from Gamma to the corners and edge centres, and a length bar in 1/Angstrom. The maximum
(every equivalent point) is marked by a bold red + (default; --no-max: none).

Every view keeps the shape of the zone. View 'cartesian' adds axes through Gamma along the cartesian
directions of the plane (k_x, k_y, ... in 1/Angstrom), view 'crystal' axes along the two reciprocal
vectors of the plane (H, K or L in r.l.u.). The axes run from border to border of the zone, without
arrows or names (they belong in the caption); the tick spacing can be set (e.g. 0.1), and the axes
are drawn in a colour that stands out from the colour map. Labels of points on the zone border sit
outside the zone, off the axes; every labelled point sits at its own H, K on the plane (Gamma where the
axes cross, also on a plane L = c /= 0), and for monoclinic and triclinic lattices chi_bz.x makes the
zone the cell -1/2 ... 1/2, so that B, D, ... lie on its border. The colour bar is a small plain bar labelled min and max (or the two
values, or a full scale), and the PDF is cropped close to the figure.

2D settings (the options give their starting values):
  part                 re or im
  --region R           first: the first Brillouin zone; extended: the whole map of chi_bz.x with all zone
                       boundaries (the first zone drawn heavier). Default: as chi_bz.x was run (zone = ...)
  --planes P           with the planes +c / -c: together (default; their average in one map), panels (side by
                       side, one colour scale), + or - (one plane)
  --panels L           labels of the two panels: plane (default; e.g. L = +1/6), letters ((a), (b)), none
  --view V             zone (default): just the zone with chi0; cartesian: k axes through Gamma (1/Angstrom);
                       crystal: axes along the reciprocal vectors (r.l.u.)
  --tick T             tick spacing on the axes, e.g. 0.1 (1/Angstrom) or 1/6 (r.l.u.); default auto
  --axis-colour C      colour of the axes; default auto: cyan for inferno/magma/hot, red for gray and
                       cividis, lime for plasma, magenta for viridis, white for turbo
  --cbar C             colour bar: minmax (default; labelled min and max), values (the bottom and top
                       values), ticks (a full scale with the quantity)
  --rings LIST         dashed outlines at these fractions of Gamma -> border, e.g. 1/3,2/3 (default none)
  --guides             guide lines from Gamma to the corners and edge centres
  --bar                a length bar in 1/Angstrom under the zone
  --no-max             do not mark the maximum (2D: a bold red +; 1D: a red dashed line; marked by default)
  --cmap NAME          inferno (default), magma, hot, gray, plasma, viridis, cividis, turbo
  --low P              bottom of the colour scale at the P-th percentile of the values (default 0)
  --high P             top of the colour scale at the P-th percentile (default 100)
  --vmin V, --vmax V   bottom / top as values, instead of --low / --high
                       (in the questions: a percentile such as 60, or a value such as =0.185)
  --gamma G            power-law colour scale; G > 1 gives more colours to the largest values (default 1)
  --levels N           levels of the filled contours and contour lines (default 12)
  --lines              contour lines over the map
  --font F             font size in points (default 12)
  --gamma-steps N      the colour scale leaves out Gamma and its equivalents up to N mesh steps away
                       (default 1); they are drawn in grey above the scale. 0: use every value
  --no-zones, --no-labels   hide the zone boundaries / the point labels
  --out FILE           the PDF (default <prefix>_<re|im>.pdf)
  --batch              no questions: write the PDF with the options (also when the input is not a
                       terminal, e.g. in a job script); needs the prefix of one plot
1D plots:
  --tick T             spacing of the x ticks of a chain, in the unit of the x axis (a cut: --cut-tick).
                       Tick labels are decimals (also for a spacing given as 1/6)
  --ytick T            spacing of the y ticks, e.g. 0.05 (default auto)
  --minor M, --yminor M  minor ticks (short, no labels) every M on the x / y axis (default none); on a
                       map --minor puts them on the axis lines of the cartesian and crystal views
  --ymin V, --ymax V   the y axis from / to (default: all values with a small margin); the terminal
                       says when the value at Gamma is far above the rest, with a --ymax that shows it
  (maximum)            a red dashed vertical line at the largest value away from Gamma (and its value
                       and H K L in the terminal). The origin (0 of the x axis) is always a black solid
                       line; labelled points on the line: solid lines with their labels on top; measured
                       q-values (experiment file) are green dash-dot lines
  --xaxis U            rlu (the crystal coordinate along the line; default of a cut) or invA (1/Angstrom)
  --planes, --panels, --font, --no-labels, --out as above (--no-zones: a chain of chi_bz.x)
Cross sections of a 2D map:
  --cut C              constant x or y of the map in r.l.u. (x: the first in-plane crystal coordinate, H of
                       an ab map; y: the second, K), any value: y=0.25, x=-1/3. On a mesh line the values
                       are those of the mesh, between two mesh lines interpolated linearly between them.
                       With --batch writes the map (with the cut as a solid line) and the 1D plot
                       <prefix>_cut_y0.25_re.pdf, using --xaxis, --ytick, --ymin, --ymax, --no-max
  --xrange FROM TO     the x axis of the cut in r.l.u. (default -0.5 0.5); any range: chi0 is periodic, so
                       the values come from the whole k-mesh (e.g. --xrange -1 1: two zones)
  --cut-tick T         the x tick spacing of the cut (--cut-minor M: its minor ticks);  --cut-only: only the 1D plot of the cut

Weak dispersion (chi0 changes by a few percent over the zone): to keep the largest values the
brightest and clearly apart from the rest,
  * use a colour map whose brightness rises steadily to the top: inferno, magma, hot or gray
    (not turbo or jet, whose top end is darker than their middle);
  * leave Gamma out of the colour scale (default): the intraband peak at q -> 0 is usually the
    largest Re chi0 and would take the top of the scale; it is drawn in grey instead. For
    Im chi0/w the Gamma peak is broader: leave out 2-3 mesh steps (--gamma-steps 3);
  * raise low to 50 ... 80: the colours span only the upper half (or fifth) of the values, and
    everything below is the darkest colour;
  * gamma 2 ... 3: of that range, most colours go to the top end, so only the maxima are bright;
  * 10 ... 15 contour levels: discrete bands, the top band is a region of its own.
  e.g.  python3 plot_bz.py chi_6.4521_300_ab re --low 60 --gamma 2.5 --levels 12 --batch
"""
import argparse
import copy
import glob
import os
import re
import sys
from fractions import Fraction

import matplotlib
matplotlib.use('Agg')                                          # PDF only: no display needed
import matplotlib.pyplot as plt                                # noqa: E402
import numpy as np                                             # noqa: E402
from matplotlib import colors                                  # noqa: E402
from matplotlib.cm import ScalarMappable                       # noqa: E402
from matplotlib.collections import Collection                  # noqa: E402
from matplotlib.contour import ContourSet                      # noqa: E402
from matplotlib.patches import Polygon                         # noqa: E402
from matplotlib.path import Path                               # noqa: E402
from matplotlib.transforms import Bbox                          # noqa: E402
from matplotlib.ticker import FuncFormatter, MaxNLocator, MultipleLocator   # noqa: E402

CMAPS = ['inferno', 'magma', 'hot', 'gray', 'plasma', 'viridis', 'cividis', 'turbo']
VIEWS = ['zone', 'cartesian', 'crystal']
CBARS = ['minmax', 'values', 'ticks']
PLANES = ['together', 'panels', '+', '-']                     # the planes +c and -c of chi_bz.x
REGIONS = ['first', 'extended']                               # the first zone, or the whole map
PANELS = ['plane', 'letters', 'none']                          # labels of the two panels
AXCOLOUR = {'inferno': 'cyan', 'magma': 'cyan', 'hot': 'cyan', 'gray': 'red', 'plasma': 'lime',
            'viridis': 'magenta', 'cividis': 'red', 'turbo': 'white'}   # axes in contrast to each colour map
CHI_LABEL = {'re': r'$\chi_{\mathrm{re}}$', 'im': r'$\chi_{\mathrm{im}}$'}   # the axis / colour bar names, no units
GAMMA_NAMES = ('g', 'gamma', 'gm')                             # point labels drawn as the letter Gamma
GUIDE = dict(color='white', lw=0.7, ls=(0, (4, 3)), alpha=0.9, zorder=4)
NICE = [0.01, 0.02, 0.05, 0.1, 0.2, 0.25, 0.5, 1.0, 2.0, 5.0]   # lengths of the scale bar (1/Angstrom)
MAXTICKS = 100                                                 # line plots: most ticks of a set spacing


def blocks(fname):
    """numeric rows of a file, split at blank lines; '#' lines skipped"""
    out, cur = [], []
    for line in open(fname):
        if line.startswith('#'):
            continue
        s = line.strip()
        if not s:
            if cur:
                out.append(np.array(cur))
                cur = []
            continue
        cur.append([float(x) for x in s.split()])
    if cur:
        out.append(np.array(cur))
    return out


class Zone:
    """the first Brillouin zone in the plane: a convex polygon around (0, 0), the projection of Gamma"""

    def __init__(self, poly):
        self.poly = np.asarray(poly)
        closed = np.vstack([self.poly, self.poly[:1]])
        self.path = Path(closed)
        self.grown = Path(closed * (1.0 + 1e-6))               # includes points on the border
        self.corners = self.poly
        self.mids = 0.5 * (self.poly + np.roll(self.poly, -1, axis=0))

    def contains(self, pts):
        return self.grown.contains_points(pts)

    def radius(self, u):
        """distance from Gamma to the border along the unit vector u"""
        best = np.inf
        for a, b in zip(self.poly, np.roll(self.poly, -1, axis=0)):
            d = b - a
            det = d[0] * u[1] - u[0] * d[1]
            if abs(det) < 1e-14:
                continue
            t = (d[0] * a[1] - a[0] * d[1]) / det
            s = (u[0] * a[1] - a[0] * u[1]) / det
            if t > 0.0 and -1e-9 <= s <= 1.0 + 1e-9:
                best = min(best, t)
        return best


def first_zone(zones):
    """the polygon of _bz.dat around Gamma"""
    best, bestd = None, np.inf
    for b in zones:
        poly = b[:-1] if len(b) > 1 and np.allclose(b[0], b[-1]) else b
        if len(poly) < 3:
            continue
        if Path(np.vstack([poly, poly[:1]])).contains_point((0.0, 0.0)):
            return Zone(poly)
        d = np.hypot(*poly.mean(axis=0))
        if d < bestd:
            best, bestd = poly, d
    return Zone(best) if best is not None else None


def periodic_fill(z, hkl, mesh):
    """NaN points (outside the first zone) get the value of the equivalent mesh point"""
    nan = ~np.isfinite(z)
    if not nan.any():
        return z
    g = np.mod(np.rint(hkl * mesh), mesh).astype(np.int64)
    flat = (g[..., 0] * mesh[1] + g[..., 1]) * mesh[2] + g[..., 2]
    table = np.full(int(np.prod(mesh)), np.nan)
    table[flat[~nan]] = z[~nan]
    out = z.copy()
    out[nan] = table[flat[nan]]
    return out


def frac_text(r):
    fr = Fraction(r).limit_denominator(12)
    return str(fr) if abs(float(fr) - r) < 1e-6 else f'{r:g}'


def bz_sections(fname):
    """the polygons of a _bz.dat. For a pair of planes chi_bz.x writes a section per plane when their zone
    outlines differ ('# plane +0.1667', '# plane -0.1667'), else one ('# plane both')"""
    secs, key, cur = {'': []}, '', []
    for line in open(fname):
        if line.startswith('# plane '):
            if cur:
                secs[key].append(np.array(cur))
                cur = []
            w = line.split()[2]
            key = '' if w == 'both' else w[0]
            secs.setdefault(key, [])
            continue
        if line.startswith('#'):
            continue
        t = line.strip()
        if not t:
            if cur:
                secs[key].append(np.array(cur))
                cur = []
            continue
        cur.append([float(x) for x in t.split()])
    if cur:
        secs[key].append(np.array(cur))
    return secs


class Data:
    """the three files of chi_bz.x: one plane, the two planes k = +c and -c on one grid (a pair), or a line.
    Each plane, and for a pair their average, is a view (a copy with its own values and zone)."""

    def __init__(self, prefix, exp):
        self.prefix, self.exp_file = prefix, exp
        self.root, self.sign = self, ''
        mapfile = prefix + '_map.dat'
        self.mapfile = mapfile
        head = [l for l in open(mapfile) if l.startswith('#')]
        self.kind = 'line' if any(l.startswith('# kind: line') for l in head) else 'plane'   # 'path': a cut
        self.extended = any('extended zones' in l for l in head[:4])     # chi_bz.x with zone = 'extended'
        self.pair = None                                       # ('L', '+0.1667', '-0.1667'): the two planes
        for l in head:
            mm = re.match(r'# pair: ([HKL]) ([+-][0-9.]+) ([+-][0-9.]+)', l)
            if mm:
                self.pair = mm.groups()
        self.axes = next(l for l in head if l.startswith('# columns:')).split(':', 1)[1].split()[:2]
        rows = blocks(mapfile)
        self.breaks, self.seg = [], None                       # a cut: no jumps, one segment
        self.axis, self.xunit, self.numeric_labels = None, 'Å$^{-1}$', False
        if self.kind == 'plane':
            g = np.array(rows)                                  # (n1, n2, columns): a structured grid
            self.X, self.Y = g[:, :, 0], g[:, :, 1]
            if self.pair:
                self.vals = {self.pair[1][0]: {'re': g[:, :, 2], 'im': g[:, :, 3]},
                             self.pair[2][0]: {'re': g[:, :, 4], 'im': g[:, :, 5]}}
                self.hkl = g[:, :, 6:9]
            else:
                self.vals = {'': {'re': g[:, :, 2], 'im': g[:, :, 3]}}
                self.hkl = g[:, :, 4:7]
        else:                                                  # a chain (1D system): x  Re  Im  H K L
            d = rows[0]
            c0 = 3
            self.X, self.Y = d[:, 0], None
            self.vals = {'': {'re': d[:, 1], 'im': d[:, 2]}}
            self.hkl = d[:, c0:c0 + 3]
        self.zonesets = bz_sections(prefix + '_bz.dat') if os.path.isfile(prefix + '_bz.dat') else {'': []}
        self.pts, self.pts_hkl = [], []                        # place in the plot, label; H K L
        for line in open(prefix + '_points.dat'):
            if line.startswith('#') or not line.strip():
                continue
            w = line.split()
            if self.kind == 'plane':
                self.pts.append(((float(w[0]), float(w[1])), w[2]))
                self.pts_hkl.append(np.array([float(x) for x in w[3:6]]))
            else:
                self.pts_hkl.append(np.array([float(x) for x in w[2:5]]))
                self.pts.append(((float(w[0]),), w[1]))      # x of the file; other units: line_units
        self.exp = np.loadtxt(exp, ndmin=2) if exp else None
        self.mesh = None
        for line in head:
            mm = re.search(r'mesh\s+(\d+)\s+(\d+)\s+(\d+)', line)
            if mm:
                self.mesh = np.array([int(x) for x in mm.groups()])
                break
        self.hkls = {k: self.hkl for k in self.vals}
        if self.kind == 'plane':
            # the crystal axes in the plane: (x, y) = H_a ba + H_b bb + o2 (o2: the point H_a = H_b = 0)
            fixed = [k for k in range(3) if np.ptp(self.hkl[..., k]) < 1e-9]
            self.ic = fixed[0] if fixed else 2
            self.ia, self.ib = [k for k in range(3) if k != self.ic]
            A = np.column_stack([self.hkl[..., self.ia].ravel(), self.hkl[..., self.ib].ravel(), np.ones(self.X.size)])
            self.ba, self.bb, self.o2 = np.linalg.lstsq(A, np.column_stack([self.X.ravel(), self.Y.ravel()]),
                                                        rcond=None)[0]
            if self.pair:                                      # the second plane: the same H_a, H_b at -c
                h2 = self.hkl.copy()
                h2[..., self.ic] = -self.hkl[..., self.ic]
                self.hkls[self.pair[2][0]] = h2
            if self.mesh is not None:                          # outside the first zone: the equivalent point
                for k, v in self.vals.items():
                    for part in ('re', 'im'):
                        v[part] = periodic_fill(v[part], self.hkls[k], self.mesh)
        if self.kind != 'plane':
            self.line_units()
        self.make_views()

    def label(self, key):
        """L = +1/6, L = -1/6, L = ±1/6 (average); '' for a single plane"""
        if not self.pair:
            return ''
        v = abs(float(self.pair[1]))
        fr = Fraction(v).limit_denominator(24)
        txt = str(fr) if abs(float(fr) - v) < 6e-5 else f'{v:g}'
        return {'+': f'{self.pair[0]} = +{txt}', '-': f'{self.pair[0]} = −{txt}',
                'avg': f'{self.pair[0]} = ±{txt} (average)'}[key]

    def line_units(self):
        """the x of a chain (1D system of chi_bz.x) in each unit: self.xs / self.pxs {unit: x of the data / of
        the labelled points}; 'invA' (q in 1/Angstrom) and 'rlu' (the crystal coordinate along the chain)"""
        xf = np.asarray(self.X, float)
        pf = np.array([x for (x,), _ in self.pts]) if self.pts else np.zeros(0)
        col = int(np.argmax(np.ptp(self.hkl, axis=0)))
        scale = np.ptp(xf) / np.ptp(self.hkl[:, col]) if np.ptp(self.hkl[:, col]) > 0 else 1.0
        self.xs = {'invA': xf, 'rlu': self.hkl[:, col]}
        self.pxs = {'invA': pf, 'rlu': pf / scale}
        self.xnames = {'invA': 'q (Å$^{-1}$)', 'rlu': f'{"HKL"[col]} (r.l.u.)'}
        self.default_unit = 'invA'
        self.units = list(self.xs)

    def make_views(self):
        self.views = {k: self._view(k) for k in list(self.vals) + (['avg'] if self.pair else [])}
        first = self.views['avg' if self.pair else '']         # the root shows the default view
        for att in ('val', 'zones', 'zone', 'inside', '_mask'):
            setattr(self, att, getattr(first, att))

    def _view(self, key):
        v = copy.copy(self)
        v.sign = key
        if key == 'avg':
            v.val = {p: 0.5 * (self.vals['+'][p] + self.vals['-'][p]) for p in ('re', 'im')}
            v.hkl = self.hkls['+']
            zones = self.zonesets.get('+') or self.zonesets.get('', [])
        else:
            v.val = self.vals[key]
            v.hkl = self.hkls[key]
            zones = self.zonesets.get(key) or self.zonesets.get('', [])
        v.zones, v._mask, v.zone, v.inside = zones, {}, None, None
        if self.kind == 'plane':
            v.zone = first_zone(zones)
            if v.zone is not None:
                pts = np.column_stack([self.X.ravel(), self.Y.ravel()])
                v.inside = v.zone.contains(pts).reshape(self.X.shape)
        return v

    def gamma_mask(self, steps):
        """Gamma and its equivalents (integer H K L) and their neighbours up to `steps` mesh steps"""
        if steps not in self._mask:
            m = np.zeros(self.hkl.shape[:-1], bool)
            if steps > 0:
                isg = np.all(np.abs(self.hkl - np.rint(self.hkl)) < 1e-6, axis=-1)
                for idx in zip(*np.nonzero(isg)):
                    m[tuple(slice(max(i - steps, 0), i + steps + 1) for i in idx)] = True
            self._mask[steps] = m
        return self._mask[steps]


def resolve(prefix):
    """the files of a prefix: <prefix>_map.dat of chi_bz.x, or the pair file <prefix>_L+-0.1667_map.dat (also
    when one plane of it, <prefix>_L+0.1667, or just <prefix> is named). None: no such files; SystemExit
    when the prefix fits several pairs"""
    if os.path.isfile(prefix + '_map.dat'):
        return prefix
    m = re.fullmatch(r'(.*_[HKL])[+-](\d+(?:\.\d*)?)', prefix)
    if m and os.path.isfile(f'{m.group(1)}+-{m.group(2)}_map.dat'):
        return f'{m.group(1)}+-{m.group(2)}'
    found = sorted(glob.glob(glob.escape(prefix) + '_[HKL]+-*_map.dat'))
    if len(found) == 1:
        return found[0][:-len('_map.dat')]
    if found:
        sys.exit(f'{prefix}: several pairs of planes: ' + ', '.join(f[:-len('_map.dat')] for f in found) +
                 '; give one of them')
    return None


def datasets(stem):
    """the files of chi_bz.x in the folder of stem: [('1d' (a chain) or '2d' (a plane), prefix)], those whose
    names start with stem (e.g. the chi file chi_6.4521_300 or chi_6.4521_300.dat; '' for all)"""
    stem = stem[:-4] if stem.endswith('.dat') else stem
    out = []
    for pt in sorted(glob.glob(os.path.join(glob.escape(os.path.dirname(stem)), '*_points.dat'))):
        base = pt[:-len('_points.dat')]
        if not base.startswith(stem):
            continue
        try:
            if os.path.isfile(base + '_map.dat'):
                with open(base + '_map.dat') as fh:
                    out.append(('1d' if fh.readline().startswith('# kind: line') else '2d', base))
        except OSError:
            continue
    return out


def load(prefix, exp):
    """the plot data of a prefix (see resolve)"""
    base = resolve(prefix)
    if base is None:
        sys.exit(f'no {prefix}_map.dat (chi_bz.x) in this folder')
    return Data(base, exp)


def choose(prefix):
    """the file, for a prefix that is not one plot (a chi file, or nothing): the maps (and 1D chains) of
    chi_bz.x in the folder"""
    sets = datasets(prefix)
    if not sets:
        sys.exit(f'no files of chi_bz.x (*_map.dat){" starting with " + prefix if prefix else ""} in this folder')
    files = [b for _, b in sets]
    note = {b: '  (a 1D chain)' if k == '1d' else '' for k, b in sets}
    if len(files) == 1:
        print(f'file: {files[0]}{note[files[0]]}')
        return files[0]
    print('files')
    for i, b in enumerate(files, 1):
        print(f'  {i:2d}  {b}{note[b]}')
    while True:
        try:
            a = input('  number [1]: ').strip()
        except EOFError:
            print()
            sys.exit(0)
        if a == '':
            return files[0]
        if a.isdigit() and 1 <= int(a) <= len(files):
            return files[int(a) - 1]
        print(f'    1 ... {len(files)}')


def parse_cut(data, a):
    """a cut through a map at constant x or y of the map, in r.l.u.: 'x=0.25', 'y=-1/3' (x: the first in-plane
    crystal coordinate, H of an ab map; y: the second, K); any value. '-' or 'none' for no cut.
    Returns ('x' or 'y', value) or None"""
    a = a.strip()
    if a.lower() in ('-', 'none', 'no', 'n', ''):
        return None
    m = re.fullmatch(r'([xyXY])\s*=\s*([-+]?[0-9./]+)', a)
    if not m:
        raise ValueError('x=value or y=value in r.l.u., e.g. y=0.25 or x=-1/3, or - for none')
    try:
        v = float(Fraction(m.group(2)))
    except (ValueError, ZeroDivisionError):
        raise ValueError(f'{m.group(2)}: a number such as 0.25 or 1/4')
    return (m.group(1).lower(), v)


def cut_text(cut):
    return f'{cut[0]}={cut[1]:g}'


def cut_axes(root, cut):
    """the crystal coordinates of a cut: (fixed, running)"""
    return (root.ia, root.ib) if cut[0] == 'x' else (root.ib, root.ia)


def parse_range(a):
    """the x range of a cut: 'from to' in r.l.u. (e.g. -1 1, or -1/2 1/2); 'auto' for -0.5 ... 0.5"""
    if a.strip().lower() in ('auto', 'a', '', '-'):
        return None
    w = a.replace(',', ' ').split()
    if len(w) != 2:
        raise ValueError('two numbers: from to (e.g. -1 1), or auto')
    lo, hi = (float(Fraction(x)) for x in w)
    if not hi > lo:
        raise ValueError('from < to')
    return (lo, hi)


def cut_ends(data, cut):
    """two far points of the cut line in the plane of the map (for drawing it over the map)"""
    big = 10.0 * (np.ptp(data.X) + np.ptp(data.Y) + 1.0)
    fixed, run = (data.ba, data.bb) if cut[0] == 'x' else (data.bb, data.ba)
    run = run / np.linalg.norm(run)
    p0 = cut[1] * fixed + data.o2
    return np.array([p0 - big * run, p0 + big * run])


def build_cut(data, s, cut, xrange=None):
    """the 1D line of a cut through the map of chi_bz.x at constant x (or y) in r.l.u.: a Data of kind 'path'.
    The values come from the whole k-mesh (chi0 is periodic), so the cut can run over any range of the other
    coordinate (xrange, r.l.u.; default -0.5 ... 0.5) at its mesh points. On a mesh line the values are
    those of the mesh; between two mesh lines they are interpolated linearly between them.
    Also writes the values to <prefix>_cut_<y0.25>.dat"""
    root = data.root
    if root.mesh is None:
        raise ValueError('the map has no mesh line in its header (chi_bz.x of an older version?)')
    L, v = cut
    i, j = cut_axes(root, cut)
    ni, nj = int(root.mesh[i]), int(root.mesh[j])
    first = list(root.vals)[0]
    base = root.hkls[first]
    gi = np.rint(base[..., i] * ni).astype(int) % ni           # the mesh point of each grid point
    gj = np.rint(base[..., j] * nj).astype(int) % nj
    tables = {}
    for k, vv in root.vals.items():                            # chi0 on the mesh of the plane (i, j)
        tables[k] = {}
        for p in ('re', 'im'):
            T = np.full((ni, nj), np.nan)
            z = vv[p]
            fin = np.isfinite(z)
            T[gi[fin], gj[fin]] = z[fin]
            tables[k][p] = T
    a, b = xrange if xrange else (-0.5, 0.5)
    m = np.arange(int(np.ceil(a * nj - 1e-9)), int(np.floor(b * nj + 1e-9)) + 1)
    if m.size < 2:
        raise ValueError(f'the range {a:g} ... {b:g} holds fewer than two mesh points (step 1/{nj})')
    g = v * ni
    lo = int(np.floor(g + 1e-9))
    w = g - lo if abs(g - round(g)) > 1e-6 else 0.0
    lo = int(round(g)) if w == 0.0 else lo
    vals, missing = {}, False
    for k in tables:
        vals[k] = {}
        for p in ('re', 'im'):
            T = tables[k][p]
            z = T[lo % ni, m % nj]
            if w:
                z = (1.0 - w) * z + w * T[(lo + 1) % ni, m % nj]
            missing |= bool(np.isnan(z).any())
            vals[k][p] = z
    names = 'HKL'
    note = ''
    if w:
        note = f'{L} = {v:g} ({names[i]}) lies between the mesh lines {names[i]} = {lo / ni:.6g} and ' \
               f'{(lo + 1) / ni:.6g}: values interpolated linearly between them (weight {w:.3f})'
    if missing:
        note += ('; ' if note else '') + 'some mesh points are not in the map (blank in the plot): ' \
            'chi_bz.x with zone = \'extended\' gives all of them'
    run = m / nj
    hkl = np.repeat(base.reshape(-1, 3)[:1], len(m), axis=0).astype(float)
    hkl[:, i], hkl[:, j] = v, run
    blen = np.linalg.norm(root.bb if j == root.ib else root.ba)
    xs = {'rlu': run, 'invA': run * blen}
    xnames = {'rlu': f'{"y" if L == "x" else "x"} = {names[j]} (r.l.u.)',
              'invA': f'q along {"abc"[j]}* (Å$^{{-1}}$)'}
    xnames['rlu'] = f'{names[j]} (r.l.u.)'
    pts, pts_hkl = [], []                                      # labelled points on the line, and their images
    for (_, name), q in zip(root.pts, root.pts_hkl):
        if abs(q[i] - v) > 1e-6 and abs((q[i] - v) - round(q[i] - v)) > 1e-6:
            continue
        for n in range(int(np.floor(a - q[j])) - 1, int(np.ceil(b - q[j])) + 2):
            x = q[j] + n
            if a - 1e-9 <= x <= b + 1e-9:
                h = q.copy()
                h[j], h[i] = x, v
                pts.append((x, name))
                pts_hkl.append(h)
    cd = Data.__new__(Data)
    tag = f'{L}{v:g}'
    cd.prefix = f'{root.prefix}_cut_{tag}'
    cd.map_prefix, cd.cut, cd.cut_range = root.prefix, (L, v), xrange
    cd.mapfile = cd.prefix + '.dat'
    cd.exp_file, cd.exp = None, None
    cd.root, cd.sign, cd.kind, cd.pair, cd.extended = cd, '', 'path', root.pair, False
    cd.axis, cd.xunit, cd.numeric_labels = None, '', False
    cd.breaks, cd.seg = [], np.zeros(len(hkl), int)
    cd.vals, cd.hkl, cd.mesh, cd.axes = vals, hkl, root.mesh, root.axes
    cd.hkls = {k: root_line_hkl(root, k, hkl) for k in vals}
    cd.pts = [((x,), name) for x, name in pts]
    cd.pts_hkl = pts_hkl
    cd.xs, cd.xnames, cd.units, cd.default_unit = xs, xnames, list(xs), 'rlu'
    cd.pxs = {'rlu': np.array([x for x, _ in pts]), 'invA': np.array([x * blen for x, _ in pts])}
    cd.zonesets, cd.X, cd.Y = {'': []}, xs['rlu'], None
    cd.note, cd.cut_axes, cd.cut_source = note, (i, j), data
    cd.make_views()
    with open(cd.mapfile, 'w') as fh:                          # the values of the cut, for other plots
        fh.write(f'# cut {L} = {v:g} ({names[i]} = {v:g}) through {root.mapfile}, {names[j]} from {run[0]:g} '
                 f'to {run[-1]:g}' + (f' ({note})' if note else '') + '\n')
        if root.pair:
            fh.write(f'# pair: {" ".join(root.pair)}\n')
        vnames = ['Re chi0', 'Im chi0/w'] if not root.pair else \
            sum([[f'Re chi0 ({root.label(k)})', f'Im chi0/w ({root.label(k)})'] for k in vals], [])
        fh.write('# columns: ' + '  '.join(['x [r.l.u.]', 'x [1/Angstrom]'] + vnames + ['H', 'K', 'L [r.l.u.]']) +
                 '\n')
        for r in range(len(hkl)):
            row = [xs['rlu'][r], xs['invA'][r]] + sum([[vals[k]['re'][r], vals[k]['im'][r]] for k in vals], []) + \
                list(hkl[r])
            fh.write('  '.join(f'{x: .8e}' if abs(x) and (abs(x) < 1e-3 or abs(x) > 1e4) else f'{x: .8f}'
                               for x in row) + '\n')
    return cd


def root_line_hkl(root, key, hkl):
    """H K L of the points of a cut on the plane of key (the second plane of a pair: the fixed coordinate at -c)"""
    if root.pair and key == root.pair[2][0]:
        h2 = hkl.copy()
        h2[:, root.ic] = -h2[:, root.ic]
        return h2
    return hkl


def chosen(data, s):
    """the views to draw: the plane, or for a pair the average, one plane, or both side by side"""
    root = data.root
    if root.pair is None:
        return [root.views['']]
    return {'together': [root.views['avg']], '+': [root.views['+']], '-': [root.views['-']],
            'panels': [root.views['+'], root.views['-']]}[s.planes]


def extended(data, s):
    """the whole map with all zones (region 'extended', or no zone in the _bz.dat), not the first zone"""
    return data.kind == 'plane' and (s.region == 'extended' or data.zone is None)


def window(data):
    """the corners of the map: a parallelogram in the plane"""
    X, Y = data.X, data.Y
    return np.array([[X[0, 0], Y[0, 0]], [X[-1, 0], Y[-1, 0]], [X[-1, -1], Y[-1, -1]], [X[0, -1], Y[0, -1]]])


def shown(data, s):
    """the points inside the drawn region: the first zone, or the whole map"""
    if data.kind == 'plane' and data.inside is not None and not extended(data, s):
        return data.inside
    return None


def limits(data, s):
    """bottom and top of the colour scale (common to the planes drawn together), and the values they
    are taken from"""
    vals = []
    for dd in chosen(data, s):
        z = dd.val[s.part]
        ok = np.isfinite(z)
        sel = shown(dd, s)
        if sel is not None:
            ok &= sel
        if s.ignore:
            ok2 = ok & ~dd.gamma_mask(max(s.gsteps, 1))
            if ok2.any():
                ok = ok2
        vals.append(z[ok])
    d = np.concatenate(vals)
    lo = s.vmin if s.vmin is not None else float(np.percentile(d, s.low))
    hi = s.vmax if s.vmax is not None else float(np.percentile(d, s.high))
    if not hi > lo:
        hi = lo + 1e-9 * max(1.0, abs(lo))
    return lo, hi, d


def yrange(data, s):
    """bottom and top of the value axis of a line plot (common to the lines drawn): all the data with a
    small margin, or the values given"""
    z = np.concatenate([d.val[s.part] for d in chosen(data, s)])
    z = z[np.isfinite(z)]
    lo, hi = float(z.min()), float(z.max())
    pad = 0.04 * (hi - lo) if hi > lo else 0.04 * max(abs(hi), 1e-9)
    lo = s.ymin if s.ymin is not None else lo - pad
    hi = s.ymax if s.ymax is not None else hi + pad
    if not hi > lo:
        hi = lo + max(abs(lo), 1e-9)
    return lo, hi


def step_format(step, frac, every=1):
    """the labels of the ticks at the multiples of step; every n-th tick is labelled"""
    def text(v, _):
        j = int(round(v / step))
        if j % every:
            return ''
        return f'{round(j * step, 4):g}'.replace('-', '\u2212')                   # decimals, also for 1/6
    return FuncFormatter(text)


def spaced(axis, step, lo, hi, frac=False):
    """ticks at the multiples of step on an axis from lo to hi; False (the ticks stay automatic) without
    a step, or when there would be more than MAXTICKS of them"""
    if not step or (hi - lo) / step > MAXTICKS:
        return False
    axis.set_major_locator(MultipleLocator(step))
    axis.set_major_formatter(step_format(step, frac))
    return True


def thin(fig, axis, step, frac):
    """an axis with a set tick spacing: label every n-th tick when the labels would touch (the tick marks stay)"""
    if not step or not isinstance(axis.get_major_locator(), MultipleLocator):
        return
    r = fig.canvas.get_renderer()
    lo, hi = sorted(axis.get_view_interval())
    labs = sorted((int(round(tk.get_loc() / step)), tk.label1.get_window_extent(r)) for tk in axis.get_major_ticks()
                  if lo - 1e-9 <= tk.get_loc() <= hi + 1e-9 and tk.label1.get_text() and tk.label1.get_visible())
    if len(labs) < 2:
        return
    gap = 0.25 * axis.get_major_ticks()[0].label1.get_fontsize() * fig.dpi / 72.0
    for n in ((1, 2, 3, 4, 6, 12, 24) if frac else (1, 2, 5, 10, 20, 50)):
        on = [b.padded(gap) for j, b in labs if j % n == 0]
        if not any(a.overlaps(b) for a, b in zip(on, on[1:])):
            break
    if n > 1:
        axis.set_major_formatter(step_format(step, frac, n))


def make_cmap(name):
    cm = (matplotlib.colormaps[name] if hasattr(matplotlib, 'colormaps') else plt.get_cmap(name)).copy()
    cm.set_over('tab:red' if name == 'gray' else '0.6')        # above the scale (Gamma): grey
    cm.set_under(cm(0.0))
    cm.set_bad(alpha=0.0)                                      # NaN (outside the first zone): blank
    return cm


def parts(artist):
    """the collections of a plot (a ContourSet is one in matplotlib >= 3.8)"""
    return [artist] if isinstance(artist, Collection) else list(getattr(artist, 'collections', [artist]))


def zone_labels(data):
    """labels of the input points that sit on a corner or an edge centre of the zone"""
    zn, lc, lm = data.zone, '', ''
    size = np.ptp(zn.poly, axis=0).max()
    for (x, y), name in data.pts:
        q = np.array([x, y])
        if np.min(np.hypot(*(zn.corners - q).T)) < 1e-3 * size:
            lc = lc or name
        if np.min(np.hypot(*(zn.mids - q).T)) < 1e-3 * size:
            lm = lm or name
    return lc, lm


def place_text(data, p):
    """where p lies in the zone: fraction of the way from Gamma to the border, direction"""
    r = float(np.hypot(*p))
    if data.zone is None or r < 1e-12:
        return f'q = ({p[0]:.3f}, {p[1]:.3f}) Å⁻¹'
    zn = data.zone
    u = p / r
    frac = r / zn.radius(u)

    def angle(v):
        return np.degrees(np.arccos(np.clip(np.dot(u, v / np.hypot(*v)), -1.0, 1.0)))
    lc, lm = zone_labels(data)
    toward = ''
    if min(angle(v) for v in zn.mids) < 2.0:
        toward = ' toward an edge centre' + (f' ({lm})' if lm else '')
    elif min(angle(v) for v in zn.corners) < 2.0:
        toward = ' toward a corner' + (f' ({lc})' if lc else '')
    return f'{frac:.2f} of the way from Γ to the zone border{toward}, |q| = {r:.3f} Å⁻¹'


def mark_max(ax, data, s, z, clip=None):
    """the largest value (away from Gamma): its description; with s.mark every equivalent point is marked
    by a bold red +. The text describes the one nearest Gamma (then top, right)"""
    cand = np.isfinite(z)
    sel = shown(data, s)
    if sel is not None:
        cand &= sel
    if s.ignore:
        c2 = cand & ~data.gamma_mask(max(s.gsteps, 1))
        cand = c2 if c2.any() else cand
    if not cand.any():
        return ''
    v = z[cand].max()
    tie = cand & (z >= v - 1e-9 * abs(v))
    xs, ys, hk = data.X[tie], data.Y[tie], data.hkl[tie]
    k = np.lexsort((-xs, -ys, np.round(np.hypot(xs, ys), 9)))[0]   # nearest Gamma, then top, then right
    nl = 0
    if s.mark:
        ax.plot(xs, ys, '+', ms=16, mew=3.0, color='red', zorder=9)   # a bold + on every equivalent point
        nl = len(xs)
    name = 'Re χ₀' if s.part == 're' else 'Im χ₀/ω'
    marked = '; marked by a red +' if nl else ''
    return (f'max {name} = {v:.4g} ({tie.sum()} equivalent point{"s" if tie.sum() > 1 else ""}{marked})'
            f'{" away from Γ" if s.ignore else ""}\n'
            f'at {place_text(data, np.array([xs[k], ys[k]]))};  H K L = '
            + ' '.join(f'{x:.4f}' for x in hk[k]))


def axis_colour(s):
    """colour of the axes (a contrast to the colour map, or the user's choice) and of its outline"""
    c = AXCOLOUR.get(s.cmap, 'white') if s.axcolour == 'auto' else s.axcolour
    r, g, b = colors.to_rgb(c)
    return c, ('black' if 0.2126 * r + 0.7152 * g + 0.0722 * b > 0.4 else 'white')


def ray(zn, p0, u):
    """distance from p0 (inside the zone) to the border along the unit vector u"""
    best = np.inf
    for a, b in zip(zn.poly, np.roll(zn.poly, -1, axis=0)):
        d, a = b - a, a - p0
        det = d[0] * u[1] - u[0] * d[1]
        if abs(det) < 1e-14:
            continue
        t = (d[0] * a[1] - a[0] * d[1]) / det
        sp = (u[0] * a[1] - a[0] * u[1]) / det
        if t > 0.0 and -1e-9 <= sp <= 1.0 + 1e-9:
            best = min(best, t)
    return best if np.isfinite(best) else 0.0


def auto_tick(length):
    """a round tick spacing with at least three ticks on the shortest half axis"""
    cands = [m * 10.0 ** e for e in range(-4, 3) for m in (1.0, 2.0, 2.5, 5.0)]
    good = [c for c in cands if c <= length / 3.0 + 1e-12]
    return max(good) if good else min(cands)


def tick_text(v, frac):
    if abs(v) < 1e-12:
        return '0'
    return f'{round(v, 4):g}'.replace('-', '\u2212')                          # decimals, also for r.l.u.


def anchor(v):
    """alignment of a text placed away from its point in the direction v"""
    ha = 'left' if v[0] > 0.35 else 'right' if v[0] < -0.35 else 'center'
    va = 'bottom' if v[1] > 0.35 else 'top' if v[1] < -0.35 else 'center'
    return ha, va


def axis_name(name):
    """qx -> k_x, q1 -> k_1"""
    return f'$k_{name[1:]}$' if name.startswith('q') and len(name) > 1 else name


def cbar_ticks(lo, hi, gamma):
    """colour-bar ticks: round values for a linear scale; with gamma, evenly spaced along the bar"""
    if abs(gamma - 1.0) < 1e-9:
        v = MaxNLocator(6).tick_values(lo, hi)
    else:
        v = lo + (hi - lo) * np.linspace(0.0, 1.0, 6) ** (1.0 / gamma)
        dec = max(0, int(-np.floor(np.log10((hi - lo) / 5.0))) + 1)
        v = np.round(v, dec)
        v[0] = np.ceil(lo * 10 ** dec) / 10 ** dec
        v[-1] = np.floor(hi * 10 ** dec) / 10 ** dec
    v = np.unique(v[(v >= lo - 1e-12 * abs(lo)) & (v <= hi + 1e-12 * abs(hi))])
    for dec in range(0, 8):                                    # as many decimals as the ticks need
        if np.all(np.abs(np.round(v, dec) - v) < 1e-9 * np.maximum(1.0, np.abs(v))):
            break
    return v, [f'{x:.{dec}f}' for x in v]


def zone_extras(ax, data, s, f):
    """the border; on request fractions of Gamma -> border, guide lines, a length bar.
    Returns the points that must stay in view."""
    zn = data.zone
    P = np.vstack([zn.poly, zn.poly[:1]])
    ax.plot(P[:, 0], P[:, 1], color='black', lw=1.3, zorder=5)
    if s.guides:
        for q in np.vstack([zn.corners, zn.mids]):
            ax.plot([0.0, q[0]], [0.0, q[1]], **GUIDE)
    if s.rings:
        top = zn.mids[np.lexsort((zn.mids[:, 0], zn.mids[:, 1]))[-1]]    # edge centre at the top
        for r in s.rings:
            ax.plot(r * P[:, 0], r * P[:, 1], **GUIDE)
            ax.annotate(frac_text(r), r * top, xytext=(3, 2), textcoords='offset points', fontsize=0.75 * f,
                        color='white', zorder=6)
    box = [zn.poly]
    (x0, y0), (x1, y1) = zn.poly.min(axis=0), zn.poly.max(axis=0)
    w, h = x1 - x0, y1 - y0
    if s.bar or s.rings:
        yb = y0 - 0.14 * h
        if s.bar:
            bar = max([v for v in NICE if v <= 0.3 * w] or NICE[:1])
            ax.plot([x0, x0 + bar], [yb, yb], color='black', lw=2.0, solid_capstyle='butt')
            ax.text(x0 + 0.5 * bar, yb - 0.02 * h, f'{bar:g} Å⁻¹', ha='center', va='top', fontsize=0.8 * f)
        if s.rings:
            ax.text(x1, yb, 'dashed: ' + ', '.join(frac_text(r) for r in s.rings) + ' of the way Γ → zone border',
                    ha='right', va='center', fontsize=0.7 * f)
        box.append(np.array([[x0, yb - 0.06 * h]]))
    return box


def draw_axes(ax, data, s, f, zn=None):
    """axes through Gamma from border to border of the zone (or of the region zn), with ticks: cartesian
    (k_x, k_y, ... in 1/A) or along the two reciprocal axes of the plane (H, K, L in r.l.u.).
    Returns the points to keep in view, the lines (point, direction, side of the tick labels), the tick
    labels of each axis as (index, text) and the axis names."""
    zn = zn if zn is not None else data.zone
    col, _ = axis_colour(s)
    size = np.ptp(zn.poly, axis=0).max()
    if s.view == 'cartesian':
        p0 = np.zeros(2)                                       # Gamma
        per = [1.0, 1.0]                                       # 1/A per unit of the tick labels
        dirs = [np.array([1.0, 0.0]), np.array([0.0, 1.0])]
        names = [axis_name(data.axes[0]) + ' (Å$^{-1}$)', axis_name(data.axes[1]) + ' (Å$^{-1}$)']
    else:
        p0 = np.asarray(data.o2)                               # H_a = H_b = 0 (Gamma for a plane through it)
        per = [float(np.hypot(*data.ba)), float(np.hypot(*data.bb))]
        dirs = [np.asarray(data.ba) / per[0], np.asarray(data.bb) / per[1]]
        names = ['HKL'[data.ia] + ' (r.l.u.)', 'HKL'[data.ib] + ' (r.l.u.)']
    ext = [(ray(zn, p0, u), ray(zn, p0, -u)) for u in dirs]
    tick = s.tick or auto_tick(min(min(e) / q for e, q in zip(ext, per)))
    pts, lines, groups, titles = [], [], [], []
    for k, (u, q, (rp, rm), name) in enumerate(zip(dirs, per, ext, names)):
        group = []
        a, b = p0 - rm * u, p0 + rp * u                        # from border to border
        side = np.array([u[1], -u[0]]) if k == 0 else np.array([-u[1], u[0]])   # tick labels below / left
        ax.plot([a[0], b[0]], [a[1], b[1]], color=col, lw=1.1, zorder=6,
                solid_capstyle='butt')
        hs = 0.012 * size * side
        ha, va = anchor(side)
        edge = 1e-6 * size                                     # no tick on the border itself
        for j in range(-int(np.floor((rm - edge) / (tick * q))), int(np.floor((rp - edge) / (tick * q))) + 1):
            if j == 0:
                continue
            t = p0 + j * tick * q * u
            ax.plot([t[0] - hs[0], t[0] + hs[0]], [t[1] - hs[1], t[1] + hs[1]], color=col, lw=1.1, zorder=6)
            group.append((j, ax.annotate(tick_text(j * tick, s.tick_frac), t, xytext=tuple((0.4 * f + 2) * side),
                                         textcoords='offset points', ha=ha, va=va, fontsize=0.7 * f, color=col, zorder=7,
                                         fontweight='bold')))
        if s.minor and s.minor < tick and (rp + rm) / (s.minor * q) <= 4 * MAXTICKS:   # minor: short, no label
            for j in range(-int(np.floor((rm - edge) / (s.minor * q))), int(np.floor((rp - edge) / (s.minor * q))) + 1):
                if j == 0 or abs(j * s.minor / tick - round(j * s.minor / tick)) < 1e-6:
                    continue                                   # Gamma and the major ticks
                t = p0 + j * s.minor * q * u
                ax.plot([t[0] - 0.5 * hs[0], t[0] + 0.5 * hs[0]], [t[1] - 0.5 * hs[1], t[1] + 0.5 * hs[1]],
                        color=col, lw=0.9, zorder=6)
        pts.append(np.array([a, b]))
        lines.append((p0, u, side))
        groups.append(group)
    return pts, lines, groups, titles


def point_places(data, s):
    """where the labelled points are drawn: in the crystal view at their own H, K on the plane (so that
    Gamma is where the H, K axes cross, also on a plane L = c /= 0 of an oblique cell); otherwise at the
    projection of the point onto the plane, as in _points.dat"""
    if s.view == 'crystal' and data.zone is not None and len(data.pts_hkl) == len(data.pts):
        return [(h[data.ia] * np.asarray(data.ba) + h[data.ib] * np.asarray(data.bb) + np.asarray(data.o2), name)
                for h, (_, name) in zip(data.pts_hkl, data.pts)]
    return [(np.array(xy, dtype=float), name) for xy, name in data.pts]


def put_label(ax, data, s, f, p, name, lines, on_map=False):
    """a labelled point: the label outside the zone for points on the border, off the axis lines"""
    zn = data.zone
    text = 'Γ' if name.lower() in GAMMA_NAMES else name
    col, _ = axis_colour(s)
    if zn is None:
        ax.plot(*p, 'o', ms=4, mfc='white', mew=0, zorder=8)
        return ax.annotate(text, p, xytext=(4, 4), textcoords='offset points', color='white', fontsize=f,
                           zorder=8, fontweight='bold')
    size = np.ptp(zn.poly, axis=0).max()
    r = float(np.hypot(*p))
    u = p / r if r > 1e-6 * size else np.array([1.0, 1.0]) / np.sqrt(2.0)
    inner = on_map or r < 1e-6 * size or r < 0.97 * zn.radius(u)
    for p0, w, side in lines:                                  # on an axis: to the side away from the ticks
        d = p - p0
        if abs(d[0] * w[1] - d[1] * w[0]) < 1e-3 * size:
            v = (w if np.dot(d, w) >= -1e-9 * size else -w) - side
            u = v / np.hypot(*v)
            break
    ax.plot(*p, 'o', ms=4, mfc=col if inner else 'black', mew=0, zorder=8)
    style = dict(color=col) if inner else dict(color='black')
    ha, va = anchor(u)
    return ax.annotate(text, p, xytext=tuple((0.3 * f + 3) * u), textcoords='offset points', ha=ha, va=va,
                       fontsize=f, zorder=8, fontweight='bold', **style)


def draw(ax, cax, data, s):
    """the map of a plane for the settings s; returns the colour range and the place of the maximum"""
    ax.cla()
    f = s.font
    z = data.val[s.part]
    lo, hi = limits(data, s)[:2]
    zlabel = CHI_LABEL[s.part]
    info = ''

    if cax is not None:
        cax.cla()
    zn = data.zone
    X, Y, Z = data.X, data.Y, np.ma.masked_invalid(z)
    cmap = make_cmap(s.cmap)
    norm = colors.PowerNorm(gamma=s.gamma, vmin=lo, vmax=hi)
    levels = lo + (hi - lo) * np.linspace(0.0, 1.0, s.levels + 1) ** (1.0 / s.gamma)   # equal colour steps
    m = ax.contourf(X, Y, Z, levels=levels, cmap=cmap, norm=norm, extend='both')   # always filled contours
    drawn = parts(m)
    if s.lines:
        drawn += parts(ax.contour(X, Y, Z, levels=levels[1:-1], colors='black', linewidths=0.5, alpha=0.6))
    ext = extended(data, s)
    win = window(data)
    clip = Polygon(win if ext else zn.poly, closed=True, transform=ax.transData, fc='none', ec='none')
    ax.add_patch(clip)                                     # the map ends at the zone border / the window
    for c in drawn:
        c.set_clip_path(clip)
    for c in drawn:                                        # the map as an image in the PDF; contours
        if not isinstance(c, ContourSet):                  # stay vector shapes (matplotlib cannot
            c.set_rasterized(True)                         # rasterize them)
    if cax is not None:
        sm = ScalarMappable(norm=norm, cmap=cmap)           # a plain bar: no pointed ends
        cb = ax.figure.colorbar(sm, cax=cax, extend='neither',
                                boundaries=levels)
        cb.outline.set_linewidth(0.6)
        cb.minorticks_off()
        if s.cbar == 'ticks':
            tv, tl = cbar_ticks(lo, hi, s.gamma)
            cb.set_ticks(tv)
            cb.set_ticklabels(tl)
            cb.set_label(zlabel, size=f)
            cb.ax.tick_params(labelsize=0.85 * f)
        else:
            cb.set_ticks([lo, hi])
            cb.set_ticklabels(['min', 'max'] if s.cbar == 'minmax' else [f'{lo:.4g}', f'{hi:.4g}'])
            cb.ax.tick_params(labelsize=0.85 * f, length=0, pad=3)
    lines, ax.tick_groups, ax.keep_texts = [], [], []
    places = point_places(data, s)
    ax.keep_points = [xy for xy, _ in places] if s.labels else []
    if not ext:                                            # the first zone
        box = zone_extras(ax, data, s, f)
        if s.view in ('cartesian', 'crystal'):
            pts, lines, ax.tick_groups, ax.keep_texts = draw_axes(ax, data, s, f)
            box += pts
        size = np.ptp(zn.poly, axis=0).max()
    else:                                                  # the whole map with every zone in it
        box = [win]
        Pw = np.vstack([win, win[:1]])
        ax.plot(Pw[:, 0], Pw[:, 1], color='black', lw=1.0, zorder=5)
        if s.zones:
            for b in data.zones:
                first = len(b) > 2 and Path(b).contains_point((0.0, 0.0))   # the zone around Gamma
                line, = ax.plot(b[:, 0], b[:, 1], color='white', lw=1.6 if first else 0.8, zorder=5)
                line.set_clip_path(clip)
        if s.view in ('cartesian', 'crystal'):
            pts, lines, ax.tick_groups, ax.keep_texts = draw_axes(ax, data, s, f, Zone(win))
            box += pts
        size = np.ptp(win, axis=0).max()
    if s.cut:                                                  # the line of the 1D cut, over the map:
        e = cut_ends(data, s.cut)                              # solid, in the colour of the axes; not drawn
        d = (e[1] - e[0]) / np.linalg.norm(e[1] - e[0])        # where it lies on an axis line
        size0 = np.ptp(data.X) + np.ptp(data.Y)
        on_axis = any(abs(d[0] * w[1] - d[1] * w[0]) < 1e-6 and
                      abs((p0 - e[0])[0] * d[1] - (p0 - e[0])[1] * d[0]) < 1e-6 * size0 for p0, w, _ in lines)
        if not on_axis:
            ln, = ax.plot(e[:, 0], e[:, 1], color=axis_colour(s)[0], lw=1.2, zorder=6)
            ln.set_clip_path(clip)
    P = np.vstack(box)
    margin = 0.02 * size                                   # labels outside are kept by the tight crop
    ax.set_xlim(P[:, 0].min() - margin, P[:, 0].max() + margin)
    ax.set_ylim(P[:, 1].min() - margin, P[:, 1].max() + margin)
    ax.set_aspect('equal')
    ax.axis('off')
    if s.labels:
        for xy, name in places:
            ax.keep_texts.append(put_label(ax, data, s, f, xy, name, lines, on_map=ext))
    if data.exp is not None:
        ax.plot(data.exp[:, 0], data.exp[:, 1], 'x', color='red', ms=9, mew=2, label='experiment', zorder=9)
        ax.legend(loc='upper left', bbox_to_anchor=(0.0, 1.0), fontsize=0.8 * f, frameon=False)
    info = mark_max(ax, data, s, z, clip)
    return lo, hi, info


def xunit(data, s):
    """the unit of the x axis of a line: the setting if this line has it, else its default"""
    return s.xaxis if s.xaxis in data.units else data.default_unit


def origins(data, unit):
    """x of the origin on a line: q = 0 of a chain, 0 of the running coordinate of a cut"""
    return [0.0]


def line_max(ax, d, s, X):
    """the largest value of a line away from Gamma (more than one mesh step), marked by red dashed lines at
    each place it occurs; returns the note for the terminal"""
    z = d.val[s.part]
    fin = np.isfinite(z)
    ok = fin & ~d.gamma_mask(1)
    if not ok.any():
        ok = fin
    m = z[ok].max()
    hit = np.flatnonzero(ok & (z >= m - 1e-6 * max(abs(m), 1e-30)))
    for x in np.unique(np.round(X[hit], 6)):
        ax.axvline(x, color='red', ls='--', lw=1.2, zorder=3)
    where = '; '.join('H K L = ' + ' '.join(f'{v:.4f}' for v in d.hkl[i]) for i in hit[:4])
    lab = d.root.label(d.sign) if d.root.pair else ''
    return (f'{lab}: ' if lab else '') + f'max {"Re χ₀" if s.part == "re" else "Im χ₀/ω"} = {m:.5g} away from Γ ' \
        f'(red dashed line{"s" if len(hit) > 1 else ""}) at {where}'


def draw_line(ax, curves, s, lo, hi, first=True, tag=''):
    """chi0 along a line: a cross section of a map, or the line of a 1D system of chi_bz.x. curves: the line drawn (one: a line, +c or -c, or their average); first: the left panel (the y
    axis name); tag: the line of the panel, above it. Returns notes for the terminal"""
    f = s.font
    data = curves[0]
    root = data.root
    unit = xunit(data, s)
    X = data.xs[unit]
    for d in curves:
        z = d.val[s.part]
        if data.kind == 'line':
            ax.plot(X, z, '-o', ms=3, lw=1.2, color='tab:blue')
            continue
        for idx in np.split(np.arange(len(X)), data.breaks):   # the pieces between jumps
            ax.plot(X[idx], z[idx], '-', lw=1.4, color='tab:blue')
        for n in np.unique(data.seg):                           # markers where the mesh gives few points
            sel = np.flatnonzero(data.seg == n)
            if len(sel) <= 16:
                ax.plot(X[sel], z[sel], 'o', ms=3.5, color='tab:blue')
    orig = origins(data, unit)
    for x in orig:                                             # the origin: always a solid line
        ax.axvline(x, color='black', lw=0.8, zorder=1)
    ticks = {}                                                 # x -> label; two labels at one x: a jump
    for x, ((_,), name) in zip(data.pxs[unit], data.pts):
        if name.lower() in GAMMA_NAMES:
            name = 'Γ'
        key = round(float(x), 6)
        ticks[key] = f'{ticks[key]}|{name}' if key in ticks and ticks[key].split('|')[-1] != name else name
    xs = sorted(ticks)
    x0, x1 = min(xs + [X.min()]), max(xs + [X.max()])
    above = False                                              # the point labels above the plot
    spaced(ax.xaxis, s.tick, x0, x1, s.tick_frac)              # numbers (automatic without a spacing)
    if data.kind == 'line':                                    # a chain: zone boundaries, labels in the plot
        if s.zones:
            for b in data.zones:
                xb = b[0, 0] if unit == 'invA' else np.interp(b[0, 0], data.xs['invA'], data.xs['rlu'])
                ax.axvline(xb, color='black', lw=0.8)
        if s.labels:
            for x in xs:
                if x not in orig:
                    ax.axvline(x, color='gray', ls=':', lw=0.8)
                ax.annotate(ticks[x], (x, 1.0), xycoords=('data', 'axes fraction'), xytext=(2, -1.2 * f),
                            textcoords='offset points', fontsize=f, fontweight='bold')
    elif s.labels and xs:                                      # a cut: lines at the labelled points, labels above
        for x in xs:
            if x not in orig:
                ax.axvline(x, color='black', lw=0.6)
        above = True
        sec = ax.secondary_xaxis('top')
        sec.set_xticks(xs)
        sec.set_xticklabels([ticks[x] for x in xs], fontsize=f, fontweight='bold')
        sec.tick_params(axis='x', length=0)
    if root.exp is not None:
        for q in root.exp[:, 0]:
            ax.axvline(q, color='tab:green', ls='-.', lw=1.2)     # red dashed: the maximum
    ax.set_xlim(x0, x1)
    ax.set_ylim(lo, hi)
    if first:
        ax.set_ylabel(CHI_LABEL[s.part], fontsize=f)
    ax.tick_params(axis='both', labelsize=0.85 * f)
    if data.xnames.get(unit):
        ax.set_xlabel(data.xnames[unit], fontsize=f)
    notes = []
    if s.tick and not isinstance(ax.xaxis.get_major_locator(), MultipleLocator):
        notes.append(f'x tick spacing {s.tick:g}: more than {MAXTICKS} ticks, the automatic ticks are used')
    if not spaced(ax.yaxis, s.ytick, lo, hi) and s.ytick:
        notes.append(f'y tick spacing {s.ytick:g}: more than {MAXTICKS} ticks, the automatic ticks are used')
    for axis, step, a0, a1, name in ((ax.xaxis, s.minor, x0, x1, 'x'), (ax.yaxis, s.yminor, lo, hi, 'y')):
        if not step:
            continue
        if (a1 - a0) / step > 4 * MAXTICKS:
            notes.append(f'{name} minor tick spacing {step:g}: more than {4 * MAXTICKS} ticks, none drawn')
            continue
        axis.set_minor_locator(MultipleLocator(step))          # short ticks without labels
    ax.tick_params(which='minor', length=2.5, width=0.6)
    ax.tick_params(which='both', direction='in')               # ticks point into the plot
    if tag:                                                    # the line of the panel, above it
        ax.set_title(tag, fontsize=f, pad=0.4 * f + (1.6 * f if above else 0))
    if s.mark:
        for d in curves:
            notes.append(line_max(ax, d, s, X))
    for d in curves:                                           # a Gamma spike that flattens the rest
        z = d.val[s.part]
        fin = np.isfinite(z)
        isg = np.all(np.abs(d.hkl - np.rint(d.hkl)) < 1e-6, axis=-1) & fin
        rest = z[fin & ~isg]
        if s.ymax is None and isg.any() and rest.size and z[isg].max() - rest.max() > np.ptp(rest) > 0.0:
            notes.append(f'the value at Γ ({z[isg].max():.5g}) is far above the rest (at most {rest.max():.5g}): '
                         f'"y axis to" (--ymax) at about {rest.max() + 0.05 * np.ptp(rest):.5g} shows the rest')
            break
    return notes


def declutter(fig, ax):
    """label every n-th tick of an axis when the labels would touch (all tick marks stay), and hide a tick
    label that runs into a point label or an axis name"""
    r = fig.canvas.get_renderer()
    keep = [t.get_window_extent(r).padded(2.0) for t in ax.keep_texts if t is not None]
    for x, y in ax.keep_points:                                # the markers of the labelled points
        cx, cy = ax.transData.transform((x, y))
        h = 4.0 * fig.dpi / 72.0
        keep.append(Bbox([[cx - h, cy - h], [cx + h, cy + h]]))
    for group in ax.tick_groups:
        if not group:
            continue
        boxes = {j: t.get_window_extent(r) for j, t in group}
        gap = 0.45 * group[0][1].get_fontsize() * fig.dpi / 72.0   # at least half a letter between labels
        for n in range(1, 11):                                 # the smallest stride without overlaps
            shown_j = sorted(j for j in boxes if j % n == 0)
            if not any(boxes[a].padded(gap / 2).overlaps(boxes[b].padded(gap / 2))
                       for a, b in zip(shown_j, shown_j[1:])):
                break
        for j, t in group:
            t.set_visible(j % n == 0 and not any(boxes[j].overlaps(k) for k in keep))


def export(data, s, out):
    """the PDF"""
    if data.kind != 'plane':
        curves = chosen(data, s)
        lo, hi = yrange(data, s)
        groups = [[c] for c in curves]                         # panels: +c and -c side by side, as in 2D
        n = len(groups)                                        # (one: a line, or the average of +c and -c)
        fig, axes = plt.subplots(1, n, figsize=(8 if n == 1 else 13, 4.5), sharey=True, squeeze=False)
        notes = []
        for i, (ax, g) in enumerate(zip(axes[0], groups)):
            tag = ''
            if n > 1 and s.panels == 'plane':
                tag = data.root.label(g[0].sign)
            elif n > 1 and s.panels == 'letters':
                tag = f'({"ab"[i]})'
            notes += [m for m in draw_line(ax, g, s, lo, hi, first=i == 0, tag=tag) if m not in notes]
        if n > 1:
            fig.subplots_adjust(wspace=0.06)
        fig.canvas.draw()
        for ax in axes[0]:
            thin(fig, ax.xaxis, s.tick, s.tick_frac)
        thin(fig, axes[0, 0].yaxis, s.ytick, False)           # shared by the panels
        for ax in axes[0]:                                     # bold tick labels
            for lab in ax.get_xticklabels() + ax.get_yticklabels():
                lab.set_fontweight('bold')
        fig.savefig(out, bbox_inches='tight', pad_inches=0.03, dpi=300)
        plt.close(fig)
        return '\n'.join(notes)
    datas = chosen(data, s)
    n = len(datas)
    fig = plt.figure(figsize=(7.5 * n, 6.5))
    axes = [fig.add_axes([0.02 + 0.9 * i / n, 0.04, 0.9 / n - 0.04, 0.92]) for i in range(n)]
    cax = fig.add_axes([0.94, 0.2, 0.02, 0.6])
    infos, tags = [], []
    for i, (ax, d) in enumerate(zip(axes, datas)):
        info = draw(ax, cax if i == n - 1 else None, d, s)[2]
        tag = ''
        if n > 1 and s.panels == 'plane':
            tag = data.root.label(d.sign)
        elif n > 1 and s.panels == 'letters':
            tag = f'({"abcdefgh"[i]})'
        tags.append(tag)
        lab = data.root.label(d.sign)
        infos.append((f'{lab}: ' if lab else '') + info if info else '')
    fig.canvas.draw()
    r = fig.canvas.get_renderer()
    W, H = fig.bbox.width, fig.bbox.height
    right = None
    for ax, tag in zip(axes, tags):                            # left to right, each clear of the last
        declutter(fig, ax)
        pos = ax.get_position()                                # the box after the equal aspect
        vis = [t.get_window_extent(r) for t in ax.texts if t.get_visible()]
        left = min([pos.x0] + [b.x0 / W for b in vis])
        rgt = max([pos.x1] + [b.x1 / W for b in vis])
        top = max([pos.y1] + [b.y1 / H for b in vis])
        shift = 0.0 if right is None else right + 0.03 - left
        if shift:
            ax.set_position([pos.x0 + shift, pos.y0, pos.width, pos.height])
        right = rgt + shift
        if tag:                                                # the plane of the panel, above it
            fig.text(0.5 * (pos.x0 + pos.x1) + shift, top + 0.012, tag, ha='center', va='bottom', fontsize=s.font)
    pos = axes[-1].get_position()
    if s.cbar == 'ticks':
        cax.set_position([right + 0.03 / n, pos.y0 + 0.05 * pos.height, 0.03 / n, 0.9 * pos.height])
    else:                                                      # small, beside the middle of the zone
        cax.set_position([right + 0.02 / n, pos.y0 + 0.3 * pos.height, 0.022 / n, 0.4 * pos.height])
    fig.savefig(out, bbox_inches='tight', pad_inches=0.03, dpi=300)
    plt.close(fig)
    return '\n'.join(i for i in infos if i)


def command(data, s, out):
    """options that give the same figure without the questions"""
    plane = data.kind == 'plane'
    root = data.root
    cut = getattr(root, 'cut', None)
    if not plane:
        c = f'python3 plot_bz.py {root.map_prefix if cut else root.prefix} {s.part}' + \
            (f' {root.exp_file}' if root.exp_file else '') + (f' --cut {cut_text(cut)} --cut-only' if cut else '')
        if cut and root.cut_range:
            c += f' --xrange {root.cut_range[0]:g} {root.cut_range[1]:g}'
        if root.pair:
            c += f' --planes {s.planes}' + (f' --panels {s.panels}' if s.planes == 'panels' and s.panels != 'plane'
                                            else '')
        c += '' if s.mark else ' --no-max'
        c += f' --xaxis {xunit(data, s)}' if xunit(data, s) != data.default_unit else ''
        if s.tick:
            c += (' --cut-tick ' if cut else ' --tick ') + (frac_text(s.tick) if s.tick_frac else f'{s.tick:g}')
        c += (f' --cut-minor {s.minor:g}' if cut else f' --minor {s.minor:g}') if s.minor else ''
        c += f' --ytick {s.ytick:g}' if s.ytick else ''
        c += f' --yminor {s.yminor:g}' if s.yminor else ''
        c += f' --ymin {s.ymin:.6g}' if s.ymin is not None else ''
        c += f' --ymax {s.ymax:.6g}' if s.ymax is not None else ''
        c += f' --font {s.font:g}' + ('' if s.labels else ' --no-labels') + \
            ('' if s.zones or data.kind != 'line' else ' --no-zones')
        return c + f' --batch --out {out}'
    c = f'python3 plot_bz.py {root.prefix} {s.part}' + (f' {root.exp_file}' if root.exp_file else '')
    if root.pair:
        c += f' --planes {s.planes}' + (f' --panels {s.panels}' if s.planes == 'panels' and s.panels != 'plane' else '')
    c += f' --region {s.region} --view {s.view} --cmap {s.cmap}' if plane else ''
    if s.tick and (not plane or s.view != 'zone'):
        c += ' --tick ' + (frac_text(s.tick) if s.tick_frac else f'{s.tick:g}')
    if s.minor and plane and s.view != 'zone':
        c += f' --minor {s.minor:g}'
    c += f' --axis-colour {s.axcolour}' if plane and s.view != 'zone' and s.axcolour != 'auto' else ''
    c += f' --cbar {s.cbar}' if plane and s.cbar != 'minmax' else ''
    if plane and s.view == 'zone':
        c += (' --rings ' + ','.join(frac_text(r) for r in s.rings)) if s.rings else ''
        c += (' --guides' if s.guides else '') + (' --bar' if s.bar else '')
    c += ' --no-max' if not s.mark and plane else ''
    c += f' --cut {cut_text(s.cut)}' if s.cut else ''
    c += f' --vmin {s.vmin:.6g}' if s.vmin is not None else f' --low {s.low:g}'
    c += f' --vmax {s.vmax:.6g}' if s.vmax is not None else f' --high {s.high:g}'
    c += f' --gamma {s.gamma:g} --levels {s.levels}' if plane else ''
    c += f' --font {s.font:g} --gamma-steps {s.gsteps if s.ignore else 0}'
    c += (' --lines' if s.lines and plane else '') + ('' if s.zones else ' --no-zones') + \
        ('' if s.labels else ' --no-labels')
    return c + f' --batch --out {out}'
HINT = ('  for a weak dispersion: colour map inferno/magma/hot/gray, ignore Gamma, low 50-80, gamma 2-3,\n'
        '  10-15 contour levels')


def summary(data):
    """the range of the values, to choose the colour scale"""
    root = data.root
    n = 'x'.join(str(k) for k in root.val['re'].shape)
    what = {'plane': '2D map (chi_bz.x)', 'path': '1D cross section', 'line': '1D line of a chain (chi_bz.x)'}
    if getattr(root, 'cut', None):
        what['path'] = f'1D cut {cut_text(root.cut)} through {root.map_prefix}'
    print(f'{os.path.basename(root.mapfile)}: {what[root.kind]}{", two planes" if root.pair and root.kind == "plane" else ", two lines" if root.pair else ""}, {n} points')
    for key, v in root.views.items():
        if key == 'avg' and root.kind != 'plane':
            continue
        if root.pair:
            print(f' {root.label(key)}:')
        one_summary(v)


def one_summary(data):
    away = ~data.gamma_mask(1)
    if data.kind != 'plane':                                   # a line: the range of the value axis
        for part, name in (('re', 'Re chi0  '), ('im', 'Im chi0/w')):
            z = data.val[part]
            fin = np.isfinite(z)
            f = z[fin & away]
            print(f'  {name}  min {z[fin].min():.5g}  max {z[fin].max():.5g}' +
                  (f'  (max away from Gamma {f.max():.5g})' if f.size and f.max() < z[fin].max() else ''))
        return
    for part, name in (('re', 'Re chi0  '), ('im', 'Im chi0/w')):
        z = data.val[part]
        fin = np.isfinite(z)
        d, f = z[fin], z[fin & away]
        line = f'  {name}  min {d.min():.5g}  max {d.max():.5g}'
        if f.size:
            p = np.percentile(f, [50, 80, 95])
            line += f'  max away from Gamma {f.max():.5g};  percentiles 50/80/95 (without Gamma): ' + \
                    ' / '.join(f'{x:.5g}' for x in p)
        print(line)


def items(data):
    """the settings that apply to this kind of plot: key, description"""
    plane = data.kind == 'plane'
    it = [('part', 'part (re, im)')]
    root = data.root
    if not plane:
        return line_items(data)
    if root.pair:
        L, v = root.pair[0], root.label('+').split('+', 1)[1]
        it += [('planes', f'planes (together: {L} = +{v} and -{v} averaged in one map; panels: side by side; '
                          f'+: {L} = +{v}; -: {L} = -{v})'),
               ('panels', 'panel labels (plane: ' + root.label('+') + '; letters: (a), (b); none)')]
    if plane:
        it += [('region', 'region (first: the first Brillouin zone; extended: the whole map of chi_bz.x '
                          'with all zone boundaries)'),
               ('view', 'view (zone: just the zone; cartesian: k axes through Γ in 1/Å; '
                        'crystal: axes along the reciprocal vectors in r.l.u.)'),
               ('tick', 'major tick spacing on the axes (labelled; 1/Å or r.l.u., e.g. 0.1 or 1/6; auto)'),
               ('minor', 'minor tick spacing on the axes (short, no label, e.g. 0.05; none)'),
               ('axcolour', 'axis colour (auto: a contrast to the colour map; or white, cyan, ...)'),
               ('cbar', 'colour bar (minmax: the words min and max; values: the two values; ticks: a full scale)'),
               ('cmap', 'colour map (' + ', '.join(CMAPS) + ')')]
    it += [('low', 'colour scale bottom (percentile 0-100, or a value as =0.185)'),
           ('high', 'colour scale top (percentile 0-100, or a value as =0.2)'),
           ('gamma', 'gamma (> 1: more colours for the largest values)'),
           ('ignore', 'leave Γ out of the colour scale and the maximum (y/n, or mesh steps around Γ, e.g. 3)'),
           ('lines', 'contour lines (y/n)'), ('levels', 'contour levels'),
           ('rings', 'dashed outlines at fractions of Γ → border (e.g. 1/3 2/3; - for none)'),
           ('guides', 'guide lines Γ → corners and edge centres (y/n)'),
           ('bar', 'length bar in 1/Angstrom (y/n)'),
           ('mark', 'mark the maximum with a bold red + (y/n)'),
           ('font', 'font size (points)'), ('zones', 'zone boundaries (y/n)'), ('labels', 'point labels (y/n)'),
           ('exp', 'experiment file (qx qy in 1/Angstrom; - for none)'), ('out', 'output PDF'),
           ('cut', f'1D cross section at constant x or y of the map, in r.l.u. (x: {"HKL"[data.root.ia]}, y: '
                   f'{"HKL"[data.root.ib]}; any value, e.g. y=0.25 or x=-1/3; - for none): a solid line on the '
                   f'map, then a 1D plot of it')]
    return it


def line_items(data):
    """the settings of a 1D plot (chi0 against the position on a line)"""
    root = data.root
    it = [('part', 'part (re, im)')]
    if root.pair:
        L, v = root.pair[0], root.label('+').split('+', 1)[1]
        it += [('planes', f'lines (together: {L} = +{v} and -{v} merged into one curve, their average; panels: '
                          f'side by side, as in the 2D plot; +: {L} = +{v} only; -: {L} = -{v} only)'),
               ('panels', 'panel labels (plane: ' + root.label('+') + '; letters: (a), (b); none)')]
    if getattr(root, 'cut', None):
        i, j = 'HKL'[root.cut_axes[0]], 'HKL'[root.cut_axes[1]]
        it += [('xrange', f'x axis range ({j} in r.l.u. along the cut {i} = {root.cut[1]:g}: from to, e.g. -1 1; '
                          f'auto: -0.5 0.5)')]
    names = {'rlu': 'rlu: r.l.u.', 'invA': 'invA: 1/Å'}
    it += [('xaxis', 'x axis unit (' + '; '.join(names[u] for u in data.units) + ')'),
           ('tick', 'x major tick spacing (labelled; in the unit of the x axis, e.g. 0.1 or 0.125; auto)'),
           ('minor', 'x minor tick spacing (short, no label, e.g. 0.05; none)')]
    it += [('ytick', 'y major tick spacing (labelled, e.g. 0.05; auto)'),
           ('yminor', 'y minor tick spacing (short, no label, e.g. 0.01; none)'),
           ('ymin', 'y axis from (a value; auto: the smallest value)'),
           ('ymax', 'y axis to (a value; auto: the largest value)'),
           ('mark', 'mark the largest value (away from Γ) with a red dashed line (y/n)'),
           ('font', 'font size (points)'), ('labels', 'point labels (y/n)')]
    if data.kind == 'line':
        it += [('zones', 'zone boundaries (y/n)')]
    it += [('exp', 'experiment file (q in the units of the x axis, one per line; - for none)'),
           ('out', 'output PDF')]
    return it


def parse_rings(a):
    """fractions such as  1/3 2/3  or  0.25,0.5,0.75 ; '-' for none"""
    if a.strip() in ('-', 'none', ''):
        return []
    out = []
    for w in a.replace(',', ' ').split():
        try:
            r = float(Fraction(w))
        except (ValueError, ZeroDivisionError):
            raise ValueError(f'{w}: fractions such as 1/3 2/3 or 0.25 0.5, or - for none')
        if not 0.0 < r < 1.0:
            raise ValueError(f'{w}: between 0 and 1')
        out.append(r)
    return sorted(set(out))


def parse_tick(a):
    """a tick spacing such as 0.1 or 1/6, or auto: (value or None, written as a fraction)"""
    if a.strip().lower() in ('auto', ''):
        return (None, False)
    try:
        v = float(Fraction(a.strip()))
    except (ValueError, ZeroDivisionError):
        raise ValueError('a spacing such as 0.1 or 1/6, or auto')
    if v <= 0.0:
        raise ValueError('must be > 0')
    return (v, '/' in a)


def outname(data, s):
    root = data.root
    if s.out:
        return s.out
    if root.pair is None or s.planes == 'together':
        return f'{root.prefix}_{s.part}.pdf'
    if s.planes == 'panels':
        return f'{root.prefix}_panels_{s.part}.pdf'
    return re.sub(r'_([HKL])\+-(?=[0-9.]+(_|$))', rf'_\g<1>{s.planes}', root.prefix, count=1) + f'_{s.part}.pdf'


def current(data, s, key):
    """the value as shown in the questions and the list"""
    if key in ('low', 'high'):
        lo, hi, _ = limits(data, s)
        if (s.vmin if key == 'low' else s.vmax) is not None:
            return f'={lo if key == "low" else hi:.6g}'
        return f'{s.low if key == "low" else s.high:g} % (= {lo if key == "low" else hi:.5g})'
    if key == 'ignore':
        return f'y, {s.gsteps} step{"s" if s.gsteps > 1 else ""}' if s.ignore else 'n'
    if key in ('lines', 'zones', 'labels', 'guides', 'bar', 'mark'):
        return 'y' if getattr(s, key) else 'n'
    if key == 'rings':
        return ' '.join(frac_text(r) for r in s.rings) or '-'
    if key == 'xaxis':
        return xunit(data, s)
    if key == 'cut':
        return cut_text(s.cut) if s.cut else '-'
    if key == 'xrange':
        r = data.cut_range or (-0.5, 0.5)
        return f'{r[0]:g} {r[1]:g}' + ('' if data.cut_range else ' (auto)')
    if key == 'tick':
        if not s.tick:
            return 'auto'
        return frac_text(s.tick) if s.tick_frac else f'{s.tick:g}'
    if key == 'ytick':
        return f'{s.ytick:g}' if s.ytick else 'auto'
    if key in ('minor', 'yminor'):
        v = getattr(s, key)
        return f'{v:g}' if v else 'none'
    if key in ('ymin', 'ymax'):
        v = getattr(s, key)
        if v is not None:
            return f'{v:.6g}'
        lo, hi = yrange(data, s)
        return f'auto (= {lo if key == "ymin" else hi:.5g})'
    if key == 'exp':
        return data.exp_file or '-'
    if key == 'out':
        return outname(data, s)
    v = getattr(s, key)
    return f'{v:g}' if isinstance(v, float) else str(v)


def parse(data, s, key, ans):
    """the new value of a setting from an answer; ValueError with a message if it is not valid"""
    a = ans.strip()
    if key == 'view' and a.lower() == 'axes':                # the name in earlier versions
        return 'cartesian'
    if key == 'cut':
        return parse_cut(data.root, a)
    if key == 'xrange':
        r = parse_range(a)
        try:                                                   # the cut again, over the new range
            new = build_cut(data.cut_source, None, data.cut, r)
        except ValueError as e:
            raise ValueError(str(e))
        data.__dict__.update(new.__dict__)
        data.root = data
        data.make_views()
        return r
    if key == 'xaxis':
        u = {'rlu': 'rlu', 'r.l.u.': 'rlu', 'r': 'rlu',
             'crystal': 'rlu', 'inva': 'invA', '1/a': 'invA', '1/å': 'invA', 'a': 'invA', 'ang': 'invA',
             'angstrom': 'invA', 'i': 'invA'}.get(a.lower())
        if u not in data.units:
            raise ValueError('one of: ' + ', '.join(data.units))
        return u
    if key == 'tick':
        return parse_tick(a)
    if key == 'ytick':
        return parse_tick(a)[0]
    if key in ('minor', 'yminor'):
        if a.strip().lower() in ('none', 'no', 'n', '-', '0', 'off'):
            return None
        v = parse_tick(a)[0]
        major = s.tick if key == 'minor' else s.ytick
        if v and major and not v < major:
            raise ValueError(f'smaller than the major spacing {major:g}')
        return v
    if key in ('ymin', 'ymax'):
        if a.lower() in ('auto', 'a'):
            return None
        v = float(a)
        other = s.ymax if key == 'ymin' else s.ymin
        if other is not None and (v >= other if key == 'ymin' else v <= other):
            raise ValueError(f'must be {"below" if key == "ymin" else "above"} {other:g}')
        return v
    if key == 'axcolour':
        if a.lower() == 'auto':
            return 'auto'
        if not colors.is_color_like(a):
            raise ValueError('auto, or a colour such as white, cyan, black, #00ffff')
        return a
    if key == 'planes':
        a2 = {'t': 'together', 'together': 'together', 'average': 'together', 'avg': 'together', '1': 'together',
              'p': 'panels', 'panels': 'panels', 'both': 'panels', '2': 'panels', '+': '+', 'plus': '+', '3': '+',
              '-': '-', 'minus': '-', '4': '-'}.get(a.lower())
        if a2 is None:
            raise ValueError('together, panels, + or -')
        return a2
    if key in ('part', 'cmap', 'view', 'cbar', 'panels', 'region'):
        opts = {'part': ['re', 'im'], 'cmap': CMAPS, 'view': VIEWS, 'cbar': CBARS,
                'panels': PANELS, 'region': REGIONS}[key]
        if a.isdigit() and 1 <= int(a) <= len(opts):
            return opts[int(a) - 1]
        hit = [o for o in opts if o.startswith(a.lower())]
        if len(hit) != 1:
            raise ValueError('one of: ' + ', '.join(opts))
        return hit[0]
    if key == 'rings':
        return parse_rings(a)
    if key == 'ignore' and a.isdigit():
        if int(a) == 0:
            return False
        if int(a) > 50:
            raise ValueError('at most 50 mesh steps')
        return int(a)
    if key in ('ignore', 'lines', 'zones', 'labels', 'guides', 'bar', 'mark'):
        if a.lower() in ('y', 'yes', '1', 'true'):
            return True
        if a.lower() in ('n', 'no', '0', 'false'):
            return False
        raise ValueError('y or n')
    if key in ('low', 'high'):
        if a.startswith('='):
            return ('value', float(a[1:]))
        v = float(a.rstrip('%'))
        if not 0.0 <= v <= 100.0:
            raise ValueError('a percentile 0 ... 100, or a value as =0.185')
        if key == 'low' and s.vmax is None and v >= s.high:
            raise ValueError(f'must be below the top ({s.high:g} %)')
        if key == 'high' and s.vmin is None and v <= s.low:
            raise ValueError(f'must be above the bottom ({s.low:g} %)')
        return ('pct', v)
    if key == 'gamma':
        v = float(a)
        if not 0.05 <= v <= 10.0:
            raise ValueError('0.05 ... 10')
        return v
    if key == 'levels':
        v = int(a)
        if not 2 <= v <= 100:
            raise ValueError('2 ... 100')
        return v
    if key == 'font':
        v = float(a)
        if not 4.0 <= v <= 48.0:
            raise ValueError('4 ... 48')
        return v
    if key == 'exp':
        if a in ('-', 'none'):
            return None
        if not os.path.isfile(a):
            raise ValueError(f'no file {a}')
        np.loadtxt(a, ndmin=2)                                 # readable?
        return a
    if key == 'out':
        return a if a.lower().endswith('.pdf') else a + '.pdf'
    raise ValueError(key)


def apply(data, s, key, v):
    if key in ('low', 'high'):
        kind, x = v
        if key == 'low':
            s.low, s.vmin = (x, None) if kind == 'pct' else (s.low, x)
        else:
            s.high, s.vmax = (x, None) if kind == 'pct' else (s.high, x)
    elif key == 'ignore' and not isinstance(v, bool):
        s.ignore, s.gsteps = True, v
    elif key == 'tick':
        s.tick, s.tick_frac = v
    elif key == 'exp':
        root = data.root
        root.exp_file = v
        root.exp = np.loadtxt(v, ndmin=2) if v else None
        root.make_views()
    else:
        setattr(s, key, v)


def ask(data, s, key, text):
    """one question; Enter keeps the value. False at the end of the input"""
    while True:
        try:
            ans = input(f'  {text} [{current(data, s, key)}]: ')
        except EOFError:
            print()
            return False
        if not ans.strip():
            return True
        try:
            v = parse(data, s, key, ans)
        except ValueError as e:
            print(f'    not valid: {e}')
            continue
        apply(data, s, key, v)
        return True


def write(data, s):
    out = outname(data, s)
    info = export(data, s, out)
    print(f'written {out}')
    if info:
        print('  ' + info.replace('\n', '\n  '))
    print('  the same without questions:  ' + command(data, s, out))


def interactive(data, s):
    """the questions, the PDF, then changes until the user is done"""
    it = items(data)
    print('settings (Enter keeps the value in brackets)')
    if data.kind == 'plane':
        print(HINT)
    for key, text in it:
        if key in ('tick', 'minor', 'axcolour') and data.kind == 'plane' and s.view == 'zone':
            continue                                           # no axes in view 'zone'
        if key == 'zones' and (data.kind == 'path' or (data.kind == 'plane' and data.zone is not None)):
            continue                                           # the zone border is always drawn; none on a cut
        if key == 'panels' and s.planes != 'panels':
            continue                                           # panel labels only with two panels
        if not ask(data, s, key, text):
            break
    pending = data.kind == 'plane' and bool(s.cut)             # a 1D plot of the cut is due
    while True:
        write(data, s)
        if pending:
            cut_session(data, s)
            pending = False
        print('settings:')
        for i, (key, text) in enumerate(it, 1):
            print(f'  {i:2d}  {text.split(" (")[0]:<36s} {current(data, s, key)}')
        try:
            ans = input('numbers to change (e.g. 4 6), Enter or q to finish: ').strip().lower()
        except EOFError:
            print()
            return
        if ans in ('', 'q', 'quit', 'exit'):
            return
        picks = []
        for w in ans.replace(',', ' ').split():
            if w.isdigit() and 1 <= int(w) <= len(it):
                picks.append(int(w) - 1)
            else:
                print(f'  {w}: not a number from the list')
        for i in picks:
            if not ask(data, s, *it[i]):
                return
            if it[i][0] == 'cut' and s.cut:
                pending = True


def line_settings(s, tick='auto', minor=None):
    """the settings of a 1D plot of a cut, from those of the map: the 1D ones at their defaults"""
    s1 = copy.copy(s)
    s1.tick, s1.tick_frac = parse_tick(tick)
    s1.minor, s1.yminor = minor, None
    s1.ytick, s1.ymin, s1.ymax, s1.xaxis, s1.out, s1.cut, s1.xrange = None, None, None, None, None, None, None
    return s1


def cut_session(data, s):
    """the 1D plot of the cut of the map: its own questions, then back to the map"""
    s1 = getattr(s, 'line_s', None) or line_settings(s)       # the 1D settings of the last cut are kept
    try:
        cd = build_cut(data, s, s.cut, s1.xrange)
    except ValueError as e:
        print(f'  cut {cut_text(s.cut)}: {e}')
        return
    s1.part, s1.planes, s1.font = s.part, s.planes, s.font
    print(f'\n--- 1D plot of the cut {cut_text(cd.cut)} (values in {cd.mapfile})')
    if cd.note:
        print(f'  {cd.note}')
    summary(cd)
    interactive(cd, s1)
    s.line_s = s1
    print('--- back to the 2D plot\n')


def main():
    ap = argparse.ArgumentParser(description=__doc__.split('\n')[0], epilog=__doc__.split('\n', 2)[2],
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument('prefix', nargs='?', default='')
    ap.add_argument('part', nargs='?', default='re', choices=['re', 'im'])
    ap.add_argument('exp', nargs='?')
    ap.add_argument('--cmap', default='inferno', choices=CMAPS)
    ap.add_argument('--low', type=float, default=0.0)
    ap.add_argument('--high', type=float, default=100.0)
    ap.add_argument('--vmin', type=float)
    ap.add_argument('--vmax', type=float)
    ap.add_argument('--gamma', type=float, default=1.0)
    ap.add_argument('--levels', type=int, default=12)
    ap.add_argument('--lines', action='store_true')
    ap.add_argument('--font', type=float, default=12.0)
    ap.add_argument('--gamma-steps', type=int, default=1)
    ap.add_argument('--no-zones', action='store_true')
    ap.add_argument('--no-labels', action='store_true')
    ap.add_argument('--out')
    ap.add_argument('--batch', action='store_true')
    ap.add_argument('--ytick', default='auto')
    ap.add_argument('--xaxis', choices=['rlu', 'invA'])
    ap.add_argument('--cut')
    ap.add_argument('--cut-only', action='store_true')
    ap.add_argument('--cut-tick', default='auto')
    ap.add_argument('--minor', type=float)
    ap.add_argument('--yminor', type=float)
    ap.add_argument('--cut-minor', type=float)
    ap.add_argument('--xrange', nargs=2, type=str, metavar=('FROM', 'TO'))
    ap.add_argument('--ymin', type=float)
    ap.add_argument('--ymax', type=float)
    ap.add_argument('--view', default='zone', choices=VIEWS + ['axes'])
    ap.add_argument('--tick', default='auto')
    ap.add_argument('--axis-colour', default='auto')
    ap.add_argument('--cbar', default='minmax', choices=CBARS)
    ap.add_argument('--planes', default='together', choices=PLANES + ['both'])
    ap.add_argument('--region', default='auto', choices=['auto'] + REGIONS)
    ap.add_argument('--panels', default='plane', choices=PANELS)
    ap.add_argument('--rings', default='-')
    ap.add_argument('--guides', action='store_true')
    ap.add_argument('--bar', action='store_true')
    ap.add_argument('--max', action='store_true', help=argparse.SUPPRESS)   # the default; kept for old commands
    ap.add_argument('--no-guides', action='store_true', help=argparse.SUPPRESS)   # earlier versions: now the default
    ap.add_argument('--no-max', action='store_true')
    a = ap.parse_args()
    if not (0.0 <= a.low < a.high <= 100.0):
        ap.error('0 <= --low < --high <= 100')
    if a.gamma <= 0.0 or a.levels < 2:
        ap.error('--gamma > 0 and --levels >= 2')
    np.seterr(over='ignore')          # the power-law scale on the open-ended contour levels

    batch = a.batch or not sys.stdin.isatty()
    base = resolve(a.prefix) if a.prefix else None
    if base is None:                                           # a chi file name, or none: ask 1D / 2D
        if batch:
            sets = datasets(a.prefix)
            sys.exit(f'{a.prefix or "no prefix"}: give the prefix of one plot' +
                     (': ' + ', '.join(b for _, b in sets) if sets else
                      ' (no files of chi_bz.x here)'))
        base = choose(a.prefix)
    data = Data(base, a.exp)
    s = argparse.Namespace(part=a.part, cmap=a.cmap, low=a.low, high=a.high, vmin=a.vmin,
                           vmax=a.vmax, gamma=a.gamma, levels=a.levels, lines=a.lines, font=a.font,
                           gsteps=max(a.gamma_steps, 1), ignore=a.gamma_steps > 0, zones=not a.no_zones,
                           labels=not a.no_labels, out=a.out, view=a.view, guides=a.guides and not a.no_guides,
                           bar=a.bar, mark=not a.no_max)   # the maximum is marked by default
    if s.view == 'axes':
        s.view = 'cartesian'
    try:
        s.rings = parse_rings(a.rings)
        s.tick, s.tick_frac = parse_tick(a.tick)
        s.ytick = parse_tick(a.ytick)[0]
    except ValueError as e:
        ap.error(f'--rings / --tick / --ytick: {e}')
    s.minor, s.yminor = a.minor, a.yminor
    for name, v in (('--minor', a.minor), ('--yminor', a.yminor), ('--cut-minor', a.cut_minor)):
        if v is not None and not v > 0:
            ap.error(f'{name}: a spacing > 0')
    s.ymin, s.ymax = a.ymin, a.ymax
    s.xaxis = a.xaxis
    s.cut = None
    if a.cut:
        if data.kind != 'plane':
            ap.error('--cut: a cross section of a 2D map (a prefix of chi_bz.x)')
        try:
            s.cut = parse_cut(data, a.cut)
        except ValueError as e:
            ap.error(f'--cut: {e}')
    if a.xaxis and data.kind != 'plane' and a.xaxis not in data.units:
        ap.error(f'--xaxis for {base}: one of ' + ', '.join(data.units))
    if s.ymin is not None and s.ymax is not None and not s.ymin < s.ymax:
        ap.error('--ymin < --ymax')
    if a.axis_colour != 'auto' and not colors.is_color_like(a.axis_colour):
        ap.error('--axis-colour: auto or a colour such as white, cyan, #00ffff')
    s.axcolour = a.axis_colour
    s.cbar = a.cbar
    s.planes, s.panels = ('panels' if a.planes == 'both' else a.planes), a.panels
    s.region = a.region if a.region != 'auto' else ('extended' if data.extended else 'first')
    if a.out and not a.out.lower().endswith('.pdf'):
        s.out = a.out + '.pdf'
    if batch:
        if not (a.cut_only and s.cut):
            write(data, s)
        if s.cut:                                              # the 1D plot of the cut
            try:
                xr = parse_range(' '.join(a.xrange)) if a.xrange else None
                cd = build_cut(data, s, s.cut, xr)
                s1 = line_settings(s, a.cut_tick, a.cut_minor)
                s1.xrange = xr
                if a.cut_only:                                  # only the cut: --out is its name
                    s1.out = s.out
                s1.ytick, s1.yminor, s1.ymin, s1.ymax, s1.xaxis = s.ytick, s.yminor, s.ymin, s.ymax, a.xaxis
            except ValueError as e:
                sys.exit(f'--cut {a.cut}: {e}')
            if cd.note:
                print(f'cut {cut_text(cd.cut)}: {cd.note}')
            write(cd, s1)
        return
    summary(data)
    interactive(data, s)


if __name__ == '__main__':
    main()
