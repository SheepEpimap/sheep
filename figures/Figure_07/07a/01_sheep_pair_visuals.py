#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
Sheep paired-category visualization suite
Outputs multiple figure styles to:
    /vol2/wulingyun/test/test_sheep/
"""

import textwrap
from pathlib import Path
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import Circle, FancyBboxPatch, PathPatch, Rectangle
from matplotlib.path import Path as MplPath
from matplotlib.gridspec import GridSpec
import numpy as np
from PIL import Image, ImageDraw

OUTPUT_DIR = Path("/vol2/wulingyun/test/test_sheep/")
OUTPUT_DIR.mkdir(parents=True, exist_ok=True)

DATA = [
    {"left_label":"modern_CEA","left_count":1104,"left_source":"Gong, M. et al.",
     "right_label":"modern_EUR","right_count":710,"right_source":"Gong, M. et al."},
    {"left_label":"ancient_CEA","left_count":21,"left_source":"Daly et al.",
     "right_label":"modern_CEA","right_count":1104,"right_source":"Gong, M. et al."},
    {"left_label":"ancient_EUR","left_count":31,"left_source":"Daly et al.",
     "right_label":"modern_EUR","right_count":710,"right_source":"Gong, M. et al."},
    {"left_label":"Demestication_domestic","left_count":738,"left_source":"Yan, Z. et al.",
     "right_label":"Wild","right_count":72,"right_source":"Yan, Z. et al."},
]

BG = "#ffffff"
INK = "#18212a"
MUTED = "#6f7d8c"
LEFT_COL = "#5B8FF9"
RIGHT_COL = "#5AD8A6"
SOURCE_COLORS = {
    "Gong, M. et al.": "#3A8FB7",
    "Daly et al.": "#8E6CBE",
    "Yan, Z. et al.": "#E08D3C",
}

plt.rcParams.update({
    "font.family": "DejaVu Sans",
    "axes.facecolor": BG,
    "figure.facecolor": BG,
    "savefig.facecolor": BG,
})

def fmt(n): return f"{n:,}"
def clean_label(s): return s.replace("_", " ")
def wrap_label(s, w=16): return "\n".join(textwrap.wrap(clean_label(s), width=w, break_long_words=False))

def save_figure(fig, stem):
    fig.savefig(OUTPUT_DIR / f"{stem}.png", dpi=320, bbox_inches="tight")
    fig.savefig(OUTPUT_DIR / f"{stem}.pdf", bbox_inches="tight")
    plt.close(fig)

def plot_mirror_lollipop(data):
    max_count = max(max(d["left_count"], d["right_count"]) for d in data)
    fig, ax = plt.subplots(figsize=(12.5, 7.8))
    fig.subplots_adjust(top=0.84, left=0.06, right=0.98, bottom=0.08)
    fig.text(0.5, 0.965, "Pair-linked categories", ha="center", va="top", fontsize=23, fontweight="bold", color=INK)
    fig.text(
        0.5, 0.91,
        "Mirror lollipop layout: counts stay quantitative, pairing stays explicit, references remain attached to each side",
        ha="center", va="top", fontsize=10.8, color=MUTED
    )
    ax.set_xlim(-1.08*max_count, 1.08*max_count)
    ax.set_ylim(-0.7, len(data)-0.3)

    for i in range(len(data)):
        if i % 2 == 0:
            ax.axhspan(i-0.45, i+0.45, color="#f8fafc", zorder=0)
    ax.axvline(0, color="#ccd6e3", lw=1.2, zorder=1)

    for yi, d in enumerate(reversed(data)):
        i = yi
        lc, rc = d["left_count"], d["right_count"]
        ax.plot([0, -lc], [i, i], color=LEFT_COL, lw=3.5, solid_capstyle="round", zorder=3)
        ax.scatter([-lc], [i], s=160, color=LEFT_COL, edgecolor="white", linewidth=1.4, zorder=4)
        ax.plot([0, rc], [i, i], color=RIGHT_COL, lw=3.5, solid_capstyle="round", zorder=3)
        ax.scatter([rc], [i], s=160, color=RIGHT_COL, edgecolor="white", linewidth=1.4, zorder=4)
        ax.plot([-lc*0.07, rc*0.07], [i, i], color="#d9e2ec", lw=8, alpha=0.65, zorder=2, solid_capstyle="round")
        ax.text(-max_count*0.03, i+0.18, wrap_label(d["left_label"], 18), ha="right", va="center", fontsize=12.5, fontweight="bold")
        ax.text(max_count*0.03, i+0.18, wrap_label(d["right_label"], 18), ha="left", va="center", fontsize=12.5, fontweight="bold")
        ax.text(-lc, i+0.26, fmt(lc), ha="center", va="bottom", fontsize=10.5, fontweight="bold", color=LEFT_COL)
        ax.text(rc, i+0.26, fmt(rc), ha="center", va="bottom", fontsize=10.5, fontweight="bold", color="#249f73")
        ax.text(-max_count*0.03, i-0.20, d["left_source"], ha="right", va="center", fontsize=9.6, color=SOURCE_COLORS[d["left_source"]])
        ax.text(max_count*0.03, i-0.20, d["right_source"], ha="left", va="center", fontsize=9.6, color=SOURCE_COLORS[d["right_source"]])

    ax.text(-max_count*0.78, len(data)-0.47, "Left set", ha="center", fontsize=11, color=MUTED)
    ax.text(max_count*0.78, len(data)-0.47, "Right set", ha="center", fontsize=11, color=MUTED)
    ax.set_xticks([]); ax.set_yticks([])
    for s in ax.spines.values():
        s.set_visible(False)
    return fig

def bezier_patch(x0, y0, x1, y1, curv=0.24, color="#999", lw=2.0, alpha=0.5):
    mx = (x0 + x1) / 2
    path = MplPath(
        [(x0, y0), (mx-curv, y0), (mx+curv, y1), (x1, y1)],
        [MplPath.MOVETO, MplPath.CURVE4, MplPath.CURVE4, MplPath.CURVE4]
    )
    return PathPatch(path, facecolor="none", edgecolor=color, lw=lw, alpha=alpha, capstyle="round")

def plot_bipartite(data):
    def dedup(nodes):
        out, seen = [], {}
        for a, b, c in nodes:
            if a not in seen:
                seen[a] = (b, c)
                out.append((a, b, c))
        return out

    left_nodes = sorted(dedup([(d["left_label"], d["left_count"], d["left_source"]) for d in data]), key=lambda x: (-x[1], x[0]))
    right_nodes = sorted(dedup([(d["right_label"], d["right_count"], d["right_source"]) for d in data]), key=lambda x: (-x[1], x[0]))
    left_y = np.linspace(0.80, 0.18, len(left_nodes))
    right_y = np.linspace(0.80, 0.18, len(right_nodes))
    left_pos = {n[0]: (0.20, y) for n, y in zip(left_nodes, left_y)}
    right_pos = {n[0]: (0.80, y) for n, y in zip(right_nodes, right_y)}
    max_count = max([n[1] for n in left_nodes + right_nodes])

    def size(c): return 0.035 + 0.07*np.sqrt(c/max_count)

    fig, ax = plt.subplots(figsize=(11.5, 8.8))
    fig.subplots_adjust(top=0.90, left=0.03, right=0.97, bottom=0.06)
    ax.set_xlim(0, 1); ax.set_ylim(0, 1)

    ax.text(0.06, 0.975, "Bipartite relation map", fontsize=22, fontweight="bold", ha="left", va="top")
    ax.text(0.06, 0.942, "Two disjoint sets connected only by the observed matched pairs", fontsize=10.5, color=MUTED, ha="left", va="top")
    ax.add_patch(FancyBboxPatch((0.06,0.08),0.30,0.78, boxstyle="round,pad=0.02,rounding_size=0.04", fc="#f8fbff", ec="#e5edf5", lw=1))
    ax.add_patch(FancyBboxPatch((0.64,0.08),0.30,0.78, boxstyle="round,pad=0.02,rounding_size=0.04", fc="#fbfffc", ec="#e5edf5", lw=1))
    ax.text(0.21, 0.855, "Left set", ha="center", fontsize=11, color=MUTED)
    ax.text(0.79, 0.855, "Right set", ha="center", fontsize=11, color=MUTED)

    for d in data:
        x0, y0 = left_pos[d["left_label"]]
        x1, y1 = right_pos[d["right_label"]]
        ax.add_patch(bezier_patch(x0+size(d["left_count"]), y0, x1-size(d["right_count"]), y1, color="#93a4b8", lw=2.3, alpha=0.55))

    def draw(nodes, pos, align):
        for label, count, source in nodes:
            x, y = pos[label]
            r = size(count)
            ax.add_patch(Circle((x, y), r, facecolor=SOURCE_COLORS[source], edgecolor="white", lw=2.2, alpha=0.95))
            ax.text(x, y, fmt(count), ha="center", va="center", fontsize=10.5, fontweight="bold", color="white")
            if align == "left":
                ax.text(x+r+0.02, y+0.018, clean_label(label), ha="left", va="center", fontsize=12.6, fontweight="bold")
                ax.text(x+r+0.02, y-0.020, source, ha="left", va="center", fontsize=9.6, color=SOURCE_COLORS[source])
            else:
                ax.text(x-r-0.02, y+0.018, clean_label(label), ha="right", va="center", fontsize=12.6, fontweight="bold")
                ax.text(x-r-0.02, y-0.020, source, ha="right", va="center", fontsize=9.6, color=SOURCE_COLORS[source])

    draw(left_nodes, left_pos, "left")
    draw(right_nodes, right_pos, "right")

    lx, ly = 0.39, 0.14
    ax.text(lx, ly+0.11, "Source", fontsize=10.5, color=MUTED, ha="left")
    for i, (src, col) in enumerate(SOURCE_COLORS.items()):
        yy = ly+0.07 - i*0.045
        ax.scatter([lx], [yy], s=110, color=col, edgecolor="white", linewidth=1.2)
        ax.text(lx+0.025, yy, src, ha="left", va="center", fontsize=9.6)

    ax.axis("off")
    return fig

def plot_story_cards(data):
    fig = plt.figure(figsize=(14.5, 8.6))
    gs = GridSpec(2, 2, figure=fig, wspace=0.16, hspace=0.18)
    fig.text(0.055, 0.955, "Pair cards for PPT / figure assembly", fontsize=22, fontweight="bold", ha="left")
    fig.text(0.055, 0.925, "Each card keeps category, count, literature, and the pair linkage in one compact object", fontsize=10.5, color=MUTED, ha="left")
    max_count = max(max(x["left_count"], x["right_count"]) for x in data)

    for idx, d in enumerate(data):
        ax = fig.add_subplot(gs[idx//2, idx%2])
        ax.set_xlim(0, 1); ax.set_ylim(0, 1)
        ax.add_patch(FancyBboxPatch((0.02,0.06),0.96,0.88, boxstyle="round,pad=0.02,rounding_size=0.05", fc="#fcfdff", ec="#e8edf3", lw=1.2))
        ax.plot([0.28,0.72],[0.50,0.50],color="#d4dde8",lw=7,solid_capstyle="round",zorder=1)
        ax.plot([0.28,0.72],[0.50,0.50],color="#b7c7d8",lw=1.8,solid_capstyle="round",zorder=2)
        r1 = 0.08 + 0.08*np.sqrt(d["left_count"]/max_count)
        r2 = 0.08 + 0.08*np.sqrt(d["right_count"]/max_count)
        ax.add_patch(Circle((0.24,0.50),r1,facecolor=LEFT_COL,edgecolor="white",lw=2.5))
        ax.add_patch(Circle((0.76,0.50),r2,facecolor=RIGHT_COL,edgecolor="white",lw=2.5))
        ax.text(0.24,0.50,fmt(d["left_count"]),ha="center",va="center",fontsize=12,fontweight="bold",color="white")
        ax.text(0.76,0.50,fmt(d["right_count"]),ha="center",va="center",fontsize=12,fontweight="bold",color="white")
        ax.text(0.24,0.76,wrap_label(d["left_label"],18),ha="center",va="center",fontsize=13,fontweight="bold")
        ax.text(0.76,0.76,wrap_label(d["right_label"],18),ha="center",va="center",fontsize=13,fontweight="bold")
        for x, src in [(0.24,d["left_source"]), (0.76,d["right_source"])]:
            w = 0.22 if len(src) < 15 else 0.27
            ax.add_patch(FancyBboxPatch((x-w/2,0.17),w,0.10, boxstyle="round,pad=0.02,rounding_size=0.06", fc=SOURCE_COLORS[src], ec="none", alpha=0.14))
            ax.text(x,0.22,src,ha="center",va="center",fontsize=9.3,color=SOURCE_COLORS[src],fontweight="bold")
        ax.text(0.50,0.50,"matched pair",ha="center",va="center",fontsize=10,color=MUTED,fontweight="bold",
                bbox=dict(boxstyle="round,pad=0.18,rounding_size=0.1",fc="white",ec="#e6edf5",lw=1))
        ax.axis("off")
    return fig

def plot_dashboard(data):
    fig, ax = plt.subplots(figsize=(13.5, 8.3))
    ax.set_xlim(0, 1); ax.set_ylim(0, 1); ax.axis("off")
    ax.text(0.05,0.95,"Structured pair dashboard",fontsize=22,fontweight="bold",ha="left")
    ax.text(0.05,0.918,"A publication-style layout balancing figure and table aesthetics",fontsize=10.5,color=MUTED,ha="left")
    cols = {"left_type":0.09,"left_n":0.31,"link":0.50,"right_type":0.69,"right_n":0.91}
    ax.text(cols["left_type"],0.84,"Left category",fontsize=11,color=MUTED,ha="center")
    ax.text(cols["left_n"],0.84,"Count",fontsize=11,color=MUTED,ha="center")
    ax.text(cols["link"],0.84,"Relation",fontsize=11,color=MUTED,ha="center")
    ax.text(cols["right_type"],0.84,"Right category",fontsize=11,color=MUTED,ha="center")
    ax.text(cols["right_n"],0.84,"Count",fontsize=11,color=MUTED,ha="center")
    ax.plot([0.05,0.95],[0.81,0.81],color="#d8e0ea",lw=1.1)

    ys = np.linspace(0.70,0.16,len(data))
    max_count = max(max(d["left_count"], d["right_count"]) for d in data)
    for i, (y, d) in enumerate(zip(ys, data)):
        if i % 2 == 0:
            ax.add_patch(Rectangle((0.05,y-0.085),0.90,0.14,facecolor="#fafcff",edgecolor="none"))
        ax.text(cols["left_type"],y+0.016,wrap_label(d["left_label"],17),ha="center",va="center",fontsize=12.6,fontweight="bold")
        ax.text(cols["left_type"],y-0.040,d["left_source"],ha="center",va="center",fontsize=9.4,color=SOURCE_COLORS[d["left_source"]])
        ax.text(cols["right_type"],y+0.016,wrap_label(d["right_label"],17),ha="center",va="center",fontsize=12.6,fontweight="bold")
        ax.text(cols["right_type"],y-0.040,d["right_source"],ha="center",va="center",fontsize=9.4,color=SOURCE_COLORS[d["right_source"]])
        bl = 0.13 * d["left_count"]/max_count
        br = 0.13 * d["right_count"]/max_count
        ax.plot([cols["left_n"]-bl, cols["left_n"]],[y,y],color=LEFT_COL,lw=9,solid_capstyle="round")
        ax.scatter([cols["left_n"]-bl],[y],s=120,color=LEFT_COL,edgecolor="white",linewidth=1.2,zorder=3)
        ax.text(cols["left_n"]+0.018,y,fmt(d["left_count"]),ha="left",va="center",fontsize=11.2,fontweight="bold")
        ax.plot([cols["right_n"], cols["right_n"]+br],[y,y],color=RIGHT_COL,lw=9,solid_capstyle="round")
        ax.scatter([cols["right_n"]+br],[y],s=120,color=RIGHT_COL,edgecolor="white",linewidth=1.2,zorder=3)
        ax.text(cols["right_n"]-0.018,y,fmt(d["right_count"]),ha="right",va="center",fontsize=11.2,fontweight="bold")
        ax.add_patch(FancyBboxPatch((cols["link"]-0.07,y-0.028),0.14,0.056,boxstyle="round,pad=0.02,rounding_size=0.03",fc="white",ec="#dfe7ef",lw=1))
        ax.text(cols["link"],y,"paired",ha="center",va="center",fontsize=10.5,color=MUTED,fontweight="bold")
        ax.plot([cols["left_n"]+0.055, cols["link"]-0.075],[y,y],color="#d8e1eb",lw=2)
        ax.plot([cols["link"]+0.075, cols["right_n"]-0.055],[y,y],color="#d8e1eb",lw=2)
    return fig

def make_overview():
    thumbs = [
        "01_mirror_lollipop_pairs.png",
        "02_bipartite_relation_map.png",
        "03_pair_story_cards.png",
        "04_structured_pair_dashboard.png",
    ]
    labels = [
        "Option A  Mirror lollipop  (best overall)",
        "Option B  Bipartite relation map",
        "Option C  Pair story cards",
        "Option D  Structured pair dashboard",
    ]
    imgs = [Image.open(OUTPUT_DIR / t).convert("RGB") for t in thumbs]
    canvas = Image.new("RGB", (2400, 1700), "white")
    positions = [(60,90), (1230,90), (60,880), (1230,880)]
    draw = ImageDraw.Draw(canvas)
    draw.text((60,28), "Sheep pair visualization overview", fill=(24,33,42))
    draw.text((60,56), "Same data rendered in four presentation styles; all exported as PNG and PDF", fill=(110,125,140))
    for img, pos, label in zip(imgs, positions, labels):
        box_w, box_h = 1040, 680
        img.thumbnail((box_w, box_h))
        x = pos[0] + (box_w - img.width)//2
        y = pos[1] + 40 + (box_h - img.height)//2
        canvas.paste(img, (x, y))
        draw.text((pos[0], pos[1]), label, fill=(24,33,42))
    canvas.save(OUTPUT_DIR / "00_overview_grid.png")
    canvas.convert("RGB").save(OUTPUT_DIR / "00_overview_grid.pdf")

def main():
    save_figure(plot_mirror_lollipop(DATA), "01_mirror_lollipop_pairs")
    save_figure(plot_bipartite(DATA), "02_bipartite_relation_map")
    save_figure(plot_story_cards(DATA), "03_pair_story_cards")
    save_figure(plot_dashboard(DATA), "04_structured_pair_dashboard")
    make_overview()
    print(f"Done. Files written to: {OUTPUT_DIR}")

if __name__ == "__main__":
    main()
