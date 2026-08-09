import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib import ticker
import scipy.stats as stats
from matplotlib.patches import Patch
from matplotlib.font_manager import FontProperties


# Processing note.
file_path = r"/vol2/mengzhu/SheepFANNG/04_ChromHMM_noblacklist_modif/Merge_chromatin_state/state_variability/AA_super_enhancer/AA_super_result_1/all_super_enhancer_summary_1.txt"

df = pd.read_csv(file_path, sep="\t", header=0)

df_selected = df[["Tissue", "E5", "E6", "E7", "E8"]].copy()

df_melted = pd.melt(
    df_selected,
    id_vars=["Tissue"],
    var_name="state",
    value_name="V2"
)

state_mapping = {
    "E5": "5 EnhA",
    "E6": "6 EnhAMe",
    "E7": "7 EnhAHet",
    "E8": "8 EnhPois"
}

df_melted["state"] = df_melted["state"].map(state_mapping)

state_order = ["5 EnhA", "6 EnhAMe", "7 EnhAHet", "8 EnhPois"]
df_melted["state"] = pd.Categorical(
    df_melted["state"],
    categories=state_order,
    ordered=True
)

categories = state_order

data_by_state = {
    state: df_melted[df_melted["state"] == state]["V2"].dropna().values
    for state in categories
}


# Processing note.
plt.rcParams["font.family"] = "Times New Roman"
plt.rcParams["axes.unicode_minus"] = False

np.random.seed(123)


# Processing note.
# Processing note.
# Processing note.
# Processing note.
# Processing note.
MM_PER_INCH = 25.4

AX_WIDTH_MM = 85.533
AX_HEIGHT_MM = 63.584

# Processing note.
# Processing note.
LEFT_MARGIN_MM = 22.0
RIGHT_MARGIN_MM = 48.0
BOTTOM_MARGIN_MM = 30.0
TOP_MARGIN_MM = 16.0

FIG_WIDTH_MM = LEFT_MARGIN_MM + AX_WIDTH_MM + RIGHT_MARGIN_MM
FIG_HEIGHT_MM = BOTTOM_MARGIN_MM + AX_HEIGHT_MM + TOP_MARGIN_MM

fig = plt.figure(
    figsize=(FIG_WIDTH_MM / MM_PER_INCH, FIG_HEIGHT_MM / MM_PER_INCH)
)

ax = fig.add_axes([
    LEFT_MARGIN_MM / FIG_WIDTH_MM,
    BOTTOM_MARGIN_MM / FIG_HEIGHT_MM,
    AX_WIDTH_MM / FIG_WIDTH_MM,
    AX_HEIGHT_MM / FIG_HEIGHT_MM
])


color_mapping = {
    "5 EnhA": "#FFFF00",
    "6 EnhAMe": "#D6A32B",
    "7 EnhAHet": "#F0C36B",
    "8 EnhPois": "#F5E1A1"
}

box_colors = [color_mapping[state] for state in categories]
violin_colors = box_colors

positions = np.arange(len(categories))
state_to_x = dict(zip(categories, positions))

box_width = 0.15
violin_width = 0.5


# Processing note.
def get_signif_mark(p_val):
    if pd.isna(p_val):
        return "ns"
    elif p_val < 0.001:
        return "***"
    elif p_val < 0.01:
        return "**"
    elif p_val < 0.05:
        return "*"
    else:
        return "ns"


ref_group = "8 EnhPois"
ref_data = data_by_state[ref_group]

# Processing note.
# Processing note.
# Processing note.
SHOW_NS_BRACKETS = True

comparison_results = {}

for state in categories:
    if state == ref_group:
        continue

    group_data = data_by_state[state]

    # Processing note.
    if len(group_data) < 2 or len(ref_data) < 2:
        p_val = np.nan
    else:
        t_stat, p_val = stats.ttest_ind(
            group_data,
            ref_data,
            equal_var=False,
            nan_policy="omit"
        )

    comparison_results[(state, ref_group)] = {
        "p_value": p_val,
        "mark": get_signif_mark(p_val)
    }


# Processing note.
for i, state in enumerate(categories):
    data_points = data_by_state[state]

    if len(data_points) == 0:
        continue

    # Processing note.
    box_pos = positions[i] - box_width / 100

    ax.boxplot(
        data_points,
        positions=[box_pos],
        widths=box_width,
        patch_artist=True,
        showfliers=False,
        notch=True,
        medianprops={
            "color": "black",
            "linewidth": 1.5
        },
        boxprops={
            "facecolor": box_colors[i],
            "edgecolor": violin_colors[i],
            "linewidth": 1.5
        },
        whiskerprops={
            "color": violin_colors[i],
            "linewidth": 1.5
        },
        capprops={
            "color": violin_colors[i],
            "linewidth": 1.5
        }
    )

    # Processing note.
    violin_pos = positions[i] + box_width / 50

    violin = ax.violinplot(
        data_points,
        positions=[violin_pos],
        widths=violin_width,
        showmeans=False,
        showmedians=False,
        showextrema=False
    )

    for pc in violin["bodies"]:
        pc.set_facecolor(violin_colors[i])
        pc.set_edgecolor(violin_colors[i])
        pc.set_alpha(0.35)

        # Processing note.
        vertices = pc.get_paths()[0].vertices
        vertices[:, 0] = np.where(
            vertices[:, 0] > violin_pos,
            vertices[:, 0],
            violin_pos
        )

    # Processing note.
    ax.scatter(
        np.random.normal(positions[i] - box_width, 0.04, len(data_points)),
        data_points,
        color=violin_colors[i],
        alpha=0.8,
        s=20,
        edgecolor="white",
        linewidth=0.5,
        zorder=3
    )


# Processing note.
ax.yaxis.set_major_formatter(ticker.FormatStrFormatter("%.3f"))

ax.set_xticks(positions)
ax.set_xticklabels(
    categories,
    fontsize=15,
    rotation=60,
    ha="right",
    color="black"
)

ax.set_xlabel("", fontsize=15)
ax.set_ylabel("Enrichment", fontsize=15, color="black")

ax.tick_params(axis="y", labelsize=15, color="black")
ax.tick_params(axis="x", color="black")

ax.grid(False)

for spine in ["top", "right", "bottom", "left"]:
    ax.spines[spine].set_visible(True)
    ax.spines[spine].set_linewidth(1)
    ax.spines[spine].set_color("black")


# Processing note.
def add_sig_bracket(
    ax,
    x1,
    x2,
    y,
    h,
    text,
    text_offset,
    line_width=1.2,
    fontsize=13
):
    """
    Add a significance bracket and label between two x positions.

    x1, x2: x-axis positions of the two groups
    y: height of the bracket base
    h: height of the bracket stems
    text: significance label, such as ns, *, **, or ***
    text_offset: distance between the text and bracket top
    """

    ax.plot(
        [x1, x1, x2, x2],
        [y, y + h, y + h, y],
        color="black",
        linewidth=line_width,
        clip_on=False,
        zorder=6
    )

    ax.text(
        (x1 + x2) / 2,
        y + h + text_offset,
        text,
        ha="center",
        va="bottom",
        fontsize=fontsize,
        weight="bold",
        color="black",
        zorder=7
    )


# Processing note.
comparisons_to_plot = [
    ("5 EnhA", "8 EnhPois"),
    ("6 EnhAMe", "8 EnhPois"),
    ("7 EnhAHet", "8 EnhPois")
]

# Processing note.
comparisons_to_plot = sorted(
    comparisons_to_plot,
    key=lambda pair: abs(state_to_x[pair[0]] - state_to_x[pair[1]])
)

# Processing note.
y_min, y_max = ax.get_ylim()
y_range = y_max - y_min

# Processing note.
bracket_start_y = y_max + y_range * 0.04
bracket_height = y_range * 0.018
bracket_gap = y_range * 0.090
text_offset = y_range * 0.008

used_bracket_count = 0

for group1, group2 in comparisons_to_plot:
    result = comparison_results.get((group1, group2), None)

    if result is None:
        continue

    mark = result["mark"]

    if (not SHOW_NS_BRACKETS) and mark == "ns":
        continue

    x1 = state_to_x[group1]
    x2 = state_to_x[group2]

    y = bracket_start_y + used_bracket_count * bracket_gap

    add_sig_bracket(
        ax=ax,
        x1=x1,
        x2=x2,
        y=y,
        h=bracket_height,
        text=mark,
        text_offset=text_offset,
        line_width=1.2,
        fontsize=13
    )

    used_bracket_count += 1

# Processing note.
if used_bracket_count > 0:
    final_top = (
        bracket_start_y
        + (used_bracket_count - 1) * bracket_gap
        + bracket_height
        + text_offset
        + y_range * 0.06
    )
    ax.set_ylim(bottom=y_min, top=final_top)


# Processing note.
legend_elements = [
    Patch(
        facecolor=color_mapping[state],
        edgecolor=color_mapping[state],
        label=state
    )
    for state in state_order
]

title_font = FontProperties(
    family="Times New Roman",
    size=16,
    weight="normal"
)

legend = ax.legend(
    handles=legend_elements,
    title="state",
    loc="center left",
    bbox_to_anchor=(1.02, 0.5),
    fontsize=16,
    prop={"family": "Times New Roman"},
    labelcolor="black",
    title_fontproperties=title_font,
    frameon=False,
    framealpha=0,
    labelspacing=1.2,
    handlelength=1.5,
    handletextpad=0.8,
    borderaxespad=0.5
)

plt.setp(legend.get_title(), color="black")


# Processing note.
# Processing note.
# Processing note.
# Processing note.
save_path = r"/vol2/wulingyun/superenhancer/Graph_beautification/05_yunyu_super_enhancer_E_boxplot.pdf"

plt.savefig(
    save_path,
    dpi=300
)

plt.close()

print(f"informationFigure saved to: {save_path}")