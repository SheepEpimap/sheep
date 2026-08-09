# Processing note.
# Processing note.
# Processing note.

# Processing note.
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.ticker as mticker  # Processing note.
from scipy import stats  # Processing note.
from matplotlib.lines import Line2D  # Processing note.
import warnings
warnings.filterwarnings('ignore')

# Processing note.
plt.rcParams['font.family'] = ['DejaVu Serif', 'serif']
plt.rcParams['axes.unicode_minus'] = False

# Processing note.
file_path = r'/vol2/mengzhu/SheepFANNG/04_ChromHMM_noblacklist_modif/Merge_chromatin_state/state_variability/AA_super_enhancer/AA_super_result_1/AA_normal_enhancer/all_super_enhancer_summary_2.txt'
data = pd.read_csv(file_path, sep="\t", header=0)

# Processing note.
print("==== Reproduce the R script data logic ====")
print("Data column names (matching R colnames(data)):")
print(data.columns.tolist())

# Processing note.
data1 = data[["Tissue", "Normal_percent", "Super_percent"]].copy()

# Processing note.
data2 = pd.melt(
    data1,
    id_vars=["Tissue"],
    value_vars=["Normal_percent", "Super_percent"],
    var_name="state",
    value_name="V2"
)

# Processing note.
categories = ["Super_percent", "Normal_percent"]
super_data = data2[data2["state"] == "Super_percent"]["V2"].dropna().values
normal_data = data2[data2["state"] == "Normal_percent"]["V2"].dropna().values

# Processing note.
t_stat, p_val = stats.ttest_ind(
    super_data, 
    normal_data, 
    equal_var=True,  # Processing note.
    nan_policy='omit'
)

# Processing note.
machine_eps = np.finfo(np.float64).eps  # Processing note.
if p_val < machine_eps:
    p_label = f'T-test, p < 2.2e-16'
else:
    p_label = f'T-test, p = {p_val:.3g}'

# Processing note.
print(f"\n[T-test results matching R]")
print(f"t statistic: {t_stat:.4f}")
print(f"raw p-value: {p_val:.2e}")
print(f"display label: {p_label}")

# Processing note.
fig, ax = plt.subplots(figsize=(4*1.5*1.5, 5*2))  # Processing note.

# Processing note.
colors = {
    "Super_percent": "#0570b0",
    "Normal_percent": "grey"
}
box_colors = [colors["Super_percent"], colors["Normal_percent"]]
violin_colors = [colors["Super_percent"], colors["Normal_percent"]]
positions = np.arange(len(categories))
box_width = 0.15
violin_width = 0.5

# Processing note.
for i, category in enumerate(categories):
    data_points = data2[data2["state"] == category]["V2"].dropna().values
    
    # Processing note.
    box_pos = positions[i] - box_width / 100
    bp = ax.boxplot(
        data_points,
        positions=[box_pos], widths=box_width,
        patch_artist=True, showfliers=False, 
        notch=True,  # Processing note.
        medianprops={'color': 'black', 'linewidth': 1.5},
        boxprops={'facecolor': box_colors[i], 'edgecolor': violin_colors[i], 'linewidth': 1.5},
        whiskerprops={'color': violin_colors[i], 'linewidth': 1.5},
        capprops={'color': violin_colors[i], 'linewidth': 1.5}
    )
    
    # Processing note.
    violin_pos = positions[i] + box_width / 50
    violin = ax.violinplot(
        data_points,
        positions=[violin_pos], widths=violin_width,
        showmeans=False, showmedians=False, showextrema=False
    )
    for pc in violin['bodies']:
        pc.set_facecolor(violin_colors[i])
        pc.set_edgecolor(violin_colors[i])
        pc.set_alpha(0.35)
        # Processing note.
        vertices = pc.get_paths()[0].vertices
        vertices[:, 0] = np.where(vertices[:, 0] > violin_pos, vertices[:, 0], violin_pos)
    
    # Processing note.
    ax.scatter(
        np.random.normal(positions[i] - box_width, 0.04, len(data_points)),
        data_points,
        color=violin_colors[i], alpha=0.6, s=12,
        edgecolor='none', zorder=3
    )

# Processing note.
# Processing note.
for spine in ['top', 'right', 'bottom', 'left']:
    ax.spines[spine].set_visible(True)
    ax.spines[spine].set_linewidth(1)

# Processing note.
ax.set_xticks(positions)
ax.set_xticklabels(categories, fontsize=13, ha='right')
ax.tick_params(axis='x', rotation=60, pad=8)

# Processing note.
ax.set_ylim(-0.05, 1.05)  # Processing note.
ax.yaxis.set_major_locator(mticker.MultipleLocator(0.1))  # Processing note.
ax.set_ylabel('Percent of enhancer', fontsize=15)
ax.set_yticklabels([f'{y:.1f}' for y in ax.get_yticks()], fontsize=15)

# Processing note.
ax.grid(False)

# Processing note.
handles = [
    Line2D([0], [0], color="#0570b0", marker='o', linestyle='None', markersize=8, label='Super_percent'),
    Line2D([0], [0], color="grey", marker='o', linestyle='None', markersize=8, label='Normal_percent')
]
ax.legend(
    handles=handles,
    loc='center left',  # Processing note.
    bbox_to_anchor=(1.02, 0.5),  # Processing note.
    fontsize=16, 
    title='state',
    borderaxespad=0.0,  # Processing note.
    frameon=False,  # Processing note.
    labelspacing=0.8,  # Processing note.
    handlelength=1.5  # Processing note.
)

# Processing note.
ax.text(
    x=0.5, y=0.98,  # Processing note.
    s=p_label,
    ha='center', va='bottom',  # Processing note.
    fontsize=12,
    bbox=dict(boxstyle='round,pad=0.3', facecolor='white', edgecolor='none', alpha=0.8)
)

# Processing note.
plt.tight_layout(rect=[0, 0, 0.85, 1])  # Processing note.
output_path = r'/vol2/wulingyun/superenhancer/Graph_beautification/08_yunyu_all_super_enhancer_summary_normal.pdf'
plt.savefig(output_path, dpi=300, bbox_inches='tight')
plt.close()

print(f"\n✅ Figure saved to: {output_path}")
print(f"✅ 1. yinformation(-0.05~1.05)；2. T-testinformation1.0information；3. information1.5information；4. information！")