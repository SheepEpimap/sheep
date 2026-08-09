# Processing note.
# Processing note.
# Processing note.
# Processing note.

import os
import pandas as pd
import plotly.graph_objects as go
import plotly.io as pio
from scipy.stats import ttest_ind

# =========================
# Processing note.
# =========================
try:
    pio.kaleido.scope.mathjax = None
except Exception:
    pass

try:
    pio.defaults.mathjax = None
except Exception:
    pass


# =========================
# Processing note.
# =========================
input_dir = (
    '/vol2/mengzhu/SheepFANNG/'
    '04_ChromHMM_noblacklist_modif/'
    'Merge_chromatin_state/'
    'state_variability/'
    'AA_super_enhancer/'
    'AA_super_result_1/'
)

normal_file = os.path.join(
    input_dir,
    'all_super_enhancer_summary_normal.txt'
)

super_file = os.path.join(
    input_dir,
    'all_super_enhancer_summary_1.txt'
)

output_dir = (
    '/vol2/wulingyun/'
    'superenhancer/'
    'Graph_beautification/'
)

os.makedirs(output_dir, exist_ok=True)

output_image = os.path.join(
    output_dir,
    '10_fenzuyunyu_normal_enhancer_E5_E8_statistics.pdf'
)


# =========================
# Processing note.
# =========================
df_normal = pd.read_csv(
    normal_file,
    sep='\t',
    header=0
)

df_normal_sub = df_normal[
    ['Tissue', 'E5', 'E6', 'E7', 'E8']
].copy()

df_normal_melt = pd.melt(
    df_normal_sub,
    id_vars=['Tissue'],
    value_vars=['E5', 'E6', 'E7', 'E8'],
    var_name='state',
    value_name='V2'
)

df_normal_melt['group'] = 'normal'


# =========================
# Processing note.
# =========================
df_super = pd.read_csv(
    super_file,
    sep='\t',
    header=0
)

df_super_sub = df_super[
    ['Tissue', 'E5', 'E6', 'E7', 'E8']
].copy()

df_super_melt = pd.melt(
    df_super_sub,
    id_vars=['Tissue'],
    value_vars=['E5', 'E6', 'E7', 'E8'],
    var_name='state',
    value_name='V2'
)

df_super_melt['group'] = 'super-enhancer'


# =========================
# Processing note.
# =========================
df_combined = pd.concat(
    [df_normal_melt, df_super_melt],
    ignore_index=True
)

state_mapping = {
    'E5': '5 EnhA',
    'E6': '6 EnhAMe',
    'E7': '7 EnhAHet',
    'E8': '8 EnhPois'
}

df_combined['state'] = df_combined['state'].map(
    state_mapping
)

state_order = [
    '5 EnhA',
    '6 EnhAMe',
    '7 EnhAHet',
    '8 EnhPois'
]

df_combined['state'] = pd.Categorical(
    df_combined['state'],
    categories=state_order,
    ordered=True
)

df_combined['V2'] = pd.to_numeric(
    df_combined['V2'],
    errors='coerce'
)


# =========================
# Processing note.
# Processing note.
# normal vs super-enhancer
# Welch two-sample t-test
# =========================
def p_to_signif(p_value):
    """Convert a p-value to significance stars."""

    if p_value is None or pd.isna(p_value):
        return 'NA'

    if p_value < 1e-4:
        return '****'
    elif p_value < 1e-3:
        return '***'
    elif p_value < 1e-2:
        return '**'
    elif p_value < 0.05:
        return '*'
    else:
        return 'ns'


signif_map = {}

for state_name in state_order:

    normal_values = df_combined.loc[
        (df_combined['state'] == state_name) &
        (df_combined['group'] == 'normal'),
        'V2'
    ].dropna()

    super_values = df_combined.loc[
        (df_combined['state'] == state_name) &
        (df_combined['group'] == 'super-enhancer'),
        'V2'
    ].dropna()

    result = ttest_ind(
        normal_values,
        super_values,
        equal_var=False,
        nan_policy='omit'
    )

    p_value = (
        float(result.pvalue)
        if result.pvalue is not None
        else None
    )

    signif_map[state_name] = (
        p_value,
        p_to_signif(p_value)
    )


# =========================
# Processing note.
# Processing note.
# =========================

# Processing note.
color_map = {
    'super-enhancer': '#0570b0',
    'normal': '#808080'
}

# Processing note.
# Processing note.
fill_color_map = {
    'super-enhancer': 'rgba(5, 112, 176, 0.35)',
    'normal': 'rgba(128, 128, 128, 0.35)'
}

# Processing note.
group_order = [
    'normal',
    'super-enhancer'
]


# =========================
# Processing note.
# =========================
fig = go.Figure()

for group_name in group_order:

    group_data = df_combined[
        df_combined['group'] == group_name
    ]

    fig.add_trace(
        go.Violin(
            x=group_data['state'],
            y=group_data['V2'],

            name=group_name,
            legendgroup=group_name,
            offsetgroup=group_name,

            # Processing note.
            side='positive',

            # Processing note.
            points='all',
            pointpos=-0.6,
            jitter=0.5,

            width=0.8,

            # Processing note.
            box_visible=False,

            # Processing note.
            meanline_visible=True,

            # Processing note.
            line=dict(
                color=color_map[group_name],
                width=1.5
            ),

            # Processing note.
            fillcolor=fill_color_map[group_name],

            # Processing note.
            opacity=1.0,

            # Processing note.
            marker=dict(
                size=4,
                color=color_map[group_name],
                opacity=0.60,
                line=dict(
                    color=color_map[group_name],
                    width=0.8
                )
            ),

            hovertemplate=(
                'State: %{x}<br>'
                'Percent: %{y:.4f}'
                '<extra>' + group_name + '</extra>'
            )
        )
    )


# =========================
# Processing note.
# =========================
if df_combined['V2'].notna().any():
    y_max = float(
        df_combined['V2'].max(skipna=True)
    )
else:
    y_max = 1.0

# Processing note.
y_top = max(
    1.10,
    y_max * 1.05
)

# Processing note.
y_bottom = -0.05

# Processing note.
y_sig = 1.02


# =========================
# Processing note.
# =========================
for state_name in state_order:

    p_value, significance_label = signif_map[state_name]

    fig.add_annotation(
        x=state_name,
        xref='x',

        y=y_sig,
        yref='y',

        text=significance_label,
        showarrow=False,

        xanchor='center',
        yanchor='bottom',

        font=dict(
            family='Times New Roman',
            size=16,
            color='black'
        )
    )


# =========================
# Processing note.
# =========================
fig.update_layout(

    width=900,
    height=650,

    # Processing note.
    violinmode='group',
    violingroupgap=0.20,

    margin=dict(
        l=90,
        r=160,
        t=60,
        b=140
    ),

    # Processing note.
    xaxis=dict(
        title='state',

        categoryorder='array',
        categoryarray=state_order,

        tickangle=60,

        tickfont=dict(
            family='Times New Roman',
            size=15,
            color='black'
        ),

        title_font=dict(
            family='Times New Roman',
            size=18,
            color='black'
        ),

        showgrid=False,
        zeroline=False,

        showline=True,
        linewidth=1,
        linecolor='black',

        ticks='outside',
        ticklen=6,
        tickwidth=1,
        tickcolor='black',

        mirror=True
    ),

    # Processing note.
    yaxis=dict(
        title='Percent of regulator in enhancer',

        range=[
            y_bottom,
            y_top
        ],

        tickformat='.2f',
        tick0=0,
        dtick=0.25,

        tickfont=dict(
            family='Times New Roman',
            size=15,
            color='black'
        ),

        title_font=dict(
            family='Times New Roman',
            size=18,
            color='black'
        ),

        showgrid=False,
        zeroline=False,

        showline=True,
        linewidth=1,
        linecolor='black',

        ticks='outside',
        ticklen=6,
        tickwidth=1,
        tickcolor='black',

        mirror=True
    ),

    # Processing note.
    legend=dict(
        title='',

        x=1.02,
        xanchor='left',

        y=0.5,
        yanchor='middle',

        traceorder='normal',

        font=dict(
            family='Times New Roman',
            size=16,
            color='black'
        ),

        borderwidth=0,

        itemsizing='constant'
    ),

    plot_bgcolor='white',
    paper_bgcolor='white',

    font=dict(
        family='Times New Roman',
        color='black'
    )
)


# =========================
# Processing note.
# =========================
pio.write_image(
    fig,
    output_image,
    format='pdf'
)

print(f'Figure saved to: {output_image}')


# =========================
# Processing note.
# =========================
print('\nSignificance test results:')
print('state\tp_value\tsignificance')

for state_name in state_order:

    p_value, significance_label = signif_map[state_name]

    if p_value is None or pd.isna(p_value):
        p_text = 'NA'
    else:
        p_text = f'{p_value:.6g}'

    print(
        f'{state_name}\t'
        f'{p_text}\t'
        f'{significance_label}'
    )