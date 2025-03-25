import pandas as pd
import plotly.express as px

# 1) Read the aggregated results
df = pd.read_csv("output/simulations/aggregated_area_change.csv")

# 2) Exclude results if 'holes_file' does NOT contain 'contract'
df = df[df["holes_file"].str.contains("contract")]

# 3) Filter for rows where time == 10
df_t10 = df[df["time"] == 10].copy()

# 4) Create an interactive scatter plot
fig = px.scatter(
    df_t10,
    x="max_shrinking",
    y="overall_relative_area_change",
    color="pattern",  # or "pattern_folder" if that's your CSV column
    hover_data=["pattern", "holes_file"],
    template="plotly_white",                  # A clean white background
    color_discrete_sequence=px.colors.qualitative.Set3,  # Pleasant pastel palette
    symbol="pattern",                         # Each pattern gets a unique symbol
    symbol_sequence=[
        "circle", "diamond", "square", "x", 
        "triangle-up", "star", "triangle-down"
    ],
)

# 5) Update marker styling (size, outline, etc.)
fig.update_traces(
    marker=dict(
        size=12, 
        line=dict(width=1, color='DarkSlateGrey')
    ),
    selector=dict(mode='markers')
)

# 6) Improve axes and layout
fig.update_layout(
    title="Max Contraction vs. Relative Area Change at t=10 (Contract-Only)",
    xaxis_title="Max Contraction (t=10)",
    yaxis_title="Overall Relative Area Change (t=10)",
    font=dict(size=16),
    legend=dict(
        title="Pattern",
        orientation="h",
        yanchor="bottom",
        y=1.02,
        xanchor="right",
        x=1
    ),
    margin=dict(l=80, r=80, t=300, b=80),
    width=800,
    height=800
)

# Make axes lines and faint gridlines
fig.update_xaxes(
    showline=True, 
    linewidth=2, 
    linecolor='black',
    mirror=True,
    gridcolor="#eeeeee"
)
fig.update_yaxes(
    showline=True, 
    linewidth=2, 
    linecolor='black',
    mirror=True,
    gridcolor="#eeeeee"
)

# 7) Show interactive figure
fig.show()

# 8) Optionally save to HTML
fig.write_html("my_interactive_plot.html")
