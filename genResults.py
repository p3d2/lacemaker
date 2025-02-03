import pandas as pd
import plotly.express as px

# Read your aggregated results
df = pd.read_csv("all_results.csv")

# Make an interactive scatter plot
# - x-axis: rel_max_contraction
# - y-axis: rel_max_area
# - color by pattern_folder (so all points from the same folder have the same color)
# - show pattern_folder & source_file on hover
fig = px.scatter(
    df,
    x="rel_max_contraction",
    y="rel_max_area",
    color="pattern_folder",
    hover_data=["pattern_folder", "source_file"],
    color_discrete_sequence=px.colors.qualitative.Light24  # Good for many categories
)

fig.update_layout(
    title="Relative Contraction vs. Relative Area",
    xaxis_title="rel_max_contraction",
    yaxis_title="rel_max_area"
)

# Display the interactive figure in a browser or notebook
fig.show()
fig.write_html("my_interactive_plot.html")