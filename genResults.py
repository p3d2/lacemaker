import pandas as pd
import numpy as np
import plotly.express as px

# 1) Read the aggregated results from the CSV file
df = pd.read_csv("all_results.csv")

# 2) Create a new column 'recoverability' as the difference between t10 and t20 area change
df["rel_max_contraction"] = -df["rel_max_contraction"]
epsilon = 0.01
# df["recoverability"] = np.exp(-df["area_change_t20"]/(df["area_change_t10"] + epsilon))
df["recoverability"] = df["area_change_t10"]-df["area_change_t20"]

# 3) Create an interactive 3D scatter plot
fig = px.scatter_3d(
	df,
	x="rel_max_contraction",		# x-axis: relative contraction at t=10
	y="area_change_t10",			# y-axis: area change at t=10
	z="recoverability",				# z-axis: recoverability = t10 - t20
	color="exp",					# Group by the experiment identifier
	symbol="exp",					# Use the same grouping for symbols
	hover_data=["exp", "contract"],	# Display additional info on hover (adjust as needed)
	template="plotly_white",
	color_discrete_sequence=px.colors.qualitative.Set3,
	symbol_sequence=[
		"circle", "circle-open", "cross", "diamond",
		"diamond-open", "square", "square-open", "x"
	],
)

# 4) Update marker styling (you might want to adjust size for 3D)
fig.update_traces(
	marker=dict(
		size=8,
		line=dict(width=1, color='DarkSlateGrey')
	)
)

# 5) Update axes and layout for the 3D scene
fig.update_layout(
	title="3D Scatter Plot: Relative Contraction, Area Change (t=10), and Recoverability",
	scene=dict(
		xaxis=dict(title="Relative Contraction (t=10)", range=[-0.05, 0.25]),
		yaxis=dict(title="Area Change (t=10)", range=[-0.05, 0.75]),
		zaxis=dict(title="Recoverability (1 - t20)", range=[-0.05, 1.0]),
		aspectmode='manual',
		aspectratio=dict(x=1, y=1, z=1) 
	),
	font=dict(size=16),
	margin=dict(l=100, r=100, t=100, b=100),
	width=1500,
	height=1200
)

# 6) Show interactive figure and optionally save to HTML
fig.write_html("my_interactive_3d_plot.html")
