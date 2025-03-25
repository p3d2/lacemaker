import os
import json
import matplotlib.font_manager as fm
import matplotlib.pyplot as plt
import numpy as np
import networkx as nx
from matplotlib.collections import LineCollection
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg
from matplotlib.patches import Patch
import tkinter as tk
from tkinter import Canvas, Scrollbar, IntVar, ttk, messagebox, filedialog
from PIL import Image, ImageTk
import subprocess
import tkinter.messagebox as msgbox
import datetime
import threading

def run_simulations():
    # First, save the data and get the new file path
    fp = save_data(False)
    if not fp:
        return  # Exit if saving failed or no file was loaded

    try:
        # Show the loading screen with GIF
        loading_window = show_gif_animation()  # Start the GIF animation

        # Extract filename and prepare commands
        filename = os.path.splitext(os.path.basename(fp))[0]
        active_checkboxes = [0] * 8
        print(checkbox_vars)
        for i, var in enumerate(checkbox_vars, start=1):
            if var.get() == 1 and i <= 8:
                active_checkboxes[i - 1] = 1 
        
        command1 = ["python", "lace_maker.py", fp,
                    "--dist_particles=0.5", "--units=1.0",
                    "--mass=1.0", "--threshold=0.0"]
        
        command2 = ["srun", "--output=log.out", "--mem=8G", "--ntasks=16", "lmp", 
                    "-var", "input_file", os.path.join("output", "lammps_data", f"{filename}_0.5_1.0_0.0_30.0_5.0.data"),
                    "-var", "fname", f"{filename}",
                    "-var", "modify_bond1", f"{active_checkboxes[0]}",
                    "-var", "modify_bond2", f"{active_checkboxes[1]}", 
                    "-var", "modify_bond3", f"{active_checkboxes[2]}", 
                    "-var", "modify_bond4", f"{active_checkboxes[3]}", 
                    "-var", "modify_bond5", f"{active_checkboxes[4]}", 
                    "-var", "modify_bond6", f"{active_checkboxes[5]}", 
                    "-var", "modify_bond7", f"{active_checkboxes[6]}", 
                    "-var", "modify_bond8", f"{active_checkboxes[7]}", 
                    "-in", "sim_run.lmp"]
        
        # Execute commands
        subprocess.run(command1, check=True)
        subprocess.run(command2, check=True)

        # Stop the GIF and close the loading window
        loading_window.destroy()
        
        # Inform the user of success
        messagebox.showinfo("Simulation Status", "The simulation finished successfully.")
    except subprocess.CalledProcessError as e:
        messagebox.showerror("Simulation Error", f"An error occurred during the simulations: {e}")
        if loading_window:
            loading_window.destroy()

def show_gif_animation():
    loading_window = tk.Toplevel()
    loading_window.title("Loading...")
    
    # Set window size
    gif_width = 400
    gif_height = 400
    loading_window.geometry(f"{gif_width}x{gif_height}")

    # Center the window
    screen_width = loading_window.winfo_screenwidth()
    screen_height = loading_window.winfo_screenheight()
    x_coordinate = int((screen_width / 2) - (gif_width / 2))
    y_coordinate = int((screen_height / 2) - (gif_height / 2))
    loading_window.geometry(f"+{x_coordinate}+{y_coordinate}")

    # Minimal window decorations
    loading_window.attributes('-type', 'splash')  # This might work on some Linux window managers

    # Load and display GIF
    frames = [tk.PhotoImage(file='ico/loading.gif', format=f'gif -index {i}') for i in range(6)]
    label = tk.Label(loading_window, image=frames[0])
    label.pack()

    # Update function for GIF
    def update_frame(index):
        frame = frames[index % len(frames)]
        label.configure(image=frame)
        loading_window.after(100, update_frame, (index+1) % len(frames))

    # Start GIF animation
    update_frame(0)

    return loading_window

def open_combobox():
    # Function to simulate opening the combobox
    file_combobox.focus() 
    file_combobox.event_generate('<Down>')

def get_files(directory, extension=".json"):
    """ List files in a directory with given extension. """
    return [file for file in os.listdir(directory) if file.endswith(extension)]

def load_data(filepath):
    with open(filepath, 'r') as file:
        data = json.load(file)
    return data

def on_file_select(event):
    global data, data_loaded, data_file
    selected_file = file_combobox.get()  # Get the selected file from combobox
    if selected_file:
        filepath = os.path.join(folder, selected_file)  # Construct the full path if needed
        data_file = filepath  # Update the global variable
        data_loaded = True
        with open(filepath, 'r') as file:
            data = json.load(file)
        clear_node_buttons(button_frame)  # Assuming this exists and is your button container
        canvas = Canvas(button_frame)
        scrollbar_button = Scrollbar(button_frame, orient="vertical", command=canvas.yview)
        canvas.configure(yscrollcommand=scrollbar_button.set)
        scrollbar_button.pack(side='right', fill='y')
        canvas.pack(side='left', fill='both', expand=True)
        interior = ttk.Frame(canvas, style = 'White.TFrame')
        canvas.create_window((0, 0), window=interior, anchor='nw')
        interior.bind("<Configure>", lambda event, canvas=canvas: onFrameConfigure(canvas))

        create_node_buttons(data, interior, canvas)  # Pass the interior and canvas
        create_path_checkboxes()
        draw_graph()

def clear_node_buttons(frame):
    """Remove all child widgets from a given frame safely."""
    global active_button
    if frame.winfo_exists():  # Check if the frame still exists
        for widget in frame.winfo_children():
            widget.destroy()
    active_button = None  # Reset the active button reference to prevent errors

def create_node_buttons(data, frame, canvas):
    """Create buttons for each node in the data inside the provided frame."""
    columns = 12  # Number of columns for button grid
    button_width = 1.5  # Set a fixed width for buttons
    button_height = 1.5  # Set a fixed height for buttons

    for i, node_id in enumerate(data.get('nodes', {}).keys()):
        btn = ttk.Button(frame, text=f"{node_id}", style="TButton", width=button_width)
        btn.grid(row=i // columns, column=i % columns, sticky='ew', ipady=button_height)
        btn.config(command=lambda nid=node_id, b=btn: node_button_click(nid, b))

    # Before configuring, check if the canvas and its parent are in a valid state
    if canvas.winfo_exists():
        frame.update_idletasks()  # Ensure all widgets inside are updated
        canvas.config(scrollregion=canvas.bbox("all"))  # Update the scrollregion of the canvas
    else:
        print("Canvas does not exist anymore. Skipping configuration.")

def onFrameConfigure(canvas):
    """Reset the scroll region to encompass the inner frame."""
    canvas.configure(scrollregion=canvas.bbox("all"))
    
def load_json_file():
    global data, data_loaded, data_file  # Include data_file in the global declaration
    filepath = filedialog.askopenfilename(filetypes=[("JSON files", "*.json")])
    if filepath:
        with open(filepath, 'r') as file:
            data = json.load(file)
        data_file = filepath  # Update the global data_file variable
        data_loaded = True
        if button_frame.winfo_exists():
            clear_node_buttons(button_frame)
            create_node_buttons(data, button_frame)
        draw_graph()
        
def save_data(trigger_message):
    if data_file:
        # Generate the current timestamp
        directory, filename = os.path.split(data_file)
        timestamp = datetime.datetime.now().strftime('%y%m%d%H%M')
        new_filename = f"{os.path.splitext(filename)[0]}_{timestamp}{os.path.splitext(filename)[1]}"
        new_filepath = os.path.join(os.path.join(directory, "custom"), new_filename)
        
        with open(new_filepath, 'w') as file:
            json.dump(data, file, indent=4, separators=(',',':'))

        if trigger_message:
            messagebox.showinfo("Save Data", "Data saved successfully to " + new_filepath)
        return new_filepath  # Return the new file path for further use
    else:
        messagebox.showerror("Error", "No file loaded to save. Please load a file first.")
        return None  # Return None if no file was saved

def setup_initial_state():
    global fig, ax, figure_canvas, fig2, ax2, figure_canvas2
    ax.clear()  # Clear any previous drawing
    ax.text(0.5, 0.5, 'No data loaded', va='center', ha='center', fontsize=12, color='gray')
    figure_canvas.draw()

def update_gui_controls():
    if data_loaded:
        # Enable buttons or entry fields here
        update_button['state'] = 'normal'
        save_button['state'] = 'normal'
    else:
        # Disable buttons or entry fields here
        update_button['state'] = 'disabled'
        save_button['state'] = 'disabled'

# Function to update the entry boxes and highlight the active button
def node_button_click(node_id, button=None):
    """Handle node button click events."""
    global active_button
    current_node_var.set(node_id)
    x_var.set(data['nodes'][node_id][0])
    y_var.set(data['nodes'][node_id][1])
    if active_button and active_button.winfo_exists():
        active_button.configure(style="TButton")  # Reset the style of the previously active button
    if button and button.winfo_exists():
        button.configure(style="Active.TButton")  # Highlight the new active button
        active_button = button

# Function to save changes to the selected node
def update_node():
    node_id = current_node_var.get()
    if node_id:
        data['nodes'][str(node_id)] = [x_var.get(), y_var.get()]
        draw_graph()

def to_int(value):
    # Remove last character if it's 'l' or 'r'
    val = str(value)
    if val[-1] in 'lr':
        return int(val[:-1])
    return int(val)

def create_path_checkboxes():
    global checkbox_vars
    checkbox_vars = []

    # Clear existing widgets
    for widget in shrink_frame.winfo_children():
        widget.destroy()

    # Add the 'Shrink' label
    shrink_label = ttk.Label(shrink_frame, text="Shrink", background='white')
    shrink_label.pack(side="left", padx=(0, 10))

    # Create checkboxes
    num_paths = len(data['unit_yarns'].items())
    for i in range(num_paths):
        var = tk.IntVar(value=0)
        chk = ttk.Checkbutton(shrink_frame, text=f"{i+1}", variable=var, style='White.TCheckbutton')
        chk.pack(side="left", padx=5)
        checkbox_vars.append(var)

def draw_graph():
    if not data_loaded:
        return setup_initial_state()  # Show placeholder if no data is loaded
    
    ax.clear()  # Clear previous drawings
    ax2.clear() 
    
    G = nx.Graph()
    node_colors = {}  
    label_dict = {}  # Dictionary for custom labels

    for node_id, pos in data['nodes'].items():
        rounded_pos = tuple(round(p, 2) for p in pos)  # Round to nearest hundredth
        G.add_node(str(node_id), pos=rounded_pos[:2])
        node_colors[str(node_id)] = 'white'
        label_dict[node_id] = str(node_id)

    pos_to_node_id = {tuple(G.nodes[node]['pos']): node for node in G.nodes()}
    pos_base = pos_to_node_id.copy()

    # Define colors for different paths for clarity
    path_colors = ['#ff6666', '#6666ff', '#ffff00', '#ff66ff', '#66ff33', '#ccffb3', '#b300ff', '#33ffff']

    # Establish paths and shifts
    paths = {}
    lace = []
    for index, path_info in enumerate(data['paths']):
        path = path_info['path']
        paths[index] = path

    # Plot unit yarns with specified path starts and shifts
    for yarn_id, yarn_data in data['unit_yarns'].items():
        path_id = yarn_data[0]
        node_start_index = yarn_data[1]
        path_start_index = yarn_data[2]
        path_z0 = yarn_data[3]
        yarn_path = [to_int(item) for item in paths[path_id]]
        
        # Unit repetitions
        unit_repetition = data['unit_repetion'][str(yarn_id)]
        rep1 = unit_repetition[0]
        vector1 = np.array([unit_repetition[1], unit_repetition[2]])  # Extract x and y components
        vector2 = np.array([unit_repetition[4], unit_repetition[5]])
        rep2 = unit_repetition[3]
        
        # Adjust path for starting point within pattern
        adjusted_path = yarn_path[path_start_index:-1] + yarn_path[:path_start_index]
        current_pos = np.array(G.nodes[str(node_start_index)]['pos'], dtype=float)  # Ensure position is float
        cumulative_shift = np.array([0.0, 0.0])  # Initialize cumulative shift as float
        
        path_z0 = yarn_data[3]
        path_zorder = []
        for i in range(len(adjusted_path)):
            path_zorder.append(path_z0)
            try:
                twists = data['nodes'][str(adjusted_path[i])][2]
                pass
            except:
                twists = 0
                pass
            
            if twists % 2 == 0: path_z0 = -path_z0

        path_save = []
        path_save.append(current_pos)
        for i in range(len(adjusted_path)): # excluding last point of the path since it is cyclic
            node_start = str(adjusted_path[i % len(adjusted_path)])
            node_end = str(adjusted_path[(i + 1) % len(adjusted_path)])
            
            shift_key = f"[{node_start}, {node_end}]"
            current_shift = np.array(data['paths'][path_id]['shifts'].get(shift_key, (0, 0)), dtype=float)
            cumulative_shift += current_shift
            end_pos = np.array(G.nodes[node_end]['pos'], dtype=float) + cumulative_shift
            rounded_pos = tuple(round(p, 2) for p in end_pos)
            
            if rounded_pos not in pos_to_node_id:
                new_node = str(end_pos)+'_'+str(rounded_pos[0])+'_'+str(rounded_pos[1])
                G.add_node(new_node, pos=tuple(rounded_pos))  # Ensure pos is tuple as your nodes expect
                node_colors[new_node] = 'white' 
                label_dict[new_node] = node_end
                pos_to_node_id[rounded_pos] = new_node
                
            G.add_edge(pos_to_node_id.get(tuple(current_pos)), pos_to_node_id.get(tuple(end_pos)), color=path_colors[int(yarn_id) % len(path_colors)], linewidth=2)
            
            # Change color of the trivial network
            if (tuple(current_pos) in pos_base) and path_zorder[i % len(adjusted_path)] > 0:    
                node_colors[node_start] = path_colors[int(yarn_id) % len(path_colors)]
                try:
                    if abs(data['nodes'][node_start][2]) > 0:
                        node_colors[node_start] = 'gray'
                        pass
                except:
                    pass
            
            path_save.append(end_pos)
                    
            # Update current point
            current_pos = end_pos

        for k2 in range(rep2):
            for k1 in range(rep1):
                replicated_path = [(x + k1*vector1[0] + k2*vector2[0], y + k1*vector1[1] + k2*vector2[1]) for x, y in path_save]
                # Assume each path uses the same color indexed by `yarn_id`
                lace.append((replicated_path, path_colors[int(yarn_id) % len(path_colors)]))  # Store paths with colors
            lace.append((None, None))  # Use None to separate paths; adjust handling in plotting

    # Draw the graph
    colors = [node_colors[node] for node in G.nodes()]
    edge_colors = [G[u][v]['color'] for u, v in G.edges()]
    
    pos = nx.get_node_attributes(G, 'pos')
    
    highlight_nodes = [node for node in G.nodes() if node_colors[node] != 'white']
    highlight_colors = [node_colors[node] for node in highlight_nodes]  # Ensure this list matches highlight_nodes
    
    nx.draw_networkx_nodes(G, pos, ax=ax2, nodelist=highlight_nodes,
                               node_color=highlight_colors, edgecolors='black', node_size=300, linewidths=1)
    nx.draw_networkx_nodes(G, pos, ax=ax2, node_color=colors, node_size=300, alpha = 0.5)
    nx.draw_networkx_edges(G, pos, ax=ax2, edge_color='black', width=12) 
    nx.draw_networkx_edges(G, pos, ax=ax2, edge_color=edge_colors, width=8)
    nx.draw_networkx_labels(G, pos, ax=ax2, labels=label_dict,
                            font_size=12, font_weight='bold',
                            horizontalalignment='center', verticalalignment='center')

    # Axis and grid setup from your example
    # Set grid spacing
    maj_grid_spacing = 4
    min_grid_spacing = 1  
    ax2.set_xticks(np.arange(np.round((min(x for x, _ in pos.values()) - 1)/maj_grid_spacing)*maj_grid_spacing, np.round((max(x for x, _ in pos.values()) + 1)/maj_grid_spacing)*maj_grid_spacing, maj_grid_spacing))
    ax2.set_yticks(np.arange(np.round((min(y for _, y in pos.values()) - 1)/maj_grid_spacing)*maj_grid_spacing, np.round((max(y for _, y in pos.values()) + 1)/maj_grid_spacing)*maj_grid_spacing, maj_grid_spacing))
    ax2.set_xticks(np.arange(min(x for x, _ in pos.values()) - 1, max(x for x, _ in pos.values()) + 1, min_grid_spacing), minor=True)
    ax2.set_yticks(np.arange(min(y for _, y in pos.values()) - 1, max(y for _, y in pos.values()) + 1, min_grid_spacing), minor=True)
    # Enable the grid
    ax2.grid(True, which='both', linestyle='-', linewidth=0.5, color='gray', zorder = -1)
    ax2.grid(True, which='minor', linestyle='--', color='lightgray', linewidth=0.5, zorder = -1)
    
    xmin, xmax = [min(x for x, _ in pos.values()), max(x for x, _ in pos.values())]
    ymin, ymax = [min(y for _, y in pos.values()), max(y for _, y in pos.values())]
    ax2.annotate('X', xy=(0.01, 1.01),  xycoords='axes fraction',
                 xytext=(0.01, 1.01),  textcoords='axes fraction',
                 ha='center', fontsize=12, color='black')
    
    ax2.annotate('', xy=(0.1, 1.02),  xycoords='axes fraction',
                 xytext=(0.02, 1.02),  textcoords='axes fraction',
            arrowprops=dict(facecolor='black', arrowstyle='->'),
            ha='center', fontsize=12, color='black')
    
    ax2.annotate('Y', xy=(-0.02, 0.9), xycoords='axes fraction',
            xytext=(-0.02, 0.98), textcoords='axes fraction',
            arrowprops=dict(facecolor='black', arrowstyle='->'),
            ha='center', fontsize=12, color='black')
    
    # Set limits and equal aspect ratio to ensure the grid lines are square if desired
    
    ax2.set_xlim(xmin - 2, xmax + 2)
    ax2.set_ylim(ymin - 2, ymax + 2)
    ax2.set_aspect('equal', adjustable='box')
    ax2.invert_yaxis()
    #ax2.margins(0.1)

    # Figure 2
    # Draw pattern
    # Filter out the paths and their colors from lace
    lines = [np.array(path) for path, color in lace if path is not None]
    colors = [color for path, color in lace if path is not None]
    
    # Create the LineCollection with the corresponding colors
    line_collection = LineCollection(lines, colors='black', alpha=1, linewidths=6, zorder = -1)
    ax.add_collection(line_collection)
    line_collection = LineCollection(lines, colors=colors, alpha=1.0, linewidths=4)
    ax.add_collection(line_collection)
    legend_handles = [Patch(color=path_colors[i % len(path_colors)], label=f'Path {i+1}') for i in range(len(data['unit_yarns'].items()))]
    ax.legend(handles=legend_handles, title="Paths", fontsize=16)
    
    ax.set_aspect('equal', adjustable='box')
    crop = data['roi_bounds']
    ax.set_xlim(crop['x_min'], crop['x_max'])
    ax.set_ylim(crop['y_min'], crop['y_max'])
    ax.invert_yaxis()
    ax.margins(0.1)

    figure_canvas.draw()  # Redraw the canvas with the updated graph
    figure_canvas2.draw()

#folder = r'\\data.triton.aalto.fi\work\silvap1\lacemaker\input\json_patterns_small'
folder = r'input/json_patterns'
file_paths = get_files(folder)
plt.rcParams['font.family'] = 'monospace'
plt.rcParams['font.monospace'] = ['Source Code Pro']

# Initialize the main window
root = tk.Tk()
root.title("Node Editor")
root.state('normal')  # Maximize the window
current_node_var = IntVar()
data = {}
data_file = None  # Initialize as None or with a default path if applicable

# Applying a theme
style = ttk.Style(root)
style.theme_use('classic')

# Style configuration for ttk buttons
style.configure("TButton", font=('Source Code Pro', 10), background="#ffffff", borderwidth=1)
style.configure("Active.TButton", background="#c0e4c0")
style.configure('White.TFrame', background='#ffffff')
style.configure("BoldLabel.TLabel", font=('Source Code Pro', 12, 'bold'), foreground="grey", background='#ffffff')
style.configure('TCombobox', arrowsize=0)  # Adjust padding as necessary
style.configure('White.TCheckbutton', background='white')

# Define left and right panels
left_panel = ttk.Frame(root, width=400, relief=tk.RIDGE, style = 'White.TFrame')
left_panel.pack(side="left", fill="y", expand=False)

right_panel = ttk.Frame(root, relief=tk.RIDGE)
right_panel.pack(side="right", fill="both", expand=True)

# Modify left panel to contain two sub-frames
controls_frame = ttk.Frame(left_panel, style = 'White.TFrame')
controls_frame.pack(side="top", fill="both", expand=True)

figure_frame = ttk.Frame(left_panel)
figure_frame.pack(side="bottom", fill="both", expand=True)

# Setup for scrollable canvas in the right panel
tk_canvas = Canvas(right_panel)
tk_canvas.pack(side="left", fill="both", expand=True)

scrollbar = Scrollbar(right_panel, orient="vertical", command=tk_canvas.yview)
scrollbar.pack(side="right", fill='y')
tk_canvas.config(yscrollcommand=scrollbar.set)

button_frame = ttk.Frame(left_panel, style = 'White.TFrame')
button_frame.pack(fill="both", expand=True)

# Embedding the matplotlib figure within the right panel's canvas
fig, ax = plt.subplots(figsize=(12, 12))
figure_canvas = FigureCanvasTkAgg(fig, master=tk_canvas)
figure_canvas.draw()
canvas_widget = figure_canvas.get_tk_widget()
canvas_widget.pack(side="top", fill="both", expand=True)

# Setup for the second matplotlib figure
fig2, ax2 = plt.subplots(figsize=(7.5, 7.5))
figure_canvas2 = FigureCanvasTkAgg(fig2, master=figure_frame)
figure_canvas2.draw()
canvas_widget2 = figure_canvas2.get_tk_widget()
canvas_widget2.pack(side="bottom", fill="both", expand=True)

# Configure the tk.Canvas to update its scrollregion automatically on resize
canvas_widget.bind("<Configure>", lambda event: tk_canvas.configure(scrollregion=tk_canvas.bbox("all")))

# Define controls in the left panel
x_var = tk.DoubleVar()
y_var = tk.DoubleVar()

# Frame for buttons within controls_frame
buttons_frame = ttk.Frame(controls_frame, style = 'White.TFrame')
buttons_frame.pack(fill='y', padx=10, expand=True)
buttons_frame.grid_propagate(False)

# Load JSON Button
#load_button = ttk.Button(buttons_frame, text="Load JSON", command=lambda: load_json_file())
#load_button.pack(side="left", padx=5, expand=True)
#icon_path = r'\\data.triton.aalto.fi\work\silvap1\lacemaker\ico\load.png'
icon_path = r'ico/load.png'
global load_ico
image = Image.open(icon_path)
load_ico = ImageTk.PhotoImage(image)

file_combobox = ttk.Combobox(buttons_frame, values=file_paths, 
                             font=('Source Code Pro', 12), background = 'White',
                             style='TCombobox')
file_combobox.pack(side="left", padx=(0,10), fill='both', expand=True)
file_combobox.bind("<<ComboboxSelected>>", on_file_select)

button = ttk.Button(buttons_frame, image=load_ico, text="Load", compound="left", command=open_combobox)
button.pack(side="left", padx=(0, 10))

# Keep a reference to the image to prevent garbage collection
button.image = load_ico

# Update Node Button
#icon_path = r'\\data.triton.aalto.fi\work\silvap1\lacemaker\ico\update.png'
icon_path = r'ico/update.png'
global update_ico
image = Image.open(icon_path)
update_ico = ImageTk.PhotoImage(image)

update_button = ttk.Button(buttons_frame, image = update_ico, text="Update", compound="left", command=update_node)
update_button.pack(side="left", padx=5, expand=True)
update_button.photo = update_ico# keep a reference to the image to avoid garbage collection
update_button.pack()

# Run Simulations Button
#icon_path_run = r'\\data.triton.aalto.fi\work\silvap1\lacemaker\ico\run.png'
icon_path_run = r'ico/run.png'
run_ico = ImageTk.PhotoImage(Image.open(icon_path_run))

run_button = ttk.Button(buttons_frame, image=run_ico, text="Run", compound="left", command=lambda: threading.Thread(target=run_simulations).start())
run_button.pack(side="left", padx=5)  # Ensure it's packed before the save button
run_button.photo = run_ico  # Keep a reference to avoid garbage collection

# Save All Changes Button
#icon_path = r'\\data.triton.aalto.fi\work\silvap1\lacemaker\ico\save.png'
icon_path = r'ico/save.png'
global save_ico
image = Image.open(icon_path)
save_ico = ImageTk.PhotoImage(image)

data_loaded = False
save_button = ttk.Button(buttons_frame, image = save_ico, text="Save", compound="left", command=lambda: save_data(True))
save_button.pack(side="left", padx=5, expand=True)
save_button.photo = save_ico# keep a reference to the image to avoid garbage collection
save_button.pack()

# Define global active_button for reference in node_button_click
active_button = None

# Frame for X and Y coordinates within controls_frame
coordinates_frame = ttk.Frame(controls_frame, style = 'White.TFrame')
coordinates_frame.pack(fill='x', padx=10, expand=True)

# X Coordinate Entry
x_label = ttk.Label(coordinates_frame, text="X:", style = 'BoldLabel.TLabel')
x_label.pack(side="left", padx=(5, 0))
x_entry = ttk.Entry(coordinates_frame, textvariable=x_var, font=('Source Code Pro', 14), width=10)
x_entry.pack(side="left", padx=(0, 40))

# Y Coordinate Entry
y_label = ttk.Label(coordinates_frame, text="Y:", style = 'BoldLabel.TLabel')
y_label.pack(side="left", padx=(5, 0))
y_entry = ttk.Entry(coordinates_frame, textvariable=y_var, font=('Source Code Pro', 14), width=10)
y_entry.pack(side="left")

shrink_frame = ttk.Frame(coordinates_frame, style='White.TFrame')
shrink_frame.pack(side="left", fill='x', padx=(5, 0))
    
setup_initial_state()    
draw_graph()
root.mainloop()