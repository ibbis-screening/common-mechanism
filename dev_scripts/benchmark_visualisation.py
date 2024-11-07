import plotly.graph_objects as go
from plotly.colors import qualitative

def time_str_to_seconds(time_str):
    # Convert "HH:MM:SS.SSS" to seconds as a float
    h, m, s = time_str.split(':')
    return int(h) * 3600 + int(m) * 60 + float(s)

def visualize_data(filename):
    # Read data from file
    entries = []
    with open(filename, 'r') as file:
        for line in file:
            parts = line.strip().split()
            if len(parts) < 5:
                continue  # Skip malformed lines

            y_pos = int(parts[0])                  # Y-axis location
            label = parts[1]                       # Task name
            duration = parts[2]                    # Duration string (for hover info)
            start_time_str = parts[3]              # Start time as string
            end_time_str = parts[4]                # End time as string

            start_time = time_str_to_seconds(start_time_str)
            end_time = time_str_to_seconds(end_time_str)
            width = end_time - start_time          # Rectangle width in seconds

            entries.append({
                "y_pos": y_pos,
                "label": label,
                "duration": duration,
                "start_time": start_time,
                "end_time": end_time
            })

    print(entries)

    # Create a Plotly figure
    fig = go.Figure()

    # Get unique labels and assign a color to each
    unique_labels = list(set(entry['label'] for entry in entries))
    colors = {label: qualitative.Plotly[(i-3) % len(qualitative.Plotly)] for i, label in enumerate(unique_labels)}

    for entry in entries:
        # Add a rectangle for each entry
        label_color = colors[entry['label']]

        # Add a rectangle shape with the designated color
        fig.add_shape(
            type="rect",
            x0=entry['start_time'], x1=entry['end_time'],  # x-axis start and end
            y0=entry['y_pos'] - 0.45, y1=entry['y_pos'] + 0.45,  # y-axis position with some height
            line=dict(width=0),  # No border around the rectangle
            fillcolor=label_color,  # Unique color for each task label
            opacity=0.5  # Transparency
        )

        fig.update()

        # Add a hover label for each rectangle as text
        fig.add_trace(go.Scatter(
            x=[(entry['start_time'] + entry['end_time']) / 2],  # Position text at the center of the rectangle
            y=[entry['y_pos']],
            text=f"{entry['label']}<br>{entry['duration']})",
            mode="text",
            showlegend=False,
            zorder=-100,
        ))



    # Configure plot layout
    fig.update_layout(
        title="Benchmark Visualization",
        xaxis=dict(
            title="Time (seconds)",
            showgrid=True,
            zeroline=False
        ),
        yaxis=dict(
            title="Call Stack",
            tickmode="linear",
            showgrid=True,
            zeroline=False
        )
    )

    # Save the plot as an HTML file
    fig.write_html("timeline_visualization.html")
    print("Plot saved as timeline_visualization.html")

# Usage:
visualize_data("benchmark.txt")