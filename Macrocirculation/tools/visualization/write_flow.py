import os
import numpy as np
import json
import argparse
import csv
from matplotlib import pyplot as plt

# Parse command-line arguments
parser = argparse.ArgumentParser(description='Extract, save, and optionally plot q data.')
parser.add_argument('--vessels', type=int, nargs='+', help='A list of IDs of the vessels to process.', default=[36])
parser.add_argument('--position', type=float, help='Position along the vessel to extract data.', default=0.5)
parser.add_argument('--t-period', type=float, default=1000.0, help='Heartbeat period duration.')
parser.add_argument('--t-systole', type=float, default=300, help='Systole duration in the heartbeat.')
parser.add_argument('--t-end', type=float, default=10000, help='End time for extraction.')
parser.add_argument('--filepath', type=str, required=True, help='Path to the metadata file.')
parser.add_argument('--output', type=str, default='output.csv', help='Output CSV file path.')
parser.add_argument('--positive-q', help='Reorients the flow q to be always positive.', action='store_true')
parser.add_argument('--use-shifted-vessel-numbers', help='Shifts the vessel index by one.', action='store_true')
parser.add_argument('--plot', help='Plot the q data for the selected vessels.', action='store_true')

args = parser.parse_args()

directory = os.path.dirname(args.filepath)

# Adjust vessel IDs if shifted numbers are used
if args.use_shifted_vessel_numbers:
    args.vessels = [vid - 1 for vid in args.vessels]

# Calculate t-start for the last full heartbeat cycle - in [ms]
t_start = (args.t_end - args.t_period) / 1000

# Load metadata
with open(args.filepath) as f:
    meta = json.loads(f.read())

def find_vessel(vessel_id):
    """Find the vessel info by ID."""
    for vessel in meta['vessels']:
        if vessel['edge_id'] == vessel_id:
            return vessel
    raise ValueError(f'Vessel with ID {vessel_id} not found.')

# Prepare data and plot if needed
output_data = []
if args.plot:
    fig, ax = plt.subplots(len(args.vessels), 1, figsize=(8, 4 * len(args.vessels)), sharex=True)
    if len(args.vessels) == 1:
        ax = [ax]  # Ensure ax is iterable for a single vessel

for idx, vessel_id in enumerate(args.vessels):
    vessel_info = find_vessel(vessel_id)

    # Load q data
    q_path = os.path.join(directory, vessel_info['filepaths']['q'])
    print(f'Loading q data from {q_path}')
    q = np.loadtxt(q_path, delimiter=',')
    q = q[:, int((q.shape[1] - 1) * args.position)]  # Extract data at the specified position

    # Adjust q values if positive-q option is enabled
    if args.positive_q and q.mean() < 0:
        q *= -1

    # Load time data
    time_path = os.path.join(directory, meta['filepath_time'])
    print(f'Loading time data from {time_path}')
    t = np.loadtxt(time_path, delimiter=',')

    # Filter data within the calculated time range
    start_index = np.sum(t < t_start)
    end_index = np.sum(t < args.t_end)
    t = t[start_index:end_index]
    q = q[start_index:end_index]

    # Normalize q by the maximum q value
    q_min = q.min() if q.size > 0 else 1  # Avoid division by zero
    q_max = q.max() if q.size > 0 else 1  # Avoid division by zero
    q_normalized = q / q_min

    print(f'Minimum flow over a period {q_min}')
    print(f'Maximum flow over a period {q_max}')

    # Shift time values to start from zero
    t_shifted = t - t.min()

    # Store data for this vessel
    for time_val, q_val, q_norm in zip(t_shifted, q, q_normalized):
        output_data.append([vessel_id, time_val, q_val, q_norm])

    # Plot if required
    if args.plot:
        ax[idx].plot(t_shifted, q, label=f'Vessel {vessel_id}')
        ax[idx].set_title(f'Vessel {vessel_id} Flow (q)')
        ax[idx].set_ylabel('Flow q [cm³/s]')
        ax[idx].grid(True)
        ax[idx].legend()

# Save extracted data to a CSV file
output_file = args.output
print(f'Saving extracted data to {output_file}')
with open(output_file, mode='w', newline='') as csvfile:
    writer = csv.writer(csvfile)
    writer.writerow(['Vessel ID', 'Time (s)', 'Flow q (cm³/s)', 'Normalized Flow q'])  # Header
    writer.writerows(output_data)

print('Data extraction and saving completed.')

# Show the plot
if args.plot:
    ax[-1].set_xlabel('Time [s]')
    plt.tight_layout()
    plt.show()
