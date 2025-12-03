from __future__ import print_function
import sys
import os

sys.path.append(r"C:\Program Files\Vicon\DataStream SDK\Win64\Python\vicon_dssdk")

from vicon_dssdk import ViconDataStream
import argparse
import sys
import csv
import time
import signal

# -------------------
# CLI Arguments
# -------------------
parser = argparse.ArgumentParser(description=__doc__)
parser.add_argument('host', nargs='?', help="Host name, in the format of server:port", default="localhost:801")
parser.add_argument('--outfile', help="CSV file to save results", default="vicon_stream.csv")
parser.add_argument('--duration', type=int, help="Duration to stream (sec)", default=5000)
args = parser.parse_args()

# This seems to work: python record.py 192.168.1.2:801 --duration 50


# -------------------
# Connect to DataStream
# -------------------
client = ViconDataStream.Client()

try:
    client.Connect(args.host)
    client.EnableSegmentData()
    client.SetBufferSize(1)
    client.SetStreamMode(ViconDataStream.Client.StreamMode.EClientPull)
    print("Connected to Vicon DataStream at", args.host)

    # -------------------
    # Prepare storage
    # -------------------
    pose_buffer = []
    start_time = time.time()

    # Graceful exit handler
    stop_stream = False
    def signal_handler(sig, frame):
        global stop_stream
        stop_stream = True
        print("\nStopping stream...")
    signal.signal(signal.SIGINT, signal_handler)

    # -------------------
    # Main streaming loop
    # -------------------
    print("Streaming poses for", args.duration, "seconds... (Ctrl+C to stop early)")
    while not stop_stream and (time.time() - start_time) < args.duration:
        if not client.GetFrame():
            print("No frame :(")
            continue

        frame_number = client.GetFrameNumber()
        timestamp = time.time() * 1e3 # Time in ms

        # For each subject/segment, grab pose
        for subject in client.GetSubjectNames():
            for segment in client.GetSegmentNames(subject):
                trans = client.GetSegmentGlobalTranslation(subject, segment)
                rot = client.GetSegmentGlobalRotationQuaternion(subject, segment)

                # Store: [time, frame, subject, segment, x, y, z, qx, qy, qz, qw]
                pose_buffer.append([
                    timestamp, frame_number, subject, segment,
                    trans[0][0], trans[0][1], trans[0][2],
                    rot[0][0], rot[0][1], rot[0][2], rot[0][3]
                ])

        # time.sleep(0.01)  # ~100 Hz max

    # -------------------
    # Write buffer to CSV
    # -------------------
    with open(args.outfile+".csv", "w", newline="") as f:
        writer = csv.writer(f)
        writer.writerow(["timestamp", "frame", "subject", "segment",
                         "x", "y", "z", "qx", "qy", "qz", "qw"])
        writer.writerows(pose_buffer)

    print("✅ Saved", len(pose_buffer), "poses to", args.outfile+".csv")

except ViconDataStream.DataStreamException as e:
    print("Handled data stream error:", e)
