import os
import csv

BASE_DIR = "benchmarks/results/linear_radix_sort/4814006"
OUTPUT_FILE = os.path.join(BASE_DIR, "benchmark.csv")

def parse_single_log(filepath):
    """Return (header, rows) from the given benchmark.log file."""
    with open(filepath, "r") as f:
        lines = f.readlines()

    # Find separator "#####"
    try:
        idx = next(i for i, line in enumerate(lines) if line.strip() == "#####")
    except StopIteration:
        raise ValueError(f"Separator ##### not found in {filepath}")

    # CSV part starts immediately after separator
    csv_lines = [line.strip() for line in lines[idx + 1:] if line.strip()]

    if not csv_lines:
        raise ValueError(f"No CSV data found in {filepath}")

    header = csv_lines[0].split(";")
    rows = [row.split(";") for row in csv_lines[1:]]

    return header, rows


def main():
    # Collect all subfolders that are integers 0..16
    folders = sorted(
        [name for name in os.listdir(BASE_DIR) if name.isdigit()],
        key=lambda x: int(x)
    )

    merged_header = None
    merged_rows = []

    for folder in folders:
        log_path = os.path.join(BASE_DIR, folder, "benchmark.log")

        if not os.path.isfile(log_path):
            print(f"WARNING: file not found: {log_path}")
            continue

        header, rows = parse_single_log(log_path)

        if merged_header is None:
            merged_header = header
        else:
            if header != merged_header:
                raise ValueError(f"Header mismatch in folder {folder}")

        merged_rows.extend(rows)

    # Write merged CSV
    with open(OUTPUT_FILE, "w", newline="") as f:
        writer = csv.writer(f, delimiter=";")
        writer.writerow(merged_header)
        writer.writerows(merged_rows)

    print(f"Merged file written to {OUTPUT_FILE}")


if __name__ == "__main__":
    main()