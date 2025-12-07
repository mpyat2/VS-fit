import re
import io
import pathlib
import pandas as pd

## Data file spec (no header row, tab-separated)
NAMES = ['Time', 'Mag']
DTYPE = {'Time': 'float64', 'Mag': 'float64'}

def load_data(fileName):
    out_lines = []

    with open(fileName, 'r', encoding='utf-8') as f:
        for line in f:

            # strip leading/trailing whitespace
            line = line.strip()
            
            # skip blanks and comments
            if not line or line.startswith('#'):
                continue

            if '\t' in line:
                # TAB mode: split on each TAB
                parts = line.split('\t')
            else:
                # SPACE mode: collapse multiple spaces → split
                parts = re.split(r'\s+', line)

            out_lines.append(','.join(parts))  # normalize to CSV with commas

    # join preprocessed lines to a single text buffer
    text = "\n".join(out_lines)

    # read via pandas
    return pd.read_csv(
        io.StringIO(text),
        names=NAMES,
        dtype={'Time': 'float64', 'Mag': 'float64'},
        header=None
    )


def save_result(fileName, data):
    filePath = pathlib.Path(fileName)
    if filePath.suffix == "":
        fileName += ".tsv"
    with open(fileName, 'w', newline='') as f:
        # write header line with '#'
        f.write('#' + '\t'.join(data.columns) + '\n')
        data.to_csv(f, index=False, header=False, sep='\t')


def create_frame(time, mag):
    return pd.DataFrame({
        "Time": time,
        "Mag": mag
    })

