"""Import Bakta JSON and run the optional local Bakta executable."""
import json
import os
from pathlib import Path
import shutil
import subprocess
import tempfile
import pandas as pd
from .analysis import parse_fasta

COLUMNS = ['sequence', 'type', 'start', 'end', 'strand', 'locus', 'gene', 'product', 'wraps_origin']


def parse_bakta(payload, expected=None):
    try:
        data = json.loads(payload)
    except (ValueError, TypeError) as exc:
        raise ValueError("Invalid JSON file.") from exc
    if not isinstance(data, dict) or not isinstance(data.get('sequences'), list) or not isinstance(data.get('features'), list):
        raise ValueError("Expected Bakta JSON with sequences and features arrays.")
    try:
        fasta = '\n'.join(f">{s['id']}\n{s['nt']}" for s in data['sequences'])
        records = parse_fasta(fasta)
    except (KeyError, TypeError) as exc:
        raise ValueError("Bakta JSON sequences must contain id and nt fields.") from exc
    if expected is not None and records != expected:
        raise ValueError("Annotation sequences do not exactly match the uploaded FASTA (IDs and bases). Use matching files or load JSON alone.")
    rows = []
    for feature in data['features']:
        try:
            sid = feature.get('sequence', feature.get('contig'))
            start, end = feature['start'], feature['stop']
            if sid not in records or type(start) is not int or type(end) is not int:
                raise ValueError("Invalid feature sequence or coordinates.")
            if not 1 <= start <= len(records[sid]) or not 1 <= end <= len(records[sid]):
                raise ValueError("Feature coordinates fall outside the sequence.")
            if start > end and not feature.get('edge', False):
                raise ValueError("Origin-spanning feature is missing Bakta's edge flag.")
            rows.append(dict(sequence=sid, type=str(feature.get('type', 'unknown')),
                             start=start, end=end, strand=str(feature.get('strand', '?')),
                             locus=str(feature.get('locus', '')), gene=str(feature.get('gene', '')),
                             product=str(feature.get('product', '')), wraps_origin=start > end))
        except (KeyError, TypeError, AttributeError) as exc:
            raise ValueError("Malformed Bakta feature record.") from exc
    provenance = {key: data[key] for key in ['version', 'run', 'genome'] if key in data}
    return records, pd.DataFrame(rows, columns=COLUMNS), provenance


def bakta_available():
    return shutil.which('bakta') is not None


def run_bakta(records, database, threads=2, timeout=3600):
    """Run without a shell; isolate all files and remove them after completion."""
    executable = shutil.which('bakta')
    if not executable:
        raise ValueError("Bakta is not installed or is not on PATH. See README setup instructions.")
    db = Path(database).expanduser().resolve()
    if not (db / 'version.json').is_file():
        raise ValueError("Choose an installed Bakta database directory containing version.json.")
    if not 1 <= threads <= 32:
        raise ValueError("Threads must be between 1 and 32.")
    with tempfile.TemporaryDirectory(prefix='baklens-') as work:
        root = Path(work)
        fasta = root / 'input.fna'
        fasta.write_text(''.join(f'>{sid}\n{seq}\n' for sid, seq in records.items()))
        command = [executable, '--db', str(db), '--output', str(root / 'result'),
                   '--prefix', 'baklens', '--keep-contig-headers', '--skip-plot',
                   '--threads', str(threads), str(fasta)]
        try:
            # New process group permits cleanup of annotation subprocesses on timeout.
            with (root / 'run.log').open('w+') as log:
                process = subprocess.Popen(command, stdout=log, stderr=subprocess.STDOUT, start_new_session=True)
                try:
                    process.wait(timeout=timeout)
                except subprocess.TimeoutExpired:
                    import signal
                    os.killpg(process.pid, signal.SIGTERM)
                    try:
                        process.wait(timeout=10)
                    except subprocess.TimeoutExpired:
                        os.killpg(process.pid, signal.SIGKILL)
                        process.wait()
                    raise ValueError("Bakta exceeded the one-hour limit. Run it separately and import its JSON.")
                log.seek(0)
                logs = log.read()[-100_000:]
            if process.returncode:
                raise ValueError(f"Bakta exited with code {process.returncode}.\n{logs[-8000:]}")
            result = (root / 'result' / 'baklens.json').read_bytes()
            parse_bakta(result, expected=records)
            return result, logs
        except OSError as exc:
            raise ValueError(f"Could not execute Bakta: {exc}") from exc
