"""Validate CMC summary counts without inventing companion phase-space data.

Each PARAM file describes one cluster state. Repeated states in one history
are not independent clusters or distinct black-hole binaries. The selected
near-reference-age record gives each available history at most one row; the
remaining summaries are descriptive only, with no interpolation or imputation.
"""
from collections import defaultdict
from decimal import Decimal, InvalidOperation
import math
from pathlib import Path
from statistics import median


STELLAR_TYPES = tuple(range(16))
TYPE_COUNT_KEYS = tuple(f'N_{kind}' for kind in STELLAR_TYPES)
BINARY_COUNT_KEYS = tuple(
    f'Nbin_{primary}_{secondary}'
    for primary in STELLAR_TYPES for secondary in range(primary + 1))
REQUIRED_COUNT_KEYS = ('Nsys', 'Nobj', 'Nbin') + TYPE_COUNT_KEYS + BINARY_COUNT_KEYS


def _count(token, key):
    """Parse integers exactly, including the upstream exponential notation."""
    try:
        value = Decimal(token)
    except InvalidOperation as exc:
        raise ValueError(f'Invalid count {key}: {token}') from exc
    if not value.is_finite() or value < 0 or value != value.to_integral_value():
        raise ValueError(f'Count {key} must be a finite nonnegative integer')
    return int(value)


def _positive_number(token, key):
    try:
        value = float(token)
    except ValueError as exc:
        raise ValueError(f'Invalid quantity {key}: {token}') from exc
    if not math.isfinite(value) or value <= 0:
        raise ValueError(f'Quantity {key} must be finite and positive')
    return value


def parse_param(text, *, model, path):
    """Return one validated count row from a CMC PARAM text record.

All stellar-type and unordered binary-type counts are required, including
explicit zeros. Unknown summary quantities are retained only while checking
syntax and duplicate keys; their values are not used. Undefined dispersions
may legitimately be NaN, but no count, host mass, or age may be nonfinite.
No source file or shell expression is executed.
"""
    if not isinstance(model, str) or not model or not isinstance(path, str) or not path:
        raise ValueError('Nonempty model and source path are required')
    tokens = {}
    for line_number, line in enumerate(text.splitlines(), 1):
        fields = line.split('#', 1)[0].split()
        if not fields:
            continue
        if len(fields) != 2:
            raise ValueError(f'Expected one key/value pair at line {line_number}')
        key, token = fields
        if key in tokens:
            raise ValueError(f'Duplicate PARAM key: {key}')
        tokens[key] = token
    required = REQUIRED_COUNT_KEYS + ('Mtot', 'age')
    missing = sorted(set(required).difference(tokens))
    if missing:
        raise ValueError('Missing PARAM quantities: ' + ', '.join(missing))
    counts = {}
    for key, token in tokens.items():
        if key.startswith('N'):
            counts[key] = _count(token, key)
        if key.startswith('N_') and key not in TYPE_COUNT_KEYS:
            raise ValueError(f'Unsupported stellar-type count: {key}')
        if key.startswith('Nbin_') and key not in BINARY_COUNT_KEYS:
            raise ValueError(f'Unsupported or noncanonical binary-type count: {key}')
    age = _positive_number(tokens['age'], 'age')
    mass = _positive_number(tokens['Mtot'], 'Mtot')
    if sum(counts[key] for key in TYPE_COUNT_KEYS) != counts['Nobj']:
        raise ValueError('Stellar-type counts do not sum to Nobj')
    if sum(counts[key] for key in BINARY_COUNT_KEYS) != counts['Nbin']:
        raise ValueError('Binary-type counts do not sum to Nbin')
    if counts['Nsys'] + counts['Nbin'] != counts['Nobj']:
        raise ValueError('Count identity Nsys + Nbin = Nobj failed')
    if counts['Nbin'] > counts['Nsys']:
        raise ValueError('Binary count exceeds system count')
    bbh = counts['Nbin_14_14']
    bh_other = sum(counts[f'Nbin_{max(14, kind)}_{min(14, kind)}']
                   for kind in STELLAR_TYPES if kind != 14)
    single_bh = counts['N_14'] - 2 * bbh - bh_other
    if single_bh < 0:
        raise ValueError('BH binary membership exceeds N_14')
    if single_bh > counts['Nsys'] - counts['Nbin']:
        raise ValueError('Single BH count exceeds all single systems')
    return {
        'model': model,
        'path': path,
        'age_gyr': age,
        'system_count': counts['Nsys'],
        'object_count': counts['Nobj'],
        'binary_count': counts['Nbin'],
        'bh_count': counts['N_14'],
        'bbh_count': bbh,
        'bh_nonbbh_binary_count': bh_other,
        'single_bh_count': single_bh,
        'host_mass_msun': mass,
    }


def load_param(local_path, *, model, source_path):
    """Read a cached file and preserve its upstream path in the result."""
    return parse_param(Path(local_path).read_text(), model=model, path=source_path)


def summarize_models(rows, *, reference_age_gyr=12.0, tolerance_gyr=0.25,
                     minimum_bbh_count=30):
    """Describe each history and select at most one record near a common age.

    Selection minimizes (absolute age distance, age, source path) inside the
    inclusive tolerance. A missing eligible record is None, never zero. The
    binary-count threshold is an operational screen, not evidence that a
    continuous distribution can be calibrated or that snapshots are independent.
    """
    if not math.isfinite(reference_age_gyr) or reference_age_gyr <= 0:
        raise ValueError('Reference age must be finite and positive')
    if not math.isfinite(tolerance_gyr) or tolerance_gyr < 0:
        raise ValueError('Age tolerance must be finite and nonnegative')
    if (isinstance(minimum_bbh_count, bool) or not isinstance(minimum_bbh_count, int)
            or minimum_bbh_count < 1):
        raise ValueError('Minimum BBH count must be a positive integer')
    histories = defaultdict(list)
    paths = set()
    for row in rows:
        if row['path'] in paths:
            raise ValueError(f'Duplicate source path: {row["path"]}')
        paths.add(row['path'])
        if not math.isfinite(row['age_gyr']) or row['age_gyr'] <= 0:
            raise ValueError('Snapshot age must be finite and positive')
        if (isinstance(row['bbh_count'], bool) or not isinstance(row['bbh_count'], int)
                or row['bbh_count'] < 0):
            raise ValueError('Snapshot BBH count must be a nonnegative integer')
        histories[row['model']].append(row)
    summaries = []
    for model, history in sorted(histories.items()):
        eligible = [row for row in history
                    if abs(row['age_gyr'] - reference_age_gyr) <= tolerance_gyr]
        selected = min(eligible, key=lambda row: (
            abs(row['age_gyr'] - reference_age_gyr), row['age_gyr'], row['path']
        )) if eligible else None
        counts = [row['bbh_count'] for row in history]
        ages = [row['age_gyr'] for row in history]
        summaries.append({
            'model': model,
            'snapshot_count': len(history),
            'minimum_age_gyr': min(ages),
            'maximum_age_gyr': max(ages),
            'minimum_bbh_count': min(counts),
            'median_bbh_count': median(counts),
            'maximum_bbh_count': max(counts),
            'saved_states_meeting_count_screen': sum(
                count >= minimum_bbh_count for count in counts),
            'selected_path': selected['path'] if selected else None,
            'selected_age_gyr': selected['age_gyr'] if selected else None,
            'selected_age_offset_gyr': (
                selected['age_gyr'] - reference_age_gyr if selected else None),
            'selected_bbh_count': selected['bbh_count'] if selected else None,
            'selected_meets_count_screen': (
                selected['bbh_count'] >= minimum_bbh_count if selected else None),
        })
    return summaries
