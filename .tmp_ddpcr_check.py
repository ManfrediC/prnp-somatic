import os
import csv
import re
import xml.etree.ElementTree as ET
import zipfile
from collections import defaultdict

root = r'C:/Projects/prnp-somatic'

# ---------- raw/ddpcr CSV ----------
raw_dir = os.path.join(root, 'raw', 'ddpcr')
status_counts = defaultdict(int)
frontal_samples = defaultdict(list)
frontal_rows = []
missing_or_zero = []

fr_pattern = re.compile(r'(?i)(?:^|[\s_\-])(fr\b|frontal)(?:\s|_|\-|$)')

def is_frontal(sample: str) -> bool:
    s = (sample or '').lower().strip()
    return bool(fr_pattern.search(s) or s.endswith('_fr') or s.endswith('-fr') or 'frontal' in s)

for fn in sorted(os.listdir(raw_dir)):
    if not fn.lower().endswith('.csv'):
        continue
    with open(os.path.join(raw_dir, fn), newline='') as f:
        rows = list(csv.DictReader(f))
    for r in rows:
        sample = (r.get('Sample') or '').strip()
        if not sample:
            continue
        status = (r.get('Status') or '').strip()
        status_counts[status] += 1
        target = (r.get('Target') or '').strip()
        accepted_txt = (r.get('Accepted Droplets') or '').replace(',', '').strip()
        try:
            accepted = float(accepted_txt) if accepted_txt != '' else None
        except Exception:
            accepted = None

        if is_frontal(sample):
            rec = {'file': fn, 'sample': sample, 'target': target, 'status': status, 'accepted': accepted}
            frontal_samples[sample].append(rec)
            frontal_rows.append(rec)
            if accepted is None or accepted <= 0:
                missing_or_zero.append((sample, fn, target, accepted, status))

frontal_zero_or_na = []
for sample, rows in frontal_samples.items():
    if not any((r['accepted'] is not None and r['accepted'] > 0) for r in rows):
        frontal_zero_or_na.append((sample, rows))

# ---------- SNV_data_final.xlsx ----------
excel_path = os.path.join(root, 'results', 'ddPCR', 'SNV_data_final.xlsx')
frontal_xlsx_total = 0
frontal_xlsx_zero = []
frontal_xlsx_by_participant = {}

ns = '{http://schemas.openxmlformats.org/spreadsheetml/2006/main}'
rel = '{http://schemas.openxmlformats.org/package/2006/relationships}'

with zipfile.ZipFile(excel_path) as zf:
    # shared strings
    shared = []
    if 'xl/sharedStrings.xml' in zf.namelist():
        s_root = ET.fromstring(zf.read('xl/sharedStrings.xml'))
        for si in s_root.findall(f'{ns}si'):
            txt = ''.join((t.text or '') for t in si.findall(f'.//{ns}t'))
            shared.append(txt)

    wb_root = ET.fromstring(zf.read('xl/workbook.xml'))
    sheet = wb_root.find(f'{ns}sheets/{ns}sheet')
    sheet_rid = sheet.attrib[f'{rel}id']

    rel_root = ET.fromstring(zf.read('xl/_rels/workbook.xml.rels'))
    target = None
    for rel_node in rel_root.findall(f'{rel}Relationship'):
        if rel_node.attrib.get('Id') == sheet_rid:
            target = rel_node.attrib.get('Target')
            break
    if target is None:
        raise RuntimeError('Could not resolve worksheet target')

    ws_path = f'xl/{target}'
    ws_root = ET.fromstring(zf.read(ws_path))


def cell_val(cell):
    t = cell.attrib.get('t')
    v = cell.find(f'{ns}v')
    if v is None:
        return ''
    text = v.text or ''
    if t == 's':
        try:
            return shared[int(text)]
        except Exception:
            return ''
    return text

rows_nodes = ws_root.find(f'{ns}sheetData').findall(f'{ns}row')
if not rows_nodes:
    raise RuntimeError('No rows found in worksheet')

# header row
header_nodes = rows_nodes[0].findall(f'{ns}c')
header = []
for c in header_nodes:
    addr = c.attrib.get('r', '')
    col_letters = re.match(r'^([A-Z]+)', addr)
    if not col_letters:
        continue
    idx = 0
    for ch in col_letters.group(1):
        idx = idx * 26 + (ord(ch) - 64)
    col_idx = idx - 1
    value = cell_val(c)
    if col_idx >= len(header):
        header.extend([''] * (col_idx - len(header) + 1))
    header[col_idx] = value

header = [h.strip() for h in header]

def col(name):
    for i, h in enumerate(header):
        if h.lower() == name.lower():
            return i
    return None

idx_participant = col('participant')
idx_region = col('brain_region')
idx_total = col('n_total_droplets')
idx_mut = col('n_mut_droplets')
idx_lod = col('detected_above_LoD')
if idx_participant is None or idx_region is None:
    raise RuntimeError('Required columns participant/brain_region not found')

for r in rows_nodes[1:]:
    vals = ['' for _ in range(len(header))]
    for c in r.findall(f'{ns}c'):
        addr = c.attrib.get('r', '')
        col_letters = re.match(r'^([A-Z]+)', addr)
        if not col_letters:
            continue
        idx = 0
        for ch in col_letters.group(1):
            idx = idx * 26 + (ord(ch) - 64)
        col_idx = idx - 1
        if col_idx >= len(vals):
            vals.extend([''] * (col_idx - len(vals) + 1))
        vals[col_idx] = cell_val(c)

    region = (vals[idx_region] if idx_region < len(vals) else '').strip().lower()
    if region in {'frontal cortex', 'frontal', 'fr'}:
        frontal_xlsx_total += 1
        participant = (vals[idx_participant] if idx_participant < len(vals) else '').strip()
        total = vals[idx_total] if idx_total is not None and idx_total < len(vals) else ''
        mut = vals[idx_mut] if idx_mut is not None and idx_mut < len(vals) else ''
        lod = vals[idx_lod] if idx_lod is not None and idx_lod < len(vals) else ''

        try:
            total_n = float(total)
        except Exception:
            total_n = None
        try:
            mut_n = float(mut)
        except Exception:
            mut_n = None

        if (total_n is None or total_n <= 0) and (mut_n is None or mut_n <= 0):
            frontal_xlsx_zero.append((participant, region, total, mut, lod))
            frontal_xlsx_by_participant[participant] = frontal_xlsx_by_participant.get(participant, 0) + 1

print('RAW_CSV_FILES', len([f for f in os.listdir(raw_dir) if f.lower().endswith('.csv')]))
print('RAW_STATUS_COUNTS', dict(sorted(status_counts.items())))
print('RAW_FRONTAL_DISTINCT_SAMPLES', len(frontal_samples))
print('RAW_FRONTAL_ROWS', len(frontal_rows))
print('RAW_FRONTAL_SAMPLES_WITH_NO_POSITIVE_ACCEPTED', len(frontal_zero_or_na))
for s, rows in sorted(frontal_zero_or_na, key=lambda x: x[0]):
    print('FRONTAL_ZERO_SAMPLE', s)
    for rr in rows:
        print('  ', rr['file'], rr['target'], rr['accepted'], rr['status'])

print('XLSX_FRONTAL_ROWS', frontal_xlsx_total)
print('XLSX_FRONTAL_ROWS_WITH_ZERO_AMPLIFICATION', len(frontal_xlsx_zero))
for rec in frontal_xlsx_zero:
    print('XLSX_ZERO', rec)
