"""Download public measurements and produce provenance-preserving, grouped data.

Only candidates.json is agent-visible. outcomes.json and raw source files belong
to the trusted vault. Labels are a declared assay threshold, not patient truth.
"""

import csv
import hashlib
import json
import math
import re
import statistics
import subprocess
from collections import Counter, defaultdict
from concurrent.futures import ThreadPoolExecutor
from operator import itemgetter
from pathlib import Path

from Bio import Align
from Bio.Align import substitution_matrices
from Bio.Seq import Seq
from Bio.SeqUtils import seq1

SCORE_SETS = {
    'FGFR1': {'pem': '10', 'fut': '4'},
    'FGFR2': {'pem': '9', 'fut': '7'},
    'FGFR3': {'pem': '5', 'fut': '2'},
    'FGFR4': {'pem': '8', 'fut': '6'},
}
API = 'https://api.mavedb.org/api/v1/score-sets/urn:mavedb:00001256-a-'
SPLIT_SALT = 'escape-lab-fgfr-v1-fixed-before-models'
LABEL_RULE = 'Each drug: median of 3 finite nonnegative replicates >1.5 and >=2 replicates >1.0; all nucleotide representations agree; otherwise inconclusive.'


def validate_source_identity(metadata, gene, drug):
    expected = 'urn:mavedb:00001256-a-' + SCORE_SETS[gene][drug]
    if metadata['urn'] != expected or metadata['targetGenes'][0]['name'] != gene:
        raise ValueError('Source identity mismatch')
    drug_name = {'pem': 'pemigatinib', 'fut': 'futibatinib'}[drug]
    if metadata['title'].lower() != f'{gene} {drug_name} resistance'.lower():
        raise ValueError('Source assay identity mismatch')
    return expected


def sha256(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def parse_protein_variant(value):
    match = re.fullmatch(r'p\.([A-Z][a-z]{2})(\d+)([A-Z][a-z]{2})', value)
    if not match:
        raise ValueError(f'Unsupported protein variant {value}')
    wt, pos, mut = match.groups()
    wt, mut = seq1(wt), seq1(mut)
    if wt not in 'ACDEFGHIKLMNPQRSTVWY' or mut not in 'ACDEFGHIKLMNPQRSTVWY':
        raise ValueError(f'Non-missense variant {value}')
    return wt, int(pos), mut


def resistance_label(replicates):
    if len(replicates) != 3 or any(v is None or not math.isfinite(v) or v < 0 for v in replicates):
        return None
    return statistics.median(replicates) > 1.5 and sum(v > 1 for v in replicates) >= 2


def parse_measurement(value):
    if value in {'NA', '', None}:
        return None
    result = float(value)
    return result if math.isfinite(result) else None


def collapse_labels(labels):
    if not labels or any(label is None for label in labels):
        return None
    return labels[0] if all(label == labels[0] for label in labels) else None


def aligned_position_map(reference, sequence):
    aligner = Align.PairwiseAligner()
    aligner.mode = 'global'
    aligner.substitution_matrix = substitution_matrices.load('BLOSUM62')
    aligner.open_gap_score = -10
    aligner.extend_gap_score = -0.5
    alignment = aligner.align(reference, sequence)[0]
    mapping = {}
    previous = 1
    for ref_index, query_index in zip(*alignment.indices, strict=True):
        if ref_index >= 0:
            previous = int(ref_index) + 1
        if query_index >= 0:
            mapping[int(query_index) + 1] = previous
    return mapping


def group_id(position):
    return f'aligned-block-{(position - 1) // 12:03d}'


def assign_split(position):
    digest = hashlib.sha256(f'{SPLIT_SALT}:{group_id(position)}'.encode()).hexdigest()
    value = int(digest[:8], 16) % 100
    if value < 60:
        return 'development'
    if value < 80:
        return 'validation'
    return 'final'


def download_sources(raw_dir):
    raw_dir = Path(raw_dir)
    raw_dir.mkdir(parents=True, exist_ok=True)
    tasks = []
    for gene, drugs in SCORE_SETS.items():
        for drug, accession in drugs.items():
            for suffix, extension in [('', 'json'), ('/scores', 'csv')]:
                tasks.append((API + accession + suffix, raw_dir / f'{gene}-{drug}.{extension}'))

    def download(task):
        url, path = task
        if not path.exists():
            temp = path.with_suffix(path.suffix + '.part')
            subprocess.run(['curl', '-fsSL', '--retry', '2', '--max-time', '60', url,
                            '-o', str(temp)], check=True, capture_output=True)
            temp.replace(path)
        return {'url': url, 'file': path.name, 'sha256': sha256(path)}

    with ThreadPoolExecutor(max_workers=4) as pool:
        manifest = list(pool.map(download, tasks))
    path = raw_dir / 'source_manifest.json'
    if path.exists() and json.loads(path.read_text()) != manifest:
        raise ValueError('Source files differ from the recorded snapshot; do not silently replace data.')
    path.write_text(json.dumps(manifest, indent=2) + '\n')
    return manifest


def prepare_dataset(raw_dir, output_dir):
    raw_dir, output_dir = Path(raw_dir), Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)
    source_manifest = json.loads((raw_dir / 'source_manifest.json').read_text())
    expected_sources = {}
    for gene, drugs in SCORE_SETS.items():
        for drug, accession in drugs.items():
            expected_sources[f'{gene}-{drug}.json'] = API + accession
            expected_sources[f'{gene}-{drug}.csv'] = API + accession + '/scores'
    actual_sources = {record['file']: record['url'] for record in source_manifest}
    if actual_sources != expected_sources or len(source_manifest) != len(expected_sources):
        raise ValueError('Unexpected source manifest paths or URLs')
    source_hashes = {record['file']: record['sha256'] for record in source_manifest}
    for record in source_manifest:
        if sha256(raw_dir / record['file']) != record['sha256']:
            raise ValueError(f"Source checksum mismatch: {record['file']}")
    sequences, metadata, tables = {}, {}, {}
    for gene, drugs in SCORE_SETS.items():
        for drug in drugs:
            meta = json.loads((raw_dir / f'{gene}-{drug}.json').read_text())
            accession = validate_source_identity(meta, gene, drug)
            if meta['license']['shortName'] != 'CC BY-SA 4.0':
                raise ValueError('Unexpected source license; inspect before using.')
            target = meta['targetGenes'][0]['targetSequence']
            sequence = str(Seq(target['sequence']).translate()).rstrip('*')
            if gene in sequences and sequences[gene] != sequence:
                raise ValueError(f'Drug assays have different reference sequences for {gene}')
            sequences[gene] = sequence
            metadata[gene, drug] = meta
            with (raw_dir / f'{gene}-{drug}.csv').open() as handle:
                rows = list(csv.DictReader(handle))
            if len({row['hgvs_nt'] for row in rows}) != len(rows):
                raise ValueError('Duplicate nucleotide identity inside a score set')
            if any(not row['accession'].startswith(accession + '#') for row in rows):
                raise ValueError('Score row accession differs from metadata')
            tables[gene, drug] = {row['hgvs_nt']: row for row in rows}
    maps = {gene: aligned_position_map(sequences['FGFR1'], seq)
            for gene, seq in sequences.items()}
    candidates, outcomes = [], {}
    row_count = 0
    for gene in SCORE_SETS:
        pem, fut = tables[gene, 'pem'], tables[gene, 'fut']
        if set(pem) != set(fut):
            raise ValueError(f'Assay pairing mismatch for {gene}')
        by_protein = defaultdict(list)
        for nucleotide in sorted(pem):
            a, b = pem[nucleotide], fut[nucleotide]
            if a['hgvs_pro'] != b['hgvs_pro']:
                raise ValueError('Drug pair protein identity mismatch')
            by_protein[a['hgvs_pro']].append((nucleotide, a, b))
        for protein, variants in sorted(by_protein.items()):
            wt, position, mut = parse_protein_variant(protein)
            if sequences[gene][position - 1] != wt:
                raise ValueError(f'Reference mismatch {gene}:{protein}')
            aligned = maps[gene][position]
            cid = f'{gene}:{protein}'
            records = []
            for nt, a, b in variants:
                row_count += 1
                entry = {'nucleotide_variant': nt, 'protein_variant': protein}
                for drug, row in [('pem', a), ('fut', b)]:
                    columns = sorted(k for k in row if '_Rep' in k)
                    replicates = [parse_measurement(row[k]) for k in columns]
                    label = resistance_label(replicates)
                    score = parse_measurement(row['score'])
                    if label is not None and not math.isclose(score, statistics.median(replicates), rel_tol=1e-7, abs_tol=1e-7):
                        raise ValueError(f'Score differs from replicate median {cid}:{drug}')
                    source = raw_dir / f'{gene}-{drug}.csv'
                    entry[drug] = {
                        'replicates': replicates, 'median': score, 'resistant': label,
                        'accession': row['accession'], 'source_file': source.name,
                        'source_sha256': source_hashes[source.name],
                    }
                records.append(entry)
            a_label = collapse_labels([record['pem']['resistant'] for record in records])
            b_label = collapse_labels([record['fut']['resistant'] for record in records])
            joint = None if a_label is None or b_label is None else a_label and b_label
            candidates.append({
                'candidate_id': cid, 'gene': gene, 'protein_variant': protein,
                'position': position, 'wt': wt, 'mut': mut, 'aligned_position': aligned,
                'group_id': group_id(aligned), 'split': assign_split(aligned),
                'context': sequences[gene][max(0, position - 4):position + 3],
                'reference_sequence_sha256': hashlib.sha256(sequences[gene].encode()).hexdigest(),
            })
            outcomes[cid] = {'dual_resistant': joint, 'records': records,
                            'evidence_ids': [record[drug]['accession'] for record in records
                                             for drug in ['pem', 'fut']]}
    if row_count != 8407:
        raise ValueError(f'Expected verified 8407 paired nucleotide rows, got {row_count}')
    # Published examples belong exclusively to development; move their entire
    # aligned blocks, never just the individual visible answer.
    demo_ids = {'FGFR2:p.Val564Gly', 'FGFR2:p.Asn549Lys'}
    demo_groups = {c['group_id'] for c in candidates if c['candidate_id'] in demo_ids}
    for candidate in candidates:
        if candidate['group_id'] in demo_groups:
            candidate['split'] = 'development'
    candidates.sort(key=itemgetter('candidate_id'))
    public_path, private_path = output_dir / 'candidates.json', output_dir / 'outcomes.json'
    public_path.write_text(json.dumps(candidates, indent=2) + '\n')
    private_path.write_text(json.dumps(outcomes, indent=2, allow_nan=False) + '\n')
    private_path.chmod(0o600)
    manifest = {
        'schema_version': 1, 'label_rule': LABEL_RULE, 'paired_nucleotide_rows': row_count,
        'unique_protein_variants': len(candidates),
        'split_counts': dict(Counter(c['split'] for c in candidates)),
        'split_group_counts': {split: len({c['group_id'] for c in candidates if c['split'] == split})
                               for split in ['development', 'validation', 'final']},
        'split_salt': SPLIT_SALT, 'grouping': '12-residue aligned blocks across all FGFRs',
        'published_demo_groups_forced_to_development': sorted(demo_groups),
        'candidates_sha256': sha256(public_path), 'outcomes_sha256': sha256(private_path),
        'source_manifest_sha256': sha256(raw_dir / 'source_manifest.json'),
        'license': 'CC BY-SA 4.0 for source and derived measurements; code is separate',
        'source_publication': 'https://doi.org/10.1038/s41588-025-02431-8',
        'final_test_status': 'sealed',
    }
    (output_dir / 'manifest.json').write_text(json.dumps(manifest, indent=2) + '\n')
    return manifest
