#!/usr/bin/env python3
"""
For each cluster in the joint SVCluster sites table that has ≥1 BioHeart and
≥1 TOB member, pick the (bioheart, tob) pair with minimum breakpoint distance
and write a TSV of the selected pairs.

Breakpoint distance = |POS_bh - POS_tb| + |END_bh - END_tb|,
tie-broken by |SVLEN_bh - SVLEN_tb|. Pairs must have compatible SVTYPEs:
identical types are allowed, and CNV is treated as compatible with DEL and
DUP (GATK-SV emits CNV when the caller can't commit to a direction). Only
directly opposing pairs (DEL <-> DUP) are rejected.

Inputs (CSVs with columns: CHROM,POS,ID,REF,ALT,QUAL,FILTER,INFO[,MEMBERS]):
  --cluster-sites   joint clustered sites table (MEMBERS column required)
  --bioheart-sites  per-cohort sites table (unprefixed IDs)
  --tob-sites       per-cohort sites table (unprefixed IDs)

Output:
  --out             TSV: cluster_id,chrom,cluster_pos,svtype,
                    bioheart_id,bioheart_pos,bioheart_end,bioheart_svlen,
                    tob_id,tob_pos,tob_end,tob_svlen,
                    bp_distance,n_bh_members,n_tb_members
"""
import click
import pandas as pd

# Pairs of SVTYPEs that should NOT be matched. CNV is intentionally omitted
# — it clusters with DEL or DUP because GATK-SV emits CNV when the caller
# cannot resolve direction.
INCOMPATIBLE_TYPE_PAIRS = frozenset({
    frozenset({'DEL', 'DUP'}),
})


def compatible_svtypes(a: str, b: str) -> bool:
    return frozenset({a, b}) not in INCOMPATIBLE_TYPE_PAIRS


def parse_info(info: str) -> dict:
    if not isinstance(info, str):
        return {}
    out = {}
    for kv in info.split(';'):
        if '=' in kv:
            k, _, v = kv.partition('=')
            out[k] = v
    return out


def load_cohort(path: str, prefix: str) -> dict:
    df = pd.read_csv(path)
    lookup = {}
    for _, r in df.iterrows():
        info = parse_info(r['INFO'])
        end = info.get('END')
        svlen = info.get('SVLEN')
        svtype = info.get('SVTYPE')
        if end is None or svtype is None:
            continue
        lookup[f'{prefix}_{r["ID"]}'] = {
            'chrom': r['CHROM'],
            'pos': int(r['POS']),
            'end': int(end),
            'svlen': int(svlen) if svlen not in (None, '.', '') else None,
            'svtype': svtype,
        }
    return lookup


@click.command()
@click.option('--cluster-sites', required=True)
@click.option('--bioheart-sites', required=True)
@click.option('--tob-sites', required=True)
@click.option('--out', required=True, help='Output TSV path')
def main(cluster_sites, bioheart_sites, tob_sites, out):
    bh = load_cohort(bioheart_sites, 'bioheart')
    tb = load_cohort(tob_sites, 'tob')

    clusters = pd.read_csv(cluster_sites)
    rows = []
    for _, c in clusters.iterrows():
        info = parse_info(c['INFO'])
        members_field = info.get('MEMBERS', '')
        if not members_field:
            continue
        members = [m for m in members_field.split(',') if m]
        bh_members = [m for m in members if m.startswith('bioheart_')]
        tb_members = [m for m in members if m.startswith('tob_')]
        if not bh_members or not tb_members:
            continue

        best = None
        for a in bh_members:
            va = bh.get(a)
            if va is None:
                continue
            for b in tb_members:
                vb = tb.get(b)
                if vb is None:
                    continue
                if not compatible_svtypes(va['svtype'], vb['svtype']):
                    continue
                bp = abs(va['pos'] - vb['pos']) + abs(va['end'] - vb['end'])
                sv_tb = (
                    abs((va['svlen'] or 0) - (vb['svlen'] or 0))
                    if va['svlen'] is not None and vb['svlen'] is not None
                    else 0
                )
                key = (bp, sv_tb)
                if best is None or key < best[0]:
                    best = (key, a, va, b, vb)
        if best is None:
            continue
        (bp, _), a_id, va, b_id, vb = best
        rows.append({
            'cluster_id': c['ID'],
            'chrom': c['CHROM'],
            'cluster_pos': c['POS'],
            'bioheart_id': a_id,
            'bioheart_svtype': va['svtype'],
            'bioheart_pos': va['pos'],
            'bioheart_end': va['end'],
            'bioheart_svlen': va['svlen'],
            'tob_id': b_id,
            'tob_svtype': vb['svtype'],
            'tob_pos': vb['pos'],
            'tob_end': vb['end'],
            'tob_svlen': vb['svlen'],
            'svtype_match': va['svtype'] == vb['svtype'],
            'bp_distance': bp,
            'n_bh_members': len(bh_members),
            'n_tb_members': len(tb_members),
        })

    out_df = pd.DataFrame(rows)
    out_df.to_csv(out, sep='\t', index=False)
    print(
        f'clusters={len(clusters)} common={len(out_df)} '
        f'multi_member={(out_df["n_bh_members"] + out_df["n_tb_members"] > 2).sum()}',
    )


if __name__ == '__main__':
    main()
