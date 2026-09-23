"""Descriptive comparisons of reported HepG2 MPRA and liver eQTL effects.

Rows are MPRA-record × eQTL-association matches, not independent variants.
Allele metadata are audited; source effect-direction definitions remain unverified.
"""
from pathlib import Path
import hashlib
import importlib.metadata
import json
import platform

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy.stats import pearsonr, spearmanr

FILES = {
    'gtex': 'GTEx_Liver_eQTLs.signif_pairs.tsv.gz',
    'finemap': 'gtex_finemap_eqtls_pip01_with_spdi.csv',
    'mohlke': 'Mohlke_HepG2_MPRA_LipidTraits.tsv',
    'tewhey': 'hepg2_combined_final_processed.tsv.gz',
}
PROVENANCE = {
    'gtex': {'source': 'https://www.synapse.org/Synapse:syn69058503.1', 'version': '1'},
    'mohlke': {'source': 'https://www.synapse.org/Synapse:syn68541428.2', 'version': '2'},
    'finemap': {'source': None, 'note': 'Original download and SPDI conversion provenance not recovered.'},
    'tewhey': {'source': None, 'note': 'Original download, processing and emVar definition not recovered.'},
}


def require_columns(df, columns):
    missing = sorted(set(columns) - set(df.columns))
    if missing:
        raise ValueError(f'Missing required columns: {missing}')


def parse_bool(values):
    """Do not treat the nonempty string 'False' as True."""
    normalized = values.astype('string').str.strip().str.lower()
    bad = normalized.isna() | ~normalized.isin(['true', 'false', '1', '0'])
    if bad.any():
        raise ValueError(f'Invalid or missing boolean labels: {values[bad].unique().tolist()}')
    return normalized.map({'true': True, 'false': False, '1': True, '0': False}).astype(bool)


def finite_rows(df, columns):
    result = df.copy()
    for col in columns:
        result[col] = pd.to_numeric(result[col], errors='raise')
    keep = np.isfinite(result[columns].to_numpy(dtype=float)).all(axis=1)
    return result.loc[keep].copy(), int((~keep).sum())


def metadata_audit(df, first, second):
    """Compare literal allele metadata with SPDI; do not infer effect direction.

    Non-SNV representations can require normalization, so report them separately.
    """
    parts = df['SPDI'].astype('string').str.split(':', expand=True)
    if parts.shape[1] != 4:
        raise ValueError('Expected four-component SPDI identifiers')
    ref, alt = parts[2], parts[3]
    a, b = df[first].astype('string'), df[second].astype('string')
    snv = ref.str.fullmatch('[ACGT]') & alt.str.fullmatch('[ACGT]')
    same = ((a == ref) & (b == alt)).fillna(False)
    reversed_order = ((a == alt) & (b == ref)).fillna(False)
    status = pd.Series('non_snv_or_unresolved', index=df.index)
    status.loc[snv.fillna(False)] = 'snv_mismatch'
    status.loc[same] = 'matches_spdi_order'
    status.loc[reversed_order & ~same] = 'reversed_spdi_order'
    return status


def load_eqtls(data_dir, pip_threshold=0.1):
    if not 0 <= pip_threshold < 1:
        raise ValueError('PIP threshold must be in [0, 1)')
    prefix = 'GTEx_Liver_eQTLs.'
    g = pd.read_csv(Path(data_dir) / FILES['gtex'], sep='\t', low_memory=False)
    require_columns(g, ['SPDI', 'ENSEMBL_ID', 'A1', 'A2', prefix+'slope',
                        prefix+'pval_nominal', prefix+'pval_nominal_threshold'])
    g['eqtl_record_id'] = np.arange(len(g))
    g = g.rename(columns={'ENSEMBL_ID': 'gene', prefix+'slope': 'eqtl_effect',
                          prefix+'pval_nominal': 'rank_score',
                          prefix+'pval_nominal_threshold': 'p_threshold'})
    g, g_bad = finite_rows(g, ['eqtl_effect', 'rank_score', 'p_threshold'])
    if not g[['rank_score', 'p_threshold']].ge(0).all().all() or not g[['rank_score', 'p_threshold']].le(1).all().all():
        raise ValueError('GTEx p-values must be in [0, 1]')
    g = g.loc[g['rank_score'] < g['p_threshold']].copy()
    g['allele_metadata'] = metadata_audit(g, 'A1', 'A2')
    # Stream the multi-tissue CSV and retain only Liver, including when absent.
    liver_chunks = []
    for chunk in pd.read_csv(Path(data_dir) / FILES['finemap'], chunksize=100000):
        require_columns(chunk, ['tissue', 'method', 'spdi', 'gene', 'pip', 'beta_marginal', 'allele1', 'allele2'])
        liver_chunks.append(chunk.loc[chunk['tissue'].eq('Liver')].copy())
    f = pd.concat(liver_chunks, ignore_index=True)
    if f.empty:
        raise ValueError('No Liver records in the fine-mapping input')
    f['eqtl_record_id'] = np.arange(len(f))
    f = f.rename(columns={'spdi': 'SPDI', 'pip': 'rank_score', 'beta_marginal': 'eqtl_effect'})
    f, f_bad = finite_rows(f, ['rank_score', 'eqtl_effect'])
    if not f['rank_score'].between(0, 1).all():
        raise ValueError('PIP must be in [0, 1]')
    f = f.loc[f['rank_score'] > pip_threshold].copy()
    f['allele_metadata'] = metadata_audit(f, 'allele1', 'allele2')
    frames = {'GTEx': g, 'FINEMAP': f.loc[f['method'].eq('FINEMAP')].copy(),
              'SUSIE': f.loc[f['method'].eq('SUSIE')].copy()}
    audit = []
    for name, frame in frames.items():
        if frame.empty:
            raise ValueError(f'No eligible records for {name}')
        if frame[['SPDI', 'gene']].isna().any().any():
            raise ValueError(f'Missing variant/gene identifiers in {name}')
        audit.append({'dataset': name, 'rows': len(frame), 'unique_variants': frame.SPDI.nunique(),
                      'unique_variant_gene_pairs': len(frame[['SPDI','gene']].drop_duplicates()),
                      'allele_metadata_counts': json.dumps(frame.allele_metadata.value_counts().to_dict(), sort_keys=True)})
    return frames, {'eqtl_inputs': audit, 'excluded_nonfinite_gtex': g_bad,
                    'excluded_nonfinite_liver_finemap': f_bad}


def load_mpra(data_dir, dataset, q_threshold=0.05):
    if dataset not in ('mohlke', 'tewhey'):
        raise ValueError('Unknown MPRA dataset')
    if not 0 < q_threshold < 1:
        raise ValueError('MPRA q threshold must be in (0, 1)')
    d = pd.read_csv(Path(data_dir) / FILES[dataset], sep='\t', low_memory=False)
    require_columns(d, ['SPDI', 'A1', 'A2'])
    d['mpra_record_id'] = np.arange(len(d))
    if dataset == 'mohlke':
        prefix = 'Mohlke_HepG2_MPRA.'
        require_columns(d, [prefix+'log2FC', prefix+'minusLog10QValue', prefix+'variantPos'])
        d = d.rename(columns={prefix+'log2FC':'mpra_effect', prefix+'minusLog10QValue':'neglog10q',
                              prefix+'variantPos':'variant_position'})
        d, dropped = finite_rows(d, ['mpra_effect', 'neglog10q', 'variant_position'])
        if d.neglog10q.lt(0).any():
            raise ValueError('Negative -log10(q) in Mohlke input')
        d['highlight'] = d.neglog10q.gt(-np.log10(q_threshold))
    else:
        require_columns(d, ['log2Skew', 'emVar'])
        d['highlight'] = parse_bool(d.emVar)
        d = d.rename(columns={'log2Skew':'mpra_effect'})
        d, dropped = finite_rows(d, ['mpra_effect'])
    if d.SPDI.isna().any():
        raise ValueError('Missing MPRA SPDI identifiers')
    d['allele_metadata'] = metadata_audit(d, 'A1', 'A2')
    audit = {'dataset': dataset, 'rows':len(d), 'unique_variants':d.SPDI.nunique(),
             'duplicated_variant_keys':int(d.SPDI.value_counts().gt(1).sum()),
             'highlight_rows':int(d.highlight.sum()), 'excluded_nonfinite':dropped,
             'allele_metadata_counts': d.allele_metadata.value_counts().to_dict()}
    return d, audit


def match_records(eqtl, mpra, method, top_only=True):
    # Ties are resolved by gene ID, effect, then source row ID (not merge order).
    e = eqtl.sort_values(['SPDI','rank_score','gene','eqtl_effect','eqtl_record_id'],
                        ascending=[True, method == 'GTEx', True, True, True], kind='stable')
    if top_only:
        e = e.drop_duplicates('SPDI')
    e = e[['SPDI','gene','eqtl_effect','eqtl_record_id','allele_metadata']].rename(
        columns={'allele_metadata':'eqtl_allele_metadata'})
    m = mpra[['SPDI','mpra_effect','highlight','mpra_record_id','allele_metadata']].rename(
        columns={'allele_metadata':'mpra_allele_metadata'})
    # Retain MPRA source row IDs; distinct constructs must not be silently collapsed.
    result = e.merge(m, on='SPDI', how='inner', validate='one_to_many' if top_only else 'many_to_many')
    result['metadata_consistent'] = (result.eqtl_allele_metadata.eq('matches_spdi_order') &
                                     result.mpra_allele_metadata.eq('matches_spdi_order'))
    if not result.columns.is_unique:
        raise ValueError('Duplicate output column labels')
    return result


def summarize(df):
    x, y = df.mpra_effect, df.eqtl_effect
    r = rho = np.nan
    if len(df) >= 2 and x.nunique() > 1 and y.nunique() > 1:
        r, rho = float(pearsonr(x, y).statistic), float(spearmanr(x, y).statistic)
    # Zero effects remain in the denominator but never count as concordant.
    concordance = float(((x > 0) & (y > 0) | (x < 0) & (y < 0)).mean()) if len(df) else np.nan
    return {'n_rows':len(df), 'n_unique_variants':df.SPDI.nunique(),
            'n_variant_gene_pairs':len(df[['SPDI','gene']].drop_duplicates()),
            'n_mpra_records':df.mpra_record_id.nunique(), 'pearson_r':r, 'spearman_rho':rho,
            'concordance':concordance, 'n_zero_effect_rows':int(((x == 0) | (y == 0)).sum())}


def compare(frames, mpra, dataset):
    versions = {'position125':mpra.loc[mpra.variant_position.eq(125)], 'all_positions':mpra} if dataset == 'mohlke' else {'all_records':mpra}
    panels, rows = {}, []
    for version, m in versions.items():
        for method, e in frames.items():
            for top_only in (True, False):
                selection = 'top' if top_only else 'all'
                joined = match_records(e, m, method, top_only)
                panels[(version, method, selection)] = joined
                for label, subset in [('all', joined), ('highlighted', joined.loc[joined.highlight])]:
                    rows.append({'dataset':dataset, 'version':version, 'eqtl':method,
                                 'selection':selection, 'subset':label, **summarize(subset)})
                # Sensitivity analysis only: matching metadata do not prove the
                # effect-allele convention, and this also excludes unresolved indels.
                for label, mask in [('metadata_consistent', joined.metadata_consistent),
                                    ('highlighted_metadata_consistent', joined.highlight & joined.metadata_consistent)]:
                    rows.append({'dataset':dataset, 'version':version, 'eqtl':method,
                                 'selection':selection, 'subset':label, **summarize(joined.loc[mask])})
    return panels, pd.DataFrame(rows)


def plot_panels(panels, dataset, out_dir, q_threshold=0.05):
    out_dir = Path(out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)
    highlighted_label = f'MPRA q < {q_threshold:g}' if dataset == 'mohlke' else 'emVar'
    xlabel = 'Reported MPRA log2FC' if dataset == 'mohlke' else 'Reported MPRA log2Skew'
    saved = []
    for version in dict.fromkeys(k[0] for k in panels):
        for methods in [('GTEx',), ('FINEMAP','SUSIE')]:
            fig, axes = plt.subplots(len(methods), 2, figsize=(14, 6*len(methods)), squeeze=False)
            for i, method in enumerate(methods):
                for j, selection in enumerate(['top','all']):
                    ax = axes[i,j]
                    df = panels[(version,method,selection)]
                    for mask, color, label in [(~df.highlight, '0.65', 'Other MPRA records'),
                                                (df.highlight, '#c93439', highlighted_label)]:
                        d = df.loc[mask]
                        ax.scatter(d.mpra_effect, d.eqtl_effect, s=12, alpha=.55, c=color,
                                   edgecolors='none', rasterized=len(df) > 10000, label=label)
                    ax.axhline(0, ls='--', color='0.5', lw=.8)
                    ax.axvline(0, ls='--', color='0.5', lw=.8)
                    text = []
                    for label, d in [('All matches', df), ('Highlighted', df.loc[df.highlight])]:
                        s = summarize(d)
                        text.append(f"{label}: {s['n_rows']:,} rows / {s['n_unique_variants']:,} variants\n"
                                    f"r={s['pearson_r']:.2f}; rho={s['spearman_rho']:.2f}; same sign={s['concordance']:.1%}")
                    ax.text(.02,.98,'\n\n'.join(text), transform=ax.transAxes, va='top', fontsize=8,
                            bbox={'facecolor':'white','alpha':.9,'edgecolor':'0.75'})
                    ax.set(xlabel=xlabel, ylabel='Reported eQTL slope' if method=='GTEx' else 'Reported marginal beta',
                           title=f"{method} | {'One eGene per variant' if selection=='top' else 'All eligible eGene associations'}")
                    ax.legend(loc='lower right', fontsize=8)
            fig.suptitle(f'{dataset.title()} HepG2 vs liver eQTL | {version.replace("_", " ")}', fontsize=14)
            fig.tight_layout(rect=(0,0,1,.97))
            stem = f'{dataset}_{version}_{"gtex" if methods == ("GTEx",) else "finemap_susie"}'
            target = out_dir / f'{stem}.png'
            fig.savefig(target, dpi=300, bbox_inches='tight')
            saved.append(target.name)
            plt.show()
            plt.close(fig)
    return saved


def describe_run(data_dir, dataset, audit, pip_threshold, q_threshold):
    """Return provenance for display in the notebook; no separate result files."""
    inputs = []
    for key in ['gtex','finemap',dataset]:
        path = Path(data_dir) / FILES[key]
        digest = hashlib.sha256()
        with path.open('rb') as f:
            for chunk in iter(lambda:f.read(2**20), b''):
                digest.update(chunk)
        inputs.append({'file':path.name, 'sha256':digest.hexdigest(), 'bytes':path.stat().st_size, **PROVENANCE[key]})
    metadata = {'dataset':dataset, 'inputs':inputs, 'audit':audit,
                'parameters':{'pip_strictly_greater_than':pip_threshold, 'mohlke_q_strictly_less_than':q_threshold},
                'python':platform.python_version(),
                'packages':{p:importlib.metadata.version(p) for p in ['numpy','pandas','scipy','matplotlib','nbformat','nbclient']},
                'analysis_unit':'MPRA record x eQTL association; repeated/linked variants are not independent',
                'direction_status':'Reported effects; source direction definitions not independently verified',
                'inference':'Descriptive correlations only; no independence-based p-values reported'}
    return metadata
