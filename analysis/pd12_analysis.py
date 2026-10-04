"""pd12_analysis.py  [CL draft, 4 October 2026]  NOT YET FROZEN.
Confirmatory analysis for P-D12 (does radius inflation by factor 3 change mantle-stripping, compared with true size,
up to a matched stage of growth?). The rule is fixed in pd12_config.json and in the preregistration; this script only
applies it. Once frozen, its SHA-256 goes into the preregistration's freeze block and it must not change.

Inputs: a run directory per run, each holding progress.csv, collision_log.csv and collision_report.txt from
formation_run v2, and pd12_config.json naming every run's directory, seed, factor and role.
Usage: python pd12_analysis.py pd12_config.json [--out verdict.json]

Definitions (all applied up to each run's matched-stage time t*):
  t*  first progress row at which the configured pool (loose or small) is at or below stage_fraction of its start.
  stripping event  a collision_log row with fragments made AND either the target lost mass (partial erosion,
      grazing partial erosion, super-catastrophic) or the projectile survived (hit-and-run).
      Rows with fragments where the projectile was absorbed are partial accretion: not stripping.
  S   total fragment mass from stripping events / starting disk mass (2.604 Earth masses).
  F   number of stripping events / number of resolved collisions (rows of the module's collision_report.txt).
Statistics (unit = run pair; bootstrap over pairs, fixed seed, percentile 90% intervals):
  R   = sum of S at factor 3 / sum of S at factor 1 (ratio of sums, defined when pairs contain zeros).
  dF  = mean over pairs of F(factor 3) - F(factor 1).
Verdict: ACCEPTED iff R interval inside [r_lo, r_hi] and dF interval inside [-f_tol, +f_tol];
  otherwise NOT ACCEPTED, unless the same-condition control (factor 1 vs factor 1) also fails, then INCONCLUSIVE.
  The reproducibility check must pass first or the verdict is INSTRUMENT FAILURE.
  rule_version v2 (approved by CBH, 4 October 2026): if either interval excludes no difference (R = 1 or dF = 0), the verdict
  is NOT ACCEPTED (difference shown) whatever the control says; v1 was the rule before that amendment.
"""
import json, sys, random

MEARTH = 3.003e-6
DISK0 = 2.604  # Earth masses

def read_progress(path):
    head, rows, cols = {}, [], None
    for line in open(path):
        if line.startswith('# seed='):
            for kv in line[2:].split():
                k, v = kv.split('='); head[k] = v
        elif line.startswith('t_yr'):
            cols = line.strip().split(',')
        elif line[:1].isdigit():
            vals = [float(x) for x in line.strip().split(',')]
            rows.append(dict(zip(cols, vals)))
    return head, rows

def stage_time(rows, measure, frac):
    col = 'loose_Mearth' if measure == 'loose' else 'small_Mearth'
    start = rows[0][col]
    for r in rows:
        if r[col] <= frac * start:
            return r['t_yr']
    return None

def stats_for_run(d, measure, frac):
    _, rows = read_progress(f"{d}/progress.csv")
    t_star = stage_time(rows, measure, frac)
    if t_star is None:
        return None
    s_mass, n_strip = 0.0, 0
    with open(f"{d}/collision_log.csv") as f:
        next(f)
        for line in f:
            t, th, tm0, tm1, ph, pm0, pm1, nfrag, mfrag, ret = line.strip().split(',')
            if float(t) > t_star: break
            if int(nfrag) > 0 and (float(tm1) < float(tm0) or float(pm1) > 0):
                s_mass += float(mfrag); n_strip += 1
    n_coll = 0
    for line in open(f"{d}/collision_report.txt"):
        try:
            if float(line.split()[0]) <= t_star: n_coll += 1
        except (ValueError, IndexError):
            continue
    return dict(t_star=t_star, S=s_mass / MEARTH / DISK0, F=(n_strip / n_coll) if n_coll else 0.0,
                n_strip=n_strip, n_coll=n_coll)

def pair_stats(pairs, B, seed):
    def R_of(ps):
        den = sum(p[0]['S'] for p in ps)
        return (sum(p[1]['S'] for p in ps) / den) if den > 0 else float('inf')
    def dF_of(ps):
        return sum(p[1]['F'] - p[0]['F'] for p in ps) / len(ps)
    rng = random.Random(seed)
    Rs, dFs = [], []
    for _ in range(B):
        ps = [pairs[rng.randrange(len(pairs))] for _ in pairs]
        Rs.append(R_of(ps)); dFs.append(dF_of(ps))
    Rs.sort(); dFs.sort()
    lo, hi = int(0.05 * B), int(0.95 * B) - 1
    return dict(R=R_of(pairs), R_ci=[Rs[lo], Rs[hi]], dF=dF_of(pairs), dF_ci=[dFs[lo], dFs[hi]])

def passes(st, cfg):
    return (cfg['r_lo'] <= st['R_ci'][0] and st['R_ci'][1] <= cfg['r_hi']
            and -cfg['f_tol'] <= st['dF_ci'][0] and st['dF_ci'][1] <= cfg['f_tol'])

def difference_shown(st, cfg):
    # rule v2 (approved by CBH): an interval that excludes 'no difference' (R = 1, dF = 0) shows a real
    # difference, which the control's noise cannot explain away
    return st['R_ci'][0] > 1 or st['R_ci'][1] < 1 or st['dF_ci'][0] > 0 or st['dF_ci'][1] < 0

def main():
    cfg = json.load(open(sys.argv[1]))
    out = sys.argv[sys.argv.index('--out') + 1] if '--out' in sys.argv else None
    m, fr = cfg['stage_measure'], cfg['stage_fraction']
    dirs = {r['id']: r['dir'] for r in cfg['runs']}
    res = dict(config=cfg)
    # reproducibility ceiling first: the repeat run must reproduce its original's collision report exactly
    rep = cfg['repeat']
    res['reproducible'] = (open(f"{dirs[rep['original']]}/collision_report.txt").read()
                           == open(f"{dirs[rep['repeat']]}/collision_report.txt").read())
    runs = {r['id']: dict(r, stats=stats_for_run(r['dir'], m, fr)) for r in cfg['runs']}
    missing = [k for k, r in runs.items() if r['stats'] is None]
    res['missing_stage'] = missing
    if missing or not res['reproducible']:
        res['verdict'] = 'INSTRUMENT FAILURE' if not res['reproducible'] else 'INCOMPLETE: stage not reached by ' + ', '.join(missing)
    else:
        main_pairs = [(runs[p['f1']]['stats'], runs[p['f3']]['stats']) for p in cfg['pairs']]
        ctrl_pairs = [(runs[p['a']]['stats'], runs[p['b']]['stats']) for p in cfg['controls']]
        res['main'] = pair_stats(main_pairs, cfg['bootstrap_B'], cfg['bootstrap_seed'])
        res['control'] = pair_stats(ctrl_pairs, cfg['bootstrap_B'], cfg['bootstrap_seed'] + 1)
        rule = cfg.get('rule_version', 'v1')
        if passes(res['main'], cfg):
            res['verdict'] = 'ACCEPTED'
        elif rule == 'v2' and difference_shown(res['main'], cfg):
            res['verdict'] = 'NOT ACCEPTED (difference shown)'
        elif not passes(res['control'], cfg):
            res['verdict'] = 'INCONCLUSIVE'
        else:
            res['verdict'] = 'NOT ACCEPTED'
        res['rule_version'] = rule
    res['per_run'] = {k: r['stats'] for k, r in runs.items()}
    txt = json.dumps(res, indent=1)
    if out: open(out, 'w').write(txt)
    print('VERDICT:', res['verdict'])
    if 'main' in res:
        print('main   R=%.3f  90%% [%.3f, %.3f]   dF=%+.4f  90%% [%+.4f, %+.4f]' % (res['main']['R'], *res['main']['R_ci'], res['main']['dF'], *res['main']['dF_ci']))
        print('control R=%.3f  90%% [%.3f, %.3f]   dF=%+.4f  90%% [%+.4f, %+.4f]' % (res['control']['R'], *res['control']['R_ci'], res['control']['dF'], *res['control']['dF_ci']))

if __name__ == '__main__':
    main()
