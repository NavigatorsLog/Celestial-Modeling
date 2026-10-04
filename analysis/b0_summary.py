"""b0_summary.py  [CL draft, 4 October 2026]  EXPLORATORY ONLY.
Reads progress.csv files from B0 and reports speed and how fast loose material falls.
Its numbers set the P-D12 matched-stage fraction and the cost estimate; they are never evidence about inflation.
Usage: python b0_summary.py run_f1/progress.csv run_f3/progress.csv
"""
import sys, csv
def load(fn):
    head, rows = {}, []
    for line in open(fn):
        if line.startswith('# seed='):
            for kv in line[2:].split():
                k, v = kv.split('='); head[k] = v
        elif line[0].isdigit():
            t, n, loose, calls, wall = line.strip().split(',')[:5]
            rows.append((float(t), int(n), float(loose), int(calls), float(wall)))
    return head, rows
for fn in sys.argv[1:]:
    head, rows = load(fn)
    if len(rows) < 2:
        print(fn, ': not enough progress rows yet'); continue
    m0 = float(head.get('mloose0_Mearth', rows[0][2]))
    t, n, loose, calls, wall = rows[-1]
    rate = t / max(wall, 1e-9) * 3600
    print(f"{fn}: factor {head.get('factor')}, reached {t:.3e} yr in {wall/3600:.2f} h ({rate:.3e} yr per wall hour), N={n}, loose={loose/m0:.3f} of start")
    for frac in (0.9, 0.75, 0.5, 0.25):
        hit = next((r for r in rows if r[2] <= frac*m0), None)
        print(f"   loose material at {frac:.2f} of start: " + (f"t={hit[0]:.3e} yr, wall={hit[4]/3600:.2f} h" if hit else "not reached"))
