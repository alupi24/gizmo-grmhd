#!/usr/bin/env python3
"""Figures of the GR hydro tests of Lupi (2023, MNRAS 519, 1115) for the Kokkos port, in the layout of the paper's own
plotting scripts (old_cluster/testIC/*/solution.py, test_particles/solutions/get_*.py, TOV/rho0_plot.py).

usage: make_figures.py <run_base> [figure ...]
  run_base: directory holding the runs (see RUNS below for the expected sub-directories)
  figure:   any of shocktubes, orbits, epicyclic, apsidal, freefall, bondi, tov (default: all that have their runs)
Figures are written next to this script as fig_*.png (and .pdf).
"""
import sys, os, glob, h5py, numpy as np, matplotlib
matplotlib.use('Agg'); import matplotlib.pyplot as plt
from scipy.stats import binned_statistic
from scipy.interpolate import UnivariateSpline
from scipy.optimize import curve_fit

OUT = os.path.dirname(os.path.abspath(__file__))
plt.rcParams.update({'lines.markersize': 3, 'markers.fillstyle': 'none', 'axes.labelsize': 16, 'xtick.labelsize': 12, 'ytick.labelsize': 12,
                     'mathtext.fontset': 'cm', 'font.family': 'serif'})
CLR = {'MFM': 'tab:blue', 'MFV': 'tab:purple'}; SYM = {'MFM': 's', 'MFV': 'o'}

RUNS = {  # sub-directories of run_base
    'shocktubes': {'marti1': 'sr/marti1_{}', 'marti2': 'sr/marti2_{}', 'marti2vt': 'sr/marti2vt_{}', 'marti2vt16k': 'sr/marti2vt16k_{}'},  # 4000 cells; marti2vt16k: the paper's 16000
    'orbits': [('k_schw_orbit_newdrift', 0.0, 10.0, 'tab:green'), ('k_kerr_orbit_newdrift', 1.0, 2.0, 'tab:orange')],
    'epicyclic': [('tp/epi_a-1.0', -1.0, 'tab:blue', 0, 99), ('tp/epi_a0.0', 0.0, 'tab:red', 0, 99), ('tp/epi_a1.0_dt001', 1.0, 'tab:green', 6, 99), ('tp/epi_a1.0_inner', 1.0, 'tab:green', 0, 6)],  # (run, spin, colour, r0 min, r0 max): the a=1 inner orbits need dt 0.002
    'apsidal': [('tp/aps_a0.1', 0.1), ('tp/aps_a-0.1', -0.1), ('tp/aps_a0.0', 0.0)],
    'freefall': 'tp/freefall',
    'bondi': {'MFM': ('k_bondi_BL105', 30, 'bl'), 'MFM (Kerr-Schild)': ('k_bondi_KS105', 30, 'ks')},  # (run, time, coordinates): the paper compares at t=30 in Boyer-Lindquist time; the KS run uses the converted IC
    'tov': {'MFM': 'tov/tov_30k_r025'},
    'polytrope': {'range factor 1.0': 'poly/n15_fac1.0', 'range factor 0.5': 'poly/n15_fac0.5', 'range factor 0.25': 'poly/n15_fac0.25', 'positivity step off': 'poly/n15_fac1e-10'},  # Newtonian n=1.5 polytrope: the positivity-limiter range test
}

def snapshots(run):
    return sorted(glob.glob(os.path.join(run, 'output', 'snapshot_*.hdf5')), key=lambda f: int(os.path.basename(f)[9:-5]))

def save(fig, name):
    for ax in fig.axes:  # dense particle scatters as raster in the PDF (vector PDFs of 1e5 points are tens of MB)
        for line in ax.lines:
            if len(line.get_xdata()) > 2000: line.set_rasterized(True)
    for ext in ('png', 'pdf'): fig.savefig(os.path.join(OUT, name + '.' + ext), dpi=200 if ext == 'png' else None, bbox_inches='tight')
    plt.close(fig); print('  saved', name)

def r_BL(x, a):
    R2 = np.sum(x*x, axis=-1); z = x[..., 2]; return np.sqrt((R2 - a*a + np.sqrt((R2 - a*a)**2 + 4*a*a*z*z))/2)

# ------------------------------------------------------------------ special-relativistic shock tubes (paper Figs 2-4, SR/solution.py)
def fig_shocktubes(base):
    import srrp
    gamma = 5./3; xs = np.linspace(0, 1, 5000); x0 = 0.5
    cases = {'marti1': (srrp.State(rho=10, pressure=40./3, vx=0, vt=0.), srrp.State(rho=1.0, pressure=1e-6, vx=0, vt=0.), 0.3, [[-1, 11], [-0.1, 1.1], [-1, 15.0], [-0.1, 2.1]], False),
             'marti2': (srrp.State(rho=1., pressure=1000., vx=0, vt=0.), srrp.State(rho=1., pressure=0.01, vx=0, vt=0.), 0.2, [[-0.5, 15], [-0.1, 1.1], [-100, 1100.0], [-100, 1700]], True),
             'marti2vt': (srrp.State(rho=1., pressure=1000., vx=0, vt=0.9), srrp.State(rho=1., pressure=0.01, vx=0, vt=0.9), 0.6, [[-0.5, 5.5], [-0.05, 0.45], [-100, 1100.0], [0.0, 1.1]], False)}
    cases['marti2vt16k'] = cases['marti2vt']
    def L2(xsim, ysim, yth): f = UnivariateSpline(xs, yth, s=0); return np.sqrt(np.mean((ysim - f(xsim))**2))/np.max(np.abs(yth))
    for name, (sL, sR, t, ylim, inset) in cases.items():
        runs = {h: os.path.join(base, RUNS['shocktubes'][name].format(h)) for h in ('MFM', 'MFV')}
        runs = {h: r for h, r in runs.items() if os.path.exists(os.path.join(r, 'output', 'snapshot_%03d.hdf5' % int(round(t*100))))}
        if not runs: print('  shock tube %s: no runs' % name); continue
        st = srrp.Solver().solve(sL, sR, gamma).getState((xs - x0)/t)
        vt = name.startswith('marti2vt')
        panels = [('RestMassDensity', st.rho, r'$\rho$'), ('vx', st.vx, r'$v_x$'), ('vy' if vt else 'u', st.vt if vt else st.pressure/(gamma-1)/st.rho, r'$v_y$' if vt else r'$u$'), ('Pressure', st.pressure, r'$P$')]
        fig, ax = plt.subplots(2, 2, figsize=[14, 8], sharex=True); plt.subplots_adjust(hspace=0); ax = [ax[0, 0], ax[0, 1], ax[1, 0], ax[1, 1]]
        for h, run in runs.items():
            g = h5py.File(os.path.join(run, 'output', 'snapshot_%03d.hdf5' % int(round(t*100))))['PartType0']
            x = g['Coordinates'][:, 0]; sel = (x > 0.5) & (x < 1.5); x = x[sel] - 0.5; v = g['Velocities'][:][sel]; rho = g['RestMassDensity'][:][sel]; u = g['InternalEnergy'][:][sel]
            data = {'RestMassDensity': rho, 'vx': v[:, 0], 'vy': v[:, 1], 'u': u, 'Pressure': (gamma-1)*rho*u}
            msg = '  %s %s (N=%d in the window):' % (name, h, sel.sum())
            for k, (q, yth, lbl) in enumerate(panels):
                ax[k].plot(x, data[q], SYM[h], color=CLR[h], label=h, ls='none'); msg += '  L2(%s)=%.2e' % (q, L2(x, data[q], yth))
            print(msg)
        for k, (q, yth, lbl) in enumerate(panels):
            ax[k].plot(xs, yth, 'k--'); ax[k].set_ylim(ylim[[0, 1, 3, 2][k]]); ax[k].set_ylabel(lbl); ax[k].yaxis.set_major_locator(plt.MaxNLocator(5))
        ax[3].legend(fontsize=12); ax[1].text(0.05, 0.4, r'$t=%g$' % t, fontsize=14, transform=ax[1].transAxes)
        ax[2].set_xlabel(r'$x$'); ax[3].set_xlabel(r'$x$'); ax[3].set_xlim(0, 1)
        if inset:
            ai = fig.add_axes([0.16, 0.6, 0.15, 0.25])
            for h, run in runs.items():
                g = h5py.File(os.path.join(run, 'output', 'snapshot_%03d.hdf5' % int(round(t*100))))['PartType0']; x = g['Coordinates'][:, 0]; sel = (x > 0.5) & (x < 1.5)
                ai.plot(x[sel] - 0.5, g['RestMassDensity'][:][sel], SYM[h], color=CLR[h], ls='none')
            ai.plot(xs, st.rho, 'k--'); ai.set_xlim(0.687, 0.702); ai.set_ylim(ylim[0]); ai.yaxis.set_major_locator(plt.MaxNLocator(5))
        save(fig, 'fig_shocktube_' + name)

# ------------------------------------------------------------------ circular orbits (paper Fig. 7, get_orbit_paper.py)
def fig_orbits(base):
    fig = plt.figure(figsize=[7, 7]); theta = np.linspace(0, 2*np.pi, 10000); found = False
    for run, a, rini, cl in RUNS['orbits']:
        files = snapshots(os.path.join(base, run))
        if not files: continue
        found = True
        for f in files:
            x = h5py.File(f)['PartType0/Coordinates'][:]; R = np.linalg.norm(x, axis=1); rp = np.sqrt(R*R - a*a)
            plt.plot(x[:, 0]/R*rp, x[:, 1]/R*rp, 'o', color=cl, ls='none', markersize=5)
        plt.plot(rini*np.cos(theta), rini*np.sin(theta), 'k', lw=3)
        X = np.array([h5py.File(f)['PartType0/Coordinates'][:] for f in files]); print('  orbit a=%g r=%g: max |r_BL - r_ini| = %.2e over %d snapshots' % (a, rini, np.abs(r_BL(X, a) - rini).max(), len(files)))
    if not found: print('  orbits: no runs'); plt.close(fig); return
    plt.gca().yaxis.set_major_locator(plt.MaxNLocator(5)); plt.xlim(-11, 11); plt.ylim(-11, 11); plt.xlabel(r'$x$'); plt.ylabel(r'$y$')
    save(fig, 'fig_orbits')

# ------------------------------------------------------------------ epicyclic frequencies (paper Fig. 8, get_epi.py; the FFT peak refined by a sinusoid fit)
def fig_epicyclic(base):
    fig = plt.figure(figsize=[7, 5]); found = False
    labelled = set()
    for run, a, cl, r0min, r0max in RUNS['epicyclic']:
        files = snapshots(os.path.join(base, run))
        if not files: continue
        found = True
        r = np.full((len(files), 32), np.nan); t = np.zeros(len(files))
        for i, f in enumerate(files):
            h = h5py.File(f); x = h['PartType0/Coordinates'][:]; ids = h['PartType0/ParticleIDs'][:].astype(int) - 1; r[i, ids] = r_BL(x, a); t[i] = h['Header'].attrs['Time']
        freq = np.fft.fftfreq(len(t), d=(t[-1]-t[0])/(len(t)-1)); kth = lambda rr: 1./(rr**1.5 + a)*np.sqrt(np.maximum(1 - (6/rr - 8*a/rr**1.5 + 3*a*a/rr/rr), 0))
        errs = []
        for i in range(32):
            if np.isnan(r[:, i]).any() or not (r0min <= r[0, i] < r0max): continue
            y = r[:, i] - r[:, i].mean(); k = np.fft.fft(y); k[np.abs(freq) < 1e-4] = 0; w = abs(freq[1:][np.argmax(np.abs(k[1:]))])*2*np.pi
            model = lambda tt, A, w, ph, c: A*np.cos(w*tt + ph) + c; best = None
            for w0 in (w, max(kth(r[0, i]), 1e-4)):
                for ph0 in (0., np.pi/2):
                    try:
                        p, _ = curve_fit(model, t, y, p0=[(y.max()-y.min())/2, w0, ph0, 0.], maxfev=20000); res = np.sum((y - model(t, *p))**2)
                        if best is None or res < best[0]: best = (res, abs(p[1]))
                    except Exception: pass
            if best: w = best[1]
            plt.plot(r[0, i], w, 'o', c=cl, markersize=6, fillstyle='full')
            if kth(r[0, i]) > 0: errs.append(abs(w/kth(r[0, i]) - 1))
        rr = np.linspace(1, 20, 1000); plt.plot(rr, kth(rr), zorder=0, lw=3, color=cl, label=None if a in labelled else r'$a=%g$' % a); labelled.add(a)
        print('  epicyclic a=%g (%s, r0 in [%g,%g)): median relative error %.1e, max %.1e (%d particles)' % (a, run, r0min, r0max, np.median(errs), np.max(errs), len(errs)))
    if not found: print('  epicyclic: no runs'); plt.close(fig); return
    plt.xlabel(r'$r$'); plt.ylabel(r'$\kappa$'); plt.xlim(0, 20); plt.ylim(0, 0.08); plt.gca().xaxis.set_major_locator(plt.MaxNLocator(5)); plt.legend(fontsize=14)
    save(fig, 'fig_epicyclic')

# ------------------------------------------------------------------ apsidal precession (paper Fig. 9, get_apsidal.py)
def fig_apsidal(base):
    fig = plt.figure(figsize=[6, 6]); found = False
    for run, a in RUNS['apsidal']:
        files = snapshots(os.path.join(base, run))[:500]
        if not files: continue
        found = True
        X = []
        for f in files:
            h = h5py.File(f); ids = np.argsort(h['PartType0/ParticleIDs'][:]); X.append(h['PartType0/Coordinates'][:][ids][0])
        X = np.array(X); plt.plot(X[:, 0], X[:, 1], lw=2, ls='-', label=r'$a=%g$' % a)
        r = r_BL(X, a); phi = np.unwrap(np.arctan2(X[:, 1], X[:, 0])); i = np.where((r[1:-1] < r[:-2]) & (r[1:-1] < r[2:]))[0] + 1
        if len(i) > 1: print('  apsidal a=%g: periapsis advance per orbit %s deg' % (a, np.round(np.degrees(np.diff(phi[i]) - 2*np.pi), 2)))
    if not found: print('  apsidal: no runs'); plt.close(fig); return
    angle = np.radians(82.4); plt.plot(90*np.cos(angle), 90*np.sin(angle), 'x', color='none', markeredgecolor='k', markersize=10)
    plt.legend(fontsize=12); plt.xlim(-95, 95); plt.ylim(-95, 95); plt.xlabel(r'$x$'); plt.ylabel(r'$y$')
    save(fig, 'fig_apsidal')

# ------------------------------------------------------------------ radial free fall (paper Fig. 6, get_vr.py)
def fig_freefall(base):
    files = snapshots(os.path.join(base, RUNS['freefall']))
    if not files: print('  freefall: no run'); return
    r = np.full((len(files), 32), np.nan); vr = np.zeros_like(r)
    for i, f in enumerate(files):
        h = h5py.File(f); x = h['PartType0/Coordinates'][:]; v = h['PartType0/Velocities'][:]; ids = h['PartType0/ParticleIDs'][:].astype(int) - 1
        rp = np.linalg.norm(x, axis=1); r[i, ids] = rp; vr[i, ids] = -np.sum(x*v, axis=1)/rp
    fig = plt.figure(figsize=[7, 5]); errs = []
    for i in range(31, -1, -1):
        ok = ~np.isnan(r[:, i])
        if r[ok, i].min() > 3: continue
        vex = (1 - 2./r[ok, i])/np.sqrt(1 - 2./r[0, i])*np.sqrt(np.maximum(2.*(1./r[ok, i] - 1./r[0, i]), 0))
        errs.append(np.sqrt(np.mean((vex - vr[ok, i])**2))/np.amax(vex))
        plt.plot(r[ok, i], vex, lw=3, c='tab:orange'); plt.plot(r[ok, i], vr[ok, i], c='k', lw=1)
    print('  free fall: rms error of |v_r| relative to its maximum: median %.1e, max %.1e (%d particles)' % (np.median(errs), np.max(errs), len(errs)))
    plt.xlabel(r'$r$'); plt.ylabel(r'$|v_r|$'); save(fig, 'fig_freefall')

# ------------------------------------------------------------------ Bondi accretion (paper Fig. 12, bondi/solution.py)
def fig_bondi(base):
    rows = [(h, os.path.join(base, run), snap, coords) for h, (run, snap, coords) in RUNS['bondi'].items() if os.path.exists(os.path.join(base, run, 'output', 'snapshot_%03d.hdf5' % (snap//10)))]
    if not rows: print('  bondi: no runs'); return
    sol = np.loadtxt(os.path.join(rows[0][1], 'bondi_sol.txt'), unpack=True)  # r, D, rho0, v_r (Schwarzschild), u
    M = 1.0; vr_ks = sol[3]/(1 + 2*M*sol[3]/(sol[0] - 2*M))               # dr/dt in Kerr-Schild time
    ylim = [[-3, -1], [-0.4, 0.], [-1.6, -0.1]]; ylbl = [r'$\log\rho$', r'$v_r$', r'$\log u$']
    fig, ax = plt.subplots(len(rows), 3, figsize=[13.5, 3.2*len(rows)], sharex=True, squeeze=False); plt.subplots_adjust(wspace=0.32, hspace=0.07)
    for j, (h, run, snap, coords) in enumerate(rows):
        f = h5py.File(os.path.join(run, 'output', 'snapshot_%03d.hdf5' % (snap//10))); g = f['PartType0']; t = float(f['Header'].attrs['Time'])
        x = g['Coordinates'][:]; r = np.linalg.norm(x, axis=1); m = g['Masses'][:]; ok = m > 0
        vr = np.sum(g['Velocities'][:]*x, axis=1)/r; rho = g['RestMassDensity'][:]; u = g['InternalEnergy'][:]
        for k, (y, yth) in enumerate(((np.log10(rho), np.log10(sol[2])), (vr, vr_ks if coords == 'ks' else sol[3]), (np.log10(u), np.log10(sol[4])))):
            ax[j, k].plot(r[ok], y[ok], '.', color=CLR.get(h, 'tab:green'), markersize=2, label=h if k == 0 else None); ax[j, k].plot(sol[0], yth, 'k--', lw=2)
            ax[j, k].set_ylim(ylim[k]); ax[j, k].set_xlim(2.1, 30.1); ax[j, k].yaxis.set_major_locator(plt.MaxNLocator(4))
            if j == 0: ax[j, k].set_ylabel(ylbl[k])
        ax[j, 0].legend(fontsize=12); ax[j, 1].text(0.82, 0.88, r'$t=%g$' % t, fontsize=14, transform=ax[j, 1].transAxes)
        bins = np.logspace(np.log10(4.2), np.log10(30), 30); rc = 0.5*(bins[1:] + bins[:-1]); med = binned_statistic(r[ok], rho[ok], 'median', bins=bins)[0]
        print('  bondi %s t=%g: median |log rho - Michel| in 4.2<r<30 = %.3f dex' % (h, t, np.nanmedian(np.abs(np.log10(med) - np.interp(rc, sol[0], np.log10(sol[2]))))))
    for k in range(3): ax[-1, k].set_xlabel(r'$r$')
    save(fig, 'fig_bondi')

# ------------------------------------------------------------------ TOV star (paper Figs 13-14, TOV/rho0_plot.py and profile_plot.py)
def fig_tov(base):
    rhoc0 = 0.129285; t_dyn = np.sqrt(1./rhoc0)
    runs = {h: os.path.join(base, run) for h, run in RUNS['tov'].items() if snapshots(os.path.join(base, run))}
    if not runs: print('  tov: no runs'); return
    fig = plt.figure(figsize=[7, 5])
    for h, run in runs.items():
        files = snapshots(run); t = []; rc = []
        for f in files:
            hh = h5py.File(f); x = hh['PartType0/Coordinates'][:]; r2 = np.sum(x*x, axis=1); t.append(float(hh['Header'].attrs['Time'])/t_dyn); rc.append(np.mean(hh['PartType0/RestMassDensity'][:][r2 < 1e-2]))
        plt.plot(t, np.array(rc)/rhoc0, label=h, color=CLR[h], ls='-' if h == 'MFM' else '--')
        print('  tov %s: rho_c/rho_c0 in [%.4f, %.4f] up to t = %.1f t_dyn' % (h, min(rc)/rhoc0, max(rc)/rhoc0, t[-1]))
    plt.xlabel(r'$t/t_{\rm dyn}$'); plt.ylabel(r'$\rho_{\rm c}/\rho_{\rm c,0}$'); plt.legend(fontsize=14); plt.ylim(0.9, 1.1); save(fig, 'fig_tov_rhoc')
    # density profiles at three times against the TOV solution (the paper shows images + profiles; here profiles only)
    for h, run in runs.items():
        tab = np.loadtxt(os.path.join(run, 'metric_table.dat'), unpack=True, comments='#'); files = snapshots(run); picks = [files[0], files[len(files)//2], files[-1]]
        fig, ax = plt.subplots(1, 3, figsize=[12, 3.6], sharey=True); plt.subplots_adjust(wspace=0.03)
        for j, f in enumerate(picks):
            hh = h5py.File(f); x = hh['PartType0/Coordinates'][:]; r = np.linalg.norm(x, axis=1); tt = float(hh['Header'].attrs['Time'])/t_dyn
            ax[j].plot(r, hh['PartType0/RestMassDensity'][:], '.', color='tab:orange', markersize=2); ax[j].plot(tab[0], tab[3], 'k', lw=1, ls='--')
            ax[j].set_xlim(-0.1, 2.0); ax[j].set_ylim(-0.01, 0.15); ax[j].text(1, 0.135, r'$t=%5.2f\,t_{\rm dyn}$' % tt, fontsize=14); ax[j].set_xlabel(r'$r$')
        ax[0].set_ylabel(r'$\rho$'); save(fig, 'fig_tov_profiles_' + h)

# ------------------------------------------------------------------ Newtonian polytrope: the positivity-limiter range (not in the paper)
def fig_polytrope(base):
    runs = {h: os.path.join(base, run) for h, run in RUNS['polytrope'].items() if snapshots(os.path.join(base, run))}
    if not runs: print('  polytrope: no runs'); return
    fig, ax = plt.subplots(1, 2, figsize=[12, 4.5]); cl = {'range factor 1.0': 'tab:red', 'range factor 0.5': 'tab:blue', 'range factor 0.25': 'tab:green', 'positivity step off': 'tab:gray'}
    for h, run in runs.items():
        r_t, rho_t = np.loadtxt(os.path.join(run, 'profile.txt'), unpack=True)[:2]; tdyn = 1/np.sqrt(rho_t[0]); files = snapshots(run)
        bins = np.linspace(0, 1.0, 26); rc = 0.5*(bins[1:]+bins[:-1]); rho_th = np.interp(rc, r_t, rho_t); t = []; rc_t = []
        for f in files:
            hh = h5py.File(f); x = hh['PartType0/Coordinates'][:]; r = np.linalg.norm(x, axis=1); rho = hh['PartType0/Density'][:]; t.append(float(hh['Header'].attrs['Time'])/tdyn); rc_t.append(np.median(rho[r < 0.1])/rho_t[0])
        hh = h5py.File(files[-1]); x = hh['PartType0/Coordinates'][:]; r = np.linalg.norm(x, axis=1); rho = hh['PartType0/Density'][:]
        ax[0].plot(rc, binned_statistic(r, rho, 'median', bins=bins)[0]/rho_th, '-', color=cl[h], label='%s, t=%.1f t$_{dyn}$' % (h, t[-1])); ax[1].plot(t, rc_t, '-', color=cl[h], label=h)
        print('  polytrope %s: t_end %.1f t_dyn, rho/rho_LE at r=0.6-0.8 R: %s' % (h, t[-1], np.round(binned_statistic(r, rho, 'median', bins=bins)[0][15:20]/rho_th[15:20], 2)))
    ax[0].axhline(1, color='k', ls='--', lw=1); ax[0].set_xlabel(r'$r/R$'); ax[0].set_ylabel(r'$\rho/\rho_{\rm Lane-Emden}$'); ax[0].set_ylim(0.5, 1.5); ax[0].legend(fontsize=10); ax[0].set_title('n=1.5 polytrope, positivity-limiter range', fontsize=11)
    ax[1].set_xlabel(r'$t/t_{\rm dyn}$'); ax[1].set_ylabel(r'$\rho_{\rm c}/\rho_{\rm c,0}$'); ax[1].legend(fontsize=10)
    save(fig, 'fig_polytrope')
    # density profiles at three times against Lane-Emden, as the TOV profile figure (rows: the two range factors)
    fig, ax = plt.subplots(len(runs), 3, figsize=[12, 3.6*len(runs)], sharey=True, squeeze=False); plt.subplots_adjust(wspace=0.03, hspace=0.25)
    for j, (h, run) in enumerate(runs.items()):
        r_t, rho_t = np.loadtxt(os.path.join(run, 'profile.txt'), unpack=True)[:2]; tdyn = 1/np.sqrt(rho_t[0]); files = snapshots(run); picks = [files[0], files[len(files)//2], files[-1]]
        for k, f in enumerate(picks):
            hh = h5py.File(f); x = hh['PartType0/Coordinates'][:]; r = np.linalg.norm(x, axis=1); tt = float(hh['Header'].attrs['Time'])/tdyn
            ax[j, k].plot(r, hh['PartType0/Density'][:], '.', color=cl[h], markersize=2); ax[j, k].plot(r_t, rho_t, 'k', lw=1, ls='--')
            ax[j, k].set_xlim(-0.05, 1.3); ax[j, k].set_ylim(-0.05, 1.1*rho_t[0]); ax[j, k].text(0.72, 0.95*rho_t[0], r'$t=%5.2f\,t_{\rm dyn}$' % tt, fontsize=13); ax[j, k].xaxis.set_major_locator(plt.FixedLocator([0, 0.25, 0.5, 0.75, 1.0]))
            if j == len(runs)-1: ax[j, k].set_xlabel(r'$r$')
        ax[j, 0].set_ylabel(r'$\rho$'); ax[j, 0].text(0.03, 0.08*rho_t[0], h, fontsize=12, color=cl[h])
    save(fig, 'fig_polytrope_profiles')

if __name__ == '__main__':
    if len(sys.argv) < 2: print(__doc__); sys.exit(1)
    base = sys.argv[1]; which = sys.argv[2:] or ['shocktubes', 'orbits', 'epicyclic', 'apsidal', 'freefall', 'bondi', 'tov', 'polytrope']
    for w in which: print(w); globals()['fig_' + w](base)
