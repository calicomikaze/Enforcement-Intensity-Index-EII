"""
================================================================================
ENFORCEMENT INTENSITY INDEX (EII)
================================================================================
Bayesian Zero-Inflated Negative Binomial Estimation
of Enforcement Intensity in Strategically Produced Data

Sean Rogers | University of South Carolina | 2026

INPUT:   EII_Panel_Definitive.csv (942 country-years, 44 states, 1945-2025)
OUTPUT:  bayesian_results_all_specs.csv, weighting_sensitivity.csv, fig1-fig4
RUNTIME: ~3 min (Windows, single core)

The EII is constructed from ten institutional primary sources across four
layers. This script takes the panel as given and estimates the model.
See codebook for construction details.

REGIME VARIABLES:
  Primary:    V-Dem Electoral Democracy Index (continuous, 0-1)
  Additional: V-Dem Free and Fair Elections Index (continuous, 0-1)
  Robustness: Polity5 (integer, -10 to +10)

WHY BAYESIAN?
  The institutions producing enforcement data have survival incentives to
  misreport. MLE assumes the DGP is honest. When the DGP is strategically
  corrupted, the likelihood surface is distorted and frequentist estimation
  may fit to artifacts. The Bayesian framework regularizes via priors,
  produces full posterior distributions, and enables explicit prior
  sensitivity analysis.

WHY ZINB?
  The DV has two sources of zeros:
    Sampling zeros:   the count process generates zero (low intensity)
    Structural zeros: the state is shielded from enforcement (sanctuary)
  These are qualitatively different. The ZINB separates them with two
  equations: the count equation models intensity; the zero-inflation
  equation models sanctuary probability.

MEASUREMENT SENSITIVITY:
  The EII triangulates ten institutional sources on the DV side because
  each source lies in a different direction for different institutional
  reasons. The fiscal pressure IV uses a single source (CBO). The
  interaction (friction x fiscal) is not supported in Bayesian estimation
  while the fiscal main effect is. This asymmetry is consistent with the
  paper's core argument: triangulation recovers signal from corrupted
  data; single-source measurement does not.
================================================================================
"""

import subprocess, sys, os
import numpy as np
import pandas as pd
import matplotlib
import matplotlib.pyplot as plt
import warnings

# ---- Install if needed ----
for pkg in ["pymc", "arviz", "statsmodels"]:
    try:
        __import__(pkg)
    except ImportError:
        subprocess.check_call([sys.executable, "-m", "pip", "install", "-q", pkg])

import pymc as pm
import arviz as az
import statsmodels.api as sm
from statsmodels.discrete.count_model import ZeroInflatedNegativeBinomialP
from scipy import stats

matplotlib.rcParams['font.family'] = 'serif'
matplotlib.rcParams['font.size'] = 11
warnings.filterwarnings('ignore')

SEED = 42
np.random.seed(SEED)


def z_score(s):
    """Standardize to mean=0, sd=1."""
    return (s - s.mean()) / s.std()


def hdi(samples, prob=0.95):
    """Highest density interval via ArviZ 1.x API."""
    return az.hdi(samples, prob=prob)


def run_freq_zinb(y, X_count, X_inflate, count_names, inflate_names, label):
    """Run frequentist ZINB and return results dict."""
    print(f"\n{'='*70}")
    print(f"FREQUENTIST ZINB: {label}")
    print(f"{'='*70}")

    try:
        result = ZeroInflatedNegativeBinomialP(
            y, X_count, exog_infl=X_inflate, p=2
        ).fit(method='bfgs', maxiter=5000, disp=0)
    except Exception:
        try:
            result = ZeroInflatedNegativeBinomialP(
                y, X_count, exog_infl=X_inflate, p=2
            ).fit(method='nm', maxiter=10000, disp=0)
        except Exception as e:
            print(f"  FAILED: {e}")
            return None

    has_se = not np.any(np.isnan(result.bse))
    converged = result.mle_retvals.get('converged', False)
    alpha_idx = len(count_names)

    # Validate parameter count
    expected = len(count_names) + 1 + len(inflate_names)
    if len(result.params) != expected:
        print(f"  Parameter count mismatch: {len(result.params)} vs {expected}")
        return None

    print(f"Converged: {converged} | Valid SEs: {has_se}")
    print(f"Log-likelihood: {result.llf:.2f} | AIC: {result.aic:.2f}")

    print(f"\nCOUNT EQUATION:")
    print(f"{'Variable':<20} {'Coef':>9} {'SE':>9} {'z':>8} {'p':>9}")
    print("-" * 60)
    for i, nm in enumerate(count_names):
        c = result.params[i]
        if has_se:
            sg = "***" if result.pvalues[i] < .001 else "**" if result.pvalues[i] < .01 else "*" if result.pvalues[i] < .05 else ""
            print(f"{nm:<20} {c:>9.4f} {result.bse[i]:>9.4f} {result.tvalues[i]:>8.3f} {result.pvalues[i]:>9.5f} {sg}")
        else:
            print(f"{nm:<20} {c:>9.4f}    [SEs unavailable]")

    print(f"\nAlpha: {result.params[alpha_idx]:.4f}")

    print(f"\nZERO-INFLATION EQUATION:")
    print(f"{'Variable':<20} {'Coef':>9} {'SE':>9} {'z':>8} {'p':>9}")
    print("-" * 60)
    for i, nm in enumerate(inflate_names):
        idx = alpha_idx + 1 + i
        c = result.params[idx]
        if has_se:
            sg = "***" if result.pvalues[idx] < .001 else "**" if result.pvalues[idx] < .01 else "*" if result.pvalues[idx] < .05 else ""
            print(f"{nm:<20} {c:>9.4f} {result.bse[idx]:>9.4f} {result.tvalues[idx]:>8.3f} {result.pvalues[idx]:>9.5f} {sg}")
        else:
            print(f"{nm:<20} {c:>9.4f}    [SEs unavailable]")

    return {'result': result, 'has_se': has_se, 'converged': converged}


def run_bayesian_zinb(label, Xc, Xz, y_obs, spec, count_names, inflate_names):
    """
    Estimate a Bayesian ZINB.

    Model:
      Count:   log(mu) = b0 + Xc @ beta
      Inflate: logit(psi) = g0 + Xz @ gamma
      DV ~ ZI-NegBin(psi, mu, alpha)

    psi = P(count process) per PyMC convention (verified: test_psi_convention.py).
    Inflate-equation coefficients are therefore interpreted on the ENGAGEMENT side:
    positive gamma = higher P(engagement) = LESS sanctuary.
    Positive beta = higher enforcement given non-sanctuary.
    """
    n_c, n_z = Xc.shape[1], Xz.shape[1]
    print(f"\n--- {label}: {spec.get('description', '')} ---")

    with pm.Model():
        b0 = pm.Normal('beta_intercept', mu=0, sigma=2.5)
        b_main = pm.Normal('beta_main', mu=0, sigma=spec['beta_sigma'], shape=n_c - 1)
        b_int = pm.Normal('beta_interaction', mu=spec['interaction_mu'], sigma=spec['interaction_sigma'])
        b_full = pm.math.concatenate([b_main, b_int.reshape((1,))])

        g0 = pm.Normal('gamma_intercept', mu=spec['gamma0_mu'], sigma=spec['gamma0_sigma'])
        g = pm.Normal('gamma', mu=0, sigma=spec['gamma_sigma'], shape=n_z)

        alpha = pm.Exponential('alpha', lam=1)

        psi = pm.math.sigmoid(g0 + pm.math.dot(Xz, g))
        mu = pm.math.exp(b0 + pm.math.dot(Xc, b_full))

        pm.ZeroInflatedNegativeBinomial('y_obs', psi=psi, mu=mu, alpha=alpha, observed=y_obs)

        trace = pm.sample(
            draws=1000, tune=1000, chains=4, cores=1,
            target_accept=0.90, random_seed=SEED,
            return_inferencedata=True, progressbar=True,
            idata_kwargs={"log_likelihood": True}
        )

    # Extract
    beta_m = trace.posterior['beta_main'].values.reshape(-1, n_c - 1)
    beta_i = trace.posterior['beta_interaction'].values.flatten()
    gamma_s = trace.posterior['gamma'].values.reshape(-1, n_z)
    beta_all = np.column_stack([beta_m, beta_i])

    # Diagnostics
    diag = az.summary(trace, var_names=['beta_main', 'beta_interaction', 'gamma', 'gamma_intercept', 'alpha'])
    rhat_max = diag['r_hat'].max()
    ess_min = min(diag['ess_bulk'].min(), diag['ess_tail'].min())
    div = int(trace.sample_stats['diverging'].sum().values)

    print(f"  R-hat max: {rhat_max:.4f} | ESS min: {ess_min:.0f} | Divergences: {div}")
    if rhat_max > 1.01:
        print(f"  *** WARNING: R-hat > 1.01 ***")
    print(f"  PLR: mean={beta_i.mean():.4f}, P(<0)={np.mean(beta_i<0):.4f}")

    return {
        'trace': trace, 'beta_all': beta_all, 'gamma': gamma_s,
        'plr': beta_i, 'rhat': rhat_max, 'ess': ess_min, 'div': div,
    }


def print_bayesian_table(r, count_names, inflate_names, title):
    """Print posterior summary with proper HDI."""
    print(f"\n{'='*70}")
    print(title)
    print(f"{'='*70}")

    print(f"\nCOUNT EQUATION (Enforcement Intensity):")
    print(f"{'Variable':<20} {'Mean':>8} {'SD':>8} {'HDI Lo':>8} {'HDI Hi':>8} {'P(dir.)':>8}")
    print("-" * 70)
    for i, name in enumerate(count_names):
        s = r['beta_all'][:, i]
        h = hdi(s)
        pd_val = max(np.mean(s > 0), np.mean(s < 0))
        sg = "***" if pd_val > 0.999 else "**" if pd_val > 0.99 else "*" if pd_val > 0.95 else ""
        print(f"{name:<20} {s.mean():>8.3f} {s.std():>8.3f} {h[0]:>8.3f} {h[1]:>8.3f} {pd_val:>8.3f} {sg}")

    print(f"\nZERO-INFLATION EQUATION (Sanctuary):")
    print(f"{'Variable':<20} {'Mean':>8} {'SD':>8} {'HDI Lo':>8} {'HDI Hi':>8} {'P(dir.)':>8}")
    print("-" * 70)
    g0 = r['trace'].posterior['gamma_intercept'].values.flatten()
    h0 = hdi(g0)
    print(f"{'intercept':<20} {g0.mean():>8.3f} {g0.std():>8.3f} {h0[0]:>8.3f} {h0[1]:>8.3f} {max(np.mean(g0>0),np.mean(g0<0)):>8.3f}")
    for i, name in enumerate(inflate_names):
        s = r['gamma'][:, i]
        h = hdi(s)
        pd_val = max(np.mean(s > 0), np.mean(s < 0))
        sg = "***" if pd_val > 0.999 else "**" if pd_val > 0.99 else "*" if pd_val > 0.95 else ""
        print(f"{name:<20} {s.mean():>8.3f} {s.std():>8.3f} {h[0]:>8.3f} {h[1]:>8.3f} {pd_val:>8.3f} {sg}")

    a = r['trace'].posterior['alpha'].values.flatten()
    ha = hdi(a)
    print(f"\nOverdispersion (alpha): {a.mean():.3f} [{ha[0]:.3f}, {ha[1]:.3f}]")


def save_results_csv(all_results, count_names, inflate_names, filename):
    """Save all posteriors to CSV."""
    rows = []
    for label, r in all_results.items():
        for i, name in enumerate(count_names):
            s = r['beta_all'][:, i]
            h = hdi(s)
            rows.append({
                'spec': label, 'equation': 'count', 'variable': name,
                'mean': s.mean(), 'sd': s.std(),
                'hdi_lo': h[0], 'hdi_hi': h[1],
                'p_dir': max(np.mean(s > 0), np.mean(s < 0)),
                'p_positive': float(np.mean(s > 0)),
                'p_negative': float(np.mean(s < 0)),
            })
        for i, name in enumerate(inflate_names):
            s = r['gamma'][:, i]
            h = hdi(s)
            rows.append({
                'spec': label, 'equation': 'inflate', 'variable': name,
                'mean': s.mean(), 'sd': s.std(),
                'hdi_lo': h[0], 'hdi_hi': h[1],
                'p_dir': max(np.mean(s > 0), np.mean(s < 0)),
                'p_positive': float(np.mean(s > 0)),
                'p_negative': float(np.mean(s < 0)),
            })
    pd.DataFrame(rows).to_csv(filename, index=False)
    print(f"Saved: {filename}")


# ======================================================================
# MAIN
# ======================================================================
def main():
    print("=" * 70)
    print("EII REPLICATION PIPELINE — v5")
    print(f"PyMC {pm.__version__} | ArviZ {az.__version__}")
    print("=" * 70)

    # ------------------------------------------------------------------
    # DATA
    # ------------------------------------------------------------------
    df = pd.read_csv('EII_Panel_Definitive.csv')
    for col in ['eiii_score', 'polity', 'friction', 'degradation',
                'resource_density', 'fiscal_pressure', 'gdp_pc_log',
                'trend', 'vdem_electoral', 'vdem_elections']:
        df[col] = pd.to_numeric(df[col], errors='coerce')
    df['eiii_score'] = df['eiii_score'].fillna(0).astype(int)
    df['degradation'] = df['degradation'].fillna(0).astype(int)

    # Estimation samples
    vdem_ivs = ['vdem_electoral', 'friction', 'resource_density', 'gdp_pc_log', 'fiscal_pressure']
    df_v = df.dropna(subset=vdem_ivs).copy()
    df_e = df.dropna(subset=['vdem_elections', 'friction', 'resource_density', 'gdp_pc_log', 'fiscal_pressure']).copy()
    df_p = df.dropna(subset=['polity', 'friction', 'resource_density', 'gdp_pc_log', 'fiscal_pressure']).copy()

    # Z-score everything
    for d in [df_v, df_e, df_p]:
        for col in ['friction', 'resource_density', 'fiscal_pressure', 'gdp_pc_log', 'trend']:
            d[f'{col}_z'] = z_score(d[col])
    df_v['regime_z'] = z_score(df_v['vdem_electoral'])
    df_e['regime_z'] = z_score(df_e['vdem_elections'])
    df_p['regime_z'] = z_score(df_p['polity'])
    for d in [df_v, df_e, df_p]:
        d['fric_x_fiscal'] = d['friction_z'] * d['fiscal_pressure_z']

    N = len(df_v)
    zero_pct = 100 * (df_v['eiii_score'] == 0).mean()
    od = df_v['eiii_score'].var() / max(df_v['eiii_score'].mean(), 1e-10)

    print(f"\nPanel: {len(df)} total | Estimation sample: {N} complete cases, {df_v['country'].nunique()} states")
    print(f"DV: {zero_pct:.1f}% zeros | Var/mean: {od:.1f}")

    # ------------------------------------------------------------------
    # VARIABLE NAMES (shared across all specs)
    # ------------------------------------------------------------------
    count_names_freq = ['const', 'resource', 'friction', 'degradation',
                        'fiscal', 'trend', 'regime', 'gdp_pc', 'friction_x_fiscal']
    inflate_names_freq = ['const', 'regime', 'friction', 'resource', 'degradation', 'trend']
    count_names = ['resource', 'friction', 'degradation', 'fiscal',
                   'trend', 'regime', 'gdp_pc', 'friction_x_fiscal']
    inflate_names = ['regime', 'friction', 'resource', 'degradation', 'trend']

    # ------------------------------------------------------------------
    # FREQUENTIST BASELINE — V-DEM
    # ------------------------------------------------------------------
    y_v = df_v['eiii_score'].values
    Xc_freq = sm.add_constant(df_v[['resource_density_z', 'friction_z', 'degradation',
                                     'fiscal_pressure_z', 'trend_z', 'regime_z',
                                     'gdp_pc_log_z', 'fric_x_fiscal']].values)
    Xi_freq = sm.add_constant(df_v[['regime_z', 'friction_z', 'resource_density_z',
                                     'degradation', 'trend_z']].values)
    freq_vdem = run_freq_zinb(y_v, Xc_freq, Xi_freq, count_names_freq, inflate_names_freq, "V-DEM ELECTORAL (PRIMARY)")

    if freq_vdem:
        fx_i = count_names_freq.index('friction_x_fiscal')
        print(f"\n>>> Frequentist PLR: b = {freq_vdem['result'].params[fx_i]:.4f}", end="")
        if freq_vdem['has_se']:
            print(f" (p = {freq_vdem['result'].pvalues[fx_i]:.6f})")
        else:
            print()

    # ------------------------------------------------------------------
    # FREQUENTIST BASELINE — POLITY (ROBUSTNESS)
    # ------------------------------------------------------------------
    Xc_pol = sm.add_constant(df_p[['resource_density_z', 'friction_z', 'degradation',
                                    'fiscal_pressure_z', 'trend_z', 'regime_z',
                                    'gdp_pc_log_z', 'fric_x_fiscal']].values)
    Xi_pol = sm.add_constant(df_p[['regime_z', 'friction_z', 'resource_density_z',
                                    'degradation', 'trend_z']].values)
    freq_pol = run_freq_zinb(df_p['eiii_score'].values, Xc_pol, Xi_pol,
                             count_names_freq, inflate_names_freq, "POLITY5 (ROBUSTNESS)")

    # ------------------------------------------------------------------
    # BAYESIAN — V-DEM PRIMARY — FOUR PRIOR SPECS
    # ------------------------------------------------------------------
    print(f"\n{'='*70}")
    print("BAYESIAN ZINB — V-DEM PRIMARY — FOUR PRIOR SPECIFICATIONS")
    print(f"{'='*70}")

    Xc = df_v[['resource_density_z', 'friction_z', 'degradation',
                'fiscal_pressure_z', 'trend_z', 'regime_z',
                'gdp_pc_log_z', 'fric_x_fiscal']].values.astype('float64')
    Xz = df_v[['regime_z', 'friction_z', 'resource_density_z',
                'degradation', 'trend_z']].values.astype('float64')
    y_obs = df_v['eiii_score'].values.astype(int)

    ZI_BASE = np.log(zero_pct / 100 / (1 - zero_pct / 100))

    specs = {
        'Spec1_Diffuse': {
            'description': 'Weakly informative (Gelman 2008)',
            'beta_sigma': 2.5, 'gamma_sigma': 2.5,
            'gamma0_mu': ZI_BASE, 'gamma0_sigma': 2.5,
            'interaction_mu': 0, 'interaction_sigma': 2.5,
        },
        'Spec2_Sequential': {
            'description': 'Sequential — proxy posteriors as priors',
            'beta_sigma': 2.5, 'gamma_sigma': 2.5,
            'gamma0_mu': ZI_BASE, 'gamma0_sigma': 2.5,
            'interaction_mu': -0.1, 'interaction_sigma': 1.0,
        },
        'Spec3_Skeptical': {
            'description': 'Skeptical — shrinks interaction toward zero',
            'beta_sigma': 1.5, 'gamma_sigma': 1.5,
            'gamma0_mu': ZI_BASE, 'gamma0_sigma': 2.5,
            'interaction_mu': 0, 'interaction_sigma': 0.5,
        },
        'Spec4_Restrictive': {
            'description': 'Restrictive — N(0,1) throughout',
            'beta_sigma': 1.0, 'gamma_sigma': 1.0,
            'gamma0_mu': ZI_BASE, 'gamma0_sigma': 2.0,
            'interaction_mu': 0, 'interaction_sigma': 1.0,
        },
    }

    results = {}
    for label, spec in specs.items():
        results[label] = run_bayesian_zinb(label, Xc, Xz, y_obs, spec, count_names, inflate_names)

    r1 = results['Spec1_Diffuse']

    # ------------------------------------------------------------------
    # BAYESIAN — V-DEM ELECTIONS
    # ------------------------------------------------------------------
    print(f"\n{'='*70}")
    print("ADDITIONAL: V-DEM FREE AND FAIR ELECTIONS")
    print(f"{'='*70}")

    Xc_e = df_e[['resource_density_z', 'friction_z', 'degradation',
                  'fiscal_pressure_z', 'trend_z', 'regime_z',
                  'gdp_pc_log_z', 'fric_x_fiscal']].values.astype('float64')
    Xz_e = df_e[['regime_z', 'friction_z', 'resource_density_z',
                  'degradation', 'trend_z']].values.astype('float64')
    y_e = df_e['eiii_score'].values.astype(int)

    r_elec = run_bayesian_zinb('VDem_Elections', Xc_e, Xz_e, y_e,
                                specs['Spec1_Diffuse'], count_names, inflate_names)

    # ------------------------------------------------------------------
    # RESULTS TABLES
    # ------------------------------------------------------------------
    print_bayesian_table(r1, count_names, inflate_names,
                         "PRIMARY RESULTS — V-DEM ELECTORAL (DIFFUSE PRIORS)")
    print_bayesian_table(r_elec, count_names, inflate_names,
                         "ADDITIONAL — V-DEM FREE AND FAIR ELECTIONS")

    # ------------------------------------------------------------------
    # PRIOR SENSITIVITY
    # ------------------------------------------------------------------
    print(f"\n{'='*70}")
    print("PRIOR SENSITIVITY — FRICTION x FISCAL INTERACTION")
    print(f"{'='*70}")
    print(f"\n{'Spec':<22} {'Mean':>8} {'SD':>7} {'HDI Lo':>7} {'HDI Hi':>7} {'P(<0)':>7}")
    print("-" * 55)
    for label, r in results.items():
        plr = r['plr']
        h = hdi(plr)
        short = label.split('_', 1)[1]
        print(f"{short:<22} {plr.mean():>8.4f} {plr.std():>7.4f} {h[0]:>7.3f} {h[1]:>7.3f} {np.mean(plr<0):>7.4f}")

    print(f"\n{'Spec':<22} {'Friction':>9} {'Fiscal':>9} {'Trend':>9} {'Regime':>9}")
    print("-" * 60)
    for label, r in results.items():
        short = label.split('_', 1)[1]
        vals = [max(np.mean(r['beta_all'][:, j] > 0), np.mean(r['beta_all'][:, j] < 0)) for j in [1, 3, 4, 5]]
        print(f"{short:<22} {vals[0]:>9.3f} {vals[1]:>9.3f} {vals[2]:>9.3f} {vals[3]:>9.3f}")

    # ------------------------------------------------------------------
    # CONVERGENT VALIDATION
    # ------------------------------------------------------------------
    print(f"\n{'='*70}")
    print("CONVERGENT VALIDATION: BAYESIAN vs. FREQUENTIST")
    print(f"{'='*70}")
    if freq_vdem:
        fr = freq_vdem['result']
        print(f"\n{'Variable':<20} {'Freq b':>9} {'Freq p':>9} {'Bayes':>9} {'P(dir.)':>8} {'Sign':>6}")
        print("-" * 65)
        for i, name in enumerate(count_names):
            fc = fr.params[i + 1]
            fp = fr.pvalues[i + 1] if freq_vdem['has_se'] else float('nan')
            bm = r1['beta_all'][:, i].mean()
            bp = max(np.mean(r1['beta_all'][:, i] > 0), np.mean(r1['beta_all'][:, i] < 0))
            ag = "Y" if np.sign(fc) == np.sign(bm) else "N"
            print(f"{name:<20} {fc:>9.4f} {fp:>9.5f} {bm:>9.4f} {bp:>8.4f} {ag:>6}")

    # ------------------------------------------------------------------
    # REGIME VARIABLE COMPARISON
    # ------------------------------------------------------------------
    print(f"\n{'='*70}")
    print("REGIME VARIABLE COMPARISON")
    print(f"{'='*70}")
    print(f"\n{'Variable':<20} {'V-Dem Elect.':>12} {'V-Dem Elec.':>12} {'Polity5':>12}")
    print("-" * 60)
    for i, name in enumerate(count_names):
        vm = r1['beta_all'][:, i].mean()
        em = r_elec['beta_all'][:, i].mean()
        pm_c = freq_pol['result'].params[i + 1] if freq_pol and i + 1 < len(freq_pol['result'].params) else float('nan')
        print(f"{name:<20} {vm:>12.3f} {em:>12.3f} {pm_c:>12.3f}")

    # ------------------------------------------------------------------
    # MCMC DIAGNOSTICS
    # ------------------------------------------------------------------
    print(f"\n{'='*70}")
    print("MCMC DIAGNOSTICS")
    print(f"{'='*70}")
    print(f"\n{'Spec':<22} {'R-hat':>8} {'ESS':>8} {'Diverg.':>8}")
    print("-" * 48)
    for label, r in list(results.items()) + [('VDem_Elections', r_elec)]:
        short = label.split('_', 1)[1] if '_' in label else label
        print(f"{short:<22} {r['rhat']:>8.4f} {r['ess']:>8.0f} {r['div']:>8}")

    # ------------------------------------------------------------------
    # DV WEIGHTING SENSITIVITY
    # ------------------------------------------------------------------
    print(f"\n{'='*70}")
    print("DV WEIGHTING SENSITIVITY — FREQUENTIST")
    print(f"{'='*70}")

    dv_specs = {
        'Weighted (primary)': 'eiii_score',
        'Equal weight': 'eiii_equal',
        'Raw count (capped)': 'eiii_rawcount',
        'Financial only': 'eiii_financial',
    }

    print(f"\n{'DV':<25} {'Friction':>9} {'Fiscal':>9} {'Trend':>9} {'FxF':>9} {'FxF p':>10}")
    print("-" * 75)
    wt_rows = []
    for spec_name, dv_col in dv_specs.items():
        y_wt = df_v[dv_col].values
        try:
            r_wt = ZeroInflatedNegativeBinomialP(
                y_wt, Xc_freq, exog_infl=Xi_freq, p=2
            ).fit(method='bfgs', maxiter=5000, disp=0)
            has_se_wt = not np.any(np.isnan(r_wt.bse))
            fc = r_wt.params[2]
            fi = r_wt.params[4]
            tr = r_wt.params[5]
            fx = r_wt.params[8]
            fp = r_wt.pvalues[8] if has_se_wt else float('nan')
            print(f"{spec_name:<25} {fc:>9.3f} {fi:>9.3f} {tr:>9.3f} {fx:>9.3f} {fp:>10.5f}")
            wt_rows.append({'dv_spec': spec_name, 'friction': fc, 'fiscal': fi, 'trend': tr, 'fxf': fx, 'fxf_p': fp})
        except Exception as e:
            print(f"{spec_name:<25} FAILED: {str(e)[:40]}")

    if wt_rows:
        pd.DataFrame(wt_rows).to_csv('weighting_sensitivity.csv', index=False)
        print("Saved: weighting_sensitivity.csv")

    # ------------------------------------------------------------------
    # FIGURES
    # ------------------------------------------------------------------
    print("\n--- Generating figures ---")

    # Fig 1: Count posteriors
    fig1, axes = plt.subplots(2, 4, figsize=(16, 8))
    fig1.suptitle('Count Equation Posteriors — V-Dem Electoral (Diffuse)', fontsize=14, fontweight='bold')
    colors = ['#2c3e50'] * 7 + ['#c0392b']
    for i, name in enumerate(count_names):
        ax = axes[i // 4, i % 4]
        s = r1['beta_all'][:, i]
        ax.hist(s, bins=50, density=True, alpha=0.8, color=colors[i], edgecolor='white', linewidth=0.5)
        ax.axvline(0, color='red', linestyle='--', alpha=0.6, linewidth=1.5)
        h = hdi(s)
        ax.axvspan(h[0], h[1], alpha=0.1, color='blue')
        pd_val = max(np.mean(s > 0), np.mean(s < 0))
        ax.set_title(f'{name}\nP(dir.)={pd_val:.3f}', fontsize=10)
    plt.tight_layout()
    plt.savefig('fig1_count_posteriors.png', dpi=300, bbox_inches='tight')
    plt.close()

    # Fig 2: Sanctuary posteriors
    fig2, axes2 = plt.subplots(1, 5, figsize=(16, 3.5))
    fig2.suptitle('Zero-Inflation (Sanctuary) Posteriors', fontsize=14, fontweight='bold')
    for i, name in enumerate(inflate_names):
        ax = axes2[i]
        s = r1['gamma'][:, i]
        ax.hist(s, bins=50, density=True, alpha=0.8, color='#2c3e50', edgecolor='white', linewidth=0.5)
        ax.axvline(0, color='red', linestyle='--', alpha=0.6)
        pd_val = max(np.mean(s > 0), np.mean(s < 0))
        ax.set_title(f'{name}\nP(dir.)={pd_val:.3f}', fontsize=10)
    plt.tight_layout()
    plt.savefig('fig2_sanctuary_posteriors.png', dpi=300, bbox_inches='tight')
    plt.close()

    # Fig 3: PLR hero
    fig3, ax3 = plt.subplots(figsize=(8, 5))
    plr = r1['plr']
    ax3.hist(plr, bins=60, density=True, alpha=0.85, color='#c0392b', edgecolor='white', linewidth=0.5)
    ax3.axvline(0, color='black', linestyle='--', linewidth=2, label='Zero')
    h_plr = hdi(plr)
    ax3.axvspan(h_plr[0], h_plr[1], alpha=0.15, color='blue', label=f'95% HDI [{h_plr[0]:.3f}, {h_plr[1]:.3f}]')
    ax3.set_title('Friction x Fiscal Pressure Interaction', fontsize=14, fontweight='bold')
    ax3.set_xlabel('Coefficient (standardized)', fontsize=12)
    ax3.set_ylabel('Posterior Density', fontsize=12)
    ax3.legend(fontsize=10)
    plt.tight_layout()
    plt.savefig('fig3_interaction.png', dpi=300, bbox_inches='tight')
    plt.close()

    # Fig 4: Prior sensitivity
    fig4, axes4 = plt.subplots(1, 4, figsize=(16, 4))
    fig4.suptitle('Prior Sensitivity — Friction x Fiscal', fontsize=13, fontweight='bold')
    spec_colors = ['#3498db', '#27ae60', '#e67e22', '#8e44ad']
    for i, (label, r) in enumerate(results.items()):
        ax = axes4[i]
        ps = r['plr']
        ax.hist(ps, bins=50, density=True, alpha=0.8, color=spec_colors[i], edgecolor='white', linewidth=0.5)
        ax.axvline(0, color='black', linestyle='--')
        short = label.split('_', 1)[1]
        ax.set_title(f'{short}\nMean={ps.mean():.3f}\nP(<0)={np.mean(ps<0):.3f}', fontsize=9)
        ax.set_xlabel('Coefficient')
    plt.tight_layout()
    plt.savefig('fig4_prior_sensitivity.png', dpi=300, bbox_inches='tight')
    plt.close()

    print("  Saved: fig1-fig4")

    # ------------------------------------------------------------------
    # SAVE
    # ------------------------------------------------------------------
    all_results = dict(results)
    all_results['VDem_Elections'] = r_elec
    save_results_csv(all_results, count_names, inflate_names, 'bayesian_results_all_specs.csv')
    df_v.to_csv('panel_estimation_sample.csv', index=False)
    print("Saved: panel_estimation_sample.csv")

    print(f"\n{'='*70}")
    print("PIPELINE COMPLETE")
    print(f"{'='*70}")
    print(f"  Panel: {len(df)} total, {N} estimation, {df_v['country'].nunique()} states")
    print(f"  Primary: V-Dem Electoral | Robustness: V-Dem Elections, Polity5")
    print(f"  Bayesian: 4 prior specs + elections spec")
    print(f"  DV sensitivity: 4 weighting schemes")
    print(f"  Seed: {SEED}")


if __name__ == '__main__':
    main()
