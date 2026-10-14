
# -*- coding: utf-8 -*-
"""
@author: SAINATH
"""


"""
Adaptive Piecewise Gaussian Fitting with Expanded Kernel Library
  – Automatic Δflux split in θ∈[0,30]°
  – Candidates: Standard, Super, Skew, Twin, Inv, DoG, GMM3, GMM4,
                Truncated, Folded, Split Gaussians
  – Model selection by AIC
  – Smooth blend ±2° at the join
  – Compute R², RMSE, MAE, SSR, AIC, AICc, BIC, CV‐RMSE, CV‐MAE

Created on Jul 20 2025
"""
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit
from scipy.stats import skewnorm
from sklearn.model_selection import KFold
import warnings


# 1) Load data
df     = pd.read_csv("C:/...../Data_From_Literature/Jeter_et_al_Data_Sun_Shape.csv") # give exact path
theta  = df["0_Angle"].values
flux   = df["0_LCR"].values
N      = len(theta)

# 2) Set split-angle manually or automatically

split_ang = 17.5 #manually
'''
#Detect split‐angle by Δflux in [0°,30°]
mask30    = (theta >= 0) & (theta <= 30)
idx30     = np.where(mask30)[0]
flux30    = flux[idx30]
dflux30   = np.diff(flux30)
split_j   = idx30[np.argmax(np.abs(dflux30)) + 1]
split_ang = theta[split_j]
print(f"Split at θ = {split_ang:.4f}° (auto set)")
'''

print(f"Split at θ = {split_ang:.4f}° (manually set)")


# 3) Partition regions
mask1 = theta <= split_ang
mask2 = theta >  split_ang
x1, y1 = theta[mask1], flux[mask1]
x2, y2 = theta[mask2], flux[mask2]

# 4) Define kernels

def gauss(x, A, mu, s):
    return A * np.exp(-((x-mu)**2)/(2*s**2))

def super_gauss(x, A, mu, s, n):
    return A * np.exp(-np.abs((x-mu)/s)**n)

def skew_gauss(x, A, a, mu, s):
    return A * skewnorm.pdf(x, a=a, loc=mu, scale=s)

def twin_gauss(x, A1,m1,s1, A2,m2,s2):
    return (A1*np.exp(-((x-m1)**2)/(2*s1**2))
          + A2*np.exp(-((x-m2)**2)/(2*s2**2)))

def inv_gauss(x, B, A, mu, s):
    return B + A*(1 - np.exp(-((x-mu)**2)/(2*s**2)))

def dog(x, A1,m1,s1, A2,m2,s2):
    return A1*np.exp(-((x-m1)**2)/(2*s1**2)) \
         - A2*np.exp(-((x-m2)**2)/(2*s2**2))

def gmm3(x, A1,m1,s1, A2,m2,s2, A3,m3,s3):
    return (A1*np.exp(-((x-m1)**2)/(2*s1**2))
          + A2*np.exp(-((x-m2)**2)/(2*s2**2))
          + A3*np.exp(-((x-m3)**2)/(2*s3**2)))

def gmm4(x, A1,m1,s1, A2,m2,s2, A3,m3,s3, A4,m4,s4):
    return gmm3(x, A1,m1,s1, A2,m2,s2, A3,m3,s3) \
         + A4*np.exp(-((x-m4)**2)/(2*s4**2))

def trunc_gauss(x, A, mu, s, r):
    y = A*np.exp(-((x-mu)**2)/(2*s**2))
    return np.where(np.abs(x-mu)<=r*s, y, 0)

def folded_gauss(x, A, mu, s):
    return A * np.exp(-((x-mu)**2)/(2*s**2)) * (x >= 0)

def split_gauss(x, A, mu, sl, sr):
    y = np.empty_like(x)
    left  = x <= mu
    right = ~left
    y[left]  = A*np.exp(-((x[left]-mu)**2)/(2*sl**2))
    y[right] = A*np.exp(-((x[right]-mu)**2)/(2*sr**2))
    return y

candidates = {
    "Gaussian":     (gauss,      3),
    "SuperGauss":   (super_gauss,4),
    "SkewGauss":    (skew_gauss, 4),
    "TwinGauss":    (twin_gauss, 6),
    "InvGauss":     (inv_gauss,  4),
    "DoG":          (dog,        6),
    "GMM3":         (gmm3,       9),
    "GMM4":         (gmm4,      12),
    "TruncGauss":   (trunc_gauss,4),
    "FoldedGauss":  (folded_gauss,3),
   "SplitGauss":   (split_gauss,4),
}

# 5) Goodness‐of‐fit metrics
def metrics(y, yhat, k):
    resid = y - yhat
    ssr   = np.sum(resid**2)
    sst   = np.sum((y - y.mean())**2)
    r2    = 1 - ssr/sst
    rmse  = np.sqrt(np.mean(resid**2))
    mae   = np.mean(np.abs(resid))
    n     = len(y)
    aic   = n*np.log(ssr/n) + 2*k
    aicc  = aic + (2*k*(k+1))/(n-k-1) if n>k+1 else np.nan
    bic   = n*np.log(ssr/n) + k*np.log(n)
    return r2, rmse, mae, ssr, aic, aicc, bic

# 6) Initial guesses
def initial_guess(name, x, y):
    span = np.ptp(x) or 1
    peak = x[np.argmax(y)]
    if name=="Gaussian":
        return [y.max(), peak, span/4]
    if name=="SuperGauss":
        return [y.max(), peak, span/4, 4.0]
    if name=="SkewGauss":
        return [y.max(), 1.0, peak, span/4]
    if name=="TwinGauss":
        return [y.max()/2, x[0], span/6, y.max()/2, x[-1], span/6]
    if name=="InvGauss":
        return [y.min(), y.max()-y.min(), peak, span/4]
    if name=="DoG":
        return [y.max(), peak, span/6, y.max()/2, peak, span/3]
    if name=="GMM3":
        mus = x[0] + np.linspace(0.2,0.8,3)*span
        A   = y.max()/3; s=span/6
        return [A,mus[0],s, A,mus[1],s, A,mus[2],s]
    if name=="GMM4":
        mus = x[0] + np.linspace(0.1,0.9,4)*span
        A   = y.max()/4; s=span/8
        return [A,mus[0],s, A,mus[1],s, A,mus[2],s, A,mus[3],s]
    if name=="TruncGauss":
        return [y.max(), peak, span/4, 2.0]
    if name=="FoldedGauss":
        return [y.max(), peak, span/4]
    if name=="SplitGauss":
        return [y.max(), peak, span/6, span/6]
    return None

def physically_plausible(name, popt, lo=0.0, hi=180.0):
    if name in ("GMM3", "GMM4", "TwinGauss", "DoG"):
        P = np.asarray(popt).reshape(-1, 3)
        mu, s = P[:, 1], np.abs(P[:, 2])
        w = hi - lo
        return bool(np.all((mu >= lo - w) & (mu <= hi + w)) and np.all(s <= w))
    return True

# 7) Fit each region
def fit_region_full(x, y, region_label="", quiet=False):
    rows, best = [], None
    for name, (fn, k) in candidates.items():
        if k > len(x):
            continue
        p0 = initial_guess(name, x, y)
        try:
            with warnings.catch_warnings():
                warnings.simplefilter("ignore")
                popt, pcov = curve_fit(fn, x, y, p0=p0, maxfev=20000)
            yhat = fn(x, *popt)
            r2, rmse, mae, ssr, aic, aicc, bic = metrics(y, yhat, len(popt))

            cov_ok = bool(np.all(np.isfinite(pcov)))
            rel_se = (np.nanmax(np.sqrt(np.abs(np.diag(pcov))) / (np.abs(popt) + 1e-12))
                      if cov_ok else np.nan)
            plausible = physically_plausible(name, popt)
            saturated = len(x) <= len(popt) + 1

            rows.append({"Model": name, "R2": r2, "RMSE": rmse, "MAE": mae,
                         "SSR": ssr, "AIC": aic, "AICc": aicc, "BIC": bic,
                         "CovOK": cov_ok, "MaxRelSE": rel_se,
                         "Plausible": plausible, "Saturated": saturated})

            cand = {"name": name, "popt": popt, "r2": r2, "rmse": rmse, "mae": mae,
                    "ssr": ssr, "aic": aic, "aicc": aicc, "bic": bic}
            if plausible and (best is None or cand["aic"] < best["aic"]):
                best = cand
        except Exception:
            continue

    df = pd.DataFrame(rows).round(4)
    if not quiet:
        print(f"\n## {region_label}: All Candidate Metrics\n")
        print(df.to_markdown(index=False))
    return best, df


res1, table1 = fit_region_full(x1, y1, region_label="Region 1")
res2, table2 = fit_region_full(x2, y2, region_label="Region 2")

for res, x, y in [(res1, x1, y1), (res2, x2, y2)]:
    fn = candidates[res["name"]][0]
    popt, pcov = curve_fit(fn, x, y, p0=res["popt"], maxfev=20000)
    se = np.sqrt(np.diag(pcov)) if np.all(np.isfinite(pcov)) else None
    print(res["name"], "params:", np.round(popt, 3))
    print("  std errors:", None if se is None else np.round(se, 3))


# 8) Evaluate fits & stitch
y1_fit = candidates[res1["name"]][0](x1, *res1["popt"])
y2_fit = candidates[res2["name"]][0](x2, *res2["popt"])

# Dense fit + blend
theta_fit = np.linspace(theta.min(), theta.max(), 600)
fit_line  = np.empty_like(theta_fit)
m1_fit    = theta_fit <= split_ang
m2_fit    = ~m1_fit

# Region baseline
fit_line[m1_fit] = candidates[res1["name"]][0](theta_fit[m1_fit], *res1["popt"])
fit_line[m2_fit] = candidates[res2["name"]][0](theta_fit[m2_fit], *res2["popt"])

# Smooth blend ±2°
bw = 2.0
mask_blend = (theta_fit>=split_ang-bw)&(theta_fit<=split_ang+bw)
w = ((split_ang+bw) - theta_fit[mask_blend])/(2*bw)

y1b = candidates[res1["name"]][0](theta_fit[mask_blend], *res1["popt"])
y2b = candidates[res2["name"]][0](theta_fit[mask_blend], *res2["popt"])
fit_line[mask_blend] = w*y1b + (1-w)*y2b

# 9) Overall metrics
y_pred = np.empty_like(flux)
y_pred[mask1] = y1_fit
y_pred[mask2] = y2_fit

resid = flux - y_pred
plt.figure(figsize=(8,3.5))
plt.axhline(0, color='k', lw=0.8); plt.axvline(split_ang, color='gray', ls='--')
plt.plot(theta, resid, 'o-', ms=4)
plt.xlabel("Receiver circumferential angle (°)"); plt.ylabel("Residual (LCR)")
plt.tight_layout(); plt.savefig("Residuals_Jeter_0deg.png", dpi=300); plt.show()

# Compute global
r2_all, rmse_all, mae_all, ssr_all, aic_all, aicc_all, bic_all = \
    metrics(flux, y_pred, len(res1["popt"])+len(res2["popt"]))

# CV‐RMSE & CV‐MAE 
kf      = KFold(n_splits=5, shuffle=True, random_state=0)
errs_rmse, errs_mae = [], []
for _,test in kf.split(theta):
    e = flux[test] - y_pred[test]
    errs_rmse.append(np.sqrt(np.mean(e**2)))
    errs_mae.append(np.mean(np.abs(e)))
cv_rmse = np.mean(errs_rmse)
cv_mae  = np.mean(errs_mae)

# 10) Plot
plt.figure(figsize=(8,5))
plt.scatter(theta, flux, s=20, label="Data")
plt.plot(theta_fit, fit_line, '-r', lw=2, label="Blended Fit")
plt.axvline(split_ang, color='k', ls='--', label=f"Split @ {split_ang:.2f}°")
plt.xlabel("θ (°)")
plt.ylabel("Flux")
plt.title(f"{res1['name']}→{res2['name']} | R²={r2_all:.4f}, RMSE={rmse_all:.4f}")
plt.legend(); plt.grid(alpha=0.3); plt.tight_layout()
plt.show()

rows = []
for sp in [12.53, 15, 17.5, 20, 22.5, 25, 27.5, 30, 37.5, 45, 67.5]:
    m1, m2 = theta <= sp, theta > sp
    r1, _ = fit_region_full(theta[m1], flux[m1], quiet=True)
    r2, _ = fit_region_full(theta[m2], flux[m2], quiet=True)
    if r1 is None or r2 is None:
        continue
    yp = np.empty_like(flux)
    yp[m1] = candidates[r1["name"]][0](theta[m1], *r1["popt"])
    yp[m2] = candidates[r2["name"]][0](theta[m2], *r2["popt"])
    k = len(r1["popt"]) + len(r2["popt"])
    R2, RM, MA, SS, AI, AC, BI = metrics(flux, yp, k)
    rows.append({"split": sp, "n1": int(m1.sum()), "Region1": r1["name"],
                 "Region2": r2["name"], "k": k, "R2": R2, "RMSE": RM, "AIC": AI, "BIC": BI})
print(pd.DataFrame(rows).round(4).to_markdown(index=False))

# 11) Summary
print("Region1:", res1["name"], "params:", np.round(res1["popt"],4))
print(" R²,RMSE,MAE,SSR,AIC,AICc,BIC =",
      f"{res1['r2']:.4f},{res1['rmse']:.4f},{res1['mae']:.4f},"
      f"{res1['ssr']:.2f},{res1['aic']:.1f},{res1['aicc']:.1f},{res1['bic']:.1f}")
print("Region2:", res2["name"], "params:", np.round(res2["popt"],4))
print(" R²,RMSE,MAE,SSR,AIC,AICc,BIC =",
      f"{res2['r2']:.4f},{res2['rmse']:.4f},{res2['mae']:.4f},"
      f"{res2['ssr']:.2f},{res2['aic']:.1f},{res2['aicc']:.1f},{res2['bic']:.1f}")
print("\nOverall: ",
      f"R²={r2_all:.4f}, RMSE={rmse_all:.4f}, MAE={mae_all:.4f}, SSR={ssr_all:.2f}")
print(f" AIC={aic_all:.1f}, AICc={aicc_all:.1f}, BIC={bic_all:.1f}")
print(f" CV-RMSE={cv_rmse:.4f}, CV-MAE={cv_mae:.4f}")


    
    