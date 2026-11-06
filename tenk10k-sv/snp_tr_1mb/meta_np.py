"""numpy re-implementation of meta::metagen(te, sete, random = TRUE, method.tau = "DL") for k studies (vectorised over loci).

Inputs are arrays of shape (n_loci, k). Returns the same fields meta_runner_sv.py writes.
"""
import numpy as np
from scipy import stats


def metagen_dl(te, se):
    te, se = np.asarray(te, float), np.asarray(se, float)
    w = 1 / se ** 2
    sw = w.sum(1)
    fixed = (w * te).sum(1) / sw
    se_fixed = np.sqrt(1 / sw)
    z = fixed / se_fixed
    k = te.shape[1]
    Q = (w * (te - fixed[:, None]) ** 2).sum(1)
    df = k - 1
    C = sw - (w ** 2).sum(1) / sw
    tau2 = np.maximum(0, (Q - df) / C)
    wr = 1 / (se ** 2 + tau2[:, None])
    random = (wr * te).sum(1) / wr.sum(1)
    se_random = np.sqrt(1 / wr.sum(1))
    zr = random / se_random
    zcrit = stats.norm.ppf(0.975)
    return dict(
        coeff_meta_random=random, se_meta_random=se_random, pval_q_meta=stats.chi2.sf(Q, df), q_meta=Q,
        tau_meta=np.sqrt(tau2), i2_meta=np.maximum(0, (Q - df) / Q), h_meta=np.sqrt(np.maximum(Q / df, 1)),
        pval_meta_random=2 * stats.norm.sf(np.abs(zr)), pval_meta_fixed=2 * stats.norm.sf(np.abs(z)),
        coeff_meta_fixed=fixed, se_meta_fixed=se_fixed,
        lowerCI_meta_random=random - zcrit * se_random, upperCI_meta_random=random + zcrit * se_random)
