# Week 3 – Fold-Change Detection and Cancer (WNT / β-catenin)  

**Course:** Physics of Molecular Diseases  
**Topic:** Fold-Change Detection and Cancer

---

## Overview

This script extends the incoherent feed-forward loop (IFFL) fold-change
detector to a model of **colon (pre-)cancer** in the WNT / β-catenin
pathway, and investigates two modifications:

- **(i)** how cancer cells "cheat" the fold-detection mechanism to escape
  apoptosis under the same WNT treatment
- **(ii)** how alternating "fasting" with a fiber-rich Mediterranean diet
  can restore apoptosis and reduce cancer risk

---

## Base Model (IFFL = fold-change detector)

```
dx/dt = S - x
dy/dt = S/x - y
S(t)  = u_ext(t) + w*(1 - alpha*fasting(t))
```

| Symbol | Meaning |
|--------|---------|
| `U` (u_ext) | external WNT signal (Mediterranean diet activates it) |
| `X` | WNT-induced inhibitor (internal controller) |
| `Y` | β-catenin level (the output) |
| `w` | autonomous / constitutive WNT-β-catenin tone (the cancer "cheat") |
| `alpha` | strength of fasting's WNT-INDEPENDENT suppression of `w` |

At steady state `x* = S` and `y* = 1` (exact adaptation); the **peak of Y**
after a step is `~ S_new / S_old`, i.e. the fold-change of the input.

---

## Interpretation of Apoptosis

The paper states apoptosis ∝ fold-change of WNT, **and** P(apoptosis) ∝ Y.
Combined, the apoptosis-driving signal is the **FCD transient of Y** — the
β-catenin level ABOVE its adapted baseline (`y - 1`), which is ∝ fold-change.

```
apoptotic dose D = integral of max(y - 1, 0) dt
P(apoptosis)     = 1 - exp(-k * D)
```

This resolves the apparent paradox: cancer cells have a **high baseline Y**
but a **small fold-response**, so they avoid the apoptotic transient.

---

## (i) How Cancer Cells Cheat

Cancer cells carry an autonomous WNT/β-catenin tone `w` (e.g. from APC
mutation). This does two things at once:

1. raises the baseline (`x* = u_base + w`) — consistent with the observed
   upregulated β-catenin in cancer cells
2. **compresses the fold-change** of any subsequent WNT input:

```
fold = (u_diet + w) / (u_base + w)  ->  1   as w increases
```

So the same WNT treatment produces a smaller Y peak → smaller apoptotic
transient → the cancer cell escapes apoptosis.

**Numerical result (single diet step):**

| Cell | w | peak Y | apoptotic response (peak−1) |
|------|---|--------|------------------------------|
| pre-cancer | 0.1 | 1.92 | 0.92 |
| cancer | 3.0 | 1.10 | 0.10 |

Panel B shows the apoptotic response decaying monotonically as the
autonomous tone `w` grows.

---

## (ii) Fasting + Mediterranean Diet

- **Mediterranean diet** (fiber-rich): activates WNT → raises `u_ext`.
- **Fasting**: downregulates β-catenin through a **WNT-INDEPENDENT**
  mechanism → modelled as suppressing the autonomous tone `w` (`alpha`),
  which resets the baseline back down to `u_base`.

Because fasting acts independently of WNT, it bypasses the cancer cell's
cheat: after a fasting period the baseline is low even in cancer cells, so
the next diet phase delivers a **large fold-change** → a big Y spike →
apoptosis.

**Numerical result — P(apoptosis) over the protocol:**

| Cell | chronic diet | alternating fast/diet |
|------|-------------|------------------------|
| pre-cancer | 0.000 | 0.974 |
| cancer | 0.000 | 0.994 |

**Key interpretation:** chronic diet gives ZERO apoptotic signal in *both*
cells — a constant signal, however high, produces no FCD response after
adaptation. The benefit comes not from the diet *level* but from the
*changes* (the fasting↔diet transitions). The alternation, combined with
fasting's WNT-independent action, defeats the cheat and culls cancer cells.

---

## Script Structure

```
week3_fcd_cancer.py
│
├── SIMULATOR
│     simulate()        — integrate the IFFL with tone w and fasting term
│     apoptotic_dose()  — integral of (y - 1)+  (FCD transient)
│     P_apoptosis()     — 1 - exp(-k * dose)
│
├── INPUT SCHEDULES
│     u_step   — single diet step           (for part i)
│     u_const  — chronic Mediterranean diet (control for part ii)
│     u_alt / fast_alt — alternating fasting <-> diet (part ii)
│
├── FIGURE (4 panels)
│     A — single diet step: pre-cancer vs cancer (the cheat)
│     B — apoptotic response vs autonomous tone w (fold compression)
│     C — cancer cell: chronic diet (flat) vs alternating (repeated spikes)
│     D — P(apoptosis) bars: both cells × both protocols
│
└── CONSOLE SUMMARY
      peak Y / apoptotic response (i) and P(apoptosis) table (ii)
```

---

## Output

**Figure (4 panels):**

| Panel | Shows |
|-------|-------|
| A | cancer cell's Y peak (1.10) ≪ pre-cancer's (1.92) for the same WNT step |
| B | apoptotic response (peak Y − 1) falls as autonomous tone w rises |
| C | chronic diet → flat Y=1; alternating → repeated spikes (~3.4) per cycle |
| D | chronic diet P≈0 for both; alternating P≈0.97–0.99 for both |

**Console:** peak-Y / apoptotic-response table for (i) and the
P(apoptosis) comparison for (ii).

---

## Requirements

```
numpy
matplotlib
scipy
```

Note: this script uses `np.trapz`. A newer NumPy (NumPy 2.0+), replace it with `np.trapezoid`.

---

## Connection to the Broader Week 3 Analysis

| Concept | Link |
|---------|------|
| Fold-change detection (IFFL) | The base detector from Shoval/Goentoro |
| WNT / β-catenin (Goentoro & Kirschner 2009) | β-catenin shows FCD; mutations link to colon cancer |
| Exact adaptation | Why chronic (constant) diet produces no apoptotic signal |
| Receptor desensitization theme | Cancer cell pre-adapts (high baseline) to dodge the signal |
| Disease state | Loss of the protective fold-change response = escape of apoptosis |

The take-home: cancer cells exploit the *relative* (fold-change) nature of
the apoptotic signal by raising their baseline β-catenin, which compresses
the fold-change of incoming WNT and lets them dodge apoptosis. A therapy
that acts on a *different* axis (fasting, WNT-independent) can reset the
baseline and restore the protective fold-change response — providing a
physics-based rationale for alternating dietary interventions.

---

## Context

- Ala Trusina, Lecture slides from the course *Physics of Molecular Diseases* (Niels Bohr Institute, 2020). The broader week covers:

- Receptor desensitization (GPCR, insulin signalling, WNT pathway)
- Adaptation and fold-change detection (Weber's law)
- Feedback algebra and conditions for perfect adaptation
