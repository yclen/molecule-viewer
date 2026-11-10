# TTA Flora 4: the math

This is the model behind `rules.py` and `engine.py` written as equations: what a state is, the 10 rules, how the rules
become one matrix acting on a vector of probabilities, and how the three ways of running it (Monte Carlo, exact, pulse)
use that matrix. Simulated time is in nanoseconds throughout, and rates are per ns.

---

## 1. The state

A molecule has $N = n + 2$ sites in a row:

| Site | Molecule |
|---|---|
| $0$ | sensitizer |
| $1 \dots n$ | DPA mediators |
| $n + 1$ | perylene |

Each site holds one letter from the alphabet $\{G, T, S\}$ (ground, triplet, singlet). The state of the molecule is
the whole string

$$
x = (x_0, x_1, \dots, x_{n+1}), \qquad x_i \in \{G, T, S\}.
$$

For example, with $n = 3$ the string `T G S G T` means the sensitizer holds a triplet, DPA 2 holds a singlet, and the
perylene holds a triplet.

### As a number

The code stores letters as $G = 0$, $T = 1$, $S = 2$ and reads the string as a base-3 number, with site 0 as the last
digit:

$$
\text{index}(x) = \sum_{i=0}^{N-1} x_i \, 3^{\,i}.
$$

So there are $3^N$ possible strings, numbered $0$ to $3^N - 1$, and index $0$ is all $G$. With the default $n = 8$
that's $3^{10} = 59{,}049$ strings.

### The probability vector

Because the rules are random, the system isn't in one string but has a probability of being in each. Collect them in a
column vector with one entry per string:

$$
\mathbf{P}(t) = \big(P_0(t),\, P_1(t),\, \dots,\, P_{3^N - 1}(t)\big)^{\!\top}, \qquad \sum_x P_x(t) = 1.
$$

A pulse run starts with every molecule in all $G$, which is $\mathbf{P}(0) = (1, 0, 0, \dots, 0)^\top$.

---

## 2. Parameters and rates

The inputs are lifetimes (ns) and ratios. The rates the rules use are worked out from them (`resolve` in `rules.py`),
unless `kHop` or `kExc` is given directly, as the live window does:

| Rate | Formula | Default |
|---|---|---|
| Excitation | $k_\text{exc} = \text{excRatio} / \tau_\text{sens}$ | $0.1 / 5 = 0.02\ \text{ns}^{-1}$ (pulse runs: $1/5 = 0.2$) |
| Intersystem crossing | $k_\text{ISC} = \dfrac{1}{\tau_\text{sens}} \cdot \dfrac{Y}{1 - Y}$, with $Y$ = tripletYield | $0.2 \cdot 0.5/0.5 = 0.2\ \text{ns}^{-1}$ |
| Sensitizer singlet decay | $1 / \tau_\text{sens}$ | $0.2\ \text{ns}^{-1}$ |
| Hopping | $k_\text{hop} = \text{hopRatio} / \tau_S$ | $100 / 5 = 20\ \text{ns}^{-1}$ |
| Singlet emission (DPA, perylene) | $1 / \tau_S$ | $0.2\ \text{ns}^{-1}$ |
| Triplet decay (DPA, perylene) | $k_\text{decay}$ | $10^{-6}\ \text{ns}^{-1}$ |
| Spin factor | $f$ (share of TTA events that make a singlet) | $0.4$ |

The ISC formula comes from ISC and singlet decay competing on the sensitizer: the share of its singlets that become
triplets is $k_\text{ISC} / (k_\text{ISC} + 1/\tau_\text{sens}) = Y$.

### FRET back transfer

A singlet on site $i$ can go back to the sensitizer instead of emitting. The share that does is

$$
E_i = \frac{1}{1 + (r_i / R_0)^6},
$$

where $r_i$ is the distance from the sensitizer to site $i$ and $R_0$ = fretRadius (default 4 nm). The distances are
measured from the real geometry (`site_distances` in `rules.py`). DPA 1 is at 1.40 nm, and each further DPA adds 1.32 nm:

| Site | DPA 1 | DPA 2 | DPA 3 | DPA 4 | DPA 5 | DPA 6 | DPA 7 | DPA 8 | perylene |
|---|---|---|---|---|---|---|---|---|---|
| $r_i$ (nm) | 1.40 | 2.72 | 4.04 | 5.36 | 6.68 | 8.00 | 9.33 | 10.65 | 11.94 |
| $E_i$ | 0.998 | 0.910 | 0.484 | 0.147 | 0.044 | 0.015 | 0.006 | 0.003 | 0.001 |

### Hop rates

Write $h_{a \to b}$ for the rate a triplet on site $a$ moves to neighbouring site $b$ (`hops` in `rules.py`):

$$
h_{0 \to 1} = k_\text{hop}, \qquad
h_{1 \to 2} = k_\text{hop}, \qquad
h_{i \to i \pm 1} = \tfrac{1}{2} k_\text{hop} \ \ \text{for } 2 \le i \le n.
$$

DPA 1 only moves forward. For the last DPA, "forward" means the perylene. Nothing hops onto the sensitizer or off the
perylene. If $n = 0$, the sensitizer hops straight to the perylene at $k_\text{hop}$.

---

## 3. The 10 rules

Each rule looks for a pattern on one site or on a pair of sites, and rewrites it at a fixed rate. In two-letter rules,
the first letter is site $a$ (where the triplet hops from) and the second is site $b$ (where it lands). Every site not
named is left unchanged.

| # | Rule | Sites | Pattern → result | Rate |
|---|---|---|---|---|
| 1 | Excitation | $0$ | $G \to S$ | $k_\text{exc}(t)$ |
| 2 | Intersystem crossing | $0$ | $S \to T$ | $k_\text{ISC}$ |
| 3 | Triplet decay | $i \ge 1$ | $T \to G$ | $k_\text{decay}$ |
| 4 | Singlet emission (a photon) | $i \ge 1$ | $S \to G$ | $(1 - E_i)/\tau_S$ if $x_0 = G$, otherwise $1/\tau_S$ |
| 5 | Hop | $a, b$ neighbours | $T\,G \to G\,T$ | $h_{a \to b}$ |
| 6 | Blocked hop | $a, b$ neighbours | $T\,S$: no change | $0$ |
| 7 | TTA makes a singlet | $a, b$ neighbours | $T\,T \to G\,S$ | $f \, h_{a \to b}$ |
| 8 | TTA fails | $a, b$ neighbours | $T\,T \to G\,T$ | $(1 - f)\, h_{a \to b}$ |
| 9 | Sensitizer singlet decay | $0$ | $S \to G$ | $1/\tau_\text{sens}$ |
| 10 | Back transfer | $0$ and $i \ge 1$ | $(x_0, x_i) = (G, S) \to (S, G)$ | $E_i / \tau_S$ |

Notes on the rules:

- **Rules 4 and 10 share one singlet.** When the sensitizer is in $G$, a singlet on site $i$ leaves at the full rate
  $1/\tau_S$, split into a photon with probability $1 - E_i$ and a back transfer with probability $E_i$. When the
  sensitizer is in $S$ or $T$ it can't accept the energy, so FRET is off and the singlet always emits. In the code,
  rule 4 is therefore written as a two-site rule on $(0, i)$ that checks the sensitizer's letter.
- **Rule 6 is the absence of a rule.** There's simply no $T\,S$ rewrite, so a triplet next to a singlet waits. (With
  `block = 0` there is one, $T\,S \to G\,T$ at $h_{a \to b}$, which is only there for comparison.)
- **Rules 7 and 8 together** fire at the hop rate, $f\,h + (1-f)\,h = h$: a triplet meeting a triplet always
  annihilates. A failed TTA leaves one triplet, on the landing site $b$.
- **The sensitizer's triplet doesn't decay.** It can only leave by hopping (rule 5) or by TTA with DPA 1 (rules 7, 8).
- **The pulse.** $k_\text{exc}(t) = k_\text{exc}$ for $0 \le t < \text{duration}$ and $0$ afterwards. In long-run
  (steady) runs, the light is always on.

Rules with rate 0 are dropped, and each pattern is written out once for every site or pair it applies to. With $n = 8$
this gives 96 concrete rules.

---

## 4. How the rules act on the vector: the master equation

The probability vector changes by the **master equation**. For each string $x$:

$$
\frac{dP_x}{dt} \;=\; \underbrace{\sum_{x' \ne x} W(x' \to x)\, P_{x'}}_{\text{flowing in}}
\;-\; \underbrace{\sum_{x' \ne x} W(x \to x')\, P_x}_{\text{flowing out}},
$$

where $W(x \to x')$ is the total rate of every rule that rewrites string $x$ into string $x'$.

For example, rule 1 contributes

$$
W(x \to x') = k_\text{exc}(t) \quad \text{if } x_0 = G,\ x'_0 = S, \text{ and } x'_i = x_i \text{ for every } i \ne 0,
$$

and $0$ for every other pair. In words: excitation takes probability from each string whose sensitizer is in $G$, and
moves it to the same string with the sensitizer in $S$.

### In matrix form

Put the rates in a $3^N \times 3^N$ matrix $Q$, the **generator** (`exact` in `engine.py` builds it):

$$
Q_{x' x} = W(x \to x') \ \ \text{for } x' \ne x, \qquad Q_{x x} = -\sum_{x' \ne x} W(x \to x').
$$

Column $x$ holds everything that can happen to string $x$. Off the diagonal are the rates it turns into each other
string. On the diagonal is minus their total, its rate of leaving. Every column adds up to zero, which is what keeps
the total probability at 1. Then the master equation is one line:

$$
\boxed{\;\frac{d\mathbf{P}}{dt} = Q(t)\, \mathbf{P}(t)\;}
$$

### Building $Q$ rule by rule

$Q$ is a sum with one term per rule, and each term only touches the sites its rule names. On a single site, use these
$3 \times 3$ matrices in the basis $(G, T, S)$:

$$
\sigma_{v \leftarrow u} = |v\rangle\langle u| \quad \text{(moves probability from letter } u \text{ to letter } v\text{)}, \qquad
n_u = |u\rangle\langle u| \quad \text{(keeps letter } u\text{)}.
$$

A one-site rule that rewrites $u \to v$ on site $i$ at rate $k$ contributes

$$
Q_\text{rule} = k \left( \sigma^{(i)}_{v \leftarrow u} - n^{(i)}_u \right),
$$

where the superscript means "acting on site $i$, and the identity on every other site". The first term moves the
probability; the second takes the same amount off the diagonal. Rule 1, written out on the sensitizer, is

$$
k_\text{exc}
\begin{pmatrix} -1 & 0 & 0 \\ 0 & 0 & 0 \\ 1 & 0 & 0 \end{pmatrix}
\quad \text{(rows and columns in the order } G, T, S\text{)},
$$

so $G$ loses probability at $k_\text{exc}$ and $S$ gains the same.

A two-site rule that rewrites $(u_a, u_b) \to (v_a, v_b)$ on sites $a, b$ at rate $k$ contributes

$$
Q_\text{rule} = k \left( \sigma^{(a)}_{v_a \leftarrow u_a}\, \sigma^{(b)}_{v_b \leftarrow u_b} - n^{(a)}_{u_a}\, n^{(b)}_{u_b} \right).
$$

When a letter doesn't change, $\sigma_{u \leftarrow u} = n_u$. That's how rule 4 checks the sensitizer without
touching it: it uses $n^{(0)}_G$ on site 0 and $\sigma^{(i)}_{G \leftarrow S}$ on site $i$.

"Acting on site $i$" is a Kronecker product. Since site 0 is the last base-3 digit of the index, it's the last factor:

$$
A^{(i)} = I \otimes \dots \otimes I \otimes \underbrace{A}_{\text{site } i} \otimes I \otimes \dots \otimes I, \qquad
\text{site order } (N-1, \dots, 1, 0).
$$

The full generator is the sum over all 96 rules:

$$
Q(t) = \sum_{\text{rules } r} Q_r, \qquad \text{with only the rule 1 term depending on } t.
$$

---

## 5. Three ways to use $Q$

### Long-run (steady) behaviour: `exact`

With the light always on, the probabilities settle to a fixed vector $\boldsymbol{\pi}$ that no longer changes:

$$
Q\, \boldsymbol{\pi} = 0, \qquad \sum_x \pi_x = 1.
$$

`exact` solves this sparse linear system directly. The rate at which any rule fires is then its rate times the share
of time spent in strings where it can fire. The quantum yield is photons per excitation:

$$
\text{QY} = \frac{\displaystyle \sum_x \pi_x \sum_{\text{photon rules } r \text{ enabled at } x} k_r}
{\displaystyle k_\text{exc} \sum_{x \,:\, x_0 = G} \pi_x}.
$$

### Pulse runs: the time-dependent solution

A pulse switches $Q$ once. Starting from all $G$:

$$
\mathbf{P}(t) =
\begin{cases}
e^{Q_\text{on}\, t}\, \mathbf{P}(0), & 0 \le t < d, \\[4pt]
e^{Q_\text{off}\, (t - d)}\, \mathbf{P}(d), & t \ge d,
\end{cases}
$$

where $d$ is the pulse duration, $Q_\text{on}$ includes rule 1 and $Q_\text{off}$ doesn't. The expected photon rate
per molecule at time $t$ is

$$
I(t) = \sum_x P_x(t) \sum_{\text{photon rules } r \text{ enabled at } x} k_r,
$$

so with $M$ molecules, the expected count in a histogram bin from $t_1$ to $t_2$ is $M \int_{t_1}^{t_2} I(t)\, dt$.
This is the exact curve the pulse histogram estimates, and `pulse_exact` in `engine.py` works it out. `live.py`
draws it over the counts when "Exact solve" is ticked.

`pulse_exact` doesn't form $e^{Qt}$. It adds three rows to $Q$ that add up perylene photons, DPA photons and
excitations, so their entries in $\mathbf{P}$ become running totals, and steps $d\mathbf{P}/dt = Q\,\mathbf{P}$
through time. Hops are far faster than emission, so the system is stiff: an ordinary step would have to be shorter than
a hop to stay stable. It uses a two-stage L-stable SDIRK instead, which stays stable with steps far longer than a hop
and damps the fast hop transients rather than tracking them. Each step is

$$
\mathbf{k}_1 = (I - \gamma h Q)^{-1}\, Q\,\mathbf{P}_n, \qquad
\mathbf{k}_2 = (I - \gamma h Q)^{-1}\, Q\,\big(\mathbf{P}_n + (1 - \gamma) h\, \mathbf{k}_1\big), \qquad
\mathbf{P}_{n+1} = \mathbf{P}_n + h\big((1 - \gamma)\,\mathbf{k}_1 + \gamma\,\mathbf{k}_2\big),
$$

with $\gamma = 1 - 1/\sqrt{2}$. Both stages solve with the same matrix $I - \gamma h Q$, by GMRES with an incomplete
LU of it worked out once per phase (pulse on, pulse off). The step is a tenth of the shorter of `tauS` and `tauSens`,
and a step ends exactly at $d$. After the window it keeps stepping until the molecules that could still emit are
under $10^{-7}$ of the excitations, so the quantum yield counts the photons after the window too. Bins are read off the
running totals with a monotone (PCHIP) interpolation, so any bin width works.

Its cost is set by the $3^N$ strings: at the defaults it takes a few seconds up to $n = 7$, about 30 s at $n = 8$ and a
few minutes at $n = 9$. The step count barely depends on `hopRatio`, which is what the implicit stepping buys. Against
Monte Carlo with millions of molecules, the histograms agree bin by bin (chi-square) and the quantum yields agree
within the Monte Carlo error.

### Totals over one pulse: `pulse_totals`, for intensity sweeps

An intensity sweep needs no histogram, only how many of each event a molecule gives over the whole pulse and
everything after it. `pulse_totals` works those out exactly, in two parts.

**During the pulse**, it adds a row to $Q$ for every event at every site, so those entries of $\mathbf{P}$ become
running totals, as in `pulse_exact`. With the light on, $Q$ doesn't change, so the answer at the end of the pulse is
just

$$
\mathbf{P}(d) = e^{Q d}\, \mathbf{P}(0),
$$

which `expm_multiply` works out directly, accurate to rounding. That matters at low power: there the photon totals
are of order $(k_\text{exc} d)^2$, far smaller than the all-$G$ entry near 1, and a stepping method's error, which is
relative to the whole vector, would swamp them.

**After the pulse**, no more light comes, so instead of following time it asks how long a molecule will spend in
each string from then on. Call the strings that can still give a photon active: two or more excited sites, or any
singlet. A string with just one triplet can't make another photon without light, so it's treated as an end, along
with all $G$. With $Q_A$ the part of $Q$ (light off) among active strings, the expected time in each active string is

$$
\boldsymbol{\tau} = (-Q_A)^{-1}\, \mathbf{P}_A(d),
$$

one sparse linear solve (an LU up to $n = 6$, GMRES with an incomplete LU above). Every later event of rule $r$ is
then $k_r$ times the time spent in the strings it fires from. A molecule that ends with one triplet is counted as one
triplet decayed, on the perylene: nothing hops off the perylene, and a hop is about $10^7$ times faster than triplet
decay, so the triplet gets there first. Hops after the pulse aren't counted.

$-Q_A$ doesn't depend on $k_\text{exc}$ or $d$, so a sweep of those factors it once and reuses it. At $n = 8$ that's
about 15 s once, then about 5 s a point, mostly the matrix exponential, which takes longer for long or bright pulses
(its cost goes as $\lVert Q \rVert\, d$). It needs no cut-off, the bookkeeping balances to rounding, and at low power
the photons rise as $k_\text{exc}^2$ to four decimal places, since a photon needs two absorptions in one pulse.

### Monte Carlo: `monte_carlo`, `pulse` and `live.py`

For large $n$, the vector is too long to work with, so instead of evolving the probabilities, the engine follows one
string at a time and makes random choices with exactly the probabilities $Q$ implies (the Gillespie method):

1. From the current string $x$, list every rule that can fire, with its rate $k_r$. Add them up:
   $R(x) = \sum_r k_r = -Q_{xx}$.
2. Wait a random time $\tau = -\ln(u_1) / R(x)$, with $u_1$ uniform on $(0, 1)$. This is the exponential waiting time
   with mean $1/R(x)$.
3. Pick which rule fires: rule $r$ with probability $k_r / R(x)$.
4. Apply it, set $t \to t + \tau$, count its event, and repeat.

Averaged over many molecules, the share of them in each string follows $d\mathbf{P}/dt = Q\,\mathbf{P}$ exactly. There's
no time step and no discretization error.

In a pulse run, if the drawn wait would cross the end of the pulse, the clock is set to the pulse end and a new wait is
drawn with $Q_\text{off}$. Exponential waits have no memory, so this is exact. A molecule is finished when it's back to
all $G$ after the pulse, where $R = 0$ and nothing else can happen.

---

## 6. The fixed time-step picture

The simplest way to think about one rule is "in each step $dt$, it fires with probability $k\,dt$". That's the master
equation with a small step (Euler's method):

$$
\mathbf{P}(t + dt) \approx (I + Q\, dt)\, \mathbf{P}(t),
$$

and for one rule on its own, $P(\text{fires in } dt) = k\, dt + O(dt^2)$. Over a finite step the exact answer is

$$
\mathbf{P}(t + \Delta t) = e^{Q\, \Delta t}\, \mathbf{P}(t), \qquad
P(\text{rule fires within } \Delta t \text{, alone}) = 1 - e^{-k \Delta t}.
$$

A fixed-step simulation only works when $R(x)\, dt \ll 1$ for the *total* rate, not just for the one rule, or two
events land in the same step. The fastest rule is hopping at $20\ \text{ns}^{-1}$, so $dt$ would have to be well under
$0.05$ ns. The Gillespie method in section 5 is the $dt \to 0$ limit of this picture, so it avoids that error and also
skips the empty steps.

---

## 7. A useful limit

If every triplet ends in TTA (no triplet decay, nothing left over), and every singlet made by TTA emits (no back
transfer), the yield has a simple form. Each TTA event uses $2$ triplets with probability $f$ and $1$ triplet with
probability $1 - f$, so on average it uses $1 + f$ triplets and makes $f$ singlets. Each absorbed photon makes a
triplet with probability $Y$, so

$$
\text{QY}_\text{ideal} = Y \cdot \frac{f}{1 + f} = 0.5 \cdot \frac{0.4}{1.4} = 14.3\% \quad \text{(defaults).}
$$

Real runs fall below this: at low excitation a lone triplet has no partner and decays, and singlets made near the
sensitizer can go back by FRET.
