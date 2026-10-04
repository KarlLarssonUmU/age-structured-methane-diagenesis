# Age-structured methane diagenesis

MATLAB implementation of a one-dimensional, steady-state
diffusion–advection–reaction model for dissolved methane in lake sediments,
with reservoir-corrected source ages, smooth bump production, and linear
loss in the upper 10 cm. The production peak is fitted subject to a minimum
depth of 10 cm. Production and loss can overlap within the shallow layer.

**The standard run optimizes all six free parameters:** upward Darcy flux
$q$, loss coefficient $k_{\mathrm{loss}}$, excess and background production
amplitudes $P_0$ and $P_{\mathrm{bg}}$, peak depth $z_{\mathrm{peak}}$, and
width of the production tail $w_{\mathrm{tail}}$.

## Model summary

Depth $z$ increases downward from the sediment–water interface over
$0\le z\le L$, with $L=3.62$ m. The steady-state methane balance and upward
dissolved flux are

$$
-\frac{dJ}{dz}+k_{\mathrm{loss}}f(z)C=P(z),
\qquad
J(z)=D_\phi(z)\frac{dC}{dz}+qC(z)
$$

Here $C$ is dissolved methane concentration per porewater volume, $q\ge0$
is the upward Darcy flux, and $D_\phi=\phi D_{\mathrm{eff}}$, where $\phi$
is porosity and $D_{\mathrm{eff}}=D_0/[1-\ln(\phi^2)]$. The molecular
diffusion coefficient is $D_0=1.0158\times10^{-9}$ m² s⁻¹. Production and
loss are expressed per bulk wet sediment volume. The loss coefficient
$k_{\mathrm{loss}}$ is expressed in s⁻¹ in this equation and reported in
day⁻¹ in the fitting bounds and model output. The mask $f(z)$ equals one in
$[0,z_{\mathrm{loss}}]$ and zero below, with $z_{\mathrm{loss}}=0.10$ m.
This is an effective removal term, without distinguishing oxidation from
other removal pathways.

The sediment–water and bottom boundary conditions are

$$
J(0)=k_{\mathrm{sw}}C(0),\qquad J(L)=0
$$

The surface condition is an empirical flux–concentration relation. At each
trial Darcy flux, $k_{\mathrm{sw}}=D_\phi(0)g/c_0+q$, where $g$ and $c_0$
are the slope and intercept fitted to the first eight measured methane
concentrations. The bottom has zero total dissolved flux.

Production is a background plus a smooth bump:

$$
P(z)=P_{\mathrm{bg}}+P_0 b(z;z_{\mathrm{peak}},w_{\mathrm{tail}}),
\qquad z_{\mathrm{peak}}\ge z_{\mathrm{loss}}
$$

The nonnegative bump has a maximum of one at $z_{\mathrm{peak}}$ and
declines to zero at $z_{\mathrm{peak}}+w_{\mathrm{tail}}$; its value and
first derivative are continuous. Thus $P_0$ is the excess peak amplitude and
$w_{\mathrm{tail}}$ is the **width of the production tail**. Both amplitudes
are nonnegative. The peak constraint still permits background production
and the rising limb within the loss layer.

For fixed parameters, the equation and boundary conditions are **linear in
concentration**, allowing methane to be decomposed into source/age classes.
Each class $C_i$ solves the same transport and loss problem for its own
production source. Total concentration and the radiocarbon signature are

$$
C(z)=\sum_i C_i(z),\qquad
F_m(z)=\frac{\sum_i F_{m,i}^{\mathrm{source}}C_i(z)}{\sum_i C_i(z)}
$$

Source signatures are assigned from the reservoir-corrected apparent
age-depth data. The same loss coefficient applies to every class. Class
surface fluxes $k_{\mathrm{sw}}C_i(0)$ give their contributions to dissolved
export; their age labels describe source carbon, not methane residence
times. The bump is implemented in [productionShape.m](src/productionShape.m).
The equation is solved using linear finite elements with a default mesh
spacing of 2.5 mm.

## Run

Requires MATLAB with Optimization Toolbox and Global Optimization Toolbox.
The code was run with MATLAB R2025b Update 2. No external workbooks or
third-party MATLAB packages are required.

Open this directory as the MATLAB current folder and run:

```matlab
result = runMainModel();
```

The code adds `src` to the MATLAB path, fits the model, and creates
`results/main/` with parameter and summary tables, concentration/Fm/production
profiles, source contributions, an age distribution, a MAT result, and
PNG/PDF figures. A repeat run replaces those generated files.
If the selected production fit or an outer search stops without convergence,
MATLAB issues a warning; the returned and saved results retain the optimizer
exit flags.

Rows in `source_contributions.csv` are already integrated contributions:
sum rows directly to combine source intervals, without multiplying by mesh
spacing. Export includes dissolved diffusion and advection. The reported
legacy fraction uses source depths strictly greater than 10 cm, with the
node at exactly 10 cm assigned to the upper interval.

## Fitting

The fitting procedure is nested. At each trial `(q,k_loss)`, all four
production parameters minimize the unweighted concentration sum of squared
residuals. The outer search selects `(q,k_loss)` by the unweighted sum of
squared residuals against the **17 individual measured methane Fm values**.
The measures minimized are

$$
E_C=\sum_{i=1}^{27}\left[C(z_i)-C_i^{\mathrm{obs}}\right]^2,
\qquad
E_F=\sum_{j=1}^{17}\left[F_m(\zeta_j)-F_{m,j}^{\mathrm{obs}}\right]^2
$$

These are separate objectives, not a weighted sum. For each trial production
shape, the two amplitudes are fitted by nonnegative least squares. Shape
fitting uses `fmincon`; the outer search uses a grid and multiple
`patternsearch` starts, with the parameter boundaries searched as well.

| Free parameter | Default constraint |
|---|---|
| Upward Darcy flux, $q$ | 0–10⁻⁸ m s⁻¹ |
| Loss coefficient, $k_{\mathrm{loss}}$ | 0–0.25 day⁻¹ |
| Excess production, $P_0$ | Nonnegative |
| Background production, $P_{\mathrm{bg}}$ | Nonnegative |
| Peak depth, $z_{\mathrm{peak}}$ | $z_{\mathrm{loss}}$–$0.60\,\mathrm{m}$ |
| Width of production tail, $w_{\mathrm{tail}}$ | 0.10–0.60 m |

## Data

| File in `data/` | Rows | Input columns |
|---|---:|---|
| `sediment_porosity.csv` | 94 | `depth_m`, `porosity_fraction` |
| `methane_concentration.csv` | 27 | `depth_m`, `CH4_mM` |
| `sediment_source_reservoir_corrected.csv` | 401 | `depth_m`, `age_BP`, `source_Fm` |
| `methane_fm_observations.csv` | 17 | `depth_cm`, `Fm`, `Fm_sd` |

Concentrations in mM equal mol m⁻³ of porewater; porosity is a volume
fraction. The Fm observations use centimetres, converted to metres on
loading. Their standard deviations are retained for plotting, not used as
fitting weights. The file also retains laboratory IDs and source-row labels.

The data comprise measured sediment porosity and porewater methane
concentration profiles, a reservoir-corrected median age-depth profile for
sediment organic carbon, and 17 individual methane radiocarbon measurements.
All 17 methane Fm measurements are used in the fit.

The sediment source ages are uncalibrated apparent radiocarbon ages. Source
Fm is calculated as `exp(-age_BP/8033)`. This adjusted source signature is a
modeling assumption, not a reconstruction of measured bulk organic-carbon
Fm. Negative ages and Fm greater than one are retained. Age and Fm are
interpolated separately; only depths 0–3.62 m are modeled.

## Configuration

```matlab
addpath('src');
cfg = defaultConfig();
cfg.outputDir = fullfile(pwd,'results','custom');
result = runMainModel(cfg);
```

`src/defaultConfig.m` defines the selected model and the search ranges.
`fitDarcy` and `fitLoss` default to `true`; fixed-value switches are available
for explicitly configured sensitivity runs. `lossDepth` controls both the
loss-layer boundary and the effective lower bound on peak depth. It must lie
on a mesh node. The default `qMagnitude` and `kLoss` placeholders are ignored
when their corresponding fit switches are true.

## Citation

If you use this code, please cite the accompanying paper:

> M. Bulínová, A. Rouillard, K. Larsson, X. Xu, J. Walker, C. Olid, J. Rydberg, C. Gudasz, S. E. Kjellman, G. Panieri, C. I. Czimczik, A. Schomacker (2026). *Legacy carbon stocks fuel methane diffusion from northern lake sediments*. Submitted manuscript. DOI: TBD

A machine-readable citation file is also included as `CITATION.cff`.

## License

The code is released under the MIT license in `LICENSE`.
