# 1D–3V Klimontovich–Maxwell Plasma Simulation (Electromagnetic PIC)

An educational Fortran implementation of a **1D–3V electromagnetic particle-in-cell (PIC) model**. It approximates discrete-particle plasma dynamics associated with the Klimontovich description: particles move in self-consistent electromagnetic fields, while the fields are advanced numerically.

## Physics model

Particle characteristics:

$$
\frac{d\mathbf{x}_p}{dt}=\mathbf{v}_p,\qquad
\frac{d\mathbf{v}_p}{dt}=\frac{q}{m}(\mathbf{E}+\mathbf{v}_p\times\mathbf{B})
$$

The corresponding discrete-particle distribution is represented conceptually by

$$
F_K(\mathbf{x},\mathbf{v},t)=\sum_p w_p\,\delta(\mathbf{x}-\mathbf{x}_p)\,\delta(\mathbf{v}-\mathbf{v}_p)
$$

In 1D–3V, the fields vary with x only, but each particle has three velocity components. The transverse fields use the 1D Maxwell equations; the longitudinal electric field is reconstructed from Gauss's law using the deposited charge density.

## What the code includes

- Fortran 2008; normalized units with `c = epsilon0 = mu0 = 1`.
- Periodic spatial domain of length `2*pi`.
- 20,000 electron macro-particles and a stationary positive-ion background.
- Cloud-in-cell charge/current deposition.
- Staggered-grid updates for transverse electromagnetic fields.
- Non-relativistic Boris particle pusher.
- A small initial transverse electromagnetic perturbation.
- Diagnostics: particle energy, electric/magnetic field energy, total energy, and a discrete Gauss-law residual.
- Output files: `diagnostics.csv`, `fields_final.dat`, and `particles_final.dat`.

## Important numerical limitations

This is a **learning baseline, not a validated production PIC solver**. The longitudinal electric field is projected from Gauss's law each step, rather than advanced with a fully charge-conserving longitudinal current scheme. The simple current deposition is not exactly charge-conserving, and the transverse particle-current deposition is basic. Therefore, small Gauss-law residuals alone do not establish that the full Maxwell system is being solved consistently. The non-relativistic pusher is intended for speeds much smaller than c. Before using this for research, implement charge-conserving current deposition, validate against electromagnetic-wave and plasma-wave benchmarks, test grid/time-step convergence, and check energy conservation.

## Reading the outputs

- `diagnostics.csv`: energy history and Gauss-law residual versus time.
- `fields_final.dat`: columns are `x Ex Ey Ez Bx By Bz`.
- `particles_final.dat`: columns are `x vx vy vz`.

The included run reaches `t = 10`. Its energy history is an initial numerical result, **not** evidence of a validated physical instability or damping process.

## Compile and run

```bash
gfortran km_pic_1d3v.f90 -o km_pic_1d3v
./km_pic_1d3v
```
