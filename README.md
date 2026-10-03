# CG50 Engineering Solver v0.4

A calculator-style fx-CG50 engineering solver source project focused on Calculus, Dynamics and Numerical Methods.

## Main menu
- RUN-MAT
- CAS SIMPLIFY / EXPAND
- CALCULUS
- EQUATION SOLVER
- NUMERICAL METHODS
- DYNAMICS
- MATRIX
- SETTINGS
- ABOUT

## Dynamics
### Position / Velocity / Acceleration
Use `x` as time `t`.
- Given position `s(t)`: evaluates position, velocity `ds/dt`, acceleration `d2s/dt2` at a selected time.
- Given velocity `v(t)`: evaluates velocity and acceleration and integrates velocity to obtain position from an initial position.
- Given acceleration `a(t)`: evaluates acceleration and integrates it to obtain velocity and position from initial conditions.

### Constant acceleration
Includes direct/rearranged forms for `u`, `v`, `a`, `t`, and displacement `s`, plus quadratic time roots when required.

### Other Dynamics
Projectile range/time/max height, normal acceleration, tangential acceleration, kinetic/potential/spring energy, momentum, impulse, restitution.

## Numerical Methods
- Bisection
- False Position
- Newton-Raphson
- Secant
- Taylor ODE method (orders 1-4; user supplies required derivatives)
- Trapezoidal Rule
- Simpson 1/3
- Euler
- Heun
- RK4

## CAS simplify / expand
A bounded polynomial simplifier intended to be dependable rather than overclaim capabilities.
Examples:
- `x+x` -> `2*x`
- `2*x+3*x-4` -> `5*x - 4`
- `(x+1)^2` -> `x^2 + 2*x + 1`

Supported polynomial syntax: constants, `x`, `+`, `-`, `*`, `^`, parentheses; maximum degree 12.

## Validation
See `VALIDATION_REPORT.txt`.

## Build
Requires fxSDK + gint + SuperH cross compiler. Build target is fx-CG50 using `fxsdk build-cg`.
