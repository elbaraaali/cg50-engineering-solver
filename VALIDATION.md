# Validation policy

This project does **not** claim that software can be mathematically error-free in every possible case.

Safety strategy:
1. Reject invalid inputs instead of manufacturing an answer.
2. Use finite-number checks after sensitive operations.
3. Validate numerical-method prerequisites.
4. Use explicit convergence status and iteration caps.
5. Use partial pivoting for determinant, inverse, and linear-system solving.
6. Keep features disabled if the UI path cannot be verified in the current environment.
7. Run host tests with `-Wall -Wextra -Werror -pedantic`.

The core test program covers parsing, trigonometry, domain errors, four root methods,
numerical differentiation, integration, ODE solvers, matrix operations, singular matrices,
and common Dynamics equations.
