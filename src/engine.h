#ifndef ENG_ENGINE_H
#define ENG_ENGINE_H
#include <stddef.h>

typedef enum { ANGLE_RAD=0, ANGLE_DEG=1 } AngleMode;

typedef enum {
    ENG_OK=0,
    ENG_ERR_PARSE,
    ENG_ERR_DOMAIN,
    ENG_ERR_DIVZERO,
    ENG_ERR_OVERFLOW,
    ENG_ERR_BRACKET,
    ENG_ERR_DERIVATIVE,
    ENG_ERR_NOCONVERGE,
    ENG_ERR_DIM,
    ENG_ERR_SINGULAR,
    ENG_ERR_BADARG
} EngError;

typedef struct {
    double x, y, ans;
    AngleMode angle;
} EvalVars;

const char *eng_error_string(EngError e);
EngError eng_eval(const char *expr, const EvalVars *vars, double *out);

typedef struct {
    double value;
    int iterations;
    double approx_error_pct;
    int converged;
    EngError error;
} SolverResult;

SolverResult eng_bisection(const char *expr, double xl, double xu, double tol_pct, int max_iter, const EvalVars *base);
SolverResult eng_false_position(const char *expr, double xl, double xu, double tol_pct, int max_iter, const EvalVars *base);
SolverResult eng_secant(const char *expr, double x0, double x1, double tol_pct, int max_iter, const EvalVars *base);
SolverResult eng_newton(const char *expr, double x0, double tol_pct, int max_iter, const EvalVars *base);

EngError eng_derivative(const char *expr, double x, double h, const EvalVars *base, double *out);
EngError eng_second_derivative(const char *expr, double x, double h, const EvalVars *base, double *out);
EngError eng_trapezoid(const char *expr, double a, double b, int n, const EvalVars *base, double *out);
EngError eng_simpson13(const char *expr, double a, double b, int n, const EvalVars *base, double *out);
EngError eng_limit_symmetric(const char *expr, double x0, double h, const EvalVars *base, double *out);
EngError eng_integral_simpson_auto(const char *expr, double a, double b, const EvalVars *base, double *out);
EngError eng_taylor_ode_step(const char *d1, const char *d2, const char *d3, const char *d4, int order, double x, double y, double h, const EvalVars *base, double *y_next);
EngError eng_simplify_polynomial(const char *expr, char *out, size_t cap);

typedef struct {
    double y_final;
    int steps;
    EngError error;
} ODEResult;
ODEResult eng_euler(const char *fxy, double x0, double y0, double h, int steps, const EvalVars *base);
ODEResult eng_heun(const char *fxy, double x0, double y0, double h, int steps, const EvalVars *base);
ODEResult eng_rk4(const char *fxy, double x0, double y0, double h, int steps, const EvalVars *base);

#define ENG_MAX_MAT 6
typedef struct {
    int rows, cols;
    double v[ENG_MAX_MAT][ENG_MAX_MAT];
} EngMatrix;
EngError eng_mat_det(const EngMatrix *a, double *det);
EngError eng_mat_inverse(const EngMatrix *a, EngMatrix *inv);
EngError eng_mat_solve(const EngMatrix *A, const double *b, double *x);
EngError eng_mat_add(const EngMatrix *a, const EngMatrix *b, EngMatrix *out);
EngError eng_mat_mul(const EngMatrix *a, const EngMatrix *b, EngMatrix *out);

double dyn_v(double u,double a,double t);
double dyn_s(double u,double a,double t);
EngError dyn_v_from_s(double u,double a,double s,double *v);
double dyn_projectile_range(double v,double theta_deg,double g);
double dyn_projectile_time(double v,double theta_deg,double g);
double dyn_projectile_hmax(double v,double theta_deg,double g);
EngError dyn_normal_accel(double v,double r,double *a_n);
double dyn_tangential_accel(double alpha,double r);
double dyn_ke(double m,double v);
double dyn_pe(double m,double g,double h);
double dyn_spring_pe(double k,double x);
double dyn_momentum(double m,double v);
double dyn_impulse(double F,double dt);
EngError dyn_restitution(double u1,double u2,double v1,double v2,double *e);

/* Constant-acceleration linear kinematics. s is displacement. */
EngError dyn_u_from_vat(double v,double a,double t,double *u);
EngError dyn_a_from_uvt(double u,double v,double t,double *a);
EngError dyn_t_from_uva(double u,double v,double a,double *t);
double dyn_s_from_uvt(double u,double v,double t);
double dyn_s_from_vat(double v,double a,double t);
EngError dyn_a_from_uvs(double u,double v,double s,double *a);
EngError dyn_t_from_us(double u,double a,double s,double *t1,double *t2,int *roots);

/* Variable-motion helpers using expressions in time variable x. */
EngError dyn_from_position(const char *s_expr,double t,const EvalVars *base,double *s,double *v,double *a);
EngError dyn_from_velocity(const char *v_expr,double t0,double t1,double s0,const EvalVars *base,double *v1,double *a1,double *s1);
EngError dyn_from_acceleration(const char *a_expr,double t0,double t1,double v0,double s0,const EvalVars *base,double *a1,double *v1,double *s1);

#endif
