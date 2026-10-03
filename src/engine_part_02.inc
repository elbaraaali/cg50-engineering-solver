    return fabs((newv-oldv)/newv)*100.0;
}
static int opposite_sign(double a,double b){
    return (a<0.0&&b>0.0)||(a>0.0&&b<0.0);
}
static SolverResult sr_error(EngError e){ SolverResult r={0,0,0,0,e}; return r; }

SolverResult eng_bisection(const char *expr,double xl,double xu,double tol_pct,int max_iter,const EvalVars *base){
    if(!expr||!finite_ok(xl)||!finite_ok(xu)||tol_pct<=0||max_iter<1||xl>=xu) return sr_error(ENG_ERR_BADARG);
    double fl,fu;
    EngError e=f_at(expr,xl,base,&fl); if(e) return sr_error(e);
    e=f_at(expr,xu,base,&fu); if(e) return sr_error(e);
    if(fl==0){ SolverResult r={xl,0,0,1,ENG_OK}; return r; }
    if(fu==0){ SolverResult r={xu,0,0,1,ENG_OK}; return r; }
    if(!opposite_sign(fl,fu)) return sr_error(ENG_ERR_BRACKET);
    double xr=xl,old=xr,fr=0,ea=HUGE_VAL;
    for(int i=1;i<=max_iter;i++){
        old=xr; xr=(xl+xu)/2.0;
        e=f_at(expr,xr,base,&fr); if(e) return sr_error(e);
        if(i>1) ea=ea_pct(xr,old);
        if(fr==0.0|| (i>1 && ea<=tol_pct)){ SolverResult r={xr,i,(fr==0?0:ea),1,ENG_OK}; return r; }
        if(opposite_sign(fl,fr)){ xu=xr; fu=fr; } else { xl=xr; fl=fr; }
    }
    SolverResult r={xr,max_iter,ea,0,ENG_ERR_NOCONVERGE}; return r;
}
SolverResult eng_false_position(const char *expr,double xl,double xu,double tol_pct,int max_iter,const EvalVars *base){
    if(!expr||!finite_ok(xl)||!finite_ok(xu)||tol_pct<=0||max_iter<1||xl>=xu) return sr_error(ENG_ERR_BADARG);
    double fl,fu;
    EngError e=f_at(expr,xl,base,&fl); if(e) return sr_error(e);
    e=f_at(expr,xu,base,&fu); if(e) return sr_error(e);
    if(fl==0){ SolverResult r={xl,0,0,1,ENG_OK}; return r; }
    if(fu==0){ SolverResult r={xu,0,0,1,ENG_OK}; return r; }
    if(!opposite_sign(fl,fu)) return sr_error(ENG_ERR_BRACKET);
    double xr=xl,old=xr,fr=0,ea=HUGE_VAL;
    for(int i=1;i<=max_iter;i++){
        double den=fl-fu;
        if(fabs(den)<=DBL_MIN) return sr_error(ENG_ERR_DIVZERO);
        old=xr; xr=xu-fu*(xl-xu)/den;
        e=f_at(expr,xr,base,&fr); if(e) return sr_error(e);
        if(i>1) ea=ea_pct(xr,old);
        if(fr==0.0||(i>1&&ea<=tol_pct)){ SolverResult r={xr,i,(fr==0?0:ea),1,ENG_OK}; return r; }
        if(opposite_sign(fl,fr)){ xu=xr; fu=fr; } else { xl=xr; fl=fr; }
    }
    SolverResult r={xr,max_iter,ea,0,ENG_ERR_NOCONVERGE}; return r;
}
SolverResult eng_secant(const char *expr,double x0,double x1,double tol_pct,int max_iter,const EvalVars *base){
    if(!expr||tol_pct<=0||max_iter<1||!finite_ok(x0)||!finite_ok(x1)||x0==x1) return sr_error(ENG_ERR_BADARG);
    double f0,f1;
    EngError e=f_at(expr,x0,base,&f0); if(e) return sr_error(e);
    e=f_at(expr,x1,base,&f1); if(e) return sr_error(e);
    double ea=HUGE_VAL;
    for(int i=1;i<=max_iter;i++){
        double den=f1-f0;
        if(fabs(den)<1e-15) return sr_error(ENG_ERR_DIVZERO);
        double x2=x1-f1*(x1-x0)/den;
        if(!finite_ok(x2)) return sr_error(ENG_ERR_OVERFLOW);
        ea=ea_pct(x2,x1);
        double f2; e=f_at(expr,x2,base,&f2); if(e) return sr_error(e);
        if(f2==0.0||ea<=tol_pct){ SolverResult r={x2,i,(f2==0?0:ea),1,ENG_OK}; return r; }
        x0=x1; f0=f1; x1=x2; f1=f2;
    }
    SolverResult r={x1,max_iter,ea,0,ENG_ERR_NOCONVERGE}; return r;
}
EngError eng_derivative(const char *expr,double x,double h,const EvalVars *base,double *out){
    if(!expr||!out||h<=0||!finite_ok(x)||!finite_ok(h)) return ENG_ERR_BADARG;
    double fm2,fm1,fp1,fp2;
    EngError e=f_at(expr,x-2*h,base,&fm2); if(e) return e;
    e=f_at(expr,x-h,base,&fm1); if(e) return e;
    e=f_at(expr,x+h,base,&fp1); if(e) return e;
    e=f_at(expr,x+2*h,base,&fp2); if(e) return e;
    *out=(fm2-8*fm1+8*fp1-fp2)/(12*h);
    return finite_ok(*out)?ENG_OK:ENG_ERR_OVERFLOW;
}
EngError eng_second_derivative(const char *expr,double x,double h,const EvalVars *base,double *out){
    if(!expr||!out||h<=0) return ENG_ERR_BADARG;
    double fm2,fm1,f0,fp1,fp2;
    EngError e=f_at(expr,x-2*h,base,&fm2); if(e) return e;
    e=f_at(expr,x-h,base,&fm1); if(e) return e;
    e=f_at(expr,x,base,&f0); if(e) return e;
    e=f_at(expr,x+h,base,&fp1); if(e) return e;
    e=f_at(expr,x+2*h,base,&fp2); if(e) return e;
    *out=(-fp2+16*fp1-30*f0+16*fm1-fm2)/(12*h*h);
    return finite_ok(*out)?ENG_OK:ENG_ERR_OVERFLOW;
}
SolverResult eng_newton(const char *expr,double x0,double tol_pct,int max_iter,const EvalVars *base){
    if(!expr||tol_pct<=0||max_iter<1||!finite_ok(x0)) return sr_error(ENG_ERR_BADARG);
    double x=x0,ea=HUGE_VAL;
    for(int i=1;i<=max_iter;i++){
        double fx; EngError e=f_at(expr,x,base,&fx); if(e) return sr_error(e);
        if(fabs(fx)<1e-14){ SolverResult r={x,i-1,0,1,ENG_OK}; return r; }
        double h=1e-5*fmax(1.0,fabs(x)), dfx;
        e=eng_derivative(expr,x,h,base,&dfx); if(e) return sr_error(e);
        if(fabs(dfx)<1e-14) return sr_error(ENG_ERR_DERIVATIVE);
        double xn=x-fx/dfx;
        if(!finite_ok(xn)) return sr_error(ENG_ERR_OVERFLOW);
        ea=ea_pct(xn,x);
        double fn; e=f_at(expr,xn,base,&fn); if(e) return sr_error(e);
        if(fabs(fn)<1e-12||ea<=tol_pct){ SolverResult r={xn,i,ea,1,ENG_OK}; return r; }
        x=xn;
    }