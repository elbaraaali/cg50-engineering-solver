    EngError e=eng_eval(s_expr,&q,s);
    if(e) return e;
    double h=1e-4*(1.0+fabs(t));e=eng_derivative(s_expr,t,h,base,v);if(e)return e;e=eng_second_derivative(s_expr,t,sqrt(h)*1e-2+1e-3,base,a);return e;
}
EngError dyn_from_velocity(const char *v_expr,double t0,double t1,double s0,const EvalVars *base,double *v1,double *a1,double *s1){
    if(!v_expr||!v1||!a1||!s1) return ENG_ERR_BADARG;
    EvalVars q={0,0,0,ANGLE_RAD};
    if(base) q=*base;
    q.x=t1;
    EngError e=eng_eval(v_expr,&q,v1);
    if(e) return e;
    double h=1e-4*(1.0+fabs(t1));e=eng_derivative(v_expr,t1,h,base,a1);if(e)return e;double ds;e=eng_integral_simpson_auto(v_expr,t0,t1,base,&ds);if(e)return e;*s1=s0+ds;return ENG_OK;
}
static EngError accel_weighted_integral(const char *a_expr,double t0,double t1,const EvalVars *base,double *out){
    int n=256; if(t0==t1){*out=0;return ENG_OK;}double h=(t1-t0)/n,sum=0;EvalVars q={0,0,0,ANGLE_RAD};if(base)q=*base;
    for(int i=0;i<=n;i++){double t=t0+i*h,av;q.x=t;EngError e=eng_eval(a_expr,&q,&av);if(e)return e;double f=(t1-t)*av;int w=(i==0||i==n)?1:(i%2?4:2);sum+=w*f;}
    *out=sum*h/3.0;return ENG_OK;
}
EngError dyn_from_acceleration(const char *a_expr,double t0,double t1,double v0,double s0,const EvalVars *base,double *a1,double *v1,double *s1){
    if(!a_expr||!a1||!v1||!s1) return ENG_ERR_BADARG;
    EvalVars q={0,0,0,ANGLE_RAD};
    if(base) q=*base;
    q.x=t1;
    EngError e=eng_eval(a_expr,&q,a1);
    if(e) return e;
    double dv;
    e=eng_integral_simpson_auto(a_expr,t0,t1,base,&dv);
    if(e) return e;
    *v1=v0+dv;
    double w;
    e=accel_weighted_integral(a_expr,t0,t1,base,&w);
    if(e) return e;
    *s1=s0+v0*(t1-t0)+w;
    return ENG_OK;
}
