    SolverResult r={x,max_iter,ea,0,ENG_ERR_NOCONVERGE}; return r;
}
EngError eng_trapezoid(const char *expr,double a,double b,int n,const EvalVars *base,double *out){
    if(!expr||!out||n<1||a==b) return ENG_ERR_BADARG;
    double h=(b-a)/n,fa,fb;
    EngError e=f_at(expr,a,base,&fa); if(e) return e;
    e=f_at(expr,b,base,&fb); if(e) return e;
    double sum=0.5*(fa+fb);
    for(int i=1;i<n;i++){ double f; e=f_at(expr,a+i*h,base,&f); if(e) return e; sum+=f; }
    *out=h*sum; return finite_ok(*out)?ENG_OK:ENG_ERR_OVERFLOW;
}
EngError eng_simpson13(const char *expr,double a,double b,int n,const EvalVars *base,double *out){
    if(!expr||!out||n<2||(n&1)||a==b) return ENG_ERR_BADARG;
    double h=(b-a)/n,fa,fb;
    EngError e=f_at(expr,a,base,&fa); if(e) return e;
    e=f_at(expr,b,base,&fb); if(e) return e;
    double sum=fa+fb;
    for(int i=1;i<n;i++){ double f; e=f_at(expr,a+i*h,base,&f); if(e) return e; sum+=(i&1?4.0:2.0)*f; }
    *out=h*sum/3.0; return finite_ok(*out)?ENG_OK:ENG_ERR_OVERFLOW;
}
EngError eng_limit_symmetric(const char *expr,double x0,double h,const EvalVars *base,double *out){
    if(!expr||!out||h<=0||!finite_ok(x0)||!finite_ok(h)) return ENG_ERR_BADARG;
    double lv[3],rv[3];
    for(int k=0;k<3;k++){
        double hk=h/pow(2.0,k);
        EngError e=f_at(expr,x0-hk,base,&lv[k]); if(e) return e;
        e=f_at(expr,x0+hk,base,&rv[k]); if(e) return e;
    }
    /* A two-sided limit must approach the same finite number from both sides.
       We also require the last two scales on each side to be stabilizing. */
    double candidate=(lv[2]+rv[2])/2.0;
    double scale=fmax(1.0,fabs(candidate));
    double side_gap=fabs(lv[2]-rv[2]);
    double left_step=fabs(lv[2]-lv[1]);
    double right_step=fabs(rv[2]-rv[1]);
    double threshold=1e-5*scale;
    if(side_gap>threshold||left_step>threshold||right_step>threshold) return ENG_ERR_NOCONVERGE;
    *out=candidate;
    return finite_ok(*out)?ENG_OK:ENG_ERR_OVERFLOW;
}
static EngError fxy(const char *expr,double x,double y,const EvalVars *base,double *out){
    EvalVars v={0,0,0,ANGLE_RAD}; if(base) v=*base; v.x=x; v.y=y; return eng_eval(expr,&v,out);
}
ODEResult eng_euler(const char *fxy_expr,double x,double y,double h,int steps,const EvalVars *base){
    ODEResult r={y,0,ENG_OK}; if(!fxy_expr||steps<1||h==0){r.error=ENG_ERR_BADARG;return r;}
    for(int i=0;i<steps;i++){ double k1; EngError e=fxy(fxy_expr,x,y,base,&k1); if(e){r.error=e;return r;} y+=h*k1; x+=h; r.steps=i+1; if(!finite_ok(y)){r.error=ENG_ERR_OVERFLOW;return r;} }
    r.y_final=y; return r;
}
ODEResult eng_heun(const char *fxy_expr,double x,double y,double h,int steps,const EvalVars *base){
    ODEResult r={y,0,ENG_OK}; if(!fxy_expr||steps<1||h==0){r.error=ENG_ERR_BADARG;return r;}
    for(int i=0;i<steps;i++){ double k1,k2; EngError e=fxy(fxy_expr,x,y,base,&k1); if(e){r.error=e;return r;} e=fxy(fxy_expr,x+h,y+h*k1,base,&k2); if(e){r.error=e;return r;} y+=h*(k1+k2)/2.0; x+=h; r.steps=i+1; if(!finite_ok(y)){r.error=ENG_ERR_OVERFLOW;return r;} }
    r.y_final=y; return r;
}
ODEResult eng_rk4(const char *fxy_expr,double x,double y,double h,int steps,const EvalVars *base){
    ODEResult r={y,0,ENG_OK}; if(!fxy_expr||steps<1||h==0){r.error=ENG_ERR_BADARG;return r;}
    for(int i=0;i<steps;i++){ double k1,k2,k3,k4; EngError e=fxy(fxy_expr,x,y,base,&k1); if(e){r.error=e;return r;} e=fxy(fxy_expr,x+h/2,y+h*k1/2,base,&k2); if(e){r.error=e;return r;} e=fxy(fxy_expr,x+h/2,y+h*k2/2,base,&k3); if(e){r.error=e;return r;} e=fxy(fxy_expr,x+h,y+h*k3,base,&k4); if(e){r.error=e;return r;} y+=h*(k1+2*k2+2*k3+k4)/6.0; x+=h; r.steps=i+1; if(!finite_ok(y)){r.error=ENG_ERR_OVERFLOW;return r;} }
    r.y_final=y; return r;
}

static int mat_ok(const EngMatrix *a){ return a&&a->rows>=1&&a->rows<=ENG_MAX_MAT&&a->cols>=1&&a->cols<=ENG_MAX_MAT; }
EngError eng_mat_det(const EngMatrix *a,double *det){
    if(!mat_ok(a)||!det||a->rows!=a->cols) return ENG_ERR_DIM;
    int n=a->rows; double m[ENG_MAX_MAT][ENG_MAX_MAT]; memcpy(m,a->v,sizeof m); double d=1.0; int sign=1;
    double maxabs=0.0; for(int i=0;i<n;i++)for(int j=0;j<n;j++)maxabs=fmax(maxabs,fabs(m[i][j]));
    double eps=1e-12*fmax(1.0,maxabs);
    for(int k=0;k<n;k++){ int p=k; double best=fabs(m[k][k]); for(int i=k+1;i<n;i++) if(fabs(m[i][k])>best){best=fabs(m[i][k]);p=i;} if(best<eps){*det=0;return ENG_OK;} if(p!=k){for(int j=0;j<n;j++){double t=m[k][j];m[k][j]=m[p][j];m[p][j]=t;}sign=-sign;} double piv=m[k][k]; d*=piv; for(int i=k+1;i<n;i++){double f=m[i][k]/piv; for(int j=k+1;j<n;j++)m[i][j]-=f*m[k][j];}}
    *det=sign*d; return finite_ok(*det)?ENG_OK:ENG_ERR_OVERFLOW;
}
EngError eng_mat_inverse(const EngMatrix *a,EngMatrix *inv){
    if(!mat_ok(a)||!inv||a->rows!=a->cols) return ENG_ERR_DIM;
    int n=a->rows; double aug[ENG_MAX_MAT][2*ENG_MAX_MAT]; double maxabs=0.0;
    for(int i=0;i<n;i++)for(int j=0;j<2*n;j++){aug[i][j]=(j<n?a->v[i][j]:(j-n==i?1.0:0.0));if(j<n)maxabs=fmax(maxabs,fabs(aug[i][j]));}
    double eps=1e-12*fmax(1.0,maxabs);
    for(int k=0;k<n;k++){ int p=k; double best=fabs(aug[k][k]); for(int i=k+1;i<n;i++) if(fabs(aug[i][k])>best){best=fabs(aug[i][k]);p=i;} if(best<eps)return ENG_ERR_SINGULAR; if(p!=k)for(int j=0;j<2*n;j++){double t=aug[k][j];aug[k][j]=aug[p][j];aug[p][j]=t;} double piv=aug[k][k]; for(int j=0;j<2*n;j++)aug[k][j]/=piv; for(int i=0;i<n;i++)if(i!=k){double f=aug[i][k]; for(int j=0;j<2*n;j++)aug[i][j]-=f*aug[k][j];}}
    inv->rows=n;inv->cols=n;for(int i=0;i<n;i++)for(int j=0;j<n;j++)inv->v[i][j]=aug[i][n+j]; return ENG_OK;
}
EngError eng_mat_solve(const EngMatrix *A,const double *b,double *x){
    if(!mat_ok(A)||!b||!x||A->rows!=A->cols) return ENG_ERR_DIM;
    int n=A->rows; double m[ENG_MAX_MAT][ENG_MAX_MAT+1]; double maxabs=0.0;
    for(int i=0;i<n;i++){for(int j=0;j<n;j++){m[i][j]=A->v[i][j];maxabs=fmax(maxabs,fabs(m[i][j]));}m[i][n]=b[i];}
    double eps=1e-12*fmax(1.0,maxabs);
    for(int k=0;k<n;k++){int p=k;double best=fabs(m[k][k]);for(int i=k+1;i<n;i++)if(fabs(m[i][k])>best){best=fabs(m[i][k]);p=i;}if(best<eps)return ENG_ERR_SINGULAR;if(p!=k)for(int j=k;j<=n;j++){double t=m[k][j];m[k][j]=m[p][j];m[p][j]=t;}for(int i=k+1;i<n;i++){double f=m[i][k]/m[k][k];for(int j=k;j<=n;j++)m[i][j]-=f*m[k][j];}}
    for(int i=n-1;i>=0;i--){double s=m[i][n];for(int j=i+1;j<n;j++)s-=m[i][j]*x[j];x[i]=s/m[i][i];}return ENG_OK;
}
EngError eng_mat_add(const EngMatrix *a,const EngMatrix *b,EngMatrix *out){
    if(!mat_ok(a)||!mat_ok(b)||!out||a->rows!=b->rows||a->cols!=b->cols)return ENG_ERR_DIM;
    out->rows=a->rows;out->cols=a->cols;for(int i=0;i<a->rows;i++)for(int j=0;j<a->cols;j++)out->v[i][j]=a->v[i][j]+b->v[i][j];return ENG_OK;
}
EngError eng_mat_mul(const EngMatrix *a,const EngMatrix *b,EngMatrix *out){
    if(!mat_ok(a)||!mat_ok(b)||!out||a->cols!=b->rows)return ENG_ERR_DIM;
    out->rows=a->rows;out->cols=b->cols;for(int i=0;i<out->rows;i++)for(int j=0;j<out->cols;j++){double s=0;for(int k=0;k<a->cols;k++)s+=a->v[i][k]*b->v[k][j];out->v[i][j]=s;}return ENG_OK;
}

double dyn_v(double u,double a,double t){return u+a*t;}
double dyn_s(double u,double a,double t){return u*t+0.5*a*t*t;}
EngError dyn_v_from_s(double u,double a,double s,double *v){if(!v)return ENG_ERR_BADARG;double q=u*u+2*a*s;if(q<0)return ENG_ERR_DOMAIN;*v=sqrt(q);return ENG_OK;}
double dyn_projectile_range(double v,double th,double g){return v*v*sin(2*th*M_PI/180.0)/g;}
double dyn_projectile_time(double v,double th,double g){return 2*v*sin(th*M_PI/180.0)/g;}
double dyn_projectile_hmax(double v,double th,double g){double vy=v*sin(th*M_PI/180.0);return vy*vy/(2*g);}
EngError dyn_normal_accel(double v,double r,double *a_n){if(!a_n||r<=0)return ENG_ERR_BADARG;*a_n=v*v/r;return ENG_OK;}