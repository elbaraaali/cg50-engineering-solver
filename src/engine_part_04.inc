double dyn_tangential_accel(double alpha,double r){return alpha*r;}
double dyn_ke(double m,double v){return 0.5*m*v*v;}
double dyn_pe(double m,double g,double h){return m*g*h;}
double dyn_spring_pe(double k,double x){return 0.5*k*x*x;}
double dyn_momentum(double m,double v){return m*v;}
double dyn_impulse(double F,double dt){return F*dt;}
EngError dyn_restitution(double u1,double u2,double v1,double v2,double *e){if(!e||fabs(u1-u2)<1e-15)return ENG_ERR_DIVZERO;*e=fabs((v2-v1)/(u1-u2));return ENG_OK;}

/* ---------- Extended numerical / symbolic / dynamics core (v0.4) ---------- */
EngError eng_integral_simpson_auto(const char *expr,double a,double b,const EvalVars *base,double *out){
    if(!expr||!out||!isfinite(a)||!isfinite(b)) return ENG_ERR_BADARG;
    if(a==b){*out=0.0;return ENG_OK;}
    double prev=0.0,cur=0.0;
    int n=32;
    EngError e=eng_simpson13(expr,a,b,n,base,&prev); if(e)return e;
    for(int k=0;k<8;k++){
        n*=2; e=eng_simpson13(expr,a,b,n,base,&cur); if(e)return e;
        if(fabs(cur-prev)<=1e-10*(1.0+fabs(cur))){*out=cur;return ENG_OK;}
        prev=cur;
    }
    *out=prev; return ENG_OK;
}

EngError eng_taylor_ode_step(const char *d1,const char *d2,const char *d3,const char *d4,int order,double x,double y,double h,const EvalVars *base,double *y_next){
    if(!d1||!y_next||order<1||order>4||!isfinite(x)||!isfinite(y)||!isfinite(h)) return ENG_ERR_BADARG;
    EvalVars v={0,0,0,ANGLE_RAD}; if(base)v=*base; v.x=x;v.y=y;
    double q1=0,q2=0,q3=0,q4=0; EngError e=eng_eval(d1,&v,&q1); if(e)return e;
    double z=y+h*q1;
    if(order>=2){ if(!d2)return ENG_ERR_BADARG; e=eng_eval(d2,&v,&q2);if(e)return e; z+=h*h*q2/2.0; }
    if(order>=3){ if(!d3)return ENG_ERR_BADARG; e=eng_eval(d3,&v,&q3);if(e)return e; z+=h*h*h*q3/6.0; }
    if(order>=4){ if(!d4)return ENG_ERR_BADARG; e=eng_eval(d4,&v,&q4);if(e)return e; z+=h*h*h*h*q4/24.0; }
    if(!isfinite(z)) return ENG_ERR_OVERFLOW;
    *y_next=z;
    return ENG_OK;
}

/* Small, deliberately bounded polynomial CAS: constants and x; + - * ^; integer powers 0..12.
   This is used for reliable simplify/expand such as x+x -> 2*x. */
#define POLY_DEG 12
typedef struct { double c[POLY_DEG+1]; int ok; } Poly;
typedef struct { const char *p; EngError err; } PParser;
static void pskip(PParser *q){while(isspace((unsigned char)*q->p))q->p++;}
static Poly pzero(void){Poly p={{0},1};return p;}
static Poly padd(Poly a,Poly b,double s){for(int i=0;i<=POLY_DEG;i++)a.c[i]+=s*b.c[i];return a;}
static Poly pmul(Poly a,Poly b,PParser *q){Poly r=pzero();for(int i=0;i<=POLY_DEG;i++)for(int j=0;j<=POLY_DEG;j++)if(a.c[i]!=0&&b.c[j]!=0){if(i+j>POLY_DEG){q->err=ENG_ERR_BADARG;return r;}r.c[i+j]+=a.c[i]*b.c[j];}return r;}
static Poly pexpr(PParser *q);
static Poly pprimary(PParser *q){
    pskip(q); Poly r=pzero();
    if(*q->p=='('){q->p++;r=pexpr(q);pskip(q);if(*q->p!=')')q->err=ENG_ERR_PARSE;else q->p++;return r;}
    if(*q->p=='x'||*q->p=='X'){r.c[1]=1;q->p++;return r;}
    if(isdigit((unsigned char)*q->p)||*q->p=='.'){
        char *e=NULL;double z=strtod(q->p,&e);if(e==q->p||!isfinite(z)){q->err=ENG_ERR_PARSE;return r;}q->p=e;r.c[0]=z;return r;
    }
    q->err=ENG_ERR_PARSE;return r;
}
static Poly ppower(PParser *q){
    Poly a=pprimary(q);if(q->err)return a;pskip(q);if(*q->p!='^')return a;q->p++;pskip(q);
    char *e=NULL;long n=strtol(q->p,&e,10);if(e==q->p||n<0||n>POLY_DEG){q->err=ENG_ERR_BADARG;return a;}q->p=e;
    Poly r=pzero();r.c[0]=1;for(long i=0;i<n;i++)r=pmul(r,a,q);return r;
}
static Poly punary(PParser *q){pskip(q);if(*q->p=='+'){q->p++;return punary(q);}if(*q->p=='-'){q->p++;Poly a=punary(q);for(int i=0;i<=POLY_DEG;i++)a.c[i]=-a.c[i];return a;}return ppower(q);}
static Poly pterm(PParser *q){Poly a=punary(q);for(;;){pskip(q);if(*q->p!='*')break;q->p++;Poly b=punary(q);a=pmul(a,b,q);if(q->err)return a;}return a;}
static Poly pexpr(PParser *q){Poly a=pterm(q);for(;;){pskip(q);char op=*q->p;if(op!='+'&&op!='-')break;q->p++;Poly b=pterm(q);a=padd(a,b,op=='+'?1.0:-1.0);}return a;}
static int app(char *out,size_t cap,size_t *n,const char *s){size_t m=strlen(s);if(*n+m+1>cap)return 0;memcpy(out+*n,s,m);*n+=m;out[*n]=0;return 1;}
EngError eng_simplify_polynomial(const char *expr,char *out,size_t cap){
    if(!expr||!out||cap<2) return ENG_ERR_BADARG;
    PParser q={expr,ENG_OK};
    Poly p=pexpr(&q);
    pskip(&q);
    if(q.err) return q.err;
    if(*q.p) return ENG_ERR_PARSE;
    size_t n=0;out[0]=0;int terms=0;
    for(int d=POLY_DEG;d>=0;d--){double c=p.c[d];if(fabs(c)<1e-12)continue;char b[48];double ac=fabs(c);
        if(terms){if(!app(out,cap,&n,c<0?" - ":" + "))return ENG_ERR_OVERFLOW;}else if(c<0){if(!app(out,cap,&n,"-"))return ENG_ERR_OVERFLOW;}
        int omit=(d>0&&fabs(ac-1.0)<1e-12);
        if(!omit){snprintf(b,sizeof b,"%.12g",ac);if(!app(out,cap,&n,b))return ENG_ERR_OVERFLOW;if(d>0&&!app(out,cap,&n,"*"))return ENG_ERR_OVERFLOW;}
        if(d==1){if(!app(out,cap,&n,"x"))return ENG_ERR_OVERFLOW;}else if(d>1){snprintf(b,sizeof b,"x^%d",d);if(!app(out,cap,&n,b))return ENG_ERR_OVERFLOW;}
        terms++;
    }
    if(!terms) app(out,cap,&n,"0");
    return ENG_OK;
}

EngError dyn_u_from_vat(double v,double a,double t,double *u){if(!u)return ENG_ERR_BADARG;*u=v-a*t;return isfinite(*u)?ENG_OK:ENG_ERR_OVERFLOW;}
EngError dyn_a_from_uvt(double u,double v,double t,double *a){if(!a)return ENG_ERR_BADARG;if(fabs(t)<=DBL_MIN)return ENG_ERR_DIVZERO;*a=(v-u)/t;return isfinite(*a)?ENG_OK:ENG_ERR_OVERFLOW;}
EngError dyn_t_from_uva(double u,double v,double a,double *t){if(!t)return ENG_ERR_BADARG;if(fabs(a)<=DBL_MIN)return ENG_ERR_DIVZERO;*t=(v-u)/a;return isfinite(*t)?ENG_OK:ENG_ERR_OVERFLOW;}
double dyn_s_from_uvt(double u,double v,double t){return 0.5*(u+v)*t;}
double dyn_s_from_vat(double v,double a,double t){return v*t-0.5*a*t*t;}
EngError dyn_a_from_uvs(double u,double v,double s,double *a){if(!a)return ENG_ERR_BADARG;if(fabs(s)<=DBL_MIN)return ENG_ERR_DIVZERO;*a=(v*v-u*u)/(2*s);return isfinite(*a)?ENG_OK:ENG_ERR_OVERFLOW;}
EngError dyn_t_from_us(double u,double a,double s,double *t1,double *t2,int *roots){
    if(!t1||!t2||!roots)return ENG_ERR_BADARG;
    if(fabs(a)<=DBL_MIN){if(fabs(u)<=DBL_MIN)return ENG_ERR_DIVZERO;*t1=s/u;*t2=*t1;*roots=1;return ENG_OK;}
    double D=u*u+2*a*s;if(D<0)return ENG_ERR_DOMAIN;double r=sqrt(D);*t1=(-u+r)/a;*t2=(-u-r)/a;*roots=(r==0?1:2);return ENG_OK;
}

EngError dyn_from_position(const char *s_expr,double t,const EvalVars *base,double *s,double *v,double *a){
    if(!s_expr||!s||!v||!a) return ENG_ERR_BADARG;
    EvalVars q={0,0,0,ANGLE_RAD};
    if(base) q=*base;
    q.x=t;