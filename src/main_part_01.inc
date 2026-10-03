}
static void result_number(const char *title,const char *label,double z){
    char b[64];snprintf(b,sizeof b,"%s = %.12g",label,z);msg(title,b,NULL);
}

static void error_box(EngError e){msg("ERROR",eng_error_string(e),"No result was accepted.");}

/* Expression builder: numbers and operators from keypad; F1..F6 add common tokens.
   x and y are available from function keys, so numerical-method expressions do not require alpha input. */
static int read_expr(const char *title,char *buf,size_t cap,int allow_y){
    int len=strlen(buf);
    for(;;){
        header(title); dtext(8,32,C_BLACK,"Expression:"); drect_border(6,52,389,86,C_NONE,1,C_BLACK); dtext(10,60,C_BLACK,buf);
        dtext(6,112,C_BLACK,"F1:x  F2:sin( F3:cos( F4:sqrt(");
        dtext(6,134,C_BLACK,allow_y?"F5:y  F6:ln(":"F5:pi F6:ln(");
        footer("EXIT: Cancel","EXE: Accept");dupdate();
        key_event_t e=getkey();int k=e.key;const char *tok=NULL;char c=0;
        if(k==KEY_EXIT)return 0;
        if(k==KEY_EXE)return len>0;
        if(k==KEY_DEL&&len>0){buf[--len]=0;continue;}
        c=digit_key(k);
        if(!c&&k==KEY_DOT)c='.';
        else if(!c&&k==KEY_ADD)c='+';
        else if(!c&&(k==KEY_SUB||k==KEY_NEG))c='-';
        else if(!c&&k==KEY_MUL)c='*';
        else if(!c&&k==KEY_DIV)c='/';
        else if(!c&&k==KEY_POWER)c='^';
        else if(!c&&k==KEY_LEFTP)c='(';
        else if(!c&&k==KEY_RIGHTP)c=')';
        else if(k==KEY_F1)tok="x";
        else if(k==KEY_F2)tok="sin(";
        else if(k==KEY_F3)tok="cos(";
        else if(k==KEY_F4)tok="sqrt(";
        else if(k==KEY_F5)tok=allow_y?"y":"pi";
        else if(k==KEY_F6)tok="ln(";
        if(c&&len<(int)cap-1){buf[len++]=c;buf[len]=0;}
        else if(tok){size_t n=strlen(tok);if(len+(int)n<(int)cap){memcpy(buf+len,tok,n+1);len+=n;}}
    }
}

static void runmat(void){
    char expr[96]=""; if(!read_expr("RUN-MAT",expr,sizeof expr,0))return;
    double z;EngError e=eng_eval(expr,&gvars,&z);if(e){error_box(e);return;}gvars.ans=z;result_number("RUN-MAT","Ans",z);
}
static void calc_menu(void){
    const char *it[]={"Derivative f'(x)","Second derivative","Definite integral","Numerical limit"};
    for(;;){int c=menu("CALCULUS",it,4);if(c<0)return;char ex[96]="";if(!read_expr("CALCULUS",ex,sizeof ex,0))continue;double x,a,b,h=1e-4,z;EngError e;
        if(c==0){if(!read_number("Derivative","x =",0,&x))continue;if(!read_number("Derivative","h =",1e-4,&h))continue;e=eng_derivative(ex,x,h,&gvars,&z);}
        else if(c==1){if(!read_number("Second Derivative","x =",0,&x))continue;if(!read_number("Second Derivative","h =",1e-3,&h))continue;e=eng_second_derivative(ex,x,h,&gvars,&z);}
        else if(c==2){double nd;if(!read_number("Integral","a =",0,&a)||!read_number("Integral","b =",1,&b)||!read_number("Integral","even n =",100,&nd))continue;e=eng_simpson13(ex,a,b,(int)nd,&gvars,&z);}
        else {if(!read_number("Limit","x0 =",0,&x)||!read_number("Limit","initial h =",1e-3,&h))continue;e=eng_limit_symmetric(ex,x,h,&gvars,&z);}
        if(e)error_box(e);else result_number("CALCULUS","Result",z);
    }
}
static void root_menu(void){
    const char *it[]={"Bisection","False Position","Newton-Raphson","Secant"};
    for(;;){int c=menu("EQUATION SOLVER",it,4);if(c<0)return;char ex[96]="";if(!read_expr("f(x)=0",ex,sizeof ex,0))continue;double a=0,b=1;SolverResult r;
        if(c<=1){if(!read_number("Root Solver","xl =",0,&a)||!read_number("Root Solver","xu =",1,&b))continue;r=(c==0)?eng_bisection(ex,a,b,g_tol,g_max_iter,&gvars):eng_false_position(ex,a,b,g_tol,g_max_iter,&gvars);}
        else if(c==2){if(!read_number("Newton","x0 =",1,&a))continue;r=eng_newton(ex,a,g_tol,g_max_iter,&gvars);}
        else {if(!read_number("Secant","x0 =",0,&a)||!read_number("Secant","x1 =",1,&b))continue;r=eng_secant(ex,a,b,g_tol,g_max_iter,&gvars);}
        if(r.error&&r.error!=ENG_ERR_NOCONVERGE){error_box(r.error);continue;} char line[64];snprintf(line,sizeof line,"x = %.12g",r.value);msg(r.converged?"ROOT":"NOT CONVERGED",line,r.converged?"Converged safely":"Increase iterations/check guess");
    }
}
static void result_three(const char *title,const char *la,double a,const char *lb,double b,const char *lc,double c){
    char x[64],y[64],z[64];
    snprintf(x,sizeof x,"%s = %.12g",la,a); snprintf(y,sizeof y,"%s = %.12g",lb,b); snprintf(z,sizeof z,"%s = %.12g",lc,c);
    header(title); dtext(12,48,C_BLACK,x); dtext(12,78,C_BLACK,y); dtext(12,108,C_BLACK,z); footer("EXIT: Back","EXE: OK"); dupdate();
    for(;;){key_event_t e=getkey();if(e.key==KEY_EXIT||e.key==KEY_EXE)break;}
}
static void cas_simplify(void){
    char ex[96]="",out[192];
    if(!read_expr("CAS SIMPLIFY",ex,sizeof ex,0))return;
    EngError e=eng_simplify_polynomial(ex,out,sizeof out);
    if(e){msg("CAS","Polynomial forms only","Use x,+,-,*,^, parentheses");return;}
    msg("SIMPLIFIED",out,NULL);
}
static void taylor_menu(void){
    double x,y,h,ord; char d1[96]="",d2[96]="",d3[96]="",d4[96]="";
    if(!read_number("TAYLOR ODE","x_n =",0,&x)||!read_number("TAYLOR ODE","y_n =",1,&y)||!read_number("TAYLOR ODE","h =",0.1,&h)||!read_number("TAYLOR ODE","order 1..4 =",2,&ord))return;
    int o=(int)ord;if(o<1||o>4||ord!=o){error_box(ENG_ERR_BADARG);return;}