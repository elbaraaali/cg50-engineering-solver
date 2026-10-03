#include "engine.h"
#include <math.h>
#include <stdlib.h>
#include <string.h>
#include <ctype.h>
#include <float.h>
#include <stdio.h>

#ifndef M_PI
#define M_PI 3.14159265358979323846264338327950288
#endif
#ifndef M_E
#define M_E 2.71828182845904523536028747135266250
#endif

typedef struct {
    const char *p;
    EvalVars vars;
    EngError err;
} Parser;

static void skip_ws(Parser *q){ while(isspace((unsigned char)*q->p)) q->p++; }
static int finite_ok(double z){ return isfinite(z); }
static double parse_expr(Parser *q);

static double angle_in(Parser *q,double z){ return q->vars.angle==ANGLE_DEG ? z*M_PI/180.0 : z; }
static double angle_out(Parser *q,double z){ return q->vars.angle==ANGLE_DEG ? z*180.0/M_PI : z; }

static int ident(Parser *q,char *buf,size_t n){
    skip_ws(q);
    if(!(isalpha((unsigned char)*q->p)||*q->p=='_')) return 0;
    size_t i=0;
    while(isalnum((unsigned char)*q->p)||*q->p=='_'){
        if(i+1<n) buf[i++]=(char)tolower((unsigned char)*q->p);
        q->p++;
    }
    buf[i]=0;
    return 1;
}
static double parse_primary(Parser *q){
    skip_ws(q);
    if(q->err) return 0.0;
    if(*q->p=='('){
        q->p++;
        double z=parse_expr(q);
        skip_ws(q);
        if(*q->p!=')'){ q->err=ENG_ERR_PARSE; return 0.0; }
        q->p++;
        return z;
    }
    if(isdigit((unsigned char)*q->p)||*q->p=='.'){
        char *end=NULL;
        double z=strtod(q->p,&end);
        if(end==q->p){ q->err=ENG_ERR_PARSE; return 0.0; }
        q->p=end;
        if(!finite_ok(z)) q->err=ENG_ERR_OVERFLOW;
        return z;
    }
    char name[24];
    if(ident(q,name,sizeof name)){
        if(strcmp(name,"x")==0) return q->vars.x;
        if(strcmp(name,"y")==0) return q->vars.y;
        if(strcmp(name,"ans")==0||strcmp(name,"last")==0) return q->vars.ans;
        if(strcmp(name,"pi")==0) return M_PI;
        if(strcmp(name,"e")==0) return M_E;
        skip_ws(q);
        if(*q->p!='('){ q->err=ENG_ERR_PARSE; return 0.0; }
        q->p++;
        double a=parse_expr(q);
        skip_ws(q);
        if(*q->p!=')'){ q->err=ENG_ERR_PARSE; return 0.0; }
        q->p++;
        if(q->err) return 0.0;
        double z=0.0;
        if(strcmp(name,"sin")==0) z=sin(angle_in(q,a));
        else if(strcmp(name,"cos")==0) z=cos(angle_in(q,a));
        else if(strcmp(name,"tan")==0) z=tan(angle_in(q,a));
        else if(strcmp(name,"asin")==0||strcmp(name,"arcsin")==0){
            if(a<-1.0||a>1.0){ q->err=ENG_ERR_DOMAIN; return 0.0; }
            z=angle_out(q,asin(a));
        }
        else if(strcmp(name,"acos")==0||strcmp(name,"arccos")==0){
            if(a<-1.0||a>1.0){ q->err=ENG_ERR_DOMAIN; return 0.0; }
            z=angle_out(q,acos(a));
        }
        else if(strcmp(name,"atan")==0||strcmp(name,"arctan")==0) z=angle_out(q,atan(a));
        else if(strcmp(name,"sqrt")==0){
            if(a<0.0){ q->err=ENG_ERR_DOMAIN; return 0.0; }
            z=sqrt(a);
        }
        else if(strcmp(name,"ln")==0||strcmp(name,"log")==0){
            if(a<=0.0){ q->err=ENG_ERR_DOMAIN; return 0.0; }
            z=log(a);
        }
        else if(strcmp(name,"log10")==0){
            if(a<=0.0){ q->err=ENG_ERR_DOMAIN; return 0.0; }
            z=log10(a);
        }
        else if(strcmp(name,"exp")==0) z=exp(a);
        else if(strcmp(name,"abs")==0) z=fabs(a);