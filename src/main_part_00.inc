#include <gint/display.h>
#include <gint/keyboard.h>
#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include "engine.h"

static EvalVars gvars={0,0,0,ANGLE_RAD};
static double g_tol=0.001;
static int g_max_iter=50;
static EngMatrix gA={0}, gB={0};

static void header(const char *title){
    dclear(C_WHITE);
    drect(0,0,395,24,C_BLACK);
    dtext(8,5,C_WHITE,title);
}
static void footer(const char *l,const char *r){
    dline(0,198,395,198,C_BLACK);
    if(l)dtext(8,205,C_BLACK,l);
    if(r)dtext(300,205,C_BLACK,r);
}
static void msg(const char *title,const char *a,const char *b){
    header(title); if(a)dtext(12,60,C_BLACK,a); if(b)dtext(12,90,C_BLACK,b); footer("EXIT: Back","EXE: OK"); dupdate();
    for(;;){key_event_t e=getkey(); if(e.key==KEY_EXIT||e.key==KEY_EXE)break;}
}
static int menu(const char *title,const char **items,int n){
    int s=0,top=0;
    for(;;){
        header(title);
        int visible=7;
        if(s<top) top=s;
        if(s>=top+visible) top=s-visible+1;
        for(int i=0;i<visible&&top+i<n;i++){
            int y=34+i*22;
            if(top+i==s){drect(6,y-2,389,y+17,C_BLACK);dtext(12,y,C_WHITE,items[top+i]);}
            else dtext(12,y,C_BLACK,items[top+i]);
        }
        footer("EXIT: Back","EXE: Select");dupdate();
        key_event_t e=getkey();
        if(e.key==KEY_UP&&s>0)s--;
        else if(e.key==KEY_DOWN&&s<n-1)s++;
        else if(e.key==KEY_EXE)return s;
        else if(e.key==KEY_EXIT)return -1;
    }
}

static char digit_key(int k){
    if(k==KEY_0) return '0';
    if(k==KEY_1) return '1';
    if(k==KEY_2) return '2';
    if(k==KEY_3) return '3';
    if(k==KEY_4) return '4';
    if(k==KEY_5) return '5';
    if(k==KEY_6) return '6';
    if(k==KEY_7) return '7';
    if(k==KEY_8) return '8';
    if(k==KEY_9) return '9';
    return 0;
}

/* Conservative numeric editor: avoids ambiguous alpha mappings.
   Expression entry uses a token palette below. */
static int read_number(const char *title,const char *label,double def,double *out){
    char buf[32]; snprintf(buf,sizeof buf,"%.10g",def); int len=strlen(buf);
    for(;;){
        header(title); dtext(12,42,C_BLACK,label); drect_border(10,70,385,100,C_NONE,1,C_BLACK); dtext(16,77,C_BLACK,buf);
        footer("DEL/EXIT","EXE: Accept"); dupdate();
        key_event_t e=getkey(); int k=e.key;
        if(k==KEY_EXIT)return 0;
        if(k==KEY_EXE){char *end;double z=strtod(buf,&end);if(end!=buf&&*end==0){*out=z;return 1;}msg("Input Error","Enter a valid number",NULL);}
        else if(k==KEY_DEL&&len>0){buf[--len]=0;}
        else {
            char c=0;
            c=digit_key(k);
            if(!c&&k==KEY_DOT)c='.';
            else if(!c&&(k==KEY_NEG||k==KEY_SUB))c='-';
            if(c&&len<(int)sizeof(buf)-1){ if(len==1&&buf[0]=='0'&&c!='.')len=0; buf[len++]=c;buf[len]=0; }
        }
    }