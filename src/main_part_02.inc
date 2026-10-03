    if(!read_expr("Enter y prime =",d1,sizeof d1,1))return;
    if(o>=2&&!read_expr("Enter y second =",d2,sizeof d2,1))return;
    if(o>=3&&!read_expr("Enter y third =",d3,sizeof d3,1))return;
    if(o>=4&&!read_expr("Enter y fourth =",d4,sizeof d4,1))return;
    double yn;EngError e=eng_taylor_ode_step(d1,o>=2?d2:NULL,o>=3?d3:NULL,o>=4?d4:NULL,o,x,y,h,&gvars,&yn);
    if(e)error_box(e);else result_number("TAYLOR ODE","y next",yn);
}
static void numeric_menu(void){
    const char *it[]={"Root methods (Bisect/Secant/...)","Taylor ODE Method","Trapezoidal Rule","Simpson 1/3","Euler ODE","Heun ODE","RK4 ODE"};
    for(;;){int c=menu("NUMERICAL METHODS",it,7);if(c<0)return;if(c==0){root_menu();continue;}if(c==1){taylor_menu();continue;}
        char ex[96]="";int allow_y=c>=4;if(!read_expr("NUMERICAL",ex,sizeof ex,allow_y))continue;double a,b,h,nn,z;
        if(c==2||c==3){if(!read_number("Integration","a =",0,&a)||!read_number("Integration","b =",1,&b)||!read_number("Integration","n =",100,&nn))continue;EngError e=(c==2)?eng_trapezoid(ex,a,b,(int)nn,&gvars,&z):eng_simpson13(ex,a,b,(int)nn,&gvars,&z);if(e)error_box(e);else result_number("NUMERICAL","Integral",z);}
        else {double x0,y0;if(!read_number("ODE","x0 =",0,&x0)||!read_number("ODE","y0 =",1,&y0)||!read_number("ODE","h =",0.1,&h)||!read_number("ODE","steps =",10,&nn))continue;ODEResult o=c==4?eng_euler(ex,x0,y0,h,(int)nn,&gvars):(c==5?eng_heun(ex,x0,y0,h,(int)nn,&gvars):eng_rk4(ex,x0,y0,h,(int)nn,&gvars));if(o.error)error_box(o.error);else result_number("ODE","y final",o.y_final);}
    }
}
static void dyn_constant_menu(void){
    const char *it[]={"v from u,a,t","u from v,a,t","a from u,v,t","t from u,v,a","s from u,a,t","s from u,v,t","s from v,a,t","v magnitude from u,a,s","a from u,v,s","t root(s) from u,a,s"};
    for(;;){int c=menu("CONST ACCEL",it,10);if(c<0)return;double u=0,v=0,a=0,t=0,ds=0,z=0,t2=0;int roots=0;EngError e=ENG_OK;
        if(c==0){if(!read_number("KINEMATICS","u =",0,&u)||!read_number("KINEMATICS","a =",0,&a)||!read_number("KINEMATICS","t =",0,&t))continue;z=dyn_v(u,a,t);}
        else if(c==1){if(!read_number("KINEMATICS","v =",0,&v)||!read_number("KINEMATICS","a =",0,&a)||!read_number("KINEMATICS","t =",0,&t))continue;e=dyn_u_from_vat(v,a,t,&z);}
        else if(c==2){if(!read_number("KINEMATICS","u =",0,&u)||!read_number("KINEMATICS","v =",0,&v)||!read_number("KINEMATICS","t =",1,&t))continue;e=dyn_a_from_uvt(u,v,t,&z);}
        else if(c==3){if(!read_number("KINEMATICS","u =",0,&u)||!read_number("KINEMATICS","v =",0,&v)||!read_number("KINEMATICS","a =",1,&a))continue;e=dyn_t_from_uva(u,v,a,&z);}
        else if(c==4){if(!read_number("KINEMATICS","u =",0,&u)||!read_number("KINEMATICS","a =",0,&a)||!read_number("KINEMATICS","t =",0,&t))continue;z=dyn_s(u,a,t);}
        else if(c==5){if(!read_number("KINEMATICS","u =",0,&u)||!read_number("KINEMATICS","v =",0,&v)||!read_number("KINEMATICS","t =",0,&t))continue;z=dyn_s_from_uvt(u,v,t);}
        else if(c==6){if(!read_number("KINEMATICS","v =",0,&v)||!read_number("KINEMATICS","a =",0,&a)||!read_number("KINEMATICS","t =",0,&t))continue;z=dyn_s_from_vat(v,a,t);}
        else if(c==7){if(!read_number("KINEMATICS","u =",0,&u)||!read_number("KINEMATICS","a =",0,&a)||!read_number("KINEMATICS","s =",0,&ds))continue;e=dyn_v_from_s(u,a,ds,&z);}
        else if(c==8){if(!read_number("KINEMATICS","u =",0,&u)||!read_number("KINEMATICS","v =",0,&v)||!read_number("KINEMATICS","s =",1,&ds))continue;e=dyn_a_from_uvs(u,v,ds,&z);}
        else {if(!read_number("KINEMATICS","u =",0,&u)||!read_number("KINEMATICS","a =",0,&a)||!read_number("KINEMATICS","s =",0,&ds))continue;e=dyn_t_from_us(u,a,ds,&z,&t2,&roots);if(!e&&roots==2){result_three("TIME ROOTS","t1",z,"t2",t2,"roots",2);continue;}}
        if(e)error_box(e);else result_number("KINEMATICS","Result",z);
    }
}
static void dyn_function_menu(void){
    const char *it[]={"Given position s(t)","Given velocity v(t)","Given acceleration a(t)"};
    for(;;){int c=menu("MOTION FUNCTIONS",it,3);if(c<0)return;char ex[96]="";if(!read_expr(c==0?"s(t), use x=t":(c==1?"v(t), use x=t":"a(t), use x=t"),ex,sizeof ex,0))continue;double t0=0,t1=0,s0=0,v0=0,p=0,v=0,a=0;EngError e;
        if(c==0){if(!read_number("s(t)","time t =",0,&t1))continue;e=dyn_from_position(ex,t1,&gvars,&p,&v,&a);if(e)error_box(e);else result_three("FROM POSITION","s",p,"v",v,"a",a);}
        else if(c==1){if(!read_number("v(t)","start t0 =",0,&t0)||!read_number("v(t)","target t =",1,&t1)||!read_number("v(t)","initial s0 =",0,&s0))continue;e=dyn_from_velocity(ex,t0,t1,s0,&gvars,&v,&a,&p);if(e)error_box(e);else result_three("FROM VELOCITY","s",p,"v",v,"a",a);}
        else {if(!read_number("a(t)","start t0 =",0,&t0)||!read_number("a(t)","target t =",1,&t1)||!read_number("a(t)","initial v0 =",0,&v0)||!read_number("a(t)","initial s0 =",0,&s0))continue;e=dyn_from_acceleration(ex,t0,t1,v0,s0,&gvars,&a,&v,&p);if(e)error_box(e);else result_three("FROM ACCEL","s",p,"v",v,"a",a);}
    }
}
static void dyn_projectile_menu(void){
    const char *it[]={"Range","Time of flight","Maximum height"};for(;;){int c=menu("PROJECTILE",it,3);if(c<0)return;double v,th,g;if(!read_number("PROJECTILE","v0 =",0,&v)||!read_number("PROJECTILE","theta deg =",45,&th)||!read_number("PROJECTILE","g =",9.81,&g))continue;if(v<0||g<=0){error_box(ENG_ERR_BADARG);continue;}double z=c==0?dyn_projectile_range(v,th,g):(c==1?dyn_projectile_time(v,th,g):dyn_projectile_hmax(v,th,g));result_number("PROJECTILE","Result",z);}
}
static void dyn_energy_menu(void){
    const char *it[]={"Kinetic energy","Grav potential energy","Spring energy","Momentum","Impulse","Coefficient restitution"};for(;;){int c=menu("ENERGY / MOMENTUM",it,6);if(c<0)return;double a,b,d,e,z;EngError er=ENG_OK;
        if(c==0){if(!read_number("KE","m =",1,&a)||!read_number("KE","v =",0,&b))continue;z=dyn_ke(a,b);}
        else if(c==1){if(!read_number("PE","m =",1,&a)||!read_number("PE","g =",9.81,&b)||!read_number("PE","h =",0,&d))continue;z=dyn_pe(a,b,d);}
        else if(c==2){if(!read_number("SPRING","k =",1,&a)||!read_number("SPRING","x =",0,&b))continue;z=dyn_spring_pe(a,b);}
        else if(c==3){if(!read_number("MOMENTUM","m =",1,&a)||!read_number("MOMENTUM","v =",0,&b))continue;z=dyn_momentum(a,b);}
        else if(c==4){if(!read_number("IMPULSE","F average =",0,&a)||!read_number("IMPULSE","dt =",0,&b))continue;z=dyn_impulse(a,b);}
        else {if(!read_number("RESTITUTION","u1 =",0,&a)||!read_number("RESTITUTION","u2 =",0,&b)||!read_number("RESTITUTION","v1 =",0,&d)||!read_number("RESTITUTION","v2 =",0,&e))continue;er=dyn_restitution(a,b,d,e,&z);}
        if(er)error_box(er);else result_number("DYNAMICS","Result",z);
    }
}
static void dynamics_menu(void){
    const char *it[]={"Position/Velocity/Acceleration","Constant Acceleration Solver","Projectile Motion","Normal acceleration v^2/rho","Tangential accel alpha*r","Energy / Momentum / Impact"};
    for(;;){int c=menu("DYNAMICS",it,6);if(c<0)return;if(c==0)dyn_function_menu();else if(c==1)dyn_constant_menu();else if(c==2)dyn_projectile_menu();else if(c==3){double v,r,z;if(read_number("NORMAL ACCEL","v =",0,&v)&&read_number("NORMAL ACCEL","rho =",1,&r)){EngError e=dyn_normal_accel(v,r,&z);if(e)error_box(e);else result_number("NORMAL ACCEL","a_n",z);}}else if(c==4){double a,r;if(read_number("TANGENTIAL","alpha =",0,&a)&&read_number("TANGENTIAL","r =",1,&r))result_number("TANGENTIAL","a_t",dyn_tangential_accel(a,r));}else dyn_energy_menu();}
}

static int edit_matrix(const char *name,EngMatrix *m){
    double rr=2,cc=2;
    if(!read_number(name,"Rows 1..6 =",m->rows?m->rows:2,&rr))return 0;
    if(!read_number(name,"Cols 1..6 =",m->cols?m->cols:2,&cc))return 0;
    int r=(int)rr,c=(int)cc;
    if(r<1||r>ENG_MAX_MAT||c<1||c>ENG_MAX_MAT||rr!=r||cc!=c){error_box(ENG_ERR_DIM);return 0;}
    EngMatrix temp={0};temp.rows=r;temp.cols=c;
    for(int i=0;i<r;i++)for(int j=0;j<c;j++){
        char label[24];snprintf(label,sizeof label,"[%d,%d] =",i+1,j+1);
        double def=(m->rows==r&&m->cols==c)?m->v[i][j]:0.0;
        if(!read_number(name,label,def,&temp.v[i][j]))return 0;
    }
    *m=temp;return 1;
}
static void show_matrix(const char *title,const EngMatrix *m){
    if(!m||m->rows<1||m->cols<1){msg(title,"Matrix is not defined.",NULL);return;}
    for(;;){
        header(title);
        for(int i=0;i<m->rows;i++){
            char line[80],num[16];line[0]=0;strcat(line,"[ ");
            for(int j=0;j<m->cols;j++){
                snprintf(num,sizeof num,"%.5g",m->v[i][j]);strncat(line,num,sizeof(line)-strlen(line)-1);