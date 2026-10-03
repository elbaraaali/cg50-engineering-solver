                if(j<m->cols-1)strncat(line,"  ",sizeof(line)-strlen(line)-1);
            }
            strncat(line," ]",sizeof(line)-strlen(line)-1);
            dtext(8,34+i*24,C_BLACK,line);
        }
        footer("EXIT: Back",NULL);dupdate();
        if(getkey().key==KEY_EXIT)return;
    }
}
static void matrix_menu(void){
    const char *it[]={"Edit Mat A","Edit Mat B","View Mat A","View Mat B","det(A)","inverse(A)","A + B","A * B","Solve A*x = b (B column)"};
    for(;;){
        int c=menu("MATRIX",it,9);if(c<0)return;EngError e;double z;EngMatrix out;
        if(c==0){edit_matrix("MAT A",&gA);continue;}
        if(c==1){edit_matrix("MAT B",&gB);continue;}
        if(c==2){show_matrix("MAT A",&gA);continue;}
        if(c==3){show_matrix("MAT B",&gB);continue;}
        if(c==4){e=eng_mat_det(&gA,&z);if(e)error_box(e);else result_number("MATRIX","det(A)",z);continue;}
        if(c==5){e=eng_mat_inverse(&gA,&out);if(e)error_box(e);else show_matrix("A^-1",&out);continue;}
        if(c==6){e=eng_mat_add(&gA,&gB,&out);if(e)error_box(e);else show_matrix("A+B",&out);continue;}
        if(c==7){e=eng_mat_mul(&gA,&gB,&out);if(e)error_box(e);else show_matrix("A*B",&out);continue;}
        if(gA.rows<1||gA.rows!=gA.cols||gB.rows!=gA.rows||gB.cols!=1){error_box(ENG_ERR_DIM);continue;}
        double b[ENG_MAX_MAT],x[ENG_MAX_MAT];for(int i=0;i<gA.rows;i++){b[i]=gB.v[i][0];x[i]=0;}
        e=eng_mat_solve(&gA,b,x);if(e){error_box(e);continue;}
        EngMatrix sol={0};sol.rows=gA.rows;sol.cols=1;for(int i=0;i<gA.rows;i++)sol.v[i][0]=x[i];show_matrix("SOLUTION x",&sol);
    }
}
static void settings(void){
    const char *it[]={"Angle: DEG/RAD","Tolerance (%)","Max iterations"};
    for(;;){int c=menu("SETTINGS",it,3);if(c<0)return;if(c==0){gvars.angle=gvars.angle==ANGLE_DEG?ANGLE_RAD:ANGLE_DEG;msg("SETTINGS",gvars.angle==ANGLE_DEG?"Angle mode: DEG":"Angle mode: RAD",NULL);}else if(c==1){double z;if(read_number("SETTINGS","Tolerance % =",g_tol,&z)&&z>0)g_tol=z;}else{double z;if(read_number("SETTINGS","Max iterations =",g_max_iter,&z)&&z>=1&&z<=1000)g_max_iter=(int)z;}}
}
int main(void){
    const char *items[]={"RUN-MAT","CAS SIMPLIFY / EXPAND","CALCULUS","EQUATION SOLVER","NUMERICAL METHODS","DYNAMICS","MATRIX","SETTINGS","ABOUT"};
    for(;;){int c=menu("ENGINEERING SOLVER",items,9);if(c<0)return 1;switch(c){case 0:runmat();break;case 1:cas_simplify();break;case 2:calc_menu();break;case 3:root_menu();break;case 4:numeric_menu();break;case 5:dynamics_menu();break;case 6:matrix_menu();break;case 7:settings();break;case 8:msg("ABOUT","CG50 Engineering Solver v0.4","Validated engineering core");break;}}
}
