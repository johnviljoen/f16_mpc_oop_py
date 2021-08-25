#include "math.h"
#include<stdio.h>

/*  Merging the nlplant.c (lofi) and nlplant_hifi.c to use
    same equations of motion, navigation equations and use 
    own look-up tables decided by a flag.                   */
    
void atmos(double,double,double*);          /* Used by both */
void accels(double*,double*,double*);       /* Used by both */

#include "lofi_F16_AeroData.c"              /* LOFI Look-up header file*/
#include "hifi_F16_AeroData.c"              /* HIFI Look-up header file*/

void Nlplant(double*,double*,int);

/*########################################*/
/*### Added for mex function in matlab ###*/
/*########################################*/

/*########################################*/
/*########################################*/

void Nlplant(double *xu, double *xdot, int fidelity){ //I modified to have fidelity as a seperate input!

  int fi_flag;

  /* #include f16_constants */
  double g    = 32.17;          /* gravity, ft/s^2 */
  double m    = 636.94;         /* mass, slugs */
  double B    = 30.0;             /* span, ft */
  double S    = 300.0;            /* planform area, ft^2 */
  double cbar = 11.32;          /* mean aero chord, ft */
  double xcgr = 0.35;      /* reference center of gravity as a fraction of cbar */
  double xcg  = 0.25;      /* center of gravity as a fraction of cbar. */

  double Heng = 0.0;              /* turbine momentum along roll axis. */
  double pi   = acos(-1);
  double r2d;                   /* radians to degrees */


/*NasaData        %translated via eq. 2.4-6 on pg 80 of Stevens and Lewis*/

  double Jy  = 55814.0;           /* slug-ft^2 */ 
  double Jxz = 982.0;             /* slug-ft^2 */     
  double Jz  = 63100.0;           /* slug-ft^2 */
  double Jx  = 9496.0;            /* slug-ft^2 */

  double *temp;

  double npos, epos, alt, phi, theta, psi, vt, alpha, beta, P, Q, R;
  double sa, ca, sb, cb, tb, st, ct, tt, sphi, cphi, spsi, cpsi;
  double T, el, ail, rud, dail, drud, lef, dlef;
  double qbar, mach, ps;
  double U, V, W, Udot,Vdot,Wdot;
  double L_tot, M_tot, N_tot, denom;

  double Cx_tot, Cx, delta_Cx_lef, dXdQ, Cxq, delta_Cxq_lef;
  double Cz_tot, Cz, delta_Cz_lef, dZdQ, Czq, delta_Czq_lef;
  double Cm_tot, Cm, eta_el, delta_Cm_lef, dMdQ, Cmq, delta_Cmq_lef, delta_Cm, delta_Cm_ds;
  double Cy_tot, Cy, delta_Cy_lef, dYdail, delta_Cy_r30, dYdR, dYdP;
  double delta_Cy_a20, delta_Cy_a20_lef, Cyr, delta_Cyr_lef, Cyp, delta_Cyp_lef;
  double Cn_tot, Cn, delta_Cn_lef, dNdail, delta_Cn_r30, dNdR, dNdP, delta_Cnbeta;
  double delta_Cn_a20, delta_Cn_a20_lef, Cnr, delta_Cnr_lef, Cnp, delta_Cnp_lef;
  double Cl_tot, Cl, delta_Cl_lef, dLdail, delta_Cl_r30, dLdR, dLdP, delta_Clbeta;
  double delta_Cl_a20, delta_Cl_a20_lef, Clr, delta_Clr_lef, Clp, delta_Clp_lef;

  temp = (double *)malloc(9*sizeof(double));  /*size of 9.1 array*/

  r2d  = 180.0/pi;     /* radians to degrees */

  /* %%%%%%%%%%%%%%%%%%%
           States
     %%%%%%%%%%%%%%%%%%% */


  npos  = xu[0];   /* north position */
  epos  = xu[1];   /* east position */
  alt   = xu[2];   /* altitude */
  phi   = xu[3];   /* orientation angles in rad. */
  theta = xu[4];
  psi   = xu[5];

  vt    = xu[6];     /* total velocity */
  alpha = xu[7]*r2d; /* angle of attack in degrees */
  beta  = xu[8]*r2d; /* sideslip angle in degrees */
  P     = xu[9];     /* Roll Rate --- rolling  moment is Lbar */
  Q     = xu[10];    /* Pitch Rate--- pitching moment is M */
  R     = xu[11];    /* Yaw Rate  --- yawing   moment is N */

  sa    = sin(xu[7]); /* sin(alpha) */
  ca    = cos(xu[7]); /* cos(alpha) */
  sb    = sin(xu[8]); /* sin(beta)  */
  cb    = cos(xu[8]); /* cos(beta)  */
  tb    = tan(xu[8]); /* tan(beta)  */

  st    = sin(theta);
  ct    = cos(theta);
  tt    = tan(theta);
  sphi  = sin(phi);
  cphi  = cos(phi);
  spsi  = sin(psi);
  cpsi  = cos(psi);

  if (vt <= 0.01) {vt = 0.01;}

  /* %%%%%%%%%%%%%%%%%%%
     Control inputs
     %%%%%%%%%%%%%%%%%%% */

  T     = xu[12];   /* thrust */
  el    = xu[13];   /* Elevator setting in degrees. */
  ail   = xu[14];   /* Ailerons mex setting in degrees. */
  rud   = xu[15];   /* Rudder setting in degrees. */
  lef   = xu[16];   /* Leading edge flap setting in degrees */
  
  //fi_flag = xu[17]/1;       /* fi_flag */
  fi_flag = fidelity; // Johns change to pass in integer seperate from numpy float array
    
  /* dail  = ail/20.0;   aileron normalized against max angle */
  /* The aileron was normalized using 20.0 but the NASA report and
     S&L both have 21.5 deg. as maximum deflection. */
  /* As a result... */
  dail  = ail/21.5;
  drud  = rud/30.0;  /* rudder normalized against max angle */
  dlef  = (1 - lef/25.0);  /* leading edge flap normalized against max angle */


  /* %%%%%%%%%%%%%%%%%%
     Atmospheric effects
     sets dynamic pressure and mach number
     %%%%%%%%%%%%%%%%%% */

atmos(alt,vt,temp);
   mach = temp[0];
   qbar = temp[1];
   ps   = temp[2];

/*
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%Dynamics%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
*/

  /* %%%%%%%%%%%%%%%%%%
     Navigation Equations
     %%%%%%%%%%%%%%%%%% */

   U = vt*ca*cb;  /* directional velocities. */
   V = vt*sb;
   W = vt*sa*cb;

/* nposdot */
xdot[0] = U*(ct*cpsi) + 
            V*(sphi*cpsi*st - cphi*spsi) + 
            W*(cphi*st*cpsi + sphi*spsi);

/* eposdot */  
xdot[1] = U*(ct*spsi) + 
            V*(sphi*spsi*st + cphi*cpsi) + 
            W*(cphi*st*spsi - sphi*cpsi);

/* altdot */
xdot[2] = U*st - V*(sphi*ct) - W*(cphi*ct);

  /* %%%%%%%%%%%%%%%%%%%
     Kinematic equations
     %%%%%%%%%%%%%%%%%%% */
/* phidot */
xdot[3] = P + tt*(Q*sphi + R*cphi);


/* theta dot */
xdot[4] = Q*cphi - R*sphi;

/* psidot */
xdot[5] = (Q*sphi + R*cphi)/ct;


/* %%%%%%%%%%%%%%%%%%
        Table lookup
     %%%%%%%%%%%%%%%%%% */

if (fi_flag == 1)          /* HIFI Table */
{
    hifi_C(alpha,beta,el,temp);
        Cx = temp[0];
        Cz = temp[1];
        Cm = temp[2];
        Cy = temp[3];
        Cn = temp[4];
        Cl = temp[5];

    hifi_damping(alpha,temp);
        Cxq = temp[0];
        Cyr = temp[1];
        Cyp = temp[2];
        Czq = temp[3];
        Clr = temp[4];
        Clp = temp[5];
        Cmq = temp[6];
        Cnr = temp[7];
        Cnp = temp[8];

    hifi_C_lef(alpha,beta,temp);
        delta_Cx_lef = temp[0];
        delta_Cz_lef = temp[1];
        delta_Cm_lef = temp[2];
        delta_Cy_lef = temp[3];
        delta_Cn_lef = temp[4];
        delta_Cl_lef = temp[5];

    hifi_damping_lef(alpha,temp);
        delta_Cxq_lef = temp[0];
        delta_Cyr_lef = temp[1];
        delta_Cyp_lef = temp[2];
        delta_Czq_lef = temp[3];
        delta_Clr_lef = temp[4];
        delta_Clp_lef = temp[5];
        delta_Cmq_lef = temp[6];
        delta_Cnr_lef = temp[7];
        delta_Cnp_lef = temp[8];

    hifi_rudder(alpha,beta,temp);
        delta_Cy_r30 = temp[0];
        delta_Cn_r30 = temp[1];
        delta_Cl_r30 = temp[2];

    hifi_ailerons(alpha,beta,temp);
        delta_Cy_a20     = temp[0];
        delta_Cy_a20_lef = temp[1];
        delta_Cn_a20     = temp[2];
        delta_Cn_a20_lef = temp[3];
        delta_Cl_a20     = temp[4];
        delta_Cl_a20_lef = temp[5];

    hifi_other_coeffs(alpha,el,temp);
        delta_Cnbeta = temp[0];
        delta_Clbeta = temp[1];
        delta_Cm     = temp[2];
        eta_el       = temp[3];
        delta_Cm_ds  = 0;        /* ignore deep-stall effect */
    
}

else if (fi_flag == 0)
{     
/* ##############################################
   ##########LOFI Table Look-up #################
   ##############################################*/

/* The lofi model does not include the
   leading edge flap.  All terms multiplied
   dlef have been set to zero but just to 
   be sure we will set it to zero. */
    
    dlef = 0.0;     

    damping(alpha,temp);
        Cxq = temp[0];
        Cyr = temp[1];
        Cyp = temp[2];
        Czq = temp[3];
        Clr = temp[4];
        Clp = temp[5];
        Cmq = temp[6];
        Cnr = temp[7];
        Cnp = temp[8];

    dmomdcon(alpha,beta, temp);
        delta_Cl_a20 = temp[0];     /* Formerly dLda in nlplant.c */
        delta_Cl_r30 = temp[1];     /* Formerly dLdr in nlplant.c */
        delta_Cn_a20 = temp[2];     /* Formerly dNda in nlplant.c */
        delta_Cn_r30 = temp[3];     /* Formerly dNdr in nlplant.c */

    clcn(alpha,beta,temp);
        Cl = temp[0];
        Cn = temp[1];

    cxcm(alpha,el,temp);
        Cx = temp[0];
        Cm = temp[1];

    Cy = -.02*beta + .021*dail + .086*drud;

    cz(alpha,beta,el,temp);
        Cz = temp[0];
/*##################################################
        
        
/*##################################################
  ##  Set all higher order terms of hifi that are ##
  ##  not applicable to lofi equal to zero. ########
  ##################################################*/
     
        delta_Cx_lef    = 0.0;
        delta_Cz_lef    = 0.0;
        delta_Cm_lef    = 0.0;
        delta_Cy_lef    = 0.0;
        delta_Cn_lef    = 0.0;
        delta_Cl_lef    = 0.0;
        delta_Cxq_lef   = 0.0;
        delta_Cyr_lef   = 0.0;
        delta_Cyp_lef   = 0.0;
        delta_Czq_lef   = 0.0;
        delta_Clr_lef   = 0.0;
        delta_Clp_lef   = 0.0;
        delta_Cmq_lef   = 0.0;
        delta_Cnr_lef   = 0.0;
        delta_Cnp_lef   = 0.0;
        delta_Cy_r30    = 0.0;
        delta_Cy_a20    = 0.0;
        delta_Cy_a20_lef= 0.0;
        delta_Cn_a20_lef= 0.0;
        delta_Cl_a20_lef= 0.0;
        delta_Cnbeta    = 0.0;
        delta_Clbeta    = 0.0;
        delta_Cm        = 0.0;
        eta_el          = 1.0;     /* Needs to be one. See equation for Cm_tot*/
        delta_Cm_ds     = 0.0;
                     
/*##################################################
  ##################################################*/ 
}


/* %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
compute Cx_tot, Cz_tot, Cm_tot, Cy_tot, Cn_tot, and Cl_tot
(as on NASA report p37-40)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% */

/* XXXXXXXX Cx_tot XXXXXXXX */

dXdQ = (cbar/(2*vt))*(Cxq + delta_Cxq_lef*dlef);

Cx_tot = Cx + delta_Cx_lef*dlef + dXdQ*Q;

   /* ZZZZZZZZ Cz_tot ZZZZZZZZ */ 

dZdQ = (cbar/(2*vt))*(Czq + delta_Cz_lef*dlef);

Cz_tot = Cz + delta_Cz_lef*dlef + dZdQ*Q;

   /* MMMMMMMM Cm_tot MMMMMMMM */ 

dMdQ = (cbar/(2*vt))*(Cmq + delta_Cmq_lef*dlef);

Cm_tot = Cm*eta_el + Cz_tot*(xcgr-xcg) + delta_Cm_lef*dlef + dMdQ*Q + delta_Cm + delta_Cm_ds;

   /* YYYYYYYY Cy_tot YYYYYYYY */

dYdail = delta_Cy_a20 + delta_Cy_a20_lef*dlef;

dYdR = (B/(2*vt))*(Cyr + delta_Cyr_lef*dlef);

dYdP = (B/(2*vt))*(Cyp + delta_Cyp_lef*dlef);

Cy_tot = Cy + delta_Cy_lef*dlef + dYdail*dail + delta_Cy_r30*drud + dYdR*R + dYdP*P;

   /* NNNNNNNN Cn_tot NNNNNNNN */ 

dNdail = delta_Cn_a20 + delta_Cn_a20_lef*dlef;

dNdR = (B/(2*vt))*(Cnr + delta_Cnr_lef*dlef);

dNdP = (B/(2*vt))*(Cnp + delta_Cnp_lef*dlef);

Cn_tot = Cn + delta_Cn_lef*dlef - Cy_tot*(xcgr-xcg)*(cbar/B) + dNdail*dail + delta_Cn_r30*drud + dNdR*R + dNdP*P + delta_Cnbeta*beta;

   /* LLLLLLLL Cl_tot LLLLLLLL */

dLdail = delta_Cl_a20 + delta_Cl_a20_lef*dlef;

dLdR = (B/(2*vt))*(Clr + delta_Clr_lef*dlef);

dLdP = (B/(2*vt))*(Clp + delta_Clp_lef*dlef);

Cl_tot = Cl + delta_Cl_lef*dlef + dLdail*dail + delta_Cl_r30*drud + dLdR*R + dLdP*P + delta_Clbeta*beta;

   /* %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      compute Udot,Vdot, Wdot,(as on NASA report p36)
      %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% */

Udot = R*V - Q*W - g*st + qbar*S*Cx_tot/m + T/m;

Vdot = P*W - R*U + g*ct*sphi + qbar*S*Cy_tot/m;

Wdot = Q*U - P*V + g*ct*cphi + qbar*S*Cz_tot/m;

   /* %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      vt_dot equation (from S&L, p82)
      %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% */

xdot[6] = (U*Udot + V*Vdot + W*Wdot)/vt;

   /* %%%%%%%%%%%%%%%%%%
      alpha_dot equation
      %%%%%%%%%%%%%%%%%% */

xdot[7] = (U*Wdot - W*Udot)/(U*U + W*W);

  /* %%%%%%%%%%%%%%%%%
     beta_dot equation
     %%%%%%%%%%%%%%%%% */

xdot[8] = (Vdot*vt - V*xdot[6])/(vt*vt*cb);



  /* %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
     compute Pdot, Qdot, and Rdot (as in Stevens and Lewis p32)
     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% */

L_tot = Cl_tot*qbar*S*B;       /* get moments from coefficients */
M_tot = Cm_tot*qbar*S*cbar;
N_tot = Cn_tot*qbar*S*B;

denom = Jx*Jz - Jxz*Jxz;

  /* %%%%%%%%%%%%%%%%%%%%%%%
     Pdot
     %%%%%%%%%%%%%%%%%%%%%%% */

xdot[9] =  (Jz*L_tot + Jxz*N_tot - (Jz*(Jz-Jy)+Jxz*Jxz)*Q*R + Jxz*(Jx-Jy+Jz)*P*Q + Jxz*Q*Heng)/denom;


  /* %%%%%%%%%%%%%%%%%%%%%%%
     Qdot
     %%%%%%%%%%%%%%%%%%%%%%% */

xdot[10] = (M_tot + (Jz-Jx)*P*R - Jxz*(P*P-R*R) - R*Heng)/Jy;

  /* %%%%%%%%%%%%%%%%%%%%%%%
     Rdot
     %%%%%%%%%%%%%%%%%%%%%%% */

xdot[11] = (Jx*N_tot + Jxz*L_tot + (Jx*(Jx-Jy)+Jxz*Jxz)*P*Q - Jxz*(Jx-Jy+Jz)*Q*R +  Jx*Q*Heng)/denom;

/*########################################*/
/*### Create accelerations anx_cg, any_cg */
/*### ans anz_cg as outputs ##############*/
/*########################################*/

accels(xu,xdot,temp);

xdot[12]  = temp[0];	/* anx_cg */
xdot[13]  = temp[1];	/* any_cg */
xdot[14]  = temp[2];	/* anz_cg */
xdot[15]  = mach;
xdot[16]  = qbar;
xdot[17]  = ps;

/*########################################*/
/*########################################*/

free(temp);

}; /*##### END of nlplant() ####*/

/*########################################*/
/*### Called Sub-Functions  ##############*/
/*########################################*/

/*########################################*/
/* Function for mach and qbar */
/*########################################*/

void atmos(double alt, double vt, double *coeff ){

    double rho0 = 2.377e-3;
    double tfac, temp, rho, mach, qbar, ps;

    tfac =1 - .703e-5*(alt);
    temp = 519.0*tfac;
    if (alt >= 35000.0) {
       temp=390;
    }

    rho=rho0*pow(tfac,4.14);
    mach = (vt)/sqrt(1.4*1716.3*temp);
    qbar = .5*rho*pow(vt,2);
    ps   = 1715.0*rho*temp;

    if (ps == 0){
      ps = 1715;
      }

    coeff[0] = mach;
    coeff[1] = qbar;
    coeff[2] = ps;
}

/*########################################*/
/*########################################*/


/*########################################*/
/*### Port from matlab fix() function ####*/
/*########################################*/

/* port from matlab sign() function */


/*########################################*/
/*########################################*/


/*########################################*/
/*### Calculate accelerations from states */
/*### and state derivatives. ############ */
/*########################################*/

void accels(double *state,
            double *xdot,
            double *y)

{

#define grav 32.174 

double sina, cosa, sinb, cosb ;
double vel_u, vel_v, vel_w ;
double u_dot, v_dot, w_dot ;
double nx_cg, ny_cg, nz_cg ;

sina = sin(state[7]) ;
cosa = cos(state[7]) ;
sinb = sin(state[8]) ;
cosb = cos(state[8]) ;
vel_u = state[6]*cosb*cosa ;
vel_v = state[6]*sinb ;
vel_w = state[6]*cosb*sina ;
u_dot =          cosb*cosa*xdot[6]
      - state[6]*sinb*cosa*xdot[8] 
      - state[6]*cosb*sina*xdot[7] ;
v_dot =          sinb*xdot[6] 
      + state[6]*cosb*xdot[8] ;
w_dot =          cosb*sina*xdot[6]
      - state[6]*sinb*sina*xdot[8] 
      + state[6]*cosb*cosa*xdot[7] ;
nx_cg = 1.0/grav*(u_dot + state[10]*vel_w - state[11]*vel_v)
      + sin(state[4]) ;
ny_cg = 1.0/grav*(v_dot + state[11]*vel_u - state[9]*vel_w)
      - cos(state[4])*sin(state[3]) ;
nz_cg = -1.0/grav*(w_dot + state[9]*vel_v - state[10]*vel_u)
      + cos(state[4])*cos(state[3]) ;

y[0] = nx_cg ;
y[1] = ny_cg ;
y[2] = nz_cg ;


} 

/*########################################*/
/*########################################*/

/*########################################*/
/*########################################*/

void Jac(double *xu, double *xdot, int fidelity, double *jac){ //I modified to have fidelity as a seperate input!

  int fi_flag;

  /* #include f16_constants */
  double g    = 32.17;          /* gravity, ft/s^2 */
  double m    = 636.94;         /* mass, slugs */
  double B    = 30.0;             /* span, ft */
  double S    = 300.0;            /* planform area, ft^2 */
  double cbar = 11.32;          /* mean aero chord, ft */
  double xcgr = 0.25;      /* reference center of gravity as a fraction of cbar */
  double xcg  = 0.25;      /* center of gravity as a fraction of cbar. */

  double Heng = 0.0;              /* turbine momentum along roll axis. */
  double pi   = acos(-1);
  double r2d;                   /* radians to degrees */


/*NasaData        %translated via eq. 2.4-6 on pg 80 of Stevens and Lewis*/

  double Jy  = 55814.0;           /* slug-ft^2 */ 
  double Jxz = 982.0;             /* slug-ft^2 */     
  double Jz  = 63100.0;           /* slug-ft^2 */
  double Jx  = 9496.0;            /* slug-ft^2 */

  double *temp;

  double npos, epos, alt, phi, theta, psi, vt, alpha, beta, P, Q, R;
  double sa, ca, sb, cb, tb, st, ct, tt, sphi, cphi, spsi, cpsi;
  double T, el, ail, rud, dail, drud, lef, dlef;
  double qbar, mach, ps;
  double U, V, W, Udot,Vdot,Wdot;
  double L_tot, M_tot, N_tot, denom;

  double Cx_tot, Cx, delta_Cx_lef, dXdQ, Cxq, delta_Cxq_lef;
  double Cz_tot, Cz, delta_Cz_lef, dZdQ, Czq, delta_Czq_lef;
  double Cm_tot, Cm, eta_el, delta_Cm_lef, dMdQ, Cmq, delta_Cmq_lef, delta_Cm, delta_Cm_ds;
  double Cy_tot, Cy, delta_Cy_lef, dYdail, delta_Cy_r30, dYdR, dYdP;
  double delta_Cy_a20, delta_Cy_a20_lef, Cyr, delta_Cyr_lef, Cyp, delta_Cyp_lef;
  double Cn_tot, Cn, delta_Cn_lef, dNdail, delta_Cn_r30, dNdR, dNdP, delta_Cnbeta;
  double delta_Cn_a20, delta_Cn_a20_lef, Cnr, delta_Cnr_lef, Cnp, delta_Cnp_lef;
  double Cl_tot, Cl, delta_Cl_lef, dLdail, delta_Cl_r30, dLdR, dLdP, delta_Clbeta;
  double delta_Cl_a20, delta_Cl_a20_lef, Clr, delta_Clr_lef, Clp, delta_Clp_lef;

  temp = (double *)malloc(9*sizeof(double));  /*size of 9.1 array*/

  r2d  = 180.0/pi;     /* radians to degrees */

  /* %%%%%%%%%%%%%%%%%%%
           States
     %%%%%%%%%%%%%%%%%%% */


  npos  = xu[0];   /* north position */
  epos  = xu[1];   /* east position */
  alt   = xu[2];   /* altitude */
  phi   = xu[3];   /* orientation angles in rad. */
  theta = xu[4];
  psi   = xu[5];

  vt    = xu[6];     /* total velocity */
  alpha = xu[7]*r2d; /* angle of attack in degrees */
  beta  = xu[8]*r2d; /* sideslip angle in degrees */
  P     = xu[9];     /* Roll Rate --- rolling  moment is Lbar */
  Q     = xu[10];    /* Pitch Rate--- pitching moment is M */
  R     = xu[11];    /* Yaw Rate  --- yawing   moment is N */

  sa    = sin(xu[7]); /* sin(alpha) */
  ca    = cos(xu[7]); /* cos(alpha) */
  sb    = sin(xu[8]); /* sin(beta)  */
  cb    = cos(xu[8]); /* cos(beta)  */
  tb    = tan(xu[8]); /* tan(beta)  */

  st    = sin(theta);
  ct    = cos(theta);
  tt    = tan(theta);
  sphi  = sin(phi);
  cphi  = cos(phi);
  spsi  = sin(psi);
  cpsi  = cos(psi);

  if (vt <= 0.01) {vt = 0.01;}

  /* %%%%%%%%%%%%%%%%%%%
     Control inputs
     %%%%%%%%%%%%%%%%%%% */

  T     = xu[12];   /* thrust */
  el    = xu[13];   /* Elevator setting in degrees. */
  ail   = xu[14];   /* Ailerons mex setting in degrees. */
  rud   = xu[15];   /* Rudder setting in degrees. */
  lef   = xu[16];   /* Leading edge flap setting in degrees */
  
  //fi_flag = xu[17]/1;       /* fi_flag */
  fi_flag = fidelity; // Johns change to pass in integer seperate from numpy float array
    
  /* dail  = ail/20.0;   aileron normalized against max angle */
  /* The aileron was normalized using 20.0 but the NASA report and
     S&L both have 21.5 deg. as maximum deflection. */
  /* As a result... */
  dail  = ail/21.5;
  drud  = rud/30.0;  /* rudder normalized against max angle */
  dlef  = (1 - lef/25.0);  /* leading edge flap normalized against max angle */


  /* %%%%%%%%%%%%%%%%%%
     Atmospheric effects
     sets dynamic pressure and mach number
     %%%%%%%%%%%%%%%%%% */

atmos(alt,vt,temp);
   mach = temp[0];
   qbar = temp[1];
   ps   = temp[2];

/*
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%Dynamics%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
*/

  /* %%%%%%%%%%%%%%%%%%
     Navigation Equations
     %%%%%%%%%%%%%%%%%% */

   U = vt*ca*cb;  /* directional velocities. */
   V = vt*sb;
   W = vt*sa*cb;

/* nposdot */
xdot[0] = U*(ct*cpsi) + 
            V*(sphi*cpsi*st - cphi*spsi) + 
            W*(cphi*st*cpsi + sphi*spsi);

/* eposdot */  
xdot[1] = U*(ct*spsi) + 
            V*(sphi*spsi*st + cphi*cpsi) + 
            W*(cphi*st*spsi - sphi*cpsi);

/* altdot */
xdot[2] = U*st - V*(sphi*ct) - W*(cphi*ct);

  /* %%%%%%%%%%%%%%%%%%%
     Kinematic equations
     %%%%%%%%%%%%%%%%%%% */
/* phidot */
xdot[3] = P + tt*(Q*sphi + R*cphi);


/* theta dot */
xdot[4] = Q*cphi - R*sphi;

/* psidot */
xdot[5] = (Q*sphi + R*cphi)/ct;


/* %%%%%%%%%%%%%%%%%%
        Table lookup
     %%%%%%%%%%%%%%%%%% */

if (fi_flag == 1)          /* HIFI Table */
{
    hifi_C(alpha,beta,el,temp);
        Cx = temp[0];
        Cz = temp[1];
        Cm = temp[2];
        Cy = temp[3];
        Cn = temp[4];
        Cl = temp[5];

    hifi_damping(alpha,temp);
        Cxq = temp[0];
        Cyr = temp[1];
        Cyp = temp[2];
        Czq = temp[3];
        Clr = temp[4];
        Clp = temp[5];
        Cmq = temp[6];
        Cnr = temp[7];
        Cnp = temp[8];

    hifi_C_lef(alpha,beta,temp);
        delta_Cx_lef = temp[0];
        delta_Cz_lef = temp[1];
        delta_Cm_lef = temp[2];
        delta_Cy_lef = temp[3];
        delta_Cn_lef = temp[4];
        delta_Cl_lef = temp[5];

    hifi_damping_lef(alpha,temp);
        delta_Cxq_lef = temp[0];
        delta_Cyr_lef = temp[1];
        delta_Cyp_lef = temp[2];
        delta_Czq_lef = temp[3];
        delta_Clr_lef = temp[4];
        delta_Clp_lef = temp[5];
        delta_Cmq_lef = temp[6];
        delta_Cnr_lef = temp[7];
        delta_Cnp_lef = temp[8];

    hifi_rudder(alpha,beta,temp);
        delta_Cy_r30 = temp[0];
        delta_Cn_r30 = temp[1];
        delta_Cl_r30 = temp[2];

    hifi_ailerons(alpha,beta,temp);
        delta_Cy_a20     = temp[0];
        delta_Cy_a20_lef = temp[1];
        delta_Cn_a20     = temp[2];
        delta_Cn_a20_lef = temp[3];
        delta_Cl_a20     = temp[4];
        delta_Cl_a20_lef = temp[5];

    hifi_other_coeffs(alpha,el,temp);
        delta_Cnbeta = temp[0];
        delta_Clbeta = temp[1];
        delta_Cm     = temp[2];
        eta_el       = temp[3];
        delta_Cm_ds  = 0;        /* ignore deep-stall effect */
    
}

else if (fi_flag == 0)
{     
/* ##############################################
   ##########LOFI Table Look-up #################
   ##############################################*/

/* The lofi model does not include the
   leading edge flap.  All terms multiplied
   dlef have been set to zero but just to 
   be sure we will set it to zero. */
    
    dlef = 0.0;     

    damping(alpha,temp);
        Cxq = temp[0];
        Cyr = temp[1];
        Cyp = temp[2];
        Czq = temp[3];
        Clr = temp[4];
        Clp = temp[5];
        Cmq = temp[6];
        Cnr = temp[7];
        Cnp = temp[8];

    dmomdcon(alpha,beta, temp);
        delta_Cl_a20 = temp[0];     /* Formerly dLda in nlplant.c */
        delta_Cl_r30 = temp[1];     /* Formerly dLdr in nlplant.c */
        delta_Cn_a20 = temp[2];     /* Formerly dNda in nlplant.c */
        delta_Cn_r30 = temp[3];     /* Formerly dNdr in nlplant.c */

    clcn(alpha,beta,temp);
        Cl = temp[0];
        Cn = temp[1];

    cxcm(alpha,el,temp);
        Cx = temp[0];
        Cm = temp[1];

    Cy = -.02*beta + .021*dail + .086*drud;

    cz(alpha,beta,el,temp);
        Cz = temp[0];
/*##################################################
        
        
/*##################################################
  ##  Set all higher order terms of hifi that are ##
  ##  not applicable to lofi equal to zero. ########
  ##################################################*/
     
        delta_Cx_lef    = 0.0;
        delta_Cz_lef    = 0.0;
        delta_Cm_lef    = 0.0;
        delta_Cy_lef    = 0.0;
        delta_Cn_lef    = 0.0;
        delta_Cl_lef    = 0.0;
        delta_Cxq_lef   = 0.0;
        delta_Cyr_lef   = 0.0;
        delta_Cyp_lef   = 0.0;
        delta_Czq_lef   = 0.0;
        delta_Clr_lef   = 0.0;
        delta_Clp_lef   = 0.0;
        delta_Cmq_lef   = 0.0;
        delta_Cnr_lef   = 0.0;
        delta_Cnp_lef   = 0.0;
        delta_Cy_r30    = 0.0;
        delta_Cy_a20    = 0.0;
        delta_Cy_a20_lef= 0.0;
        delta_Cn_a20_lef= 0.0;
        delta_Cl_a20_lef= 0.0;
        delta_Cnbeta    = 0.0;
        delta_Clbeta    = 0.0;
        delta_Cm        = 0.0;
        eta_el          = 1.0;     /* Needs to be one. See equation for Cm_tot*/
        delta_Cm_ds     = 0.0;
                     
/*##################################################
  ##################################################*/ 
}


/* %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
compute Cx_tot, Cz_tot, Cm_tot, Cy_tot, Cn_tot, and Cl_tot
(as on NASA report p37-40)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% */

/* XXXXXXXX Cx_tot XXXXXXXX */

dXdQ = (cbar/(2*vt))*(Cxq + delta_Cxq_lef*dlef);

Cx_tot = Cx + delta_Cx_lef*dlef + dXdQ*Q;

   /* ZZZZZZZZ Cz_tot ZZZZZZZZ */ 

dZdQ = (cbar/(2*vt))*(Czq + delta_Cz_lef*dlef);

Cz_tot = Cz + delta_Cz_lef*dlef + dZdQ*Q;

   /* MMMMMMMM Cm_tot MMMMMMMM */ 

dMdQ = (cbar/(2*vt))*(Cmq + delta_Cmq_lef*dlef);

Cm_tot = Cm*eta_el + Cz_tot*(xcgr-xcg) + delta_Cm_lef*dlef + dMdQ*Q + delta_Cm + delta_Cm_ds;

   /* YYYYYYYY Cy_tot YYYYYYYY */

dYdail = delta_Cy_a20 + delta_Cy_a20_lef*dlef;

dYdR = (B/(2*vt))*(Cyr + delta_Cyr_lef*dlef);

dYdP = (B/(2*vt))*(Cyp + delta_Cyp_lef*dlef);

Cy_tot = Cy + delta_Cy_lef*dlef + dYdail*dail + delta_Cy_r30*drud + dYdR*R + dYdP*P;

   /* NNNNNNNN Cn_tot NNNNNNNN */ 

dNdail = delta_Cn_a20 + delta_Cn_a20_lef*dlef;

dNdR = (B/(2*vt))*(Cnr + delta_Cnr_lef*dlef);

dNdP = (B/(2*vt))*(Cnp + delta_Cnp_lef*dlef);

Cn_tot = Cn + delta_Cn_lef*dlef - Cy_tot*(xcgr-xcg)*(cbar/B) + dNdail*dail + delta_Cn_r30*drud + dNdR*R + dNdP*P + delta_Cnbeta*beta;

   /* LLLLLLLL Cl_tot LLLLLLLL */

dLdail = delta_Cl_a20 + delta_Cl_a20_lef*dlef;

dLdR = (B/(2*vt))*(Clr + delta_Clr_lef*dlef);

dLdP = (B/(2*vt))*(Clp + delta_Clp_lef*dlef);

Cl_tot = Cl + delta_Cl_lef*dlef + dLdail*dail + delta_Cl_r30*drud + dLdR*R + dLdP*P + delta_Clbeta*beta;

   /* %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      compute Udot,Vdot, Wdot,(as on NASA report p36)
      %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% */

Udot = R*V - Q*W - g*st + qbar*S*Cx_tot/m + T/m;

Vdot = P*W - R*U + g*ct*sphi + qbar*S*Cy_tot/m;

Wdot = Q*U - P*V + g*ct*cphi + qbar*S*Cz_tot/m;

   /* %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      vt_dot equation (from S&L, p82)
      %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% */

xdot[6] = (U*Udot + V*Vdot + W*Wdot)/vt;

   /* %%%%%%%%%%%%%%%%%%
      alpha_dot equation
      %%%%%%%%%%%%%%%%%% */

xdot[7] = (U*Wdot - W*Udot)/(U*U + W*W);

  /* %%%%%%%%%%%%%%%%%
     beta_dot equation
     %%%%%%%%%%%%%%%%% */

xdot[8] = (Vdot*vt - V*xdot[6])/(vt*vt*cb);



  /* %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
     compute Pdot, Qdot, and Rdot (as in Stevens and Lewis p32)
     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% */

L_tot = Cl_tot*qbar*S*B;       /* get moments from coefficients */
M_tot = Cm_tot*qbar*S*cbar;
N_tot = Cn_tot*qbar*S*B;

denom = Jx*Jz - Jxz*Jxz;

  /* %%%%%%%%%%%%%%%%%%%%%%%
     Pdot
     %%%%%%%%%%%%%%%%%%%%%%% */

xdot[9] =  (Jz*L_tot + Jxz*N_tot - (Jz*(Jz-Jy)+Jxz*Jxz)*Q*R + Jxz*(Jx-Jy+Jz)*P*Q + Jxz*Q*Heng)/denom;


  /* %%%%%%%%%%%%%%%%%%%%%%%
     Qdot
     %%%%%%%%%%%%%%%%%%%%%%% */

xdot[10] = (M_tot + (Jz-Jx)*P*R - Jxz*(P*P-R*R) - R*Heng)/Jy;

  /* %%%%%%%%%%%%%%%%%%%%%%%
     Rdot
     %%%%%%%%%%%%%%%%%%%%%%% */

xdot[11] = (Jx*N_tot + Jxz*L_tot + (Jx*(Jx-Jy)+Jxz*Jxz)*P*Q - Jxz*(Jx-Jy+Jz)*Q*R +  Jx*Q*Heng)/denom;

/*########################################*/
/*### Create accelerations anx_cg, any_cg */
/*### ans anz_cg as outputs ##############*/
/*########################################*/

accels(xu,xdot,temp);

xdot[12]  = temp[0];	/* anx_cg */
xdot[13]  = temp[1];	/* any_cg */
xdot[14]  = temp[2];	/* anz_cg */
xdot[15]  = mach;
xdot[16]  = qbar;
xdot[17]  = ps;

/*########################################*/
/*########################################*/

jac[0] = 0; 
jac[1] = 0; 
jac[2] = 0; 
jac[3] = sb*vt*(sphi*spsi + cphi*cpsi*st) + cb*sa*vt*(cphi*spsi - cpsi*sphi*st); 
jac[4] = cpsi*ct*sb*sphi*vt - ca*cb*cpsi*st*vt + cb*cphi*cpsi*ct*sa*vt; 
jac[5] = cb*sa*vt*(cpsi*sphi - cphi*spsi*st) - sb*vt*(cphi*cpsi + sphi*spsi*st) - ca*cb*ct*spsi*vt; 
jac[6] = cb*sa*(sphi*spsi + cphi*cpsi*st) - sb*(cphi*spsi - cpsi*sphi*st) + ca*cb*cpsi*ct; 
jac[7] = ca*cb*vt*(sphi*spsi + cphi*cpsi*st) - cb*cpsi*ct*sa*vt; 
jac[8] = - cb*vt*(cphi*spsi - cpsi*sphi*st) - sa*sb*vt*(sphi*spsi + cphi*cpsi*st) - ca*cpsi*ct*sb*vt; 
jac[9] = 0; 
jac[10] = 0; 
jac[11] = 0; 
jac[12] = 0; 
jac[13] = 0; 
jac[14] = 0; 
jac[15] = - sb*vt*(cpsi*sphi - cphi*spsi*st) - cb*sa*vt*(cphi*cpsi + sphi*spsi*st); 
jac[16] = ct*sb*sphi*spsi*vt - ca*cb*spsi*st*vt + cb*cphi*ct*sa*spsi*vt; 
jac[17] = cb*sa*vt*(sphi*spsi + cphi*cpsi*st) - sb*vt*(cphi*spsi - cpsi*sphi*st) + ca*cb*cpsi*ct*vt; 
jac[18] = sb*(cphi*cpsi + sphi*spsi*st) - cb*sa*(cpsi*sphi - cphi*spsi*st) + ca*cb*ct*spsi; 
jac[19] = - ca*cb*vt*(cpsi*sphi - cphi*spsi*st) - cb*ct*sa*spsi*vt; 
jac[20] = cb*vt*(cphi*cpsi + sphi*spsi*st) + sa*sb*vt*(cpsi*sphi - cphi*spsi*st) - ca*ct*sb*spsi*vt; 
jac[21] = 0; 
jac[22] = 0; 
jac[23] = 0; 
jac[24] = 0; 
jac[25] = 0; 
jac[26] = 0; 
jac[27] = cb*ct*sa*sphi*vt - cphi*ct*sb*vt; 
jac[28] = ca*cb*ct*vt + sb*sphi*st*vt + cb*cphi*sa*st*vt; 
jac[29] = 0; 
jac[30] = ca*cb*st - ct*sb*sphi - cb*cphi*ct*sa; 
jac[31] = - cb*sa*st*vt - ca*cb*cphi*ct*vt; 
jac[32] = cphi*ct*sa*sb*vt - ca*sb*st*vt - cb*ct*sphi*vt; 
jac[33] = 0; 
jac[34] = 0; 
jac[35] = 0; 
jac[36] = 0; 
jac[37] = 0; 
jac[38] = 0; 
jac[39] = tt*(Q*cphi - R*sphi); 
jac[40] = (pow(tt,2) + 1)*(R*cphi + Q*sphi); 
jac[41] = 0; 
jac[42] = 0; 
jac[43] = 0; 
jac[44] = 0; 
jac[45] = 1; 
jac[46] = sphi*tt; 
jac[47] = cphi*tt; 
jac[48] = 0; 
jac[49] = 0; 
jac[50] = 0; 
jac[51] = - R*cphi - Q*sphi; 
jac[52] = 0; 
jac[53] = 0; 
jac[54] = 0; 
jac[55] = 0; 
jac[56] = 0; 
jac[57] = 0; 
jac[58] = cphi; 
jac[59] = -sphi; 
jac[60] = 0; 
jac[61] = 0; 
jac[62] = 0; 
jac[63] = (Q*cphi - R*sphi)/ct; 
jac[64] = (st*(R*cphi + Q*sphi))/pow(ct,2); 
jac[65] = 0; 
jac[66] = 0; 
jac[67] = 0; 
jac[68] = 0; 
jac[69] = 0; 
jac[70] = sphi/ct; 
jac[71] = cphi/ct; 
jac[72] = 0; 
jac[73] = 0; 
jac[74] = -((0.000000034590341700000001074671459439479*S*sb*pow(vt,3)*pow((1.0 - 0.0000070299999999999996077638432512291*alt),(157/50))*(Cy + (delta_Cy_r30*rud)/30 + (2*ail*(delta_Cy_a20 - delta_Cy_a20_lef*(lef/25 - 1)))/43 - delta_Cy_lef*(lef/25 - 1) + (B*P*(Cyp - delta_Cyp_lef*(lef/25 - 1)))/(2*vt) + (B*R*(Cyr - delta_Cyr_lef*(lef/25 - 1)))/(2*vt)))/m + (0.000000034590341700000001074671459439479*S*ca*cb*pow(vt,3)*pow((1.0 - 0.0000070299999999999996077638432512291*alt),(157/50))*(Cx - delta_Cx_lef*(lef/25 - 1) + (Q*cbar*(Cxq - delta_Cxq_lef*(lef/25 - 1)))/(2*vt)))/m + (0.000000034590341700000001074671459439479*S*cb*sa*pow(vt,3)*pow((1.0 - 0.0000070299999999999996077638432512291*alt),(157/50))*(Cz - delta_Cz_lef*(lef/25 - 1) + (Q*cbar*(Czq - delta_Cz_lef*(lef/25 - 1)))/(2*vt)))/m)/vt; 
jac[75] = (cphi*ct*g*sb*vt - cb*ct*g*sa*sphi*vt)/vt; 
jac[76] = -(ca*cb*ct*g*vt + g*sb*sphi*st*vt + cb*cphi*g*sa*st*vt)/vt; 
jac[77] = 0; 
jac[78] = (sb*(ct*g*sphi - R*ca*cb*vt + P*cb*sa*vt + (0.0011885000000000001032368635023317*S*pow(vt,2)*pow((1.0 - 0.0000070299999999999996077638432512291*alt),(207/50))*(Cy + (delta_Cy_r30*rud)/30 + (2*ail*(delta_Cy_a20 - delta_Cy_a20_lef*(lef/25 - 1)))/43 - delta_Cy_lef*(lef/25 - 1) + (B*P*(Cyp - delta_Cyp_lef*(lef/25 - 1)))/(2*vt) + (B*R*(Cyr - delta_Cyr_lef*(lef/25 - 1)))/(2*vt)))/m) + cb*sa*(cphi*ct*g - P*sb*vt + Q*ca*cb*vt + (0.0011885000000000001032368635023317*S*pow(vt,2)*pow((1.0 - 0.0000070299999999999996077638432512291*alt),(207/50))*(Cz - delta_Cz_lef*(lef/25 - 1) + (Q*cbar*(Czq - delta_Cz_lef*(lef/25 - 1)))/(2*vt)))/m) + ca*cb*(T/m - g*st + R*sb*vt - Q*cb*sa*vt + (0.0011885000000000001032368635023317*S*pow(vt,2)*pow((1.0 - 0.0000070299999999999996077638432512291*alt),(207/50))*(Cx - delta_Cx_lef*(lef/25 - 1) + (Q*cbar*(Cxq - delta_Cxq_lef*(lef/25 - 1)))/(2*vt)))/m) - sb*vt*(R*ca*cb - P*cb*sa - (0.0023770000000000002064737270046635*S*vt*pow((1.0 - 0.0000070299999999999996077638432512291*alt),(207/50))*(Cy + (delta_Cy_r30*rud)/30 + (2*ail*(delta_Cy_a20 - delta_Cy_a20_lef*(lef/25 - 1)))/43 - delta_Cy_lef*(lef/25 - 1) + (B*P*(Cyp - delta_Cyp_lef*(lef/25 - 1)))/(2*vt) + (B*R*(Cyr - delta_Cyr_lef*(lef/25 - 1)))/(2*vt)))/m + (0.0011885000000000001032368635023317*S*pow(vt,2)*((B*P*(Cyp - delta_Cyp_lef*(lef/25 - 1)))/(2*pow(vt,2)) + (B*R*(Cyr - delta_Cyr_lef*(lef/25 - 1)))/(2*pow(vt,2)))*pow((1.0 - 0.0000070299999999999996077638432512291*alt),(207/50)))/m) - cb*sa*vt*(P*sb - Q*ca*cb - (0.0023770000000000002064737270046635*S*vt*pow((1.0 - 0.0000070299999999999996077638432512291*alt),(207/50))*(Cz - delta_Cz_lef*(lef/25 - 1) + (Q*cbar*(Czq - delta_Cz_lef*(lef/25 - 1)))/(2*vt)))/m + (0.00059425000000000005161843175116587*Q*S*cbar*(Czq - delta_Cz_lef*(lef/25 - 1))*pow((1.0 - 0.0000070299999999999996077638432512291*alt),(207/50)))/m) + ca*cb*vt*(R*sb - Q*cb*sa + (0.0023770000000000002064737270046635*S*vt*pow((1.0 - 0.0000070299999999999996077638432512291*alt),(207/50))*(Cx - delta_Cx_lef*(lef/25 - 1) + (Q*cbar*(Cxq - delta_Cxq_lef*(lef/25 - 1)))/(2*vt)))/m - (0.00059425000000000005161843175116587*Q*S*cbar*(Cxq - delta_Cxq_lef*(lef/25 - 1))*pow((1.0 - 0.0000070299999999999996077638432512291*alt),(207/50)))/m))/vt - (sb*vt*(ct*g*sphi - R*ca*cb*vt + P*cb*sa*vt + (0.0011885000000000001032368635023317*S*pow(vt,2)*pow((1.0 - 0.0000070299999999999996077638432512291*alt),(207/50))*(Cy + (delta_Cy_r30*rud)/30 + (2*ail*(delta_Cy_a20 - delta_Cy_a20_lef*(lef/25 - 1)))/43 - delta_Cy_lef*(lef/25 - 1) + (B*P*(Cyp - delta_Cyp_lef*(lef/25 - 1)))/(2*vt) + (B*R*(Cyr - delta_Cyr_lef*(lef/25 - 1)))/(2*vt)))/m) + cb*sa*vt*(cphi*ct*g - P*sb*vt + Q*ca*cb*vt + (0.0011885000000000001032368635023317*S*pow(vt,2)*pow((1.0 - 0.0000070299999999999996077638432512291*alt),(207/50))*(Cz - delta_Cz_lef*(lef/25 - 1) + (Q*cbar*(Czq - delta_Cz_lef*(lef/25 - 1)))/(2*vt)))/m) + ca*cb*vt*(T/m - g*st + R*sb*vt - Q*cb*sa*vt + (0.0011885000000000001032368635023317*S*pow(vt,2)*pow((1.0 - 0.0000070299999999999996077638432512291*alt),(207/50))*(Cx - delta_Cx_lef*(lef/25 - 1) + (Q*cbar*(Cxq - delta_Cxq_lef*(lef/25 - 1)))/(2*vt)))/m))/pow(vt,2); 
jac[79] = -(Q*pow(ca,2)*pow(cb,2)*pow(vt,2) - sb*vt*(P*ca*cb*vt + R*cb*sa*vt) + Q*pow(cb,2)*pow(sa,2)*pow(vt,2) + cb*sa*vt*(T/m - g*st + R*sb*vt - Q*cb*sa*vt + (0.0011885000000000001032368635023317*S*pow(vt,2)*pow((1.0 - 0.0000070299999999999996077638432512291*alt),(207/50))*(Cx - delta_Cx_lef*(lef/25 - 1) + (Q*cbar*(Cxq - delta_Cxq_lef*(lef/25 - 1)))/(2*vt)))/m) - ca*cb*vt*(cphi*ct*g - P*sb*vt + Q*ca*cb*vt + (0.0011885000000000001032368635023317*S*pow(vt,2)*pow((1.0 - 0.0000070299999999999996077638432512291*alt),(207/50))*(Cz - delta_Cz_lef*(lef/25 - 1) + (Q*cbar*(Czq - delta_Cz_lef*(lef/25 - 1)))/(2*vt)))/m))/vt; 
jac[80] = -(sb*vt*(P*sa*sb*vt - R*ca*sb*vt) - cb*vt*(ct*g*sphi - R*ca*cb*vt + P*cb*sa*vt + (0.0011885000000000001032368635023317*S*pow(vt,2)*pow((1.0 - 0.0000070299999999999996077638432512291*alt),(207/50))*(Cy + (delta_Cy_r30*rud)/30 + (2*ail*(delta_Cy_a20 - delta_Cy_a20_lef*(lef/25 - 1)))/43 - delta_Cy_lef*(lef/25 - 1) + (B*P*(Cyp - delta_Cyp_lef*(lef/25 - 1)))/(2*vt) + (B*R*(Cyr - delta_Cyr_lef*(lef/25 - 1)))/(2*vt)))/m) + sa*sb*vt*(cphi*ct*g - P*sb*vt + Q*ca*cb*vt + (0.0011885000000000001032368635023317*S*pow(vt,2)*pow((1.0 - 0.0000070299999999999996077638432512291*alt),(207/50))*(Cz - delta_Cz_lef*(lef/25 - 1) + (Q*cbar*(Czq - delta_Cz_lef*(lef/25 - 1)))/(2*vt)))/m) + ca*sb*vt*(T/m - g*st + R*sb*vt - Q*cb*sa*vt + (0.0011885000000000001032368635023317*S*pow(vt,2)*pow((1.0 - 0.0000070299999999999996077638432512291*alt),(207/50))*(Cx - delta_Cx_lef*(lef/25 - 1) + (Q*cbar*(Cxq - delta_Cxq_lef*(lef/25 - 1)))/(2*vt)))/m) + cb*sa*vt*(P*cb*vt + Q*ca*sb*vt) - ca*cb*vt*(R*cb*vt + Q*sa*sb*vt))/vt; 
jac[81] = (sb*vt*(cb*sa*vt + (0.00059425000000000005161843175116587*B*S*vt*(Cyp - delta_Cyp_lef*(lef/25 - 1))*pow((1.0 - 0.0000070299999999999996077638432512291*alt),(207/50)))/m) - cb*sa*sb*pow(vt,2))/vt; 
jac[82] = -(ca*cb*vt*(cb*sa*vt - (0.00059425000000000005161843175116587*S*cbar*vt*(Cxq - delta_Cxq_lef*(lef/25 - 1))*pow((1.0 - 0.0000070299999999999996077638432512291*alt),(207/50)))/m) - cb*sa*vt*(ca*cb*vt + (0.00059425000000000005161843175116587*S*cbar*vt*(Czq - delta_Cz_lef*(lef/25 - 1))*pow((1.0 - 0.0000070299999999999996077638432512291*alt),(207/50)))/m))/vt; 
jac[83] = -(sb*vt*(ca*cb*vt - (0.00059425000000000005161843175116587*B*S*vt*(Cyr - delta_Cyr_lef*(lef/25 - 1))*pow((1.0 - 0.0000070299999999999996077638432512291*alt),(207/50)))/m) - ca*cb*sb*pow(vt,2))/vt; 
jac[84] = 0; 
jac[85] = 0; 
jac[86] = -((0.000000034590341700000001074671459439479*S*ca*cb*pow(vt,3)*pow((1.0 - 0.0000070299999999999996077638432512291*alt),(157/50))*(Cz - delta_Cz_lef*(lef/25 - 1) + (Q*cbar*(Czq - delta_Cz_lef*(lef/25 - 1)))/(2*vt)))/m - (0.000000034590341700000001074671459439479*S*cb*sa*pow(vt,3)*pow((1.0 - 0.0000070299999999999996077638432512291*alt),(157/50))*(Cx - delta_Cx_lef*(lef/25 - 1) + (Q*cbar*(Cxq - delta_Cxq_lef*(lef/25 - 1)))/(2*vt)))/m)/(pow(ca,2)*pow(cb,2)*pow(vt,2) + pow(cb,2)*pow(sa,2)*pow(vt,2)); 
jac[87] = -(ca*cb*ct*g*sphi*vt)/(pow(ca,2)*pow(cb,2)*pow(vt,2) + pow(cb,2)*pow(sa,2)*pow(vt,2)); 
jac[88] = (cb*ct*g*sa*vt - ca*cb*cphi*g*st*vt)/(pow(ca,2)*pow(cb,2)*pow(vt,2) + pow(cb,2)*pow(sa,2)*pow(vt,2)); 
jac[89] = 0; 
jac[90] = ((2*pow(ca,2)*pow(cb,2)*vt + 2*pow(cb,2)*pow(sa,2)*vt)*(cb*sa*vt*(T/m - g*st + R*sb*vt - Q*cb*sa*vt + (0.0011885000000000001032368635023317*S*pow(vt,2)*pow((1.0 - 0.0000070299999999999996077638432512291*alt),(207/50))*(Cx - delta_Cx_lef*(lef/25 - 1) + (Q*cbar*(Cxq - delta_Cxq_lef*(lef/25 - 1)))/(2*vt)))/m) - ca*cb*vt*(cphi*ct*g - P*sb*vt + Q*ca*cb*vt + (0.0011885000000000001032368635023317*S*pow(vt,2)*pow((1.0 - 0.0000070299999999999996077638432512291*alt),(207/50))*(Cz - delta_Cz_lef*(lef/25 - 1) + (Q*cbar*(Czq - delta_Cz_lef*(lef/25 - 1)))/(2*vt)))/m)))/pow((pow(ca,2)*pow(cb,2)*pow(vt,2) + pow(cb,2)*pow(sa,2)*pow(vt,2)),2) - (cb*sa*(T/m - g*st + R*sb*vt - Q*cb*sa*vt + (0.0011885000000000001032368635023317*S*pow(vt,2)*pow((1.0 - 0.0000070299999999999996077638432512291*alt),(207/50))*(Cx - delta_Cx_lef*(lef/25 - 1) + (Q*cbar*(Cxq - delta_Cxq_lef*(lef/25 - 1)))/(2*vt)))/m) - ca*cb*(cphi*ct*g - P*sb*vt + Q*ca*cb*vt + (0.0011885000000000001032368635023317*S*pow(vt,2)*pow((1.0 - 0.0000070299999999999996077638432512291*alt),(207/50))*(Cz - delta_Cz_lef*(lef/25 - 1) + (Q*cbar*(Czq - delta_Cz_lef*(lef/25 - 1)))/(2*vt)))/m) + ca*cb*vt*(P*sb - Q*ca*cb - (0.0023770000000000002064737270046635*S*vt*pow((1.0 - 0.0000070299999999999996077638432512291*alt),(207/50))*(Cz - delta_Cz_lef*(lef/25 - 1) + (Q*cbar*(Czq - delta_Cz_lef*(lef/25 - 1)))/(2*vt)))/m + (0.00059425000000000005161843175116587*Q*S*cbar*(Czq - delta_Cz_lef*(lef/25 - 1))*pow((1.0 - 0.0000070299999999999996077638432512291*alt),(207/50)))/m) + cb*sa*vt*(R*sb - Q*cb*sa + (0.0023770000000000002064737270046635*S*vt*pow((1.0 - 0.0000070299999999999996077638432512291*alt),(207/50))*(Cx - delta_Cx_lef*(lef/25 - 1) + (Q*cbar*(Cxq - delta_Cxq_lef*(lef/25 - 1)))/(2*vt)))/m - (0.00059425000000000005161843175116587*Q*S*cbar*(Cxq - delta_Cxq_lef*(lef/25 - 1))*pow((1.0 - 0.0000070299999999999996077638432512291*alt),(207/50)))/m))/(pow(ca,2)*pow(cb,2)*pow(vt,2) + pow(cb,2)*pow(sa,2)*pow(vt,2)); 
jac[91] = -(cb*sa*vt*(cphi*ct*g - P*sb*vt + Q*ca*cb*vt + (0.0011885000000000001032368635023317*S*pow(vt,2)*pow((1.0 - 0.0000070299999999999996077638432512291*alt),(207/50))*(Cz - delta_Cz_lef*(lef/25 - 1) + (Q*cbar*(Czq - delta_Cz_lef*(lef/25 - 1)))/(2*vt)))/m) + ca*cb*vt*(T/m - g*st + R*sb*vt - Q*cb*sa*vt + (0.0011885000000000001032368635023317*S*pow(vt,2)*pow((1.0 - 0.0000070299999999999996077638432512291*alt),(207/50))*(Cx - delta_Cx_lef*(lef/25 - 1) + (Q*cbar*(Cxq - delta_Cxq_lef*(lef/25 - 1)))/(2*vt)))/m))/(pow(ca,2)*pow(cb,2)*pow(vt,2) + pow(cb,2)*pow(sa,2)*pow(vt,2)); 
jac[92] = - (ca*sb*vt*(cphi*ct*g - P*sb*vt + Q*ca*cb*vt + (0.0011885000000000001032368635023317*S*pow(vt,2)*pow((1.0 - 0.0000070299999999999996077638432512291*alt),(207/50))*(Cz - delta_Cz_lef*(lef/25 - 1) + (Q*cbar*(Czq - delta_Cz_lef*(lef/25 - 1)))/(2*vt)))/m) - sa*sb*vt*(T/m - g*st + R*sb*vt - Q*cb*sa*vt + (0.0011885000000000001032368635023317*S*pow(vt,2)*pow((1.0 - 0.0000070299999999999996077638432512291*alt),(207/50))*(Cx - delta_Cx_lef*(lef/25 - 1) + (Q*cbar*(Cxq - delta_Cxq_lef*(lef/25 - 1)))/(2*vt)))/m) + ca*cb*vt*(P*cb*vt + Q*ca*sb*vt) + cb*sa*vt*(R*cb*vt + Q*sa*sb*vt))/(pow(ca,2)*pow(cb,2)*pow(vt,2) + pow(cb,2)*pow(sa,2)*pow(vt,2)) - ((2*pow(ca,2)*cb*sb*pow(vt,2) + 2*cb*pow(sa,2)*sb*pow(vt,2))*(cb*sa*vt*(T/m - g*st + R*sb*vt - Q*cb*sa*vt + (0.0011885000000000001032368635023317*S*pow(vt,2)*pow((1.0 - 0.0000070299999999999996077638432512291*alt),(207/50))*(Cx - delta_Cx_lef*(lef/25 - 1) + (Q*cbar*(Cxq - delta_Cxq_lef*(lef/25 - 1)))/(2*vt)))/m) - ca*cb*vt*(cphi*ct*g - P*sb*vt + Q*ca*cb*vt + (0.0011885000000000001032368635023317*S*pow(vt,2)*pow((1.0 - 0.0000070299999999999996077638432512291*alt),(207/50))*(Cz - delta_Cz_lef*(lef/25 - 1) + (Q*cbar*(Czq - delta_Cz_lef*(lef/25 - 1)))/(2*vt)))/m)))/pow((pow(ca,2)*pow(cb,2)*pow(vt,2) + pow(cb,2)*pow(sa,2)*pow(vt,2)),2); 
jac[93] = -(ca*cb*sb*pow(vt,2))/(pow(ca,2)*pow(cb,2)*pow(vt,2) + pow(cb,2)*pow(sa,2)*pow(vt,2)); 
jac[94] = (ca*cb*vt*(ca*cb*vt + (0.00059425000000000005161843175116587*S*cbar*vt*(Czq - delta_Cz_lef*(lef/25 - 1))*pow((1.0 - 0.0000070299999999999996077638432512291*alt),(207/50)))/m) + cb*sa*vt*(cb*sa*vt - (0.00059425000000000005161843175116587*S*cbar*vt*(Cxq - delta_Cxq_lef*(lef/25 - 1))*pow((1.0 - 0.0000070299999999999996077638432512291*alt),(207/50)))/m))/(pow(ca,2)*pow(cb,2)*pow(vt,2) + pow(cb,2)*pow(sa,2)*pow(vt,2)); 
jac[95] = -(cb*sa*sb*pow(vt,2))/(pow(ca,2)*pow(cb,2)*pow(vt,2) + pow(cb,2)*pow(sa,2)*pow(vt,2)); 
jac[96] = 0; 
jac[97] = 0; 
jac[98] = (sb*((0.000000034590341700000001074671459439479*S*sb*pow(vt,3)*pow((1.0 - 0.0000070299999999999996077638432512291*alt),(157/50))*(Cy + (delta_Cy_r30*rud)/30 + (2*ail*(delta_Cy_a20 - delta_Cy_a20_lef*(lef/25 - 1)))/43 - delta_Cy_lef*(lef/25 - 1) + (B*P*(Cyp - delta_Cyp_lef*(lef/25 - 1)))/(2*vt) + (B*R*(Cyr - delta_Cyr_lef*(lef/25 - 1)))/(2*vt)))/m + (0.000000034590341700000001074671459439479*S*ca*cb*pow(vt,3)*pow((1.0 - 0.0000070299999999999996077638432512291*alt),(157/50))*(Cx - delta_Cx_lef*(lef/25 - 1) + (Q*cbar*(Cxq - delta_Cxq_lef*(lef/25 - 1)))/(2*vt)))/m + (0.000000034590341700000001074671459439479*S*cb*sa*pow(vt,3)*pow((1.0 - 0.0000070299999999999996077638432512291*alt),(157/50))*(Cz - delta_Cz_lef*(lef/25 - 1) + (Q*cbar*(Czq - delta_Cz_lef*(lef/25 - 1)))/(2*vt)))/m) - (0.000000034590341700000001074671459439479*S*pow(vt,3)*pow((1.0 - 0.0000070299999999999996077638432512291*alt),(157/50))*(Cy + (delta_Cy_r30*rud)/30 + (2*ail*(delta_Cy_a20 - delta_Cy_a20_lef*(lef/25 - 1)))/43 - delta_Cy_lef*(lef/25 - 1) + (B*P*(Cyp - delta_Cyp_lef*(lef/25 - 1)))/(2*vt) + (B*R*(Cyr - delta_Cyr_lef*(lef/25 - 1)))/(2*vt)))/m)/(cb*pow(vt,2)); 
jac[99] = -(sb*(cphi*ct*g*sb*vt - cb*ct*g*sa*sphi*vt) - cphi*ct*g*vt)/(cb*pow(vt,2)); 
jac[100] = (sb*(ca*cb*ct*g*vt + g*sb*sphi*st*vt + cb*cphi*g*sa*st*vt) - g*sphi*st*vt)/(cb*pow(vt,2)); 
jac[101] = 0; 
jac[102] = (2*(sb*(sb*vt*(ct*g*sphi - R*ca*cb*vt + P*cb*sa*vt + (0.0011885000000000001032368635023317*S*pow(vt,2)*pow((1.0 - 0.0000070299999999999996077638432512291*alt),(207/50))*(Cy + (delta_Cy_r30*rud)/30 + (2*ail*(delta_Cy_a20 - delta_Cy_a20_lef*(lef/25 - 1)))/43 - delta_Cy_lef*(lef/25 - 1) + (B*P*(Cyp - delta_Cyp_lef*(lef/25 - 1)))/(2*vt) + (B*R*(Cyr - delta_Cyr_lef*(lef/25 - 1)))/(2*vt)))/m) + cb*sa*vt*(cphi*ct*g - P*sb*vt + Q*ca*cb*vt + (0.0011885000000000001032368635023317*S*pow(vt,2)*pow((1.0 - 0.0000070299999999999996077638432512291*alt),(207/50))*(Cz - delta_Cz_lef*(lef/25 - 1) + (Q*cbar*(Czq - delta_Cz_lef*(lef/25 - 1)))/(2*vt)))/m) + ca*cb*vt*(T/m - g*st + R*sb*vt - Q*cb*sa*vt + (0.0011885000000000001032368635023317*S*pow(vt,2)*pow((1.0 - 0.0000070299999999999996077638432512291*alt),(207/50))*(Cx - delta_Cx_lef*(lef/25 - 1) + (Q*cbar*(Cxq - delta_Cxq_lef*(lef/25 - 1)))/(2*vt)))/m)) - vt*(ct*g*sphi - R*ca*cb*vt + P*cb*sa*vt + (0.0011885000000000001032368635023317*S*pow(vt,2)*pow((1.0 - 0.0000070299999999999996077638432512291*alt),(207/50))*(Cy + (delta_Cy_r30*rud)/30 + (2*ail*(delta_Cy_a20 - delta_Cy_a20_lef*(lef/25 - 1)))/43 - delta_Cy_lef*(lef/25 - 1) + (B*P*(Cyp - delta_Cyp_lef*(lef/25 - 1)))/(2*vt) + (B*R*(Cyr - delta_Cyr_lef*(lef/25 - 1)))/(2*vt)))/m)))/(cb*pow(vt,3)) - (vt*(R*ca*cb - P*cb*sa - (0.0023770000000000002064737270046635*S*vt*pow((1.0 - 0.0000070299999999999996077638432512291*alt),(207/50))*(Cy + (delta_Cy_r30*rud)/30 + (2*ail*(delta_Cy_a20 - delta_Cy_a20_lef*(lef/25 - 1)))/43 - delta_Cy_lef*(lef/25 - 1) + (B*P*(Cyp - delta_Cyp_lef*(lef/25 - 1)))/(2*vt) + (B*R*(Cyr - delta_Cyr_lef*(lef/25 - 1)))/(2*vt)))/m + (0.0011885000000000001032368635023317*S*pow(vt,2)*((B*P*(Cyp - delta_Cyp_lef*(lef/25 - 1)))/(2*pow(vt,2)) + (B*R*(Cyr - delta_Cyr_lef*(lef/25 - 1)))/(2*pow(vt,2)))*pow((1.0 - 0.0000070299999999999996077638432512291*alt),(207/50)))/m) + sb*(sb*(ct*g*sphi - R*ca*cb*vt + P*cb*sa*vt + (0.0011885000000000001032368635023317*S*pow(vt,2)*pow((1.0 - 0.0000070299999999999996077638432512291*alt),(207/50))*(Cy + (delta_Cy_r30*rud)/30 + (2*ail*(delta_Cy_a20 - delta_Cy_a20_lef*(lef/25 - 1)))/43 - delta_Cy_lef*(lef/25 - 1) + (B*P*(Cyp - delta_Cyp_lef*(lef/25 - 1)))/(2*vt) + (B*R*(Cyr - delta_Cyr_lef*(lef/25 - 1)))/(2*vt)))/m) + cb*sa*(cphi*ct*g - P*sb*vt + Q*ca*cb*vt + (0.0011885000000000001032368635023317*S*pow(vt,2)*pow((1.0 - 0.0000070299999999999996077638432512291*alt),(207/50))*(Cz - delta_Cz_lef*(lef/25 - 1) + (Q*cbar*(Czq - delta_Cz_lef*(lef/25 - 1)))/(2*vt)))/m) + ca*cb*(T/m - g*st + R*sb*vt - Q*cb*sa*vt + (0.0011885000000000001032368635023317*S*pow(vt,2)*pow((1.0 - 0.0000070299999999999996077638432512291*alt),(207/50))*(Cx - delta_Cx_lef*(lef/25 - 1) + (Q*cbar*(Cxq - delta_Cxq_lef*(lef/25 - 1)))/(2*vt)))/m) - sb*vt*(R*ca*cb - P*cb*sa - (0.0023770000000000002064737270046635*S*vt*pow((1.0 - 0.0000070299999999999996077638432512291*alt),(207/50))*(Cy + (delta_Cy_r30*rud)/30 + (2*ail*(delta_Cy_a20 - delta_Cy_a20_lef*(lef/25 - 1)))/43 - delta_Cy_lef*(lef/25 - 1) + (B*P*(Cyp - delta_Cyp_lef*(lef/25 - 1)))/(2*vt) + (B*R*(Cyr - delta_Cyr_lef*(lef/25 - 1)))/(2*vt)))/m + (0.0011885000000000001032368635023317*S*pow(vt,2)*((B*P*(Cyp - delta_Cyp_lef*(lef/25 - 1)))/(2*pow(vt,2)) + (B*R*(Cyr - delta_Cyr_lef*(lef/25 - 1)))/(2*pow(vt,2)))*pow((1.0 - 0.0000070299999999999996077638432512291*alt),(207/50)))/m) - cb*sa*vt*(P*sb - Q*ca*cb - (0.0023770000000000002064737270046635*S*vt*pow((1.0 - 0.0000070299999999999996077638432512291*alt),(207/50))*(Cz - delta_Cz_lef*(lef/25 - 1) + (Q*cbar*(Czq - delta_Cz_lef*(lef/25 - 1)))/(2*vt)))/m + (0.00059425000000000005161843175116587*Q*S*cbar*(Czq - delta_Cz_lef*(lef/25 - 1))*pow((1.0 - 0.0000070299999999999996077638432512291*alt),(207/50)))/m) + ca*cb*vt*(R*sb - Q*cb*sa + (0.0023770000000000002064737270046635*S*vt*pow((1.0 - 0.0000070299999999999996077638432512291*alt),(207/50))*(Cx - delta_Cx_lef*(lef/25 - 1) + (Q*cbar*(Cxq - delta_Cxq_lef*(lef/25 - 1)))/(2*vt)))/m - (0.00059425000000000005161843175116587*Q*S*cbar*(Cxq - delta_Cxq_lef*(lef/25 - 1))*pow((1.0 - 0.0000070299999999999996077638432512291*alt),(207/50)))/m)) - ct*g*sphi + R*ca*cb*vt - P*cb*sa*vt - (0.0011885000000000001032368635023317*S*pow(vt,2)*pow((1.0 - 0.0000070299999999999996077638432512291*alt),(207/50))*(Cy + (delta_Cy_r30*rud)/30 + (2*ail*(delta_Cy_a20 - delta_Cy_a20_lef*(lef/25 - 1)))/43 - delta_Cy_lef*(lef/25 - 1) + (B*P*(Cyp - delta_Cyp_lef*(lef/25 - 1)))/(2*vt) + (B*R*(Cyr - delta_Cyr_lef*(lef/25 - 1)))/(2*vt)))/m)/(cb*pow(vt,2)); 
jac[103] = (vt*(P*ca*cb*vt + R*cb*sa*vt) + sb*(Q*pow(ca,2)*pow(cb,2)*pow(vt,2) - sb*vt*(P*ca*cb*vt + R*cb*sa*vt) + Q*pow(cb,2)*pow(sa,2)*pow(vt,2) + cb*sa*vt*(T/m - g*st + R*sb*vt - Q*cb*sa*vt + (0.0011885000000000001032368635023317*S*pow(vt,2)*pow((1.0 - 0.0000070299999999999996077638432512291*alt),(207/50))*(Cx - delta_Cx_lef*(lef/25 - 1) + (Q*cbar*(Cxq - delta_Cxq_lef*(lef/25 - 1)))/(2*vt)))/m) - ca*cb*vt*(cphi*ct*g - P*sb*vt + Q*ca*cb*vt + (0.0011885000000000001032368635023317*S*pow(vt,2)*pow((1.0 - 0.0000070299999999999996077638432512291*alt),(207/50))*(Cz - delta_Cz_lef*(lef/25 - 1) + (Q*cbar*(Czq - delta_Cz_lef*(lef/25 - 1)))/(2*vt)))/m)))/(cb*pow(vt,2)); 
jac[104] = - (vt*(P*sa*sb*vt - R*ca*sb*vt) - sb*(sb*vt*(P*sa*sb*vt - R*ca*sb*vt) - cb*vt*(ct*g*sphi - R*ca*cb*vt + P*cb*sa*vt + (0.0011885000000000001032368635023317*S*pow(vt,2)*pow((1.0 - 0.0000070299999999999996077638432512291*alt),(207/50))*(Cy + (delta_Cy_r30*rud)/30 + (2*ail*(delta_Cy_a20 - delta_Cy_a20_lef*(lef/25 - 1)))/43 - delta_Cy_lef*(lef/25 - 1) + (B*P*(Cyp - delta_Cyp_lef*(lef/25 - 1)))/(2*vt) + (B*R*(Cyr - delta_Cyr_lef*(lef/25 - 1)))/(2*vt)))/m) + sa*sb*vt*(cphi*ct*g - P*sb*vt + Q*ca*cb*vt + (0.0011885000000000001032368635023317*S*pow(vt,2)*pow((1.0 - 0.0000070299999999999996077638432512291*alt),(207/50))*(Cz - delta_Cz_lef*(lef/25 - 1) + (Q*cbar*(Czq - delta_Cz_lef*(lef/25 - 1)))/(2*vt)))/m) + ca*sb*vt*(T/m - g*st + R*sb*vt - Q*cb*sa*vt + (0.0011885000000000001032368635023317*S*pow(vt,2)*pow((1.0 - 0.0000070299999999999996077638432512291*alt),(207/50))*(Cx - delta_Cx_lef*(lef/25 - 1) + (Q*cbar*(Cxq - delta_Cxq_lef*(lef/25 - 1)))/(2*vt)))/m) + cb*sa*vt*(P*cb*vt + Q*ca*sb*vt) - ca*cb*vt*(R*cb*vt + Q*sa*sb*vt)) + cb*(sb*vt*(ct*g*sphi - R*ca*cb*vt + P*cb*sa*vt + (0.0011885000000000001032368635023317*S*pow(vt,2)*pow((1.0 - 0.0000070299999999999996077638432512291*alt),(207/50))*(Cy + (delta_Cy_r30*rud)/30 + (2*ail*(delta_Cy_a20 - delta_Cy_a20_lef*(lef/25 - 1)))/43 - delta_Cy_lef*(lef/25 - 1) + (B*P*(Cyp - delta_Cyp_lef*(lef/25 - 1)))/(2*vt) + (B*R*(Cyr - delta_Cyr_lef*(lef/25 - 1)))/(2*vt)))/m) + cb*sa*vt*(cphi*ct*g - P*sb*vt + Q*ca*cb*vt + (0.0011885000000000001032368635023317*S*pow(vt,2)*pow((1.0 - 0.0000070299999999999996077638432512291*alt),(207/50))*(Cz - delta_Cz_lef*(lef/25 - 1) + (Q*cbar*(Czq - delta_Cz_lef*(lef/25 - 1)))/(2*vt)))/m) + ca*cb*vt*(T/m - g*st + R*sb*vt - Q*cb*sa*vt + (0.0011885000000000001032368635023317*S*pow(vt,2)*pow((1.0 - 0.0000070299999999999996077638432512291*alt),(207/50))*(Cx - delta_Cx_lef*(lef/25 - 1) + (Q*cbar*(Cxq - delta_Cxq_lef*(lef/25 - 1)))/(2*vt)))/m)))/(cb*pow(vt,2)) - (sb*(sb*(sb*vt*(ct*g*sphi - R*ca*cb*vt + P*cb*sa*vt + (0.0011885000000000001032368635023317*S*pow(vt,2)*pow((1.0 - 0.0000070299999999999996077638432512291*alt),(207/50))*(Cy + (delta_Cy_r30*rud)/30 + (2*ail*(delta_Cy_a20 - delta_Cy_a20_lef*(lef/25 - 1)))/43 - delta_Cy_lef*(lef/25 - 1) + (B*P*(Cyp - delta_Cyp_lef*(lef/25 - 1)))/(2*vt) + (B*R*(Cyr - delta_Cyr_lef*(lef/25 - 1)))/(2*vt)))/m) + cb*sa*vt*(cphi*ct*g - P*sb*vt + Q*ca*cb*vt + (0.0011885000000000001032368635023317*S*pow(vt,2)*pow((1.0 - 0.0000070299999999999996077638432512291*alt),(207/50))*(Cz - delta_Cz_lef*(lef/25 - 1) + (Q*cbar*(Czq - delta_Cz_lef*(lef/25 - 1)))/(2*vt)))/m) + ca*cb*vt*(T/m - g*st + R*sb*vt - Q*cb*sa*vt + (0.0011885000000000001032368635023317*S*pow(vt,2)*pow((1.0 - 0.0000070299999999999996077638432512291*alt),(207/50))*(Cx - delta_Cx_lef*(lef/25 - 1) + (Q*cbar*(Cxq - delta_Cxq_lef*(lef/25 - 1)))/(2*vt)))/m)) - vt*(ct*g*sphi - R*ca*cb*vt + P*cb*sa*vt + (0.0011885000000000001032368635023317*S*pow(vt,2)*pow((1.0 - 0.0000070299999999999996077638432512291*alt),(207/50))*(Cy + (delta_Cy_r30*rud)/30 + (2*ail*(delta_Cy_a20 - delta_Cy_a20_lef*(lef/25 - 1)))/43 - delta_Cy_lef*(lef/25 - 1) + (B*P*(Cyp - delta_Cyp_lef*(lef/25 - 1)))/(2*vt) + (B*R*(Cyr - delta_Cyr_lef*(lef/25 - 1)))/(2*vt)))/m)))/(pow(cb,2)*pow(vt,2)); 
jac[105] = (vt*(cb*sa*vt + (0.00059425000000000005161843175116587*B*S*vt*(Cyp - delta_Cyp_lef*(lef/25 - 1))*pow((1.0 - 0.0000070299999999999996077638432512291*alt),(207/50)))/m) - sb*(sb*vt*(cb*sa*vt + (0.00059425000000000005161843175116587*B*S*vt*(Cyp - delta_Cyp_lef*(lef/25 - 1))*pow((1.0 - 0.0000070299999999999996077638432512291*alt),(207/50)))/m) - cb*sa*sb*pow(vt,2)))/(cb*pow(vt,2)); 
jac[106] = (sb*(ca*cb*vt*(cb*sa*vt - (0.00059425000000000005161843175116587*S*cbar*vt*(Cxq - delta_Cxq_lef*(lef/25 - 1))*pow((1.0 - 0.0000070299999999999996077638432512291*alt),(207/50)))/m) - cb*sa*vt*(ca*cb*vt + (0.00059425000000000005161843175116587*S*cbar*vt*(Czq - delta_Cz_lef*(lef/25 - 1))*pow((1.0 - 0.0000070299999999999996077638432512291*alt),(207/50)))/m)))/(cb*pow(vt,2)); 
jac[107] = -(vt*(ca*cb*vt - (0.00059425000000000005161843175116587*B*S*vt*(Cyr - delta_Cyr_lef*(lef/25 - 1))*pow((1.0 - 0.0000070299999999999996077638432512291*alt),(207/50)))/m) - sb*(sb*vt*(ca*cb*vt - (0.00059425000000000005161843175116587*B*S*vt*(Cyr - delta_Cyr_lef*(lef/25 - 1))*pow((1.0 - 0.0000070299999999999996077638432512291*alt),(207/50)))/m) - ca*cb*sb*pow(vt,2)))/(cb*pow(vt,2)); 
jac[108] = 0; 
jac[109] = 0; 
jac[110] = (0.000000034590341700000001074671459439479*B*Jz*S*pow(vt,2)*pow((1.0 - 0.0000070299999999999996077638432512291*alt),(157/50))*(Cl + beta*delta_Clbeta + (delta_Cl_r30*rud)/30 + (2*ail*(delta_Cl_a20 - delta_Cl_a20_lef*(lef/25 - 1)))/43 - delta_Cl_lef*(lef/25 - 1) + (B*P*(Clp - delta_Clp_lef*(lef/25 - 1)))/(2*vt) + (B*R*(Clr - delta_Clr_lef*(lef/25 - 1)))/(2*vt)) + 0.000000034590341700000001074671459439479*B*Jxz*S*pow(vt,2)*pow((1.0 - 0.0000070299999999999996077638432512291*alt),(157/50))*(Cn + beta*delta_Cnbeta + (delta_Cn_r30*rud)/30 + (2*ail*(delta_Cn_a20 - delta_Cn_a20_lef*(lef/25 - 1)))/43 - delta_Cn_lef*(lef/25 - 1) - (cbar*(xcgr - 1/4)*(Cy + (delta_Cy_r30*rud)/30 + (2*ail*(delta_Cy_a20 - delta_Cy_a20_lef*(lef/25 - 1)))/43 - delta_Cy_lef*(lef/25 - 1) + (B*P*(Cyp - delta_Cyp_lef*(lef/25 - 1)))/(2*vt) + (B*R*(Cyr - delta_Cyr_lef*(lef/25 - 1)))/(2*vt)))/B + (B*P*(Cnp - delta_Cnp_lef*(lef/25 - 1)))/(2*vt) + (B*R*(Cnr - delta_Cnr_lef*(lef/25 - 1)))/(2*vt)))/(pow(Jxz,2) - Jx*Jz); 
jac[111] = 0; 
jac[112] = 0; 
jac[113] = 0; 
jac[114] = -(0.0023770000000000002064737270046635*B*Jz*S*vt*pow((1.0 - 0.0000070299999999999996077638432512291*alt),(207/50))*(Cl + beta*delta_Clbeta + (delta_Cl_r30*rud)/30 + (2*ail*(delta_Cl_a20 - delta_Cl_a20_lef*(lef/25 - 1)))/43 - delta_Cl_lef*(lef/25 - 1) + (B*P*(Clp - delta_Clp_lef*(lef/25 - 1)))/(2*vt) + (B*R*(Clr - delta_Clr_lef*(lef/25 - 1)))/(2*vt)) + 0.0023770000000000002064737270046635*B*Jxz*S*vt*pow((1.0 - 0.0000070299999999999996077638432512291*alt),(207/50))*(Cn + beta*delta_Cnbeta + (delta_Cn_r30*rud)/30 + (2*ail*(delta_Cn_a20 - delta_Cn_a20_lef*(lef/25 - 1)))/43 - delta_Cn_lef*(lef/25 - 1) - (cbar*(xcgr - 1/4)*(Cy + (delta_Cy_r30*rud)/30 + (2*ail*(delta_Cy_a20 - delta_Cy_a20_lef*(lef/25 - 1)))/43 - delta_Cy_lef*(lef/25 - 1) + (B*P*(Cyp - delta_Cyp_lef*(lef/25 - 1)))/(2*vt) + (B*R*(Cyr - delta_Cyr_lef*(lef/25 - 1)))/(2*vt)))/B + (B*P*(Cnp - delta_Cnp_lef*(lef/25 - 1)))/(2*vt) + (B*R*(Cnr - delta_Cnr_lef*(lef/25 - 1)))/(2*vt)) - 0.0011885000000000001032368635023317*B*Jz*S*pow(vt,2)*((B*P*(Clp - delta_Clp_lef*(lef/25 - 1)))/(2*pow(vt,2)) + (B*R*(Clr - delta_Clr_lef*(lef/25 - 1)))/(2*pow(vt,2)))*pow((1.0 - 0.0000070299999999999996077638432512291*alt),(207/50)) - 0.0011885000000000001032368635023317*B*Jxz*S*pow(vt,2)*pow((1.0 - 0.0000070299999999999996077638432512291*alt),(207/50))*((B*P*(Cnp - delta_Cnp_lef*(lef/25 - 1)))/(2*pow(vt,2)) - (cbar*((B*P*(Cyp - delta_Cyp_lef*(lef/25 - 1)))/(2*pow(vt,2)) + (B*R*(Cyr - delta_Cyr_lef*(lef/25 - 1)))/(2*pow(vt,2)))*(xcgr - 1/4))/B + (B*R*(Cnr - delta_Cnr_lef*(lef/25 - 1)))/(2*pow(vt,2))))/(pow(Jxz,2) - Jx*Jz); 
jac[115] = 0; 
jac[116] = -(0.0011885000000000001032368635023317*B*Jxz*S*delta_Cnbeta*pow(vt,2)*pow((1.0 - 0.0000070299999999999996077638432512291*alt),(207/50)) + 0.0011885000000000001032368635023317*B*Jz*S*delta_Clbeta*pow(vt,2)*pow((1.0 - 0.0000070299999999999996077638432512291*alt),(207/50)))/(pow(Jxz,2) - Jx*Jz); 
jac[117] = -(Jxz*Q*(Jx - Jy + Jz) + 0.0011885000000000001032368635023317*B*Jxz*S*pow(vt,2)*((B*(Cnp - delta_Cnp_lef*(lef/25 - 1)))/(2*vt) - (cbar*(Cyp - delta_Cyp_lef*(lef/25 - 1))*(xcgr - 1/4))/(2*vt))*pow((1.0 - 0.0000070299999999999996077638432512291*alt),(207/50)) + 0.00059425000000000005161843175116587*pow(B,2)*Jz*S*vt*(Clp - delta_Clp_lef*(lef/25 - 1))*pow((1.0 - 0.0000070299999999999996077638432512291*alt),(207/50)))/(pow(Jxz,2) - Jx*Jz); 
jac[118] = -(Heng*Jxz - R*(pow(Jxz,2) - Jz*(Jy - Jz)) + Jxz*P*(Jx - Jy + Jz))/(pow(Jxz,2) - Jx*Jz); 
jac[119] = -(0.0011885000000000001032368635023317*B*Jxz*S*pow(vt,2)*((B*(Cnr - delta_Cnr_lef*(lef/25 - 1)))/(2*vt) - (cbar*(Cyr - delta_Cyr_lef*(lef/25 - 1))*(xcgr - 1/4))/(2*vt))*pow((1.0 - 0.0000070299999999999996077638432512291*alt),(207/50)) - Q*(pow(Jxz,2) - Jz*(Jy - Jz)) + 0.00059425000000000005161843175116587*pow(B,2)*Jz*S*vt*(Clr - delta_Clr_lef*(lef/25 - 1))*pow((1.0 - 0.0000070299999999999996077638432512291*alt),(207/50)))/(pow(Jxz,2) - Jx*Jz); 
jac[120] = 0; 
jac[121] = 0; 
jac[122] = -(0.000000034590341700000001074671459439479*S*cbar*pow(vt,2)*pow((1.0 - 0.0000070299999999999996077638432512291*alt),(157/50))*(delta_Cm + delta_Cm_ds + Cm*eta_el + (xcgr - 1/4)*(Cz - delta_Cz_lef*(lef/25 - 1) + (Q*cbar*(Czq - delta_Cz_lef*(lef/25 - 1)))/(2*vt)) - delta_Cm_lef*(lef/25 - 1) + (Q*cbar*(Cmq - delta_Cmq_lef*(lef/25 - 1)))/(2*vt)))/Jy; 
jac[123] = 0; 
jac[124] = 0; 
jac[125] = 0; 
jac[126] = -(0.0011885000000000001032368635023317*S*cbar*pow(vt,2)*((Q*cbar*(Cmq - delta_Cmq_lef*(lef/25 - 1)))/(2*pow(vt,2)) + (Q*cbar*(Czq - delta_Cz_lef*(lef/25 - 1))*(xcgr - 1/4))/(2*pow(vt,2)))*pow((1.0 - 0.0000070299999999999996077638432512291*alt),(207/50)) - 0.0023770000000000002064737270046635*S*cbar*vt*pow((1.0 - 0.0000070299999999999996077638432512291*alt),(207/50))*(delta_Cm + delta_Cm_ds + Cm*eta_el + (xcgr - 1/4)*(Cz - delta_Cz_lef*(lef/25 - 1) + (Q*cbar*(Czq - delta_Cz_lef*(lef/25 - 1)))/(2*vt)) - delta_Cm_lef*(lef/25 - 1) + (Q*cbar*(Cmq - delta_Cmq_lef*(lef/25 - 1)))/(2*vt)))/Jy; 
jac[127] = 0; 
jac[128] = 0; 
jac[129] = -(R*(Jx - Jz) + 2*Jxz*P)/Jy; 
jac[130] = (0.0011885000000000001032368635023317*S*cbar*pow(vt,2)*((cbar*(Cmq - delta_Cmq_lef*(lef/25 - 1)))/(2*vt) + (cbar*(Czq - delta_Cz_lef*(lef/25 - 1))*(xcgr - 1/4))/(2*vt))*pow((1.0 - 0.0000070299999999999996077638432512291*alt),(207/50)))/Jy; 
jac[131] = -(Heng + P*(Jx - Jz) - 2*Jxz*R)/Jy; 
jac[132] = 0; 
jac[133] = 0; 
jac[134] = (0.000000034590341700000001074671459439479*B*Jxz*S*pow(vt,2)*pow((1.0 - 0.0000070299999999999996077638432512291*alt),(157/50))*(Cl + beta*delta_Clbeta + (delta_Cl_r30*rud)/30 + (2*ail*(delta_Cl_a20 - delta_Cl_a20_lef*(lef/25 - 1)))/43 - delta_Cl_lef*(lef/25 - 1) + (B*P*(Clp - delta_Clp_lef*(lef/25 - 1)))/(2*vt) + (B*R*(Clr - delta_Clr_lef*(lef/25 - 1)))/(2*vt)) + 0.000000034590341700000001074671459439479*B*Jx*S*pow(vt,2)*pow((1.0 - 0.0000070299999999999996077638432512291*alt),(157/50))*(Cn + beta*delta_Cnbeta + (delta_Cn_r30*rud)/30 + (2*ail*(delta_Cn_a20 - delta_Cn_a20_lef*(lef/25 - 1)))/43 - delta_Cn_lef*(lef/25 - 1) - (cbar*(xcgr - 1/4)*(Cy + (delta_Cy_r30*rud)/30 + (2*ail*(delta_Cy_a20 - delta_Cy_a20_lef*(lef/25 - 1)))/43 - delta_Cy_lef*(lef/25 - 1) + (B*P*(Cyp - delta_Cyp_lef*(lef/25 - 1)))/(2*vt) + (B*R*(Cyr - delta_Cyr_lef*(lef/25 - 1)))/(2*vt)))/B + (B*P*(Cnp - delta_Cnp_lef*(lef/25 - 1)))/(2*vt) + (B*R*(Cnr - delta_Cnr_lef*(lef/25 - 1)))/(2*vt)))/(pow(Jxz,2) - Jx*Jz); 
jac[135] = 0; 
jac[136] = 0; 
jac[137] = 0; 
jac[138] = -(0.0023770000000000002064737270046635*B*Jxz*S*vt*pow((1.0 - 0.0000070299999999999996077638432512291*alt),(207/50))*(Cl + beta*delta_Clbeta + (delta_Cl_r30*rud)/30 + (2*ail*(delta_Cl_a20 - delta_Cl_a20_lef*(lef/25 - 1)))/43 - delta_Cl_lef*(lef/25 - 1) + (B*P*(Clp - delta_Clp_lef*(lef/25 - 1)))/(2*vt) + (B*R*(Clr - delta_Clr_lef*(lef/25 - 1)))/(2*vt)) + 0.0023770000000000002064737270046635*B*Jx*S*vt*pow((1.0 - 0.0000070299999999999996077638432512291*alt),(207/50))*(Cn + beta*delta_Cnbeta + (delta_Cn_r30*rud)/30 + (2*ail*(delta_Cn_a20 - delta_Cn_a20_lef*(lef/25 - 1)))/43 - delta_Cn_lef*(lef/25 - 1) - (cbar*(xcgr - 1/4)*(Cy + (delta_Cy_r30*rud)/30 + (2*ail*(delta_Cy_a20 - delta_Cy_a20_lef*(lef/25 - 1)))/43 - delta_Cy_lef*(lef/25 - 1) + (B*P*(Cyp - delta_Cyp_lef*(lef/25 - 1)))/(2*vt) + (B*R*(Cyr - delta_Cyr_lef*(lef/25 - 1)))/(2*vt)))/B + (B*P*(Cnp - delta_Cnp_lef*(lef/25 - 1)))/(2*vt) + (B*R*(Cnr - delta_Cnr_lef*(lef/25 - 1)))/(2*vt)) - 0.0011885000000000001032368635023317*B*Jxz*S*pow(vt,2)*((B*P*(Clp - delta_Clp_lef*(lef/25 - 1)))/(2*pow(vt,2)) + (B*R*(Clr - delta_Clr_lef*(lef/25 - 1)))/(2*pow(vt,2)))*pow((1.0 - 0.0000070299999999999996077638432512291*alt),(207/50)) - 0.0011885000000000001032368635023317*B*Jx*S*pow(vt,2)*pow((1.0 - 0.0000070299999999999996077638432512291*alt),(207/50))*((B*P*(Cnp - delta_Cnp_lef*(lef/25 - 1)))/(2*pow(vt,2)) - (cbar*((B*P*(Cyp - delta_Cyp_lef*(lef/25 - 1)))/(2*pow(vt,2)) + (B*R*(Cyr - delta_Cyr_lef*(lef/25 - 1)))/(2*pow(vt,2)))*(xcgr - 1/4))/B + (B*R*(Cnr - delta_Cnr_lef*(lef/25 - 1)))/(2*pow(vt,2))))/(pow(Jxz,2) - Jx*Jz); 
jac[139] = 0; 
jac[140] = -(0.0011885000000000001032368635023317*B*Jx*S*delta_Cnbeta*pow(vt,2)*pow((1.0 - 0.0000070299999999999996077638432512291*alt),(207/50)) + 0.0011885000000000001032368635023317*B*Jxz*S*delta_Clbeta*pow(vt,2)*pow((1.0 - 0.0000070299999999999996077638432512291*alt),(207/50)))/(pow(Jxz,2) - Jx*Jz); 
jac[141] = -(Q*(pow(Jxz,2) + Jx*(Jx - Jy)) + 0.0011885000000000001032368635023317*B*Jx*S*pow(vt,2)*((B*(Cnp - delta_Cnp_lef*(lef/25 - 1)))/(2*vt) - (cbar*(Cyp - delta_Cyp_lef*(lef/25 - 1))*(xcgr - 1/4))/(2*vt))*pow((1.0 - 0.0000070299999999999996077638432512291*alt),(207/50)) + 0.00059425000000000005161843175116587*pow(B,2)*Jxz*S*vt*(Clp - delta_Clp_lef*(lef/25 - 1))*pow((1.0 - 0.0000070299999999999996077638432512291*alt),(207/50)))/(pow(Jxz,2) - Jx*Jz); 
jac[142] = -(P*(pow(Jxz,2) + Jx*(Jx - Jy)) + Heng*Jx - Jxz*R*(Jx - Jy + Jz))/(pow(Jxz,2) - Jx*Jz); 
jac[143] = -(0.0011885000000000001032368635023317*B*Jx*S*pow(vt,2)*((B*(Cnr - delta_Cnr_lef*(lef/25 - 1)))/(2*vt) - (cbar*(Cyr - delta_Cyr_lef*(lef/25 - 1))*(xcgr - 1/4))/(2*vt))*pow((1.0 - 0.0000070299999999999996077638432512291*alt),(207/50)) - Jxz*Q*(Jx - Jy + Jz) + 0.00059425000000000005161843175116587*pow(B,2)*Jxz*S*vt*(Clr - delta_Clr_lef*(lef/25 - 1))*pow((1.0 - 0.0000070299999999999996077638432512291*alt),(207/50)))/(pow(Jxz,2) - Jx*Jz); 


free(temp);

}; /*##### END of nlplant() ####*/