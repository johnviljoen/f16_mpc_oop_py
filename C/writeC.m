close all
clear all

% sim variables
syms g m B S cbar xcgr xcg Heng Jy Jxz Jz Jx

% primary state variables
syms npos epos alt phi theta psi vt alpha beta P Q R

sa = sin(alpha);
ca = cos(alpha);
sb = sin(beta);
cb = cos(beta);
tb = tan(beta);

st = sin(theta);
ct = cos(theta);
tt = tan(theta);
sphi = sin(phi);
cphi = cos(phi);
spsi = sin(psi);
cpsi = cos(psi);

% actuation state variables
syms T el ail rud lef

dail = ail/21.5;
drud = rud/30.0;
dlef = (1-lef/25.0);

%% the atmos function is then called - mach and ps not used
rho0 = 2.377e-3;

tfac =1 - .703e-5*(alt);
rho=rho0*tfac^4.14;
qbar = .5*rho*vt^2;

%% navigation equations

U = vt*ca*cb;
V = vt*sb;
W = vt*sa*cb;

dot_npos = U*(ct*cpsi) + ...
            V*(sphi*cpsi*st - cphi*spsi) + ...
            W*(cphi*st*cpsi + sphi*spsi);
        
dot_epos = U*(ct*spsi) + ...
            V*(sphi*spsi*st + cphi*cpsi) + ...
            W*(cphi*st*spsi - sphi*cpsi);
        
dot_alt = U*st - V*(sphi*ct) - W*(cphi*ct);

%% kinematic equations

dot_phi = P + tt*(Q*sphi + R*cphi);

dot_theta = Q*cphi - R*sphi;

dot_psi = (Q*sphi + R*cphi)/ct;

%% table lookup occurs

% hifi C - f(alph, beta, el)
syms Cx Cz Cm Cy Cn Cl

% hifi damping - f(alpha)
syms Cxq Cyr Cyp Czq Clr Clp Cmq Cnr Cnp

% hifi C lef - f(alpha, beta)
syms delta_Cx_lef delta_Cz_lef delta_Cm_lef delta_Cy_lef delta_Cn_lef delta_Cl_lef

% hifi damping lef - f(alpha, beta)
syms delta_Cxq_lef delta_Cyr_lef delta_Cyp_lef delta_Czq_lef delta_Clr_lef delta_Clp_lef delta_Cmq_lef delta_Cnr_lef delta_Cnp_lef

% hifi rudder - f(alpha, beta)
syms delta_Cy_r30 delta_Cn_r30 delta_Cl_r30

% hifi ailerons - f(alpha, beta)
syms delta_Cy_a20 delta_Cy_a20_lef delta_Cn_a20 delta_Cn_a20_lef delta_Cl_a20 delta_Cl_a20_lef

% hifi other coeffs - f(alpha, el)
syms delta_Cnbeta delta_Clbeta delta_Cm eta_el delta_Cm_ds

%%

dXdQ = (cbar/(2*vt))*(Cxq + delta_Cxq_lef*dlef);

Cx_tot = Cx + delta_Cx_lef*dlef + dXdQ*Q;

dZdQ = (cbar/(2*vt))*(Czq + delta_Cz_lef*dlef);

Cz_tot = Cz + delta_Cz_lef*dlef + dZdQ*Q;

dMdQ = (cbar/(2*vt))*(Cmq + delta_Cmq_lef*dlef);

Cm_tot = Cm*eta_el + Cz_tot*(xcgr-xcg) + delta_Cm_lef*dlef + dMdQ*Q + delta_Cm + delta_Cm_ds;

dYdail = delta_Cy_a20 + delta_Cy_a20_lef*dlef;

dYdR = (B/(2*vt))*(Cyr + delta_Cyr_lef*dlef);

dYdP = (B/(2*vt))*(Cyp + delta_Cyp_lef*dlef);

Cy_tot = Cy + delta_Cy_lef*dlef + dYdail*dail + delta_Cy_r30*drud + dYdR*R + dYdP*P;

dNdail = delta_Cn_a20 + delta_Cn_a20_lef*dlef;

dNdR = (B/(2*vt))*(Cnr + delta_Cnr_lef*dlef);

dNdP = (B/(2*vt))*(Cnp + delta_Cnp_lef*dlef);

Cn_tot = Cn + delta_Cn_lef*dlef - Cy_tot*(xcgr-xcg)*(cbar/B) + dNdail*dail + delta_Cn_r30*drud + dNdR*R + dNdP*P + delta_Cnbeta*beta;

dLdail = delta_Cl_a20 + delta_Cl_a20_lef*dlef;

dLdR = (B/(2*vt))*(Clr + delta_Clr_lef*dlef);

dLdP = (B/(2*vt))*(Clp + delta_Clp_lef*dlef);

Cl_tot = Cl + delta_Cl_lef*dlef + dLdail*dail + delta_Cl_r30*drud + dLdR*R + dLdP*P + delta_Clbeta*beta;

Udot = R*V - Q*W - g*st + qbar*S*Cx_tot/m + T/m;

Vdot = P*W - R*U + g*ct*sphi + qbar*S*Cy_tot/m;

Wdot = Q*U - P*V + g*ct*cphi + qbar*S*Cz_tot/m;

dot_vt = (U*Udot + V*Vdot + W*Wdot)/vt;

dot_alpha = (U*Wdot - W*Udot)/(U*U + W*W);

dot_beta = (Vdot*vt - V*dot_vt)/(vt*vt*cb);

L_tot = Cl_tot*qbar*S*B; 
M_tot = Cm_tot*qbar*S*cbar;
N_tot = Cn_tot*qbar*S*B;

denom = Jx*Jz - Jxz*Jxz;

dot_P = (Jz*L_tot + Jxz*N_tot - (Jz*(Jz-Jy)+Jxz*Jxz)*Q*R + Jxz*(Jx-Jy+Jz)*P*Q + Jxz*Q*Heng)/denom;

dot_Q = (M_tot + (Jz-Jx)*P*R - Jxz*(P*P-R*R) - R*Heng)/Jy;

dot_R = (Jx*N_tot + Jxz*L_tot + (Jx*(Jx-Jy)+Jxz*Jxz)*P*Q - Jxz*(Jx-Jy+Jz)*Q*R +  Jx*Q*Heng)/denom;

syms sa ca sb cb tb st ct tt sphi cphi spsi cpsi

%%

vars = [sin(alpha), cos(alpha), sin(beta), cos(beta), tan(beta), sin(theta), cos(theta), tan(theta), sin(phi), cos(phi), sin(psi), cos(psi)];
replacements = [sa, ca, sb, cb, tb, st, ct, tt, sphi, cphi, spsi, cpsi];

vars2 = [sa^2, ca^2, sb^2, cb^2, tb^2, st^2, ct^2, tt^2, sphi^2, cphi^2, spsi^2, cpsi^2];
% replacements2 = [pow(sa,2), pow(ca,2), pow(sb,2), pow(cb,2), pow(tb,2), pow(st,2), pow(ct,2), pow(tt,2), pow(sphi,2), pow(cphi,2), pow(spsi,2), pow(cpsi,2)];

%%

% subs(diff(dot_P, R), vars, replacements)

%% list of vars to run through
dot_states = [dot_npos, dot_epos, dot_alt, dot_phi, dot_theta, dot_psi, dot_vt, dot_alpha, dot_beta, dot_P, dot_Q, dot_R];
states = [npos, epos, alt, phi, theta, psi, vt, alpha, beta, P, Q, R];

fileID = fopen('jac_C.txt','w');

k = 0;

for i = 1:12
   fprintf(fileID, '//%ss relationship with the other states \n', states(i));
   for j = 1:12
       
       formatSpec = 'dx%ddot_dx%d = %s';
       A1 = i-1;
       A2 = j-1;
       A3 = subs(diff(dot_states(i), states(j)), vars, replacements);
%        A3 = diff(dot_states(i), states(j));
%        A3 = subs(A3,vars2,replacements2);
       
%        fprintf(fileID,'dx%idot_dx%i = %s \n',A1,A2,A3);
       fprintf(fileID,'jac[%i] = %s; \n', k, A3);
       k = k + 1;
   end
   fprintf(fileID, '\n');
end

% ddot_P_da = diff(dot_P,alpha);
% 
% syms dCxda ddelta_Cnbetada
% 
% ddot_P_da = subs(ddot_P_da, diff(Cx(alpha,beta,el), alpha), dCxda)
% 
% subs(ddot_P_da, diff(delta_Cnbeta, alpha), ddelta_Cnbetada)
    
fclose(fileID);