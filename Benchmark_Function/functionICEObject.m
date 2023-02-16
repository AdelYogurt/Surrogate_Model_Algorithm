function fval=functionICEObject(x)
% ICE problem
%
%
%

% x=[75; 6.42; 26; 39; 7500];
% x=[83.33; 9.45; 30.99; 37.34; 6070];
b=x(1);
c_r=x(2);
d_E=x(3);
d_I=x(4);
omega=x(5)*1e-3;

% parameter
K_0=1/120;
rou=1.225;
gama=1.33;
V=1.859*1e6;
Q=43958;
N_c=4;
C_s=0.44;
A_f=14.6;

if omega > 5.25
    eta_vb=1.067-0.038*exp(omega-5.25);
else
    eta_vb=0.637+0.13*omega-0.014*omega^2+0.00066*omega^3;
end
eta_tad=0.8595*(1-c_r^(-0.33));
eta_V=eta_vb*(1+5.96*1e-3*omega^2)/(1+((9.428*1e-5)*(4*V/pi/N_c/C_s)*(omega/d_I^2))^2);
S_V=0.83*((8+4*c_r)+1.5*(c_r-1)*(pi*N_c/V)*b^3)/((2+c_r)*b);
eta_t=eta_tad-S_V*(1.5/omega)^0.5;
V_P=(8*V/pi/N_c)*omega*b^(-2);
FMEP=4.826*(c_r-9.2)+(7.97+0.253*V_P+9.7*(1e-6)*V_P^2);

fval=K_0*(FMEP-(rou*Q/A_f)*eta_t*eta_V)*omega;
end