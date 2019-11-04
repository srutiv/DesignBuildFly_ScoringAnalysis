rho = 1.225;
Cd0c = 0.008;
e = 0.95;
S = 0.1778;
b = 1.5;
P = 300;
T = 45;
mt = 3;
eff = 0.6;
Cdb = 0;
nu = 2;
g = 9.8;
Cl0c = 0.0025;
Cdc = Cd0c + (Cl0c^2/(pi*0.9*0.8)); %cruise Cd


% cruise velocity when propulsion power is equal to drag power

vc_matM2 = [Cdc 0 0 ((-2*P)/(rho*S)) ((4*mt^2*g^2)/(rho^2*S*b^2*pi*e))];
vb_matM2 = [Cdc 0 0 ((-2*P)/(rho*S)) ((4*nu^2*mt^2*g^2)/(rho^2*S*b^2*pi*e))];
rc = roots(vc_matM2);
rb = roots(vb_matM2);
vcM2 = max(rc(imag(rc)==0));
vbM2 = max(rb(imag(rb)==0));

%%%%%%%%%
v1 = P*eff/T;

%%%%%%%%%%
syms v2
%eqns = T*vcM2 - P == 0, 
eqns = [(0.5*rho*(v2^3)*S)*(Cd0c + (((2*mt*9.8/(S*rho*v2^2))^2)/(pi*e*(b^2/S)))) - P == 0];
J = solve(eqns, v2);
v2 = double(max(abs(J)));

%%%%%%
syms v3
eqns = [(0.5*rho*(v3^2)*S)*(Cd0c + (((2*mt*9.8/(S*rho*v3^2))^2)/(pi*e*(b^2/S)))) - T*0.6 == 0];
K = solve(eqns, v3);
v3 = double(max(abs(K)));

%%%%%%%%%%
syms v4
%eqns = T*vcM2 - P == 0, 
eqns = [(0.5*rho*(v4^3)*S)*(Cd0c + (((2*mt*9.8/(S*rho*v4^2))^2)/(pi*e*(b^2/S)))) - P*eff == 0];
J = solve(eqns, v4);
v4 = double(max(abs(J)));
