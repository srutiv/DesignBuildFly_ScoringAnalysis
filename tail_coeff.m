AR_w = 5.16; %input('Wing aspect ratio: ');
b_w = 1.5; %input('Input wing span: '); % wingspan
S_w = 0.44;  %((b_w)^2)/AR_w; % area of wing
c_w = 0.39; % wing chord length
V_v = 0.04; %input('Input vertical volume ratio for tail (0.04 for GA SE): ');
V_h = 0.5; %input('Input horizontal volume ratio for tail (0.70 for GA SE): ');
l_v = 0.5*b_w; %input('Input distance between aerodynamic center of wing and vertical tail: ');
l_h = 0.5*b_w;%input('Input distance between aerodynamic center of wing and horizontal tail: ');
downwash = 0.22; %assuming 13 deg downwash
AR_t = 1.2;
AR_h = 1.2;
AR_v = 1.2;
eta = 0.9; %assuming high effectiveness
SM = 0.15; %15%
a0h = 2*pi; %lfit curve slope of hstab; assuming small aoa
a0w = 2*pi;
alphawc = 0; %Sweep angle at quarter chord for wing
alphahc = 0; %Sweep angle at quarter chord for tail
f = (2*0.0075)/(S_w*c_w); %fuselage pitch stiffness

S_v = (V_v*S_w*b_w)/l_v; % area of vertical tail
S_h = (V_h*S_w*c_w)/l_h; % area of horizontal tail

H_span = sqrt(AR_h*S_h); % span of horizontal tail
H_chord = S_h / H_span; % average chord length of horizontal tail

H_stab_chord = H_chord*(0.75);
Elev_chord = H_chord*(0.25);

V_chord = H_chord; % average chord length of vetical tail
V_span = sqrt(S_v/V_chord); % length of vertical tail

V_stab_chord = V_chord*(0.60);
Rud_chord = V_chord*(0.40);

claw = (pi*AR_w)/(1+sqrt(1+((pi*AR_w)/(a0w*cos(alphawc)))^2));
clah = (pi*AR_h)/(1+sqrt(1+((pi*AR_h)/(a0h*cos(alphahc)))^2));
cla = claw + clah + eta*(S_h/S_w)*(1-downwash);
cmaf = (2*f)/(S_w*c_w);
clat = 0.001; %lift of tail

% syms xnp xcg
% A = solve(0.25+((eta*V_h*clat*(1-downwash)-cmaf)/cla) - xnp == 0, ...
%     0 == (xnp/c_w) - (xcg/c_w) + SM,xnp, xcg);
% xnp = double(A.xnp);
% xcg = double(A.xcg);

xnp = 0;
syms xcg
xcg = double(solve(0 == (xnp/c_w) - (xcg/c_w) + SM,xcg));

fprintf('Area of vertical tail is %f \n',S_v)
fprintf('Span of vertical tail is %f \n',V_span)
fprintf('Average chord length of vertical stabilizer is %f \n',V_stab_chord)
fprintf('Average chord length of rudder surface is %f \n',Rud_chord)

fprintf('\n')
fprintf('Area of horizontal  tail is %f \n',S_h)
fprintf('Span of horizontal tail is %f \n',H_span)
fprintf('Average chord length of horizontal stabilizer is %f \n',H_stab_chord)
fprintf('Average chord length of elevator surface is %f \n',Elev_chord)
fprintf('\n')

fprintf('nuetral point distance in front of CL %f \n',xnp)
fprintf('CG distance in front of CL %f \n',xcg)