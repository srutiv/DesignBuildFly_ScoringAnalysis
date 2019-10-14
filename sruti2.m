%aerodynamics: input lt, tt, vs; iterate: mt, b<5, S, power, thrust, m_prop
%output: banner length, vc, n laps
%use the outputs for scoring

% Q = struct;
% global Q
% 
% Q.S = []; Q.vc = []; Q.vt = [];
% Q.x = []; Q.laps = []; Q.mt = []; Q.b = []; Q.P = []; Q.T = [];

for mt = 1:1:10
    for b = 1:0.2:3
        for P = 250:50:500
            for T = 15:5:50
                [S, vc, vt, x, laps, mt, b, P, T] = calculate_values(mt,b,P,T)
                print(
            end
        end
    end
end

% S_array = array(Q.S);
% vc_array = Q.vc;
% vt_array = Q.vt;
% x_array = Q.x;
% laps_array = Q.laps;
% mt_array = Q.mt;
% b_array = Q.b;
% P_array = Q.P;
% T_array = Q.T; 
% 
% T = table(S_array, vc_array, vt_array, ...
%     x_array, laps_array, mt_array, b_array, P_array, T_array);

function [S, vc, vt, x, laps, mt, b, P, T] = calculate_values(mt,b,P,T)
%global Q
%constants/fixed parameters
%S = area, vc = cruise velocity, x = banner length, laps = number of laps

e = 0.95; %Oswald spanwise efficiency
rho = 1.180; %density in wichita kg/m^3
g = 9.81; %gravitational acceleration in m/s^2
nu = 2; %load factor
Clmax = 1.44; %max coefficient of lift
Cd0 = 0.003; %zero lift coefficient of drag
x = 1; %banner length
vs = 10; %stall speed in m/s
lt = 6.096; %takeoff distance 20ft in m
mu = 0.02; %rolling friction during taxi
m_pass = 0; %passenger mass
m_lugg = 0; %luggage mass

%prop shit

% mu_bat = 0.023;
eff_motor = 0.4; %overall propulsive effiiency
m_bat = 0.17803 ; %battery weight in kg; Venom Lipo
voltage = 9; %total battery voltage; Venom Lipo
I = 42; %current in amps; Venom Lipo
% m_bat = (mu_bat*P)/(10.8*eff_motor); %battery mass
m_mot = 0.0425 ; %motor weight
mprop = m_mot + m_bat; %total propulsion system mass

S = (2*mt*g)/(Clmax*rho*vs^2); %function of mt and vs; in m/s
AR = b^2/S; %aspect ratio

%structural shit
mo = 1.2*(-5.871+.8538*(S+.5976)+.03113*(S+13.12).^2); %poly fit for structural weight; input ft^2, output lb %mo = 0.454 * mo_imp; % structural weight; in kg
m_pay = (mt - mprop - mo); %mass ALLOCATED for banner and payload

%drag
Cdb = 0.0001; %banner drag; rando
Clc = 0.01; %idk how to do this because dep. on velocity
Cdc = Cd0 + (Clc^2/(pi*0.9*0.8)) + x*Cdb; %cruise Cd; Cdb = drag of banner per unit length
Cdmax = Cd0 + (Clmax^2/(pi*40*0.8));
%incorporate banner somehow into drag (which will drive total drag and vc)

A = sqrt((2*T)/(Cd0*rho*S));
B = sqrt((T*Cd0*rho*S)/(2*mt^2));

tk = (1/B)*acosh(exp(lt*(B/A))); %tk = takeoff time
vt = A*tanh(B*tk); %vt = takeoff velocity

if vt < 1.3*vs 
     vt = 0; %the plane cannot takeoff
end

%cruise velocity when propulsion power is equal to drag power
%Cd0*vc^4 - (2*P/(rho*S))*vc + ((4*m^2*g^2)/(rho^2*S*b^2*pi*e)) = 0
%Cd0*vb^4 - (2*P/(rho*S))*vb + ((4*nu^2*m^2*g^2)/(rho^2*S*b^2*pi*e)) = 0
%use numerical root finder
vc_mat = [1.44 0 0 ((-2*P)/(rho*S)) ((4*mt^2*g^2)/(rho^2*S*b^2*pi*e))]; %what do we use for drag here?
vb_mat = [1.44 0 0 ((-2*P)/(rho*S)) ((4*nu^2*mt^2*g^2)/(rho^2*S*b^2*pi*e))]; %what do we use for drag here?

vc = max(real(roots(vc_mat)));
vb = max(real(roots(vb_mat)));

tt = (2000/vc) + ((2*pi*mt)/(Clmax*rho*5*vb)); %total lap time function of power, total mass
laps = ceil(600/tt); %actual number laps the plane can take in 10 minutes; optimal plane x=N

%[S, vc, vt, x, laps, mt, b, P, T];

% Q.S = [Q.S S]; 
% Q.vc = [Q.vc vc];
% Q.vt = [Q.vt vt];
% Q.x  = [Q.x x];
% Q.laps = [Q.laps laps];
% Q.mt = [Q.mt mt];
% Q.b = [Q.b b];
% Q.P = [Q.P P];
% Q.T = [Q.T T];

end

% function [S,lt,tk] = TakeoffOde(S,Clmax,rho,e,g,Cd0,lt,accuracy)
%     %this script takes in motor and battery power ratings, along with masses of
%     %the plane and a previous iteration's span and surface area. it iterates
%     %over spans and surface areas to minimize span, then surface area, finding
%     %an aircraft just able to reach v=Vto within R meters - Yamaan Atiq
%     %03/18/18
%     %Smax = 0.5*(1+obj.tap)*obj.b/2; % max area based on rules
% 
%     stepsize = 0.1; % initial stepsize for area iteration
%     grounded = 0; % can the plane take off
% 
%     % Iterate down the wing area until it can't takeoff
%     while stepsize >= accuracy && grounded
%         cwto=0; %condition stands for "can we take off?", follows boolean logic
%         while cwto==0
%             if find
%                 S = S+stepsize;
%             end
%             
%             vs = sqrt(2*m*g/(rho*S*Clmax)); % stall speed
%             vt=1.3*vs; %takeoff velocity in m/s
% 
%             z0 = [0;0];
%             t_array = linspace(0,10,1000);
%             [tout,zout] = ode45(@take_off_run(t, z, e, Clmax, rho, Cd0, b, S, m,T),t_array,z0);
% 
%             x = zout(:,1);
%             v = zout(:,2);
% 
%             % takeoffindex is first index where v> vto - i.e. it has taken off
%             takeoffindex = 500;
%             for i=1:length(v)
%                 if v(i) >= vt
%                     takeoffindex = i;
%                     break
%                 end
%             end
%             %takeoffindex = find(intermediate,1,'first');
%             lt = x(takeoffindex);
%             tk = tout(takeoffindex);
%             if lt<lt && find
%                 cwto=1;
%                 S = S-stepsize;
%                 stepsize = stepsize/10;
%             end
%         end
%     end
% end
%         
% function zd = take_off_run(~,t, z, e, Clmax, rho, Cd0, b, S, m, T)
%     AR=b^2/S;
%     K=(1/(e*pi*AR));
%     Cd=Cd0+K*Clmax^2;
% 
%     x=z(1);
%     xd=z(2);
%     xdd=(0.(1/m)*8*T-0.5*rho*xd^2*S*Cd);
% 
%     zd=[xd; xdd];
% end
% 
% function total_score = calc_score(x,laps)
%     GM = 1;
%     M1 = 1;
%     M2 = 1;
%     M3 = laps*x;
%     total_score = GM + M1 + M2 + M3;
% end
              
    