%{
mission 3 scoring
we want to maximize banner length AND # of laps

equation 1: Takeoff: iterate over mtot and thrust; fixed lt, b, tt
 --> prop capacity (mass, wattage, thrust) --> vt

equation 2 Flight: 10min = (course distance)(#of laps)(vc)
iterate over #of laps --> vcruise (where banner length comes in) --> vt

--> match vt from equations 1 and 2
%}


%equation 1
D = (1/2)*Cd*rho*S*v^2; %drag is a function of velocity
Cd = Cd0 + Cdi; %Cd is a sum of airfoil (zero lift drag) and lift induced drag
Cdi = Clmax^2/(pi*(AR)*e);
AR = b^2/S;

classdef aircraftIter
    properties
        airfoil % name of airfoil
        n % max load factor; n= L/W
        m % total mass
        P % mechanic power
        T % static thrust
        S % wing planform area
        b % wing span
        e % efficiency factor
        tap % taper ratio
        mp % mass of propulsion system
        mn % mass of payloads
        m0 % mass of structure
        N % number of payloads
        x % number of laps
        vc % straight speed
        vb % turn speed
        vt % takeoff speed
        ts % straight time
        tb % turn time
        tl % lap time
        tt % takeoff time
        grounded
    end
    methods
        % constructor
        function obj = aircraftIter(m,P,T,b,airfoil,airfoilData)
            obj.m = m; % total mass
            obj.P = P; % mechanic power
            obj.T = T; % static thrust
            obj.b = b; % wing span
            obj.grounded = 0;
            
            % unpacking airfoil data
            obj.airfoil = airfoil; % airfoil name
            Clmaxt = airfoilData(1); % max Cl for takeoff
            Cd0t = airfoilData(2); % Cd0 during takeoff
            Cd0c = airfoilData(3); % Cd0 during cruise
            foilD = airfoilData(4); % arclength
            foilA = airfoilData(5); % cross sectional area
            
            % constants
            rho = 1.225; % standard air density kg/m^3
            %rho = 1.084; % density in Tucson in April: 80F, 2,389ft  kg/m^3
            Cdon = 0.05; %drag coefficient of american football
            An = 0.0456; % frontal area of dart (3 inch diameter circle)
            obj.e = 0.9; % oswald efficiency factor
            obj.tap = 1; %0.45; % taper ratio
            g = 9.81; % gravitational acceleration
            tod = 6; % takeoff distance
            obj.n = 2; % load factor
            Teff = 1.0; % thrust efficiency for takeoff
            fos = 1.2; % factor of safety on takeoff velocity
            
            % calculations
            obj=obj.Takeoff(Clmaxt,rho,obj.e,g,Cd0t,tod,Teff,fos); % get wing area
            obj=obj.mass(foilD,foilA); % get all the masses
            obj=obj.velocities(rho,g,Cd0c,Cdon,An); % get the cruise, and turn velocities
            obj=obj.times(g); % use it all to get times
        end
        % calculate things
        function obj = mass(obj,foilD,foilA)
            rootChord = 2*obj.S/obj.b/(1+obj.tap);
            tipChord = rootChord*obj.tap;
            
            % propulsion mass (Power method)
            eff_motor = 0.4;
            V_cell = 1.2; % Voltage per cell
            cap_cell = 1.5; % mAh capacity of cell
            I_cell = cap_cell/(1/6); % max current for the cell to last 10 minutes (1/6 hour)
            mu_bat = 0.023; % cell density
            P_cell = I_cell*V_cell; % Power output of one cell
            m_bat = (mu_bat*obj.P)/(P_cell*eff_motor); % total battery mass, where P/P_cell is number of cells
            obj.mp = 0.4 + m_bat; % total propulsion system mass
            
            
            % structure mass
            mechanismMass = 0.75; % the approximate mass of all the mechanisms
            lam_pole = 0.0905; % linear density of carbon fiber rod
            sig_wrap = 0.0610; % area density of the wrap (microlite/monokote)
            mu_sheet = 160; % density of the sheeting
            thich_sheet = 0.00079375; % thickness of sheeting (usually 1/32")
            sig_skin = sig_wrap+thich_sheet*mu_sheet; % area density of the skin
            mu_rib = 160; % density of the ribs
            thick_rib = 0.0015875; % thickness of the ribs (should be 1/16")
            l_pole = 2; % length of carbon fiber rod
            ribSpacing = 0.0762; % 3" rib spacing per Amil?
            %hmax = rootChord*0.089; % thickest part of the spar (airfoil thickness)
            %[moss] = obj.wing_spar(obj.m,obj.b,hmax,obj.n); % getting wing spar details
            %sparMass = moss(2); % mass of the wing spar
            %thick = 0.003175;
            sparDensity = 0.125;
            sparMass = sparDensity*obj.b;
            wettedA = obj.S*foilD+2*tipChord*foilA; % total wetted area of the wing
            ribA = obj.S/obj.b*foilA; % area of the average rib
            n_rib = obj.b/ribSpacing; % total number of ribs
            ribMass = n_rib*ribA*mu_rib*thick_rib; % mass of all the ribs
            skinMass = wettedA*sig_skin; % mass of the wing skin
            poleMass = lam_pole*l_pole; % mass of the carbon fiber rod
            obj.m0 = sparMass+skinMass+ribMass+2*poleMass+mechanismMass; % total structural mass
            
            % mass ALLOCATED for payload
            obj.mn = (obj.m - obj.mp - obj.m0);
        end
        function [Clmaxt, Cd0t, Cd0c] = airfoilData(obj)
            AR = 8;
            airfoilF = strcat(obj.airfoil,'_flap');
            % Takeoff
            Re = 200000;
            path = strcat('./polar/',airfoilF,'_',num2str(Re),'.txt');
            [alpha,clt,cdt] = parsePolar(path);
            cl0t=min(abs(clt)); % smallest cl in file (closest to zero lift)
%             alpha0=alpha(abs(clt)==cl0t); % the zero lift angle
            Cd0t=cdt(abs(clt)==cl0t); % Cd at zero lift
            Clmaxt=max(clt); % 2d wing maximum CL
%             Clmax25=clt(clt==min(abs(clt-cl0t*0.25)));
%             alpha25=clt(clt==Clmax25);
%             Clmax75=clt(clt==min(abs(clt-cl0t*0.75)));
%             alpha75=clt(clt==Clmax75);
%             Cla=(Clmax75-Clmax25)/(alpha75-alpha25)
%             CLa=Cla*(AR/(2+sqrt(4+AR^2))); % lift curve slope for 3d wing
%             Clmaxt=CLa*(alpha(clt==cl0t)-alpha0); % 3d wing coefficient of lift (AR>3)
            CLmaxt=0.9*Clmaxt; % approximation
            
            % Cruise
            Re = 500000;
            path = strcat('./polar/',obj.airfoil,'_',num2str(Re),'.txt');
            [alpha,cl,cd] = parsePolar(path);
            cl0c=min(abs(cl)); % smallest cl in file (closest to zero lift)
            Cd0c=cd(abs(cl)==cl0c); % Cd at zero lift for cruise conditions
        end
        function obj = TakeoffOde(obj,Clmax,rho,e,g,Cd0,R,accuracy)
            %this script takes in motor and battery power ratings, along with masses of
            %the plane and a previous iteration's span and surface area. it iterates
            %over spans and surface areas to minimize span, then surface area, finding
            %an aircraft just able to reach v=Vto within R meters - Yamaan Atiq
            %03/18/18
            %Smax = 0.5*(1+obj.tap)*obj.b/2; % max area based on rules
            Smax = 1.3;
            
            stepsize = 0.1; % initial stepsize for area iteration
            obj.grounded = 0; % can the plane take off
            
            if obj.S == 0
                find = 1;
            else
                find = 0;
            end
            
            % Iterate down the wing area until it can't takeoff
            while stepsize >= accuracy && ~obj.grounded
                cwto=0; %condition stands for "can we take off?", follows boolean logic
                while cwto==0
                    if find
                        obj.S = obj.S+stepsize;
                    end
                    
                    Vs = sqrt(2*obj.m*g/(rho*obj.S*Clmax)); % stall speed
                    Vto=1.2*Vs; %m/s
                    % Vto = Vs;
                    
                    z0 = [0;0];
                    tarray = linspace(0,10,1000);
                    [tout,zout] = ode45(@(t,z) obj.take_off_run(t, z, e, Clmax, rho, Cd0, obj.b, obj.S, obj.m, obj.T),tarray,z0);
                    
                    x = zout(:,1);
                    v = zout(:,2);
                    
                    % takeoffindex is first index where v> vto - i.e. it has taken off
                    takeoffindex = 500;
                    for i=1:length(v)
                        if v(i) >= Vto
                            takeoffindex = i;
                            break
                        end
                    end
                    %takeoffindex = find(intermediate,1,'first');
                    tod = x(takeoffindex);
                    tot = tout(takeoffindex);
                    if tod<R && find
                        cwto=1;
                        obj.S = obj.S-stepsize;
                        stepsize = stepsize/10;
                    end
                    if obj.S > Smax || (tod<R && ~find)
                        obj.grounded=1; % it cannot take off
                        break
                    end
                end
            end
            obj.tt = tot;
            obj.vt = Vto;
        end
        function zd = take_off_run(~,t, z, e, Clmax, rho, Cd0, b, S, m, T)
            AR=b^2/S;
            K=(1/(e*pi*AR));
            Cd=Cd0+K*Clmax^2;
            
            x=z(1);
            xd=z(2);
            xdd=(1/m)*(0.8*T-0.5*rho*xd^2*S*Cd);
            
            zd=[xd; xdd];
        end
        function obj = Takeoff(obj,Clmax,rho,e,g,Cd0,R,Teff,fos)
            F = obj.T*Teff;
            Cd = Cd0+Clmax^2/(pi*e*6);
            obj.S = -(obj.m/(R*Cd*rho))*log(1-fos^2*obj.m*g*Cd/(Clmax*F));
            
            Smax = 10; % maximum wing area based on folding and rules, 2ft chord x 10ft span
            if obj.S>Smax || real(obj.S)~=obj.S
                obj.grounded = obj.S>Smax;
            end
            
            A = sqrt(2*F/(Cd*rho*obj.S));
            B = sqrt(F*Cd*rho*obj.S/(2*obj.m^2));
            obj.vt = fos*sqrt(obj.m*g/(0.5*rho*obj.S*Clmax));
            obj.tt = 1/B*atanh(obj.vt/A);
        end
        function [d,A] = airfoilMeasures(obj)
            a=importdata(strcat(obj.airfoil,'.dat')); % Get the polar data
            x = a.data(:,1); % split it into x
            y = a.data(:,2); % and y
            np=length(x); % find number of data points
            
            % Length
            dd=zeros(np,1); % preallocate distance array
            i=1;
            % loop over the points and calculate the distance between each
            while i < np
                i=i+1;
                dx=x(i)-x(i-1);
                dy=y(i)-y(i-1);
                dd(i)=sqrt(dx^2+dy^2); % distance formula
            end
            d=sum(dd); % sum up all the small distances to get total distance
            
            % Area
            endpoint = find(x==min(x)); % find the point where the polar reverses direction
            % split the curve to integrate
            topCurvex = flipud(x(1:endpoint)); % we reverse this so it goes from 0->1
            topCurvey = flipud(abs(y(1:endpoint))); % using abs so we get a positive area
            bottomCurvex = x(endpoint+1:np);
            bottomCurvey = abs(y(endpoint+1:np)); % using abs so we get a positive area
            % integrate both halves
            topA = trapz(topCurvex,topCurvey);
            bottomA = trapz(bottomCurvex,bottomCurvey);
            A = topA+bottomA; % add them together to get total
        end
        function obj = Oswald(obj,b,cr,airfoil,g,d,m,v,divisions)
            % Iterating taper ratios to find best e %
            wingSpan = b; % max 3m due to rules
            rootChord = cr; % max 0.5m due to folding
            wingAirfoil = strcat(airfoil,'.dat'); % must include the file extension
            mass = m; % plane mass
            es = zeros(divisions-1,2); % To store all the e's we find
            
            % Looping over lambda which will have values from 0 to 1
            for lambda = 0:1/divisions:1
                tipChord = lambda*rootChord; % by the definition of taper ratio
                wingArea = (rootChord+tipChord)*wingSpan/2; % area of a trapezoid
                wingTaper = lambda;
                
                % Inputting the geometry and run case parameters
                wing = Main_Wing(wingSpan,wingArea,wingTaper,0,wingAirfoil,0,0,[],[0 0 0],'wing');
                Case = inputs(g,d,mass,0,0,0,0,v,0,false);
                
                % Building the aircraft, .avl file, and running it
                plane = Aircraft(wing,[],[],Case,[],[0,0,0],mass,[1,1,1],'tapir');
                plane.build_file;
                plane.run_avl;
                
                % Reading the header data from the sb file
                Header = parseRunCaseHeader(strcat('.\AVL_DATA\','tapir','_DATA.st'));
                e = Header.e; % Oswald efficiency factor
                
                % Printing where the script is, among other things
                fprintf('%f %f\n', lambda, e);
                
                % Adding the result to the tracker
                es(fix(lambda*divisions)+1,:) = [lambda,e];
            end
            
            % for some reason there are random (0,0) rows in es, so we delete them
            rowsToDelete = find(es(:,2)==0); % find them!
            es(rowsToDelete,:) = []; % send them to meet their maker
            
            % Find the max e value and the corresponding lambda value
            maxe = max(es(:,2));
            tapmaxs = es((es(:,2)==maxe),1);
            tapmax = tapmaxs(fix(length(tapmaxs)/2)+1); % there might be multiple taper ratios that give the same e
            obj.tap = tapmax;
            obj.e = maxe;
        end
        function obj = velocitiesRoot(obj,rho,g,e,Cd0,Cdon,An)
            obj.N = obj.mn/0.086; % number of payloads
            Cdn = (obj.N*(Cdon*An/obj.S)); % normalized payload drag coefficient
            
            % velocity is a 4th degree polynomial based on stuff
            pvc = [(Cd0+Cdn) 0 0 -2*obj.P/(rho*obj.S) 4*obj.m^2*g^2/(rho^2*obj.S*obj.b^2*pi*e)];
            pvb = [(Cd0+Cdn) 0 0 -2*obj.P/(rho*obj.S) 4*obj.n^2*obj.m^2*g^2/(rho^2*obj.S*obj.b^2*pi*e)];
            try % find the maximum real root of the polynomial
                obj.vc = max(real(roots(pvc))); % straight speed
                obj.vb = max(real(roots(pvb))); % turn speed
            catch % if it doesn't exist make your cruise speed equal to takeoff speed
                obj.vc = 1; obj.vb = 1;
            end
            
            
        end
        function obj = velocities(obj,rho,g,Cd0,Cdon,An)
            obj.N = obj.mn/0.086; % number of payloads
            if obj.N > 0
                Cdn = (obj.N*(Cdon*An/obj.S)); % normalized payload drag coefficient
                D = (Cd0+Cdn);
                F = 2*obj.P/(rho*obj.S);
                J = 4*obj.m^2*g^2/(rho^2*obj.S*obj.b^2*pi*obj.e);
                
                A = sqrt(3)*(sqrt(27*D^2*F^4-256*D^3*J^3*obj.n^3)+9*D*F^2)^(1/3);
                B = (4*(2/3)^(1/3)*J*obj.n);
                C = (2^(1/3)*3^(2/3)*D);
                G = sqrt(B/A+A/C);
                obj.vb = 1/2*G+1/2*sqrt((2*F)/(D*G)-G);
                
                A = sqrt(3)*(sqrt(27*D^2*F^4-256*D^3*J^3)+9*D*F^2)^(1/3);
                B = (4*(2/3)^(1/3)*J);
                G = sqrt(B/A+A/C);
                obj.vc = 1/2*G+1/2*sqrt((2*F)/(D*G)-G);
            end
            if obj.N <= 0 || real(obj.vc)~=obj.vc || real(obj.vb)~=obj.vb
                obj.vb = 10^-3;
                obj.vc = 10^-3;
            end
        end
        function obj = times(obj,g)
            R = obj.vb^2/(g*sqrt(obj.n^2-1)); % radius of turn
            circ = pi*R; % length of a 180deg turn
            
            obj.ts = 305/obj.vc; % straight portion time
            obj.tb = circ/obj.vb; % turn time
            obj.tl = 2*obj.ts+4*obj.tb; %lap time function of power, total mass
            if obj.grounded
                obj.x = 0;
            else
                obj.x =(600-2)/obj.tl; %actual number laps the plane can take in 10 minutes; optimal plane x=N
            end
        end
    end
end

