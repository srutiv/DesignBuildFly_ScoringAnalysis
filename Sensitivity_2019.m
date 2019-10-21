function Sensitivity_2019
% Michael Zakoworotny
% Convention: 
% ALL_CAPS: global variables and constants
% Leading_Caps: important local variables (ie. computed maximum values)
% lowercaseNames: local variables

%CONVERSIONS
mToft = 100/2.54/12;
ftToin = 12;
minTosec = 60;

resolution = 100;

% CONSTANTS
M3_TIME = 10*minTosec; %10 minutes
LAP_DISTANCE = 3570; % ft
%DETERMINE GOOD VALUES FOR MAX PARAMETERS
MAX_PASS = 12; MIN_TIME = 180;%68; %s (from 2017 script)
MAX_LAPS = 9; MAX_BAN = 0.889*mToft; %ft (equal to max wingspan)
MAX_NPASS_TIME = MAX_PASS/MIN_TIME; %M2
MAX_NLAPS_BAN = MAX_LAPS*MAX_BAN; %M3
MIN_GM = 20; %s %GM

%Baseline numbers from scoring
nPass_base = 5; 
banLen_base = 1;
speed_base = 20 * mToft; %m/s -> ft/s
GM_base = 2*minTosec;
sc_base = getScore(nPass_base, banLen_base, speed_base, GM_base);

% Parameters to vary over
range = linspace(0.5, 1.5, resolution);
nPass = floor(range*nPass_base);
banLen = range*banLen_base;
speed = range*speed_base;
GM_time = range*GM_base;

% Call functions to run
sensitivityStudy();

    function sc = getScore(nPass, banLen, speed, GM_time)
        %Parametrize with just one variable (speed)
        M2_time = 3*LAP_DISTANCE/speed;
        nLaps = floor(speed*M3_TIME/LAP_DISTANCE);
        Max_NPass_Time = max(nPass/M2_time, MAX_NPASS_TIME); %M2 criterion
        Max_NLaps_Ban = max(nLaps*banLen, MAX_NLAPS_BAN); %M3 criterion

        M1 = 1.0;
        M2 = 1 + (nPass/M2_time)/Max_NPass_Time;
        M3 = 2 + (nLaps*banLen)/Max_NLaps_Ban;
        GM = MIN_GM / GM_time;

        % Assuming successful flight each time, scores can vary from
        % 4 to 7
        sc = M1 + M2 + M3 + GM;
    end

    function sensitivityStudy()
        scores = zeros(resolution, 4);
        % Compute scores over all varied parameters
        for i = 1:resolution
           scores(i,1) = getScore(nPass(i),banLen_base,speed_base,GM_base)/sc_base; 
        end
        for i = 1:resolution
           scores(i,2) = getScore(nPass_base,banLen(i),speed_base,GM_base)/sc_base; 
        end
        for i = 1:resolution
           scores(i,3) = getScore(nPass_base,banLen_base,speed(i),GM_base)/sc_base; 
        end
        for i = 1:resolution
           scores(i,4) = getScore(nPass_base,banLen_base,speed_base,GM_time(i))/sc_base; 
        end
        
        %Plotting
        plot(range,scores);
%         figure(1); hold on;
%         plot(nPass,scores(:,1)*sc_base);
%         plot(banLen,scores(:,2)*sc_base);
%         plot(speed,scores(:,3)*sc_base);
%         plot(GM_time,scores(:,4)*sc_base);
        xlabel('Deviation of parameter from baseline');
        ylabel('Deviation of score from baseline');
        title('Sensitivity study of global competition parameters');
        legend({'Number of passengers','Banner Length','Speed','GM Time'},...
                'Location', 'eastoutside');
        
    end

end