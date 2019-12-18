function stats = iterSensitivity(mt,P,T,span,S_sens,ban_length)
    
    global max_M2 max_M3 best_config top_score

    resolution = 100;

    %Baseline numbers from scoring
    mass_base = mt; 
    power_base = P;
    thrust_base = T; %m/s -> ft/s

    % Parameters to vary over: mt, P, T; keep b, S, banner length same
    range = linspace(0.5, 1.5, resolution);
    mass = floor(range*mass_base);
    power = range*power_base;
    thrust = range*thrust_base;
    
    %control variables: span, S_sens, ban_length
    
    scores_sens = zeros(resolution, 3);
    % Compute scores over all varied parameters
    for i = 1:resolution
       scores_sens(i,1) = getScore(mass(i),power_base,thrust_base, span,S_sens,ban_length)/top_score;
    end
    for i = 1:resolution
       scores_sens(i,2) = getScore(mass_base,power(i),thrust_base, span,S_sens,ban_length)/top_score; 
    end
    for i = 1:resolution
       scores_sens(i,3) = getScore(mass_base,power_base,thrust(i), span,S_sens,ban_length)/top_score ; 
    end   
%        check1 = isnan(scores_sens(i,1)); %returns 1 if NaN, 0 otherwise
%        check2 = isnan(scores_sens(i,2)); 
%        check3 = isnan(scores_sens(i,3)); 
       


        %stats     
        %calculate stats for each varied parameter and package
        stats.massMax = max(scores_sens(1,:)); 
        stats.massMean = mean(scores_sens(1,:));
        stats.massStd = std(scores_sens(1,:));
        stats.massVar = var(scores_sens(1,:));

        stats.powerMax = max(scores_sens(2,:)); 
        stats.powerMean = mean(scores_sens(2,:));
        stats.powerStd = std(scores_sens(2,:));
        stats.powerVar = var(scores_sens(2,:));

        stats.thrustMax = max(scores_sens(3,:)); 
        stats.thrustMean = mean(scores_sens(3,:));
        stats.thrustStd = std(scores_sens(3,:));
        stats.thrustVar = var(scores_sens(3,:));
%     %Plotting
%     plot(range,scores);
%     figure(1); hold on;
%     plot(nPass,scores(:,1)*sc_base);
%     plot(banLen,scores(:,2)*sc_base);
%     plot(speed,scores(:,3)*sc_base);
%     plot(GM_time,scores(:,4)*sc_base);
%     xlabel('Deviation of parameter from baseline');
%     ylabel('Deviation of score from baseline');
%     title('Sensitivity study of global competition parameters');
%     legend({'Number of passengers','Banner Length','Speed','GM Time'},...
%             'Location', 'eastoutside');
    
end

function sc = getScore(mass, power, thrust, span, S_sens, ban_length)
    global max_M2 max_M3 best_config top_score
    
    [SM2, vcM2, vtM2, lapsM2, flyM2, peeps, prod2,tt, profileM2, mo] = calculate_valuesM2(mass,span,power,thrust,ban_length, S_sens);
    [SM3, vcM3, vtM3, lapsM3, flyM3, prod3, profileM3] = calculate_valuesM3(mass,span,power,thrust,ban_length, mo, S_sens);

    if flyM2 == 0 || flyM3 == 0
        M2 = 0;
        M3 = 0;
    else
        M2 = 1 + (prod2)/max_M2; %(peeps/tt2)/max_M2;
        M3 = 2 + (prod3)/max_M3; %(lapsM3*xl)/max_M3;
        
    end
    sc = M2 + M3;
end