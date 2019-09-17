% This is a rewriting of TheOneScript.m to reduce clutter and make it
% easier to examine specific performances

run = 1; % whether to run the full script or not
clean = 1; % to clean up the score matrix

% iteration parameters
mit = 0.1; Pit = 10; Tit = 0.5;
airfoil = 'mh32';

% To run manually run lines 9,15,16,17 in the command window. Then run line 23 with proper inputs

tic;
if run
    ac = aircraftIter(1,1,1,3,airfoil,[1,1,1,1,1]);
    [Clmax, Cd0t, Cd0c] = ac.airfoilData;
    [d,A] = ac.airfoilMeasures;
    score = zeros(80,100,200,10);
    % Iterate over total mass and mechanical power and static thrust
    for m = 5:80
        parfor P = 10:100
            for T = 20:200
                ac = aircraftIter(m*mit,P*Pit,T*Tit,3,airfoil,[Clmax,Cd0t,Cd0c,d,A]);
                score(m,P,T,:) = [ac.S, ac.N, ac.x, ac.mn/ac.m, ac.mp, ac.mp/ac.m, ac.m0/ac.m, ac.vc, ac.vb, ac.vt];
                %fprintf('M=%.2fkg, P=%.2fW, T=%.2fN, Takeoff: %d\n', m*mit, P*Pit, T*Tit, ac.x>0)
            end
        end
    end
end
toc

tic;
if clean
    i=0;
    cleanScore = zeros((80-5)*(100-10)*(200-20),11);
    for m = 5:80
        for P = 10:100
            for T = 20:200
                i = i+1;
                if score(m,P,T,3) >= 1 && score(m,P,T,2) >= score(m,P,T,3) && score(m,P,T,1) <= 2% && score(m,P,T,5) >= T*Tit*(2*0.9/30)+0.4
                %if cleanScore(i,1) == 0
                    cleanScore(i,:) = [m*mit, P*Pit, T*Tit, score(m,P,T,1), score(m,P,T,2), score(m,P,T,3), score(m,P,T,4),...
                        score(m,P,T,5), score(m,P,T,8), score(m,P,T,9), score(m,P,T,10)]; % record the successful planes
                
                end
            end
        end
    end
    rowsToDelete = find(cleanScore(:,1)==0); % find them!
    cleanScore(rowsToDelete,:) = []; % send them to meet their maker
end
toc

% figure
% mplot = scatter(cleanScore(:,1),cleanScore(:,6),5,'filled');
% xlabel('mass (kg)')
% ylabel('Number of laps')
% figure
% Pplot = scatter(cleanScore(:,2),cleanScore(:,6),5,'filled');
% xlabel('Mechanical Power (W)')
% ylabel('Number of laps')
% figure
% Tplot = scatter(cleanScore(:,3),cleanScore(:,6),5,'filled');
% xlabel('Static Thrust (N)')
% ylabel('Number of laps')