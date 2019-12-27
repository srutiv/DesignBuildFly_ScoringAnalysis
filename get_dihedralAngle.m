dihedral_angle(wingSpan,wingArea,rootChord,tipChord, wingTaper, airfoil,g,d,mass,v,range)

function [maxgamma,maxfeature] = dihedral_angle(wingSpan,wingArea,rootChord,tipChord, wingTaper, airfoil,g,d,mass,v,range)
    % Iterating taper ratios to find best e %
    wingAirfoil = strcat(airfoil,'.dat'); % must include the file extension
    
    features = zeros(range-1,2); % To store all the feature we find

    % Looping over lambda which will have values from 0 to 1
    for gamma = -(range*pi/180):(pi/180):(range*pi/180)

        % Inputting the geometry and run case parameters
        % Main_Wing(B,SB,S,Taper,Sweep,Airfoil,Dihedral,Incidence,Surfaces,Coord,Name) %% class constructor
        wing = Main_Wing(wingSpan,wingArea,wingTaper,0,wingAirfoil,gamma,0,[],[0 0 0],'wing');
        Case = inputs(g,d,mass,0,0,0,0,v,0,false);

        % Building the aircraft, .avl file, and running it
        plane = Aircraft(wing,[],[],Case,[],[0,0,0],mass,[1,1,1],'tapir');
        plane.build_file;
        plane.run_avl;

        % Reading the header data from the sb file
        Header = parseRunCaseHeader(strcat('.\AVL_DATA\','tapir','_DATA.st'));
        feature = Header.afasdfasfasf; % what header values are important for us?

        % Printing where the script is, among other things
        fprintf('%f %f\n', gamma, feature);

        % Adding the result to the tracker
        features(fix(gamma*range)+1,:) = [gamma,feature];
    end

    % for some reason there are random (0,0) rows in es, so we delete them
    rowsToDelete = find(features(:,2)==0); % find them!
    features(rowsToDelete,:) = []; % clean them up

    % Find the max feature value and the corresponding gamma value
    maxfeature = max(features(:,2));
    maxgammamas = features((features(:,2)==maxfeature),1);
    maxgamma= maxgammamas(fix(length(maxgammamas)/2)+1); % there might be multiple dihedral angle (gamma) that give the same feature
end