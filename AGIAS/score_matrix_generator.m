function score_matrix=score_matrix_generator()
% generates the score matrix for the genetic algorithm. Containts target
% values for each aircraft parameter, and weights for each parameter. 

%% targets 
CL=1.1;
CD=.1;
CL_CD=200;
e=.95;
Clb=-.1; %% IMPORTANT roll moment relative to side slip
Clp=-.1; %% IMPORTANT roll moment relative to roll rate
Clr=0;
Cnb=.2; %% IMPORTANT yaw moment side slip
Cnp=0;
Cnr=-.1; %% IMPORTANT yaw moment yaw rate
Cma=-.5; %% IMPORTANT pitch moment angle of attack
Cmq=-8;
AR=7;
Lift=1;
roll_d=.7;
dutch_d=0;
short_d=1;
spiral_d=0;
phugoid_d=.5;



%% Gains
CL_gain=1;
CD_gain=0;
CL_CD_gain=0;
e_gain=0;
Clb_gain=1; %% IMPORTANT
Clp_gain=1; %% IMPORTANT
Clr_gain=0;
Cnb_gain=1; %% IMPORTANT
Cnp_gain=0;
Cnr_gain=1; %% IMPORTANT
Cma_gain=1; %% IMPORTANT 
Cmq_gain=0;
AR_gain=1;
Lift_Gain=0;
roll_dg=0;
dutch_dg=0;
short_dg=0;
spiral_dg=0;
phugoid_dg=0;

score_matrix(:,1)=[CL;CD;CL_CD;e;Clb;Clp;Clr;Cnb;Cnp;Cnr;Cma;Cmq;AR;Lift;roll_d;dutch_d;short_d;spiral_d;phugoid_d];
score_matrix(:,2)=[CL_gain;CD_gain;CL_CD_gain;e_gain;Clb_gain;Clp_gain;Clr_gain;Cnb_gain;Cnp_gain;Cnr_gain;Cma_gain;Cmq_gain;AR_gain;Lift_Gain;roll_dg;dutch_dg;short_dg;spiral_dg;phugoid_dg];
end