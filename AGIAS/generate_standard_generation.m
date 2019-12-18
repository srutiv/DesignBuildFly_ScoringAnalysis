num=10; % number of aircraft per generation 
generations=10; % number of generations 
generation_name='test'; % base name for aircraft geometry files that are generated 
score_matrix=score_matrix_generator();

%% create the first generation of aircraft using the generate_standard_aircraft function
for n=1:num
    plane_name=sprintf('%s_%i',generation_name, n);
    plane(n)=generate_standard_aircraft(plane_name);
end

%% initialize waitbar
clear time
steps=num*generations;
time_step=3;
time_left=steps*time_step;
wait_string=sprintf('%.0f s to complete',time_left);
h=waitbar(0,wait_string);
count=0;
time=[];
for p=1:generations
    
    %% run avl analysis on the current generation and score them with the score matrix 
    for n=1:num
        tic
        plane(n).build_file
        plane(n).build_mass
        plane(n)=plane(n).run_avl;
        score(n)=score_aircraft(score_matrix,plane(n));
        time_left=(steps-count)*time_step;
        wait_string=sprintf('%.0f s to complete',time_left);
        count=count+1;
        waitbar(count/steps,h,wait_string);
        time=[time toc];
        time_step=mean(time);
    end
    [B,I]=sort(score);
    %% create a new generation, cromprised of the 50th percentile, and an evolution of every aircraft in the 50th percentile 
    for n=1:num/2
        nplane(n*2-1)=plane(I(n));
        nplane(n*2)=evolve(plane(I(n)),1);
    end
    
    %% plot the 9 best aircraft in this generation 
    for n=1:9 
       figure(p)
       subplot(3,3,n)
       plane(I(n)).plot_all
    end
    %%  
        planeold=plane;
        plane = nplane;
end


