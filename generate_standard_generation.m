num=20;
generations=8;
generation_name='test';


%%
for n=1:num
    plane_name=sprintf('%s_%i',generation_name, n);
    plane(n)=generate_standard_aircraft(plane_name);
end

%%
% steps=num*generations;
% time_left=steps*3.0875;
% wait_string=sprintf('%f s to complete',time_left);
% h=waitbar(0,wait_string);
% 
% time_left=steps*3.0875;
% count=0;

for p=1:generations
    
    %%
    for n=1:num
        plane(n).build_file
        plane(n).build_mass

    %%
        plane(n)=plane(n).run_avl;
%         count=count+1;
%         time_left=time_left-3.0875;
%         close(h)
%         wait_string=sprintf('%f s to complete',time_left);

%      %%
%     for n=1:num
%         figure(1)
%         subplot(3,num/3,n)
%         title(plane(n).name)
%         plane(n).plot_all
%     end
    %%
        score(n)=score_aircraft(score_matrix,plane(n));
    end
    %%
    [B,I]=sort(score);
%     %%
%     
%     nplane(1)=plane(I(1));
%     nplane(2)=evolve(plane(I(1)),1);
%     nplane(3)=evolve(plane(I(1)),2);
%     nplane(4)=plane(I(2));
%     nplane(5)=evolve(plane(I(2)),1);
%     nplane(6)=evolve(plane(I(2)),2);
%     nplane(7)=plane(I(3));
%     nplane(8)=evolve(plane(I(3)),1);
%     nplane(9)=evolve(plane(I(3)),2);
%     plane=nplane;
%     
    %%
    for n=1:num/2
        nplane(n*2-1)=plane(I(n));
        nplane(n*2)=evolve(plane(I(n)),1);
    end
%     for n=1:3
%       nplane(n)=mate(plane(I(1)),plane(I(2)),n)  
%     end
%     for n=4:6
%     nplane(n)=mate(plane(I(1)),plane(I(3)),n)
%     end
    
    %%
    for n=1:9 
       figure(p)
       subplot(3,3,n)
       plane(I(n)).plot_all
    end
    %%
        planeold=plane;
        plane = nplane;
end