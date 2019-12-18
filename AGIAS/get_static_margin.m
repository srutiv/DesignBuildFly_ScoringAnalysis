plane.build_file
plane=plane.run_avl ;

k=1;
while plane.st.NP >.51 || plane.st.NP < .49 && k < 100
    
    if plane.st.NP > .51
        coord=plane.htails.coord;
        plane.htails.coord=coord+[coord(1)*.95,0,0];
    elseif plane.st.NP <.49
        coord=plane.htails.coord;
        plane.htails.coord=coord+[coord(1)*1.05,0,0];
        
    end
    plane.build_file
    plane=plane.run_avl;
    k=k+1;
end
