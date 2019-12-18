classdef Surface 
    properties
        chinge % %of chord
        sign % +1 for same direction or -1 opposite direction across symmetry...
        %(i.e. +1 for elevator, -1 for aileron)
        bi % fraction of span for start of surface (i.e. .25)
        be % fraction of span for end of surface  (i.e. .75)
        name
        %the example above would take up half of the span
    end
    methods
        function obj=Surface(Chinge,Sign,Bi,Be,Name)
            % constructor
            obj.chinge=Chinge
            obj.sign=Sign
            obj.bi=Bi
            obj.be=Be
            obj.name=Name
        end
    end
end