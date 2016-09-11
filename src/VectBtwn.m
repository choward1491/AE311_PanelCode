function Vect = VectBtwn(PtA,PtB)
%Creator: Christian
%Created in October 2011
%Made to find the Vector between two points
%PtB input can have a matrix of multiple points but not PtA

[r,c]=size(PtB);
dx = zeros(r,1);
dy = zeros(r,1);
    for i=1:r

    dx(i) = PtB(i,1)-PtA(1);
    dy(i) = PtB(i,2)-PtA(2);
    end

Vect = [dx,dy];

end