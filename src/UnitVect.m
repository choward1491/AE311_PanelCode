function [ux,uy] = UnitVect(Vector)
%
%Creator: Christian Howard
%Created in October 2011
%
%Find the unit vector based on the input Vector(s)

[r,c]=size(Vector);
unitVector = zeros(r,c);
q=linspace(0,0,r);
Magnitude = zeros(1,r);


        for j=1:r
            for i=1:c
                q(j)= q(j) +Vector(j,i)^2;
            end

        Magnitude(j) = sqrt(q(j));
        end

    for i=1:r
        unitVector(i,:) = (1/Magnitude(i))*Vector(i,:);
    end    
    ux=unitVector(:,1);
    uy=unitVector(:,2);
end