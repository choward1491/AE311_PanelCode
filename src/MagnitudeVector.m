function value = MagnitudeVector(Vector)
%Author: Christian Howard
%Made in October 2011
%Takes a vector of N dimensions and finds the magnitude of it
[r,c]=size(Vector);

    for i=1:r
    value(i) = sqrt(Vector(i,1)^2+Vector(i,2)^2);
    end

end

