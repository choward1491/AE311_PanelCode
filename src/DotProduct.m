function [value anglBtwn ]= DotProduct(Vector1,Vector2)
%Creator: Christian
%Created in October 2011
%Made to find the dot product between two 2D vectors and angle between vectors
%Angle in radians

value = Vector1(1)*Vector2(1)+Vector1(2)*Vector2(2);
Mag1 = sqrt(Vector1(1)^2+Vector1(2)^2);
Mag2 = sqrt(Vector2(1)^2+Vector2(2)^2);

anglBtwn = acos(value / (Mag1*Mag2));


end