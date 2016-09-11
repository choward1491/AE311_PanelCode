function VectorField2D(fun,xRange,yRange)
%By Christian Howard
%Created November 2011
%Input Variables:
%~fun: the function for the vector field
%~xRange: the range of x values you want the plot to show the vector field
%~yRange: same as xRange but for the y values
%~dim: the dimension, 2 or 3 are the options

xDiff = diff(xRange);
yDiff = diff(yRange);

x = linspace(xRange(1),xRange(2),30);
y = linspace(yRange(1),yRange(2),30);

slope = size(length(x),length(y));

%Find slope values at each point desired in the vector field
for i=1:length(x)
    for j=1:length(y)
        
        slope(i,j) = fun(i,j);
        
    end
end

%Find the length each slope line should be
t = linspace(0,1,10)*sqrt(xDiff.^2+yDiff.^2);

%Find the line paramatrizations for each slope and plot them
for i=1:length(x)
    for j=1:length(y)
        
    end
end





end