function value = crossValue(fun,xRange)
%Author: Christian Howard
%Use recursion to find when the linear line crosses the y=0 axis
%fun : the function handle of the function you are using
%xRange: the range that houses the x value where the line crosses y=0

%Find the first minimum value of the abs(fun(x)) to see around where the
%line crosses the y=o axis
x=linspace(xRange(1),xRange(2),10000);
y=fun(x);
y=abs(y);
xmin=x(y==min(y));


xplus=x(x>xmin(1));
xminus=x(x<xmin(1));

%Check if the value of xmin gives an output for the function that is
%roughly equal to zero and if not, use recursion with the new, smaller x
%range until the value is roughly zero
    if(abs(fun(xmin(1)))<1e-8)
        value = xmin(1);
    else
        value = crossValue(fun,[xplus(1),xminus(length(xminus))]);
    end
    
end