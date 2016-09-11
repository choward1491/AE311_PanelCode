function output = Derivative(fun,X)
%Creator: Christian Howard
%Created in October 2011
%Input the function and x value you want the derivative at

h=1e-6;
output = (fun(X+h)-fun(X))/h;
end