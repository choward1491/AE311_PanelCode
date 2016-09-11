function [value] = SimpsInt(fun,xRange,n)
%Creator: Christian Howard
%Created in October 2011
%Find the integral value of any function between any x values
%Approximate this integral by n+1 terms

h=diff(xRange)/n;
f=zeros(n+1,1);
x1=xRange(1);
x2=xRange(2);


    for i=1:n+1
        
       if(i==n+1)
       f(i)=(1/3)*h*fun(x2);


       else if(i-1==0)
       f(i)=(1/3)*h*fun(x1);   


           else if(mod((i-1),2)~=0)
       f(i)=(1/3)*h*4*fun(x1+(i-1)*h);


               else if(mod((i-1),2)==0)
       f(i)=(1/3)*h*2*fun(x1+(i-1)*h);

                   end
               end
           end
       end    
    end

value = sum(f);


end
