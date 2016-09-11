function output = Ycamb(input)
output=zeros(1,length(input));
p=.4;
m=.01;
    for i=1:length(input)
        if input(i)<=p
            output(i)=(m/p^2)*(2*p*input(i)-input(i)^2);
        else
            output(i)=(m/(1-p)^2)*((1-2*p)+2*p*input(i)-input(i)^2);
        end
    end

end