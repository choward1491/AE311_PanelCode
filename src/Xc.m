function output = Xc(N)
output = zeros(1,N);
    for i=1:N
        thetai=2*pi*((i-1)/(N-1));
        output(i)=.5*(1+cos(thetai));
    end
end