function output = Xs(Ns)
output = zeros(1,Ns);
    for j=1:Ns
        thetaj=pi*((j-1)/(Ns-1));
        if(cos(thetaj)==1)
            output(j)=(1-.5*(1+cos(pi*((j+1-1)/(Ns-1)))))/2 + .5*(1+cos(pi*((j+1-1)/(Ns-1))));
        else if(cos(thetaj)==-1)
            output(j)=(0+.5*(1+cos(pi*((j+1-1)/(Ns-1)))))/2 ;
            else
            output(j)=.5*(1+cos(thetaj));
            end
        end
    end
    

end