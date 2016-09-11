function output = Xv(Nv)
xs = Xs(Nv+1);
output = zeros(1,Nv);
    for k=1:Nv

        output(k)=.5*(xs(k)+xs(k+1));

    end    

end