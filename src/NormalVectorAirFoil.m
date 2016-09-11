function n = NormalVectorAirFoil(fun,x)
%
%Creator: Christian Howard
%Created in October 2011
%
%Find the normal of any 2D function
%Input the function and the x value(s) where a normal is desired to be
%found
%
%Choose the normal pointing out by inputing 'out' for input 3
%Choose the normal pointing in by inputing 'in' for input 3
%
%For functions with negative values, the normals will have y components
%pointing downward and vise versa


df = Derivative(fun,x);
r=length(x);
        
            %Outward Normal
                  for i=1:r
                    f = fun(x(i));
                    a= df(i);
                    %Make unit vector of slope of f at that x value
                        [tx,ty]=UnitVect([1,a]);
                    
                    %If the function has positive values at the X input
                    if(f>0)
%-----------------------                        
                        %Find normal, based on value of the slope
                        if(a==Inf)
                            if(df(i+5)~=Inf)
                                if(df(i+5)>0)
                                    x1=-1;
                                    y1=0;
                                else
                                    x1=1;
                                    y1=0;
                                end
                            else if(Derivative(fun,x(i)-1e-3)~=Inf)
                                    if(Derivative(fun,x(i)-1e-3)>0)
                                    x1=-1;
                                    y1=0;
                                    else
                                    x1=1;
                                    y1=0;
                                    end
                                end
                            end
                        else if(a>0)
                            
                            x1=-1;
                            y1=-((x1)*tx)/(ty);
                            
                            else if(a<0)
                                x1=1;
                                y1=-((x1)*tx)/(ty);
                                else if(a==0)
                                    x1=0;
                                    y1=1;
                                    end
                                end
                            end
                        end
%-------------------------
                        vect(1,:)=[x1,y1];
                        [nx(i),ny(i)]=UnitVect(vect);
%--------------------------------------------------------

                    %If the function has negative values at the X input
                    else if(f<0)
%-------------------------                           
                        %Find normal, based on value of the slope
                        if(a==Inf)
                            if(df(i+5)~=Inf)
                                if(df(i+5)>0)
                                    x1=1;
                                    y1=0;
                                else
                                    x1=-1;
                                    y1=0;
                                end
                            else if(Derivative(fun,x(i)-1e-3)~=Inf)
                                    if(Derivative(fun,x(i)-1e-3)>0)
                                    x1=1;
                                    y1=0;
                                    else
                                    x1=-1;
                                    y1=0;
                                    end
                                end
                            end
                        else if(a>0)
                            
                            x1=1;
                            y1=-((x1)*tx)/(ty);
                            
                            else if(a<0)
                                x1=-1;
                                y1=-((x1)*tx)/(ty);
                                else if(a==0)
                                    x1=0;
                                    y1=-1;
                                    end
                                end
                            end
                        end
%-------------------------
                        vect(1,:)=[x1,y1];
                        [nx(i),ny(i)]=UnitVect(vect);
%-------------------------------------------------
                        %If the value of the function is 0 
                        else if(f==0)
                                if(fun(x(i)+1e-3)>0 && fun(x(i)-1e-3)<0)
                                    %Find normal, based on value of the slope
                                            if(a==Inf)
                                                if(df(i+5)~=Inf)
                                                    if(df(i+5)>0)
                                                        x1=1;
                                                        y1=0;
                                                    else
                                                        x1=-1;
                                                        y1=0;
                                                    end
                                                else if(Derivative(fun,x(i)-1e-3)~=Inf)
                                                        if(Derivative(fun,x(i)-1e-3)>0)
                                                        x1=1;
                                                        y1=0;
                                                        else
                                                        x1=-1;
                                                        y1=0;
                                                        end
                                                    end
                                                end
                                            else if(a>0)

                                                x1=1;
                                                y1=-((x1)*tx)/(ty);

                                                else if(a<0)
                                                    x1=-1;
                                                    y1=-((x1)*tx)/(ty);
                                                    else if(a==0)
                                                        x1=0;
                                                        y1=-1;
                                                        end
                                                    end
                                                end
                                            end
                    %-------------------------
                                            vect(1,:)=[x1,y1];
                                            [nx(i),ny(i)]=UnitVect(vect);
                    %-------------------------------------------------
                                else if(fun(x(i)-1e-3)>0 && fun(x(i)+1e-3)<0)
                                        if(a==Inf)
                                                if(df(i+5)~=Inf)
                                                    if(df(i+5)>0)
                                                        x1=-1;
                                                        y1=0;
                                                    else
                                                        x1=1;
                                                        y1=0;
                                                    end
                                                else if(Derivative(fun,x(i)-1e-3)~=Inf)
                                                        if(Derivative(fun,x(i)-1e-3)>0)
                                                        x1=-1;
                                                        y1=0;
                                                        else
                                                        x1=1;
                                                        y1=0;
                                                        end
                                                    end
                                                end
                                            else if(a>0)

                                                x1=-1;
                                                y1=-((x1)*tx)/(ty);

                                                else if(a<0)
                                                    x1=1;
                                                    y1=-((x1)*tx)/(ty);
                                                    else if(a==0)
                                                        x1=0;
                                                        y1=1;
                                                        end
                                                    end
                                                end
                                            end
                    %-------------------------
                                            vect(1,:)=[x1,y1];
                                            [nx(i),ny(i)]=UnitVect(vect);
                                    end
                                end
                                %-------
                                
                                            
                    %--------------------------------------------------------
                                else
                                  
                                end
                            end
                        
                        end
                    end

        
        
        
        nx=nx';
        ny=ny';
        n=[nx,ny];
       for k=1:length(x) 
            if(fun(x(k))==0)
            n(fun(x)==0,:)=[-1,0];
            break;
            end
       end
end