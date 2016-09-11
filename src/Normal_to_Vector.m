function NormVect = Normal_to_Vector(MnPt,Pts)
%Creator: Christian
%Created in October 2011
%Made to find the unit Normal vector to a Vector, with the orientation such
%that the normal would be going in the  +y direction for vector <1,0>

Vectors = VectBtwn(MnPt,Pts);
[r c] = size(Vectors);
nx = zeros(r,1);
ny = zeros(r,1);

    for i=1:r
        x=Vectors(i,1);
        y=Vectors(i,2);
       if(y>0) 
           x1 = -1;
           y1 = -((x1*x)/y);
           [nx(i),ny(i)]=UnitVect([x1,y1]);
       else if(y<0)
               x1 = 1;
               y1 = -((x1*x)/y);
               [nx(i),ny(i)]=UnitVect([x1,y1]);
           else 
               if(x>0)
                   nx(i)=0;
                   ny(i)=1;
               else if(x<0)
                   nx(i)=0;
                   ny(i)=-1;
                   end
               end
           end
       end



    end    

NormVect = [nx,ny];



end