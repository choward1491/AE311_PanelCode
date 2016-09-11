function [Data alphaZL r2 ] = AirfoilCalc(Title,Ycamb,Ythick,Xc,Xs,Xv,Ns,Nv,U,angleOfAttack,angleType,showPlot)
%Creator: Christian Howard
%Created in October 2011
%Purpose: 
%~To find and plot the coefficient of lift based on the angle of attack
%~To find the coefficient of pressure distribution based on the angle of
% attack and plot it
%~To find an approximate for the zero lift angle of attack
%~To find and plot the sources' and vortices' strengths verse the x
%  position along the airfoil's chord.
%~To plot the airfoil design
%--------------------------------------------------------------------------
%AirfoilCalc(Title,Ycamb,Ythick,Xc,Xs,Xv,Ns,Nv,U,angleOfAttack,angleType,showPlot)
%
%Inputs and descriptions:
%~Title: The title you want the graph of the wing to say, most likely the
%wings name
%~Ycamb: The function for the camber line of the wing
%~Ythick: The function for the thickness of the wing
%~Xc: The function for the x position of the control points
%~Xs: The function for the x position of the sources
%~Xv: The function for the x position of the vortices
%~Ns: The number of sources
%~Nv: The number of vorticies
%~U: The freestream velocity
%~angleOfAttack: The angle at which airflow comes in  relative to the wing
%  can be in the form of a 1xN or Nx1 matrix
%~angleType: String of either 'radians' or 'degrees' to know what the angle
%  of attack angles are measured in
%~showPlot: String of either 'yes' or 'no', those answers deciding whether
%  or not you wish to show the any of the plots, aside from the plot for
%  the coefficient of lift which will show up no matter what
%
%Outputs and Descriptions:
%~Data: It is a Mx2 matrix with column 1 containing the angle of attacks in
%  degrees and then column 2 contains the coefficient of lift values that are
%  related to their angle of attack values in their row
%~alphaZL: It is the approximated zero lift angle of attack 
%~r2: This is the measure of how good the fit if the fitted linear line for
%  Cl vs angle of attack is. This gives an idea how accurate the alphaZL
%  value might be as well. r2 is between 1 and 0, and the closer to 1 it is,
%  the better the fit it is.
%
%Example: [Data Alphazl r2]= AirfoilCalc('NACA 1408',@Ycamb,@Ythick,@Xc,@Xs,@Xv,26,25,15,-10:3:10,'degrees','yes');

%Close any open figures
close ALL

%Check if the angle of attack is in degrees or radians
if(strcmp(angleType,'degrees') || strcmp(angleType,'Degrees')|| strcmp(angleType,'degree')|| strcmp(angleType,'Degree'))
    for l=1:length(angleOfAttack)
        angle(l) = angleOfAttack(l)*pi/180;
    end
else
    for l=1:length(angleOfAttack)
        angle(l) = angleOfAttack(l);
    end
end



%Total number of building block flows used
N=Ns+Nv;

%Make equations of the top and bottom pieces to the airfoil
Yup = @(x) Ycamb(x) + Ythick(x);
Ydown = @(x) Ycamb(x) - Ythick(x);

%The x/c position of the leading edge and the x/c position of the trailing
%edge
Start = 0;
End = 1;

%Matrix of numbers between the start and end to draw the foil
XCrange =linspace(Start,End,500);

%The y/c values of the top of the airfoil and the bottom of the airfoil
FoilUp = Yup(XCrange);
FoilDown = Ydown(XCrange);


%Check if the plots are desired to be shown
if(strcmp(showPlot,'Yes') || strcmp(showPlot,'yes'))
    %Plot the y/c values for the two pieces of the wing vs x/c
    %Label the important pieces
    figure
    plot(XCrange,FoilUp,'b-',XCrange,FoilDown,'b-');
    axis([Start-.2,End+.2,min(FoilDown)-.5,max(FoilUp)+.5])
    title(Title);
    xlabel('x/c');
    ylabel('y/c');
    hold on
else
end
%We need to plot where the sources, vortices and the control points are

%Use the functions of the x/c position of the Controls, Sources and
%Vortices to find the exact x/c positions of each of these
xc = Xc(N);
for t=1:length(xc)
    if(xc(t)==0)
        xct(t)=xc(t);
        xcb=xc(t+1:end);
        break;
    else
        xct(t)=xc(t);
    end
end
xs = Xs(Ns);
xv = Xv(Nv);

%Find the y/c position of the Control points using the functions for the
%top of the airfoil and the bottom of the airfoil
Cup = Yup(xct);

Cdown = Ydown(xcb);

%Find the y/c position of the vortices and sources using the function for
%the camber
Ys = Ycamb(xs);
Yv = Ycamb(xv);

%Check if the plots are desired to be shown
if(strcmp(showPlot,'Yes') || strcmp(showPlot,'yes'))
    
    %Plot where all of the Control points, vortices and sources are
    %Use red circles for control points
    %Use green triangles pointing up for the sources
    %Use blue X's for the vortices
    plot(xct,Cup,'ro',xcb,Cdown,'ro',xs,Ys,'g^',xv,Yv,'bX')
else
end

%Gather up the positions of the Controls, Vortices and Sources
Vortices = [xv',Yv'];
Sources = [xs',Ys'];
ControlUp = [xct',Cup'];
ControlDwn = [xcb',Cdown'];

%Find the normals along the airfoil surface at the control points
nCup = NormalVectorAirFoil(Yup,xct);
nDwn = NormalVectorAirFoil(Ydown,xcb);
A=size([N,N]);

%Start making the A matrix to solve for the magnitudes of the sources and
%the circulation
for r = 1:N
    c1=(1/(2*pi));
    %If the rows go past the number of control points on top, we want
    %to have the loops look at the control points on bottom
    if(r<length(xct)+1)
        for c=1:N
            %If the columns go past the number of Sources, we want to
            %have the loop looking at setting up the circulation part
            %of the matrix A
            if(c<Ns+1)
                R = VectBtwn(Sources(c,:),ControlUp(r,:));
                MgR=MagnitudeVector(R);
                [rx,ry]=UnitVect(R);
                rHat=[rx,ry];
                n=nCup(r,:);
                vTop=DotProduct(rHat,n);
                A(r,c)=c1*vTop/MgR;
            else
                R = VectBtwn(Vortices(c-Ns,:),ControlUp(r,:));
                MgR=MagnitudeVector(R);
                tHat=Normal_to_Vector(Vortices(c-Ns,:),ControlUp(r,:));
                n=nCup(r,:);
                vTop=DotProduct(tHat,n);
                A(r,c)=c1*vTop/MgR;
            end
        end
    else
        for c=1:N
            if(c<Ns+1)
                R = VectBtwn(Sources(c,:),ControlDwn(r-length(xct),:));
                MgR=MagnitudeVector(R);
                [rx,ry]=UnitVect(R);
                rHat=[rx,ry];
                n=nDwn(r-length(xct),:);
                vTop=DotProduct(rHat,n);
                A(r,c)=c1*vTop/MgR;
            else
                R = VectBtwn(Vortices(c-Ns,:),ControlDwn(r-length(xct),:));
                MgR=MagnitudeVector(R);
                tHat=Normal_to_Vector(Vortices(c-Ns,:),ControlDwn(r-length(xct),:));
                n=nDwn(r-length(xct),:);
                vTop=DotProduct(tHat,n);
                A(r,c)=c1*vTop/MgR;
            end
        end
        
    end
end


b=zeros(N,1);
%Make the direction vectors for each angle of attack
for i=1:length(angleOfAttack)
    alphaHat(i,:) = [cos(angle(i)),sin(angle(i))];
end

%Start putting the b matrix/matrices together
for i=1:length(angleOfAttack)
    for r=1:N
        if(r<length(xct)+1)
            b(r,1)=U*DotProduct(-alphaHat(i,:),nCup(r,:));
        else
            b(r,1)=U*DotProduct(-alphaHat(i,:),nDwn(r-length(xct),:));
        end
    end
    %We want the solutions for the magnitudes for each inputed angle of
    %attack
    X(:,i)=A\b;
    
end

%Make the plots showing the sources' and vortices' strengths vs x position
for i=1:length(angleOfAttack)
    if(strcmp(showPlot,'Yes') || strcmp(showPlot,'yes'))
    figure
    subplot(2,1,1)
    plot(xs,X(1:26,i),'b*')
    hold on
    alpha = angle(i)*180/pi;
    a = num2str(alpha);
    b='\Lambda vs Position along airfoil at \alpha=';
    b2 =' degrees';
    c = [b,a,b2];
    xlabel('Position (x/c)')
    ylabel('\Lambda Strength')
    title(c)
    
    subplot(2,1,2)
    plot(xv,X(27:end,i),'b*')
    hold on
    alpha = angle(i)*180/pi;
    a = num2str(alpha);
    b='\Gamma vs Position along airfoil at \alpha=';
    b2 =' degrees';
    c = [b,a,b2];
    xlabel('Position (x/c)')
    ylabel('\Gamma Strength')
    title(c)
    else
    end
end
%Velocity magnitude at each Control Point for each angle of attack
for i=1:length(angleOfAttack)
    for r = 1:N
        c1=(1/(2*pi));
        %If the rows go past the number of control points on top, we want
        %to have the loops look at the control points on bottom
        if(r<length(xct)+1)
            for c=1:N
                %If the columns go past the number of Sources, we want to
                %have the loop looking at setting up the circulation part
                %of the matrix A
                if(c<Ns+1)
                    R = VectBtwn(Sources(c,:),ControlUp(r,:));
                    MgR=MagnitudeVector(R);
                    [rx,ry]=UnitVect(R);
                    rHat=[rx,ry];
                    Mag = X(c,i);
                    %Break up the velocities into the x and y
                    %components
                    Usx(r,c)= (c1*Mag/MgR)*rHat(1);
                    Usy(r,c)= (c1*Mag/MgR)*rHat(2);
                else
                    R = VectBtwn(Vortices(c-Ns,:),ControlUp(r,:));
                    MgR=MagnitudeVector(R);
                    tHat=Normal_to_Vector(Vortices(c-Ns,:),ControlUp(r,:));
                    Mag = X(c,i);
                    %Break up the velocities into the x and y
                    %components
                    Uvx(r,c)= (c1*Mag/MgR)*tHat(1);
                    Uvy(r,c)= (c1*Mag/MgR)*tHat(2);
                    
                end
            end
        else
            for c=1:N
                if(c<Ns+1)
                    R = VectBtwn(Sources(c,:),ControlDwn(r-length(xct),:));
                    MgR=MagnitudeVector(R);
                    [rx,ry]=UnitVect(R);
                    rHat=[rx,ry];
                    Mag = X(c,i);
                    %Break up the velocities into the x and y
                    %components
                    Usx(r,c)= (c1*Mag/MgR)*rHat(1);
                    Usy(r,c)= (c1*Mag/MgR)*rHat(2);
                else
                    R = VectBtwn(Vortices(c-Ns,:),ControlDwn(r-length(xct),:));
                    MgR=MagnitudeVector(R);
                    tHat=Normal_to_Vector(Vortices(c-Ns,:),ControlDwn(r-length(xct),:));
                    Mag = X(c,i);
                    %Break up the velocities into the x and y
                    %components
                    Uvx(r,c)= (c1*Mag/MgR)*tHat(1);
                    Uvy(r,c)= (c1*Mag/MgR)*tHat(2);
                    
                end
            end
            
        end
    end
    Uvx=Uvx(:,27:end);
    Uvy=Uvy(:,27:end);
    
    %Preallocate space for the final summed velocities caused by the
    %sources and virtices
    Usf = zeros([N,2]);
    Uvf = zeros([N,2]);
    
    %Find the sum of all the velocities caused by each source at each
    % control point
    for r = 1:N
        Usf(r,:) = [0,0];
        %r is the looping through all the control points
        for c = 1:Ns
            %We will be looking at the columns of Us, which has length Ns
            Usf(r,:) = Usf(r,:)+[Usx(r,c),Usy(r,c)];
        end
    end
    
    %Find the sum of all the velocities caused by each source at each
    % control point
    for r = 1:N
        Uvf(r,:) = [0,0];
        %r is the looping through all the control points
        for c = 1:Nv
            %We will be looking at the columns of Uv, which has length
            %Nv
            Uvf(r,:) = Uvf(r,:)+[Uvx(r,c),Uvy(r,c)];
        end
    end
    
    %Preallocate matrix for final velocity vectors for each control
    % point
    Uf = zeros([N,2]);
    
    %Find the freestream velocity based on the angle-of-attack's
    %direction vector
    Uinf = U*alphaHat(i,:);
    
    %Various things are found in the next few steps:
    %1. Total velocity vector at each control point is found
    %2. Magnitude of the velocity at each control point is found
    %3. The function for the Pressure Coefficient on top found
    %4. The function for the Pressure Coefficient on bottom found
    %5. The Coefficient of Lift value for the ith alpha is found
    
    %For finding the Coefficient of Pressure, I used:
    %~Cp = 1-(V/Uinf)^2 where V is the magnitude of the velocity at any point on the
    %   surface of the wing and Uinf is the magnitude of the velocity of the freestream
    
    %For finding the Coefficient of lift, I used:
    %~Cl = Integral(Cp,bottom(x)-Cp,top(x)) from the leading edge to the
    %   trailing edge
    
    for r = 1:N
        Uf(r,:)=Usf(r,:)+Uvf(r,:)+Uinf;
        MagU(r) = MagnitudeVector(Uf(r,:));
        Cp(r)  = 1 - (MagU(r)/U)^2;
    end
    
    %Find how many rows the control points on top take up to use to
    %find where the Cp matrix starts looking at the Cp coefficients on
    %the bottom matrix
    [Lt,Wt] = size(ControlUp);
    
    %Define the coefficient of pressure on the top and bottom
    Cpt = Cp(1:Lt);
    Cpb = Cp(Lt+1:end);
    
    %Useing the least squares method, fit the data we have to find
    %functions for Cp on top and Cp on bottom
    CptFun = LstSqrsFitCustom(xct,Cpt,'poly6');
    CpbFun = LstSqrsFitCustom(xcb,Cpb,'poly6');
    
    %Next, we will plot the Cp vs x position on the top surface on one graph with the
    %fitted line and then do the same thing with Cp vs x position on
    %the lower surface and put both these plots on the same figure
    if(strcmp(showPlot,'Yes') || strcmp(showPlot,'yes'))
    figure
    subplot(2,1,1)
    plot(xct,Cpt,'b*')
    hold on
    X1=linspace(min(xct),max(xct),100);
    Y1 = CptFun(X1);
    plot(X1,Y1,'r--')
    alpha = angle(i)*180/pi;
    a = num2str(alpha);
    b='Cp vs Position along airfoil at \alpha=';
    b2 =' degrees';
    c = [b,a,b2];
    xlabel('Position (x/c)')
    ylabel('Coefficient of Pressure on Top')
    title(c)
    
    subplot(2,1,2)
    plot(xcb,Cpb,'b*')
    hold on
    X2=linspace(min(xcb),max(xcb),100);
    Y2 = CpbFun(X2);
    plot(X2,Y2,'r--')
    alpha = angle(i)*180/pi;
    a = num2str(alpha);
    b='Cp vs Position along airfoil at \alpha=';
    b2 =' degrees';
    c = [b,a,b2];
    xlabel('Position (x/c)')
    ylabel('Coefficient of Pressure on Bottom')
    title(c)
    else
    end
    %Make a new function representing Cp,bottom(x)-Cp,top(x)
    CpNet = @(x) CpbFun(x)-CptFun(x);
    
    %Integrate Cp,bottom(x)-Cp,top(x) from the leading edge to the
    %   trailing edge, using a quadratic integration approximation with
    %   5000 terms
    Cl(i) = SimpsInt(CpNet,[0,1],5000);
    
end

%Convert the angle back to degrees for use on the plot
for l=1:length(angleOfAttack)
    angle(l) = angle(l)*180/pi;
end

%Create Another figure showing the coefficient of lift plotted vs
%angle-of-attack in degrees
figure
plot(angle,Cl,'b-')
xlabel('Angle of Attack (Degrees)')
ylabel('Coefficient of Lift')
title('Angle of Attack vs Lift Coefficient')

    %Find the angle of attack for zero lift, AlphaZL
    %Fit a line to the data with the highest possible r^2 value
    %Keep reducing number of data points till best linear fit is made
    if(min(angle)<0 && max(angle)>0)
     for m=0:length(angle)
        [Clfun T r2]= LstSqrsFitCustom(angle(1+m:end-m),Cl(1+m:end-m),'line');
            if(r2>=.99999)
                break
            else
            end
     end
    elseif(min(angle)>=0)
        for m=0:length(angle)    
        [Clfun T r2]= LstSqrsFitCustom(angle(1:end-m),Cl(1:end-m),'line');
            if(r2>=.99999)
                break
            else
            end
        end
    elseif(max(angle)<=0)
        for m=0:length(angle)    
            [Clfun T r2]= LstSqrsFitCustom(angle(1+m:end),Cl(1+m:end),'line');
            if(r2>=.99999)
                break
            else
            end
        end
    end




%Using found best fit line, find the x value where 
%this function crosses Cl=0
alphaZL=crossValue(@(x) Clfun(x),[-90,90]);

%Return the coefficient of lift with its equivilant angle of attack
[r1,c1]=size(angleOfAttack);
    if(r1==1)
    Data = [angleOfAttack',Cl'];
    else
    Data = [angleOfAttack,Cl'];
    end
    
    

end