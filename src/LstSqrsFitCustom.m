function [outFunc t r2] = LstSqrsFitCustom(x,y,choice)
%
%Creator: Christian Howard
%Created in October 2011
%
%Find an approximate function to fit the given x and y values
%Input many x and y values to get a better fit
%choice option includes:
%all: Fit a curve with polynomials from x^.5 to x^6
%line: find line that fits function of f(x) = ax + b
%poly2: Find line that fits function of f(x) = ax^2 +bx +c
%poly3: Find line that fits function of f(x) = ax^3 +bx^2 +cx +d
%poly4: Find line that fits function of f(x) = ax^4 +bx^3 +cx^2 +dx +e
%poly5: Find line that fits function of f(x) = ax^5 +bx^4 +cx^3 +dx^2 +ex + f
%poly6: Find line that fits function of f(x) = ax^6 +bx^5 +cx^4 +dx^3 +ex^2 + fx + g
%sine: Find line that fits function of f(x) = Asin(x) + b

f0 = @(x) x.^(1/3);
f1 = @(x) sqrt(x);
f2 = @(x) x;
f3 = @(x) x.^2;
f4 = @(x) x.^3;
f5 = @(x) x.^4;
f6 = @(x) x.^5;
f7 = @(x) x.^6;
f8 = @(x) 1;
f9 = @(x) sin(x);
 


    if(strcmp(choice,'all'))
        A=zeros(length(x),9);
        [r c]=size(y);

                if(r==1)
                    y=y';
                else
                end
            for i=1:length(x)
            A(i,1)=f0(x(i));
            A(i,2)=f1(x(i));
            A(i,3)=f2(x(i));
            A(i,4)=f3(x(i));
            A(i,5)=f4(x(i));
            A(i,6)=f5(x(i));
            A(i,7)=f6(x(i));
            A(i,8)=f7(x(i));
            A(i,9)=f8(x(i));
            end
            At = A';
    Af = At*A;
    yf = At*y;

    t = gauss(Af,yf);

    outFunc = @(x) t(1)*f0(x)+t(2)*f1(x)+t(3)*f2(x)+t(4)*f3(x)+t(5)*f4(x)+t(6)*f5(x)+t(7)*f6(x)+t(8)*f7(x);
    
    end

    if(strcmp(choice,'line'))
        A=zeros(length(x),2);
        [r c]=size(y);

                if(r==1)
                    y=y';
                else
                end
            for i=1:length(x)
            A(i,1)=f2(x(i));
            A(i,2)=f8(x(i));

            end
            At = A';
    Af = At*A;
    yf = At*y;

    t = gauss(Af,yf);

    outFunc = @(x) t(1)*f2(x)+t(2)*f8(x);
    end

    if(strcmp(choice,'poly2'))
        A=zeros(length(x),3);
        [r c]=size(y);

                if(r==1)
                    y=y';
                else
                end
            for i=1:length(x)
            A(i,1)=f2(x(i));
            A(i,2)=f3(x(i));
            A(i,3)=f8(x(i));
            end

            At = A';
    Af = At*A;
    yf = At*y;

    t = gauss(Af,yf);

    outFunc = @(x) t(1)*f2(x)+t(2)*f3(x)+t(3)*f8(x);
    end

    if(strcmp(choice,'poly3'))
        A=zeros(length(x),4);
        [r c]=size(y);

                if(r==1)
                    y=y';
                else
                end
            for i=1:length(x)

            A(i,1)=f2(x(i));
            A(i,2)=f3(x(i));
            A(i,3)=f4(x(i));
            A(i,4)=f8(x(i));
            end
            At = A';
    Af = At*A;
    yf = At*y;

    t = gauss(Af,yf);

    outFunc = @(x) t(1)*f2(x)+t(2)*f3(x)+t(3)*f4(x)+t(4)*f8(x);
    end

    if(strcmp(choice,'poly4'))
        A=zeros(length(x),5);
        [r c]=size(y);

                if(r==1)
                    y=y';
                else
                end
            for i=1:length(x)
            A(i,1)=f2(x(i));
            A(i,2)=f3(x(i));
            A(i,3)=f4(x(i));
            A(i,4)=f5(x(i));
            A(i,5)=f8(x(i));
            end
            At = A';
    Af = At*A;
    yf = At*y;

    t = gauss(Af,yf);

    outFunc = @(x) t(1)*f2(x)+t(2)*f3(x)+t(3)*f4(x)+t(4)*f5(x)+t(5)*f8(x);
    end

    if(strcmp(choice,'poly5'))
        A=zeros(length(x),6);
        [r c]=size(y);

                if(r==1)
                    y=y';
                else
                end
            for i=1:length(x)
            A(i,1)=f2(x(i));
            A(i,2)=f3(x(i));
            A(i,3)=f4(x(i));
            A(i,4)=f5(x(i));
            A(i,5)=f6(x(i));
            A(i,6)=f8(x(i));
            end
            At = A';
    Af = At*A;
    yf = At*y;

    t = gauss(Af,yf);

    outFunc = @(x) t(1)*f2(x)+t(2)*f3(x)+t(3)*f4(x)+t(4)*f5(x)+t(5)*f6(x)+t(6)*f8(x);
    end
    if(strcmp(choice,'poly6'))
        A=zeros(length(x),7);
        [r c]=size(y);

                if(r==1)
                    y=y';
                else
                end
            for i=1:length(x)
            A(i,1)=f2(x(i));
            A(i,2)=f3(x(i));
            A(i,3)=f4(x(i));
            A(i,4)=f5(x(i));
            A(i,5)=f6(x(i));
            A(i,6)=f7(x(i));
            A(i,7)=f8(x(i));
            end
            At = A';
    Af = At*A;
    yf = At*y;

    t = gauss(Af,yf);

    outFunc = @(x) t(1)*f2(x)+t(2)*f3(x)+t(3)*f4(x)+t(4)*f5(x)+t(5)*f6(x)+t(6)*f7(x)+t(7)*f8(x);
    end

    if(strcmp(choice,'poly.5'))
        A=zeros(length(x),2);
        [r c]=size(y);

                if(r==1)
                    y=y';
                else
                end
            for i=1:length(x)
            A(i,1)=f1(x(i));
            A(i,2)=f8(x(i));
            end
            At = A';
    Af = At*A;
    yf = At*y;

    t = gauss(Af,yf);

    outFunc = @(x) t(1)*f1(x)+t(2)*f8(x);
    end
    
    if(strcmp(choice,'sine'))
        A=zeros(length(x),2);
        [r c]=size(y);

                if(r==1)
                    y=y';
                else
                end
            for i=1:length(x)
            A(i,1)=f9(x(i));
            A(i,2)=f8(x(i));
            end
            At = A';
    Af = At*A;
    yf = At*y;

    t = gauss(Af,yf);

    outFunc = @(x) t(1)*f9(x)+t(2)*f8(x);
    end
    
    if(strcmp(choice,'sine line'))
        A=zeros(length(x),2);
        [r c]=size(y);

                if(r==1)
                    y=y';
                else
                end
            for i=1:length(x)
            A(i,1)=f9(x(i))*f2(x(i));
            A(i,2)=f8(x(i));
            end
            At = A';
    Af = At*A;
    yf = At*y;

    t = gauss(Af,yf);

    outFunc = @(x) t(1)*f9(x)*f2(x)+t(2)*f8(x);
    end
    
    %Find how good of a fit the line is with the data. 
    %1 means perfect fit,0 means no fit at all thus 0<=r2<=1
    ybar = mean(y);
    [r1 c1]=size(y);
    [r2 c2]=size(outFunc(x));
    
    SSM=sum((y-ybar).^2);
    if(r1~=r2)
    residual = y'-outFunc(x);
    residual2 = (y'-outFunc(x)).^2;
    else
    residual = y-outFunc(x);
    residual2 = (y-outFunc(x)).^2;
    end
    SSE=sum(residual2);
    
    r2 = (SSM-SSE)/SSM;
end