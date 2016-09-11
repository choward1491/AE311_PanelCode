function z = gauss(A,b)
%Gaussian Elimination  By Christian Howard
%Input any NxM matrix A and any Nx1 matrix b
%Used for solving linear equations

%Define number of rows,r, and number of columns,c, of matrix A
[r c]=size(A);

%Create matrix that will be used for Gaussian Elimination
mSolve = zeros(r,c+1);

%Fill up the created matrix with the initial values
for i=1:c
   mSolve(:,i)=A(:,i); 
end
   mSolve(:,c+1)=b(:);

%Work to get matrix to the simplest row reduced form
%----------------------


    for i=1:r
        
        %Check and implementation of row swapping if needed
        %--------------------------------
        if(mSolve(i,i) ~=0)
        mSolve(i,:)=mSolve(i,:)./mSolve(i,i);
        else
            for s=1+i:r
                if(s==i)
                    %----
                end    
                if(mSolve(s,i)~=0)
                    m2=mSolve;
                    mSolve(i,:)=m2(s,:);
                    mSolve(s,:)=m2(i,:);
                    break;
                end    
            end
        mSolve(i,:)=mSolve(i,:)./mSolve(i,i);
        end
        %--------------------------------
        
        
        for t=1:r
            
            if (t==i)
                
            else
              const = -mSolve(t,i);  
              mSolve(t,:)=mSolve(t,:)+(const.*mSolve(i,:));
            end
        

        end 
    end
%----------------------

    
%$  <Check to see if an exact solution has been obtained>  %$

%If there are free variables (when more columns than rows):
       %Will leave the row reduced matrix of the equations
%Otherwise, will give exact solution
%----------------------
    for i=2:c
        if(mSolve(i-1,i:c) == 0)
            if(mSolve(i-1,i-1)~=1)
                z=mSolve;
                break;
            end
                if(i==c)
                    z=mSolve(:,c+1);
                end
        else 
            z=mSolve;
            break;
        end
    end        
%-----------------------
end