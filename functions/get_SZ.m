function [SS,ZZ]=get_SZ(X)
    if nargin~=1
        error('Pogresan broj argumenata');
    end
    
    % put all restrictions into one matrix (impacts)         
Y=[];
restrictions=keys(X);
for i=1:size(X,1)
    s=restrictions{i};
    Y=[Y X(s)];
end
Y=Y';
    
    %when we want to impose sign and zero restrictions at several horizons 
    %we stack the IRFs into single matrix : impact pattern matrix (Y)
    %each coloumn of impact pattern matrix represents restrictions for
    %given structural shock (for all horizons that restriczions are
    %imposed)
    %number of rows and columns of impact pattern matrix
    number_of_rows=size(Y,1);
    number_of_columns=size(Y,2);
    
    %SIGN RESTRICTIONS:
    %representation of sign restrictions: cell array SS
    %each element (S) of cell array SS is selection matrix for given structural
    %shock, if there are no sign restrictions, S matrix doesn't exist, else
    %dimension of S matrix is c x(n*r) where c is number of sign
    %restrictions for given structural shock, n is number of variables and
    %r is number of types of restrictions ( for example restrictions to
    %different horizons of Irfs, structural parameters, etc)
    SS=cell(number_of_columns,1);
    for j=1:number_of_columns % for each structural shock
        %number of sign restrictions for j-th structural shock
        check_sign=(Y(:,j)==1 | Y(:,j)==-1);
        c=sum(check_sign); 
        if c~=0 %if there is at least one sign restriction for given shock
            %initialize selection matrix to zero matrix
            S=zeros(c,number_of_rows);
            d=1;
            for i=1: number_of_rows
                if (Y(i,j)==1) | (Y(i,j)==-1)
                    S(d,i)=Y(i,j); 
                    %skip to next row (exactly one non-zero entry in each
                    %row)
                    d=d+1;
                end
            end    
            SS{j,1}=S;
            clear S;
        end
        if c==0
            SS{j,1}=zeros(0,number_of_rows);
        end
            
    end 
    
    %ZERO RESTRICTIONS:
    %representation of zero restrictions: cell array ZZ
    %each element of cell array ZZ is selection matrix for given structural
    %shock
    ZZ=cell(number_of_columns,1);
    for j=1:number_of_columns % for each structural shock
        %number of zero restrictions for j-th structural shock
        c=sum(Y(:,j)==0);
        if c~=0 %if there is at least one zero restriction for given shock
            Z=zeros(c,number_of_rows);
            d=1;
            for i=1:number_of_rows
                if (Y(i,j)==0)
                    Z(d,i)=1;
                    d=d+1;
                end
            end
            ZZ{j,1}=Z;
            clear Z;
        end
        if c==0
            ZZ{j,1}=zeros(0,number_of_rows);
        end
    end

