function [J] = lin_interp_mat(xo,xi)
    %
    %     Compute the interpolation matrix from xi to xo
    %
    [xs,is]=sort(xo); m=length(xo);n=length(xi);
    J=spalloc(m,n,3*m);


    j0=1;j1=2;x0=xi(j0);x1=xi(j1);
    for i=1:m;
        while xs(i) > x1;
            j0=j1; j1=j0+1;;x0=x1;x1=xi(j1);
        end;
        J(i,j0)=(xs(i)-x1)/(x0-x1);
        J(i,j1)=1-J(i,j0);
    end;
    J(is,:)=J(:,:);
   
end;