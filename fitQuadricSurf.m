function [coefs,err]=fitQuadricSurf(x,y,z)
        %% Function to fit a quadratic surface to a set of nodes
        % Created by Thor Andreassen
        % 1/10/22
        % the function takes in x, y, z coordinates as array of nodes
        % connected on a portion of a mesh, and fits a quadratic surface to
        % these nodes and reruns the corresponding coefficients.
        % the output coeficicents following the following pattern in order:
        % Q(x,y,z)=ax^2+bxy+cy^2+dx+ey+f
        b=z;
        A=zeros(length(x),6);
        for counti=1:length(x)
                A(counti,1)=x(counti).^2;
                A(counti,2)=x(counti).*y(counti);
                A(counti,3)=y(counti).^2;
                A(counti,4)=x(counti);
                A(counti,5)=y(counti);
                A(counti,6)=1;
        end
        coefs=pinv(A)*b;
        err=norm(A*coefs-b);
end