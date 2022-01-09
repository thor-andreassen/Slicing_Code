function [coefs,err]=fitQuadricSurf(x,y,z)
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