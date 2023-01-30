function [u,res,iter]= Jacobi(Nx,Ny,xmin,xmax,ymin,ymax,itgmr,omega,tol)

    ii=[]; ismax=[]; x=[]; y =[]; u=[]; rhs=[]; qrr=[]; 
    du=[]; uold=[];
    
    string();
    
    %GRID Coordinates

    dx = (xmax - xmin)/(Nx-1);
    dy = (ymax - ymin)/(Ny-1);
    for i=2:Nx+1
        for j=2:Ny+1
            L=ii(i,j);
            x(L)=xmin+dx*(i-2);
            y(L)=ymin+dy*(j-2);
            u(L)=0;
        end
    end
   
    
    rhsDiri();
    
    

    solver();
    
    

        function string()
    
            %initialization of ii
            for i=1:Nx+2
                for j=1:Ny+2
                    ii(i,j) = 1;
                end
            end
    
            ismax=0; %1D counter of nodes
            for i=2:Nx+1
                for j=2:Ny+1
                    ismax=ismax+1;
                    ii(i,j)=ismax;
                end
            end
        end
        
        function rhsDiri()

            for i=2:Nx+1
                for j=2:Ny+1 
                    L=ii(i,j);
                        if i==2 | i==Nx+1 | j==2 | j==Ny+1
                            rhs(L)=x(L)^2 + y(L)^2;
                            u(L)=rhs(L);
                        else
                            rhs(L)=-4;
                        end
                end
            end
    
          
        end

           
    function solver()

        for L=1:ismax+1
            du(L)= 0;
            uold(L)=0;
        end

        itmmn = 1;

        for iter=itmmn:itgmr

            %Delta Formulation - Compute Residual

            uold=u;

            for i=3:Nx
                for j=3:Ny

                    L=ii(i,j);
                    ie=ii(i+1,j);
                    iw=ii(i-1,j);
                    in=ii(i,j+1);
                    is=ii(i,j-1);

                    
                    du(L)=(-1/(2*(dx^2+dy^2)/((dx^2)*(dy^2))))*(-rhs(L)-...
                        (1/dx^2*uold(ie)+1/dx^2*uold(iw)+ ...
                        1/dy^2*uold(in)+1/dy^2*uold(is)))-uold(L);
                    u(L) = omega*(uold(L)+du(L)) + (1-omega)*uold(L);
                    

                end
            end


            resit=0;

            for i=2:Nx+1
                for j=2:Ny+1

                    L=ii(i,j);
                    ie=ii(i+1,j);
                    iw=ii(i-1,j);
                    in=ii(i,j+1);
                    is=ii(i,j-1);
                    

                    if i==2 | i==Nx+1 | j==2 | j==Ny+1
                        qrr(L) = u(L)-rhs(L);
                    else
                        qrr(L) = -( ...
                            ...
                            1/dx^2*u(ie)+1/dx^2*u(iw)+1/dy^2*u(in)+1/dy^2*u(is)  ... 
                            -2*(dx^2+dy^2)/((dx^2)*(dy^2))*u(L)+rhs(L));
                    end

                    resit=resit+qrr(L)*qrr(L);
                    
                end

            end

            resit=sqrt(resit)/ismax;

            res(iter)=resit;

            if resit<tol; break; end
            
            
            
           
            
        end

    end

end


            

    



