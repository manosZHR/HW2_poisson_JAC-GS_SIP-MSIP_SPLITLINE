function [u,res,iter]= MSIP5(Nx,Ny,xmin,xmax,ymin,ymax,psi,itgmr,tol)

    ii=[]; ismax=[]; x=[]; y =[]; u=[]; rhs=[]; qrr=[]; 
    am=[]; bm=[]; cm=[]; dm=[]; em=[]; fm=[]; gm=[];  hm=[]; km=[]; aux=[]; du=[];
    
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
    
    

    msip5();
    
    solver();
    
    

        function string()
    
            %initialization of ii
            for i=1:Nx+2
                for j=1:Ny+2
                    ii(i,j) = 1;
                end
            end
    
            ismax=0; %1D counter of nodes
            %καθε φορα που θελω να λουπαρω στα nodes θα παιρνω i=2:Nx+1,
            %j=2:Ny+1 αλλα το L μεσω του ii θα βγαινει σωστο
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

           
        
        function msip5()

         
            for L=1:ismax+1
                am(L)=0;
                bm(L)=0;
                cm(L)=0;
                dm(L)=0;
                em(L)=0;
                fm(L)=0;
                gm(L)=0;
                hm(L)=0;
                km(L)=0;
            end

            %MSIP, stencil 9 nodes
         

            for i=2:Nx+1
                for j=2:Ny+1
                    L=ii(i,j);
                    ie=ii(i+1,j);
                    iw=ii(i-1,j);
                    in=ii(i,j+1);
                    is=ii(i,j-1);
                   
                    %capital letters
                    if i==2 | i==Nx+1 | j==2 | j==Ny+1
                        ec=1;
                        hc=0;
                        bc=0;
                        fc=0;
                        dc=0;
                    else
                        %change when dx,dy are inserted
                        ec=-2*(dx^2+dy^2)/((dx^2)*(dy^2));;
                        hc=1/dx^2;
                        bc=1/dx^2;
                        fc=1/dy^2;
                        dc=1/dy^2;
                    end

                    %small letters, elements of L and U

                    bm(L)=bc/(1+psi*fm(iw));
                    dm(L)=dc/(1+psi*hm(is));
                    em(L) = ec - bm(L)*hm(iw) - dm(L)*fm(is) ... 
                        + psi * ( bm(L)*fm(iw) + dm(L)*hm(is) );
                    if abs(em(L))<1e-10; break; end
                    em(L)=1/em(L); %attention keeps the INVERSE of EPSILON
                    fm(L) = ( fc - psi*bm(L)*fm(iw) ) * em(L);
                    hm(L) = ( hc - psi*dm(L)*hm(is) ) * em(L);
                end
            end
        
        
        end



    function solver()

        itmmn = 1;

        for iter=itmmn:itgmr

            %Delta Formulation - Compute Residual
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
            
            
            %Back-Front Substitution

            backfron();

            for L=1:ismax
                u(L)=u(L)+du(L);
            end
            
            
        end


    end

    function backfron()
        
        for L=1:ismax+1
            aux(L) = 0;
            du(L)= 0;
        end

        %Front Substitution

        for i=2:Nx+1
            for j=2:Ny+1

                L=ii(i,j);
                iw=ii(i-1,j);
                is=ii(i,j-1);
               
               
                comw=bm(L)*aux(iw);
              
                coms=dm(L)*aux(is);
                aux(L)=(qrr(L)-comw-coms)*em(L);
            end
        end
        
        %Back Substitution

        for i=Nx+1:-1:2
            for j=Ny+1:-1:2
                L=ii(i,j);
                ie=ii(i+1,j);
                in=ii(i,j+1);
                
                comn=fm(L)*du(in);
               
                come=hm(L)*du(ie);
                
                du(L)=aux(L)-comn-come;
            end
        end
    end
end


            

    



