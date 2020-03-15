%% periodic flow my version
clear all;
clc;

% Macroscopic parameters
NX=128;
NY=128;
NPOP=9;
NSTEPS=200;
NOUTPUT=10;

% Parameters of the binary-liquid model
ksurf=0.04;
gamma=1.0;
a=0.04;
tau_gas=0.7;
tau_liq=2.5;
radius=20;
tau_phi=1.0;

% Stencils parameters
wxx=[0.0, 1.0/3.0, -1.0/6.0, 1.0/3.0, -1.0/6.0, -1.0/24.0, -1.0/24.0, -1.0/24.0, -1.0/24.0];
wyy=[0.0, -1.0/6.0, 1.0/3.0, -1.0/6.0, 1.0/3.0, -1.0/24.0, -1.0/24.0, -1.0/24.0, -1.0/24.0];
wxy=[0.0, 0.0, 0.0, 0.0, 0.0, 1.0/4.0, -1.0/4.0, 1.0/4.0, -1.0/4.0];

gradstencilx=[0.0,4.0/12.0,0.0,-4.0/12.0,0.0,1.0/12.0,-1.0/12.0,-1.0/12.0,1.0/12.0];
gradstencily=[0.0,0.0,4.0/12.0,0.0,-4.0/12.0,1.0/12.0,1.0/12.0,-1.0/12.0,-1.0/12.0];
laplacestencil=[-20.0/6.0,4.0/6.0,4.0/6.0,4.0/6.0,4.0/6.0,1.0/6.0,1.0/6.0,1.0/6.0,1.0/6.0];

% Initialization 
rho=ones(NX,NY);
phi=zeros(NX,NY);
ux=zeros(NX,NY);
uy=zeros(NX,NY);

feq=zeros(NPOP);
geq=zeros(NPOP);
f1=zeros(NPOP,NX,NY);
f2=zeros(NPOP,NX,NY);
g1=zeros(NPOP,NX,NY);
g2=zeros(NPOP,NX,NY);

% Parameters of the lattice
weights=[4/9 1/9 1/9 1/9 1/9 1/36 1/36 1/36 1/36];
cx=[0 1 0 -1 0 1 -1 -1 1];
cy=[0 0 1 0 -1 1 1 -1 -1];

% Initialization
for y=1:NY
    for x=1:NX
       if (x-NX/2)^2+(y-NY/2)^2<=radius*radius
            phi(x,y)=1;
        else
            phi(x,y)=-1;
        end
 
    end
end
    

for y=1:NY
    for x=1:NX      
        rho(x,y)=1;
        ux(x,y)=0;
        uy(x,y)=0;
        vx=ux(x,y);
        vy=uy(x,y);
        
        gradx=0.0;
        grady=0.0;
        laplace=0.0;
        
        for k=1:NPOP
            newx=1+mod(x-1+cx(k)+NX,NX);
            newy=1+mod(y-1+cy(k)+NY,NY);
            gradx=gradx+gradstencilx(k)*phi(newx,newy);
            grady=grady+gradstencily(k)*phi(newx,newy);
            laplace=laplace+laplacestencil(k)*phi(newx,newy);
        end
        
        phase=phi(x,y);
        dense=rho(x,y);
        vx=ux(x,y);
        vy=uy(x,y);
        
        
        sum_dense=0.0;
        sum_phase=0.0;
        phase_square=phi(x,y)*phi(x,y);
        pressure_bulk=dense/3.0+a*(-0.5*phase_square+3.0/4.0*phase_square*phase_square)-ksurf*phase*laplace;
        chemical=gamma*(a*(-phase+phase*phase*phase)-ksurf*laplace);
        
        for k=2:NPOP
            
            feq(k)=weights(k)*(3.0*pressure_bulk+3.0*dense*(cx(k)*vx+cy(k)*vy) ...
                +4.5*dense*((cx(k)*cx(k)-1.0/3.0)*vx*vx+(cy(k)*cy(k)-1.0/3.0)*vy*vy ...
                +2.0*vx*vy*cx(k)*cy(k)))...
                +ksurf*(wxx(k)*gradx*gradx+wyy(k)*grady*grady+wxy(k)*gradx*grady);
            geq(k)=weights(k)*(3.0*chemical+3.0*phase*(cx(k)*vx+cy(k)*vy) ...
                +4.5*phase*((cx(k)*cx(k)-1.0/3.0)*vx*vx+(cy(k)*cy(k)-1.0/3.0)*vy*vy ...
                +2.0*vx*vy*cx(k)*cy(k)));
            sum_dense=sum_dense+feq(k);
            sum_phase=sum_phase+geq(k);
            
            f1(k,x,y)=feq(k);
            g1(k,x,y)=geq(k);
        end
        
        f1(1,x,y)=dense-sum_dense;
        g1(1,x,y)=phase-sum_phase;
        
        
    end
end

%% Main iteration loop
counter_frame=1;
for counter=1:NSTEPS
    
    % Calculation of the macroscopic quantities
    for y=1:NY
        for x=1:NX
            dense=0;
            phase=0;
            vx=0;
            vy=0;
            for k=1:NPOP
                dense=dense+f1(k,x,y);
                phase=phase+g1(k,x,y);
                vx=vx+cx(k)*f1(k,x,y);
                vy=vy+cy(k)*f1(k,x,y);
            end
            
            rho(x,y)=dense;
            ux(x,y)=vx/dense;
            uy(x,y)=vy/dense;
            
            phi(x,y)=phase;
        end
    end
    
    % main loop
    for y=1:NY
        for x=1:NX
             
            % Calculation of the laplacians and gradients and equilibrium
            % functions

         
            % Change tau_rho depending on phase and tau_liq with tau_gas 
            
            for k=1:NPOP
   
                %forcepop=weights(k)*(1-0.5*omega)*((3*(cx(k)-vx)+9*cx(k)*(cx(k)*vx+cy(k)*vy))*forcex(x,y)...
                %   +(3*(cy(k)-vy)+9*cy(k)*(cx(k)*vx+cy(k)*vy))*forcey(x,y));
               
                newx=1+mod(x-1+cx(k)+NX,NX);
                newy=1+mod(y-1+cy(k)+NY,NY);
         
                f1(k,x,y)=f1(k,x,y)*(1.0-1.0/tau_rho)+feq(k)*1.0/tau_rho;
                f2(k,newx,newy)=f1(k,x,y);
                
                g1(k,x,y)=g1(k,x,y)*(1.0-1.0/tau_phi)+geq(k)*1.0/tau_phi;
                g2(k,newx,newy)=g1(k,x,y);
            end
           
        end
    end
    
    f1=f2;
    g1=g2;
    counter

    if mod(counter,NOUTPUT)==0
        imagesc(phi);
        %surf(phi);
        %sum(sum(phi))
        %zlim([-1.1 1.1])
        F(counter_frame) = getframe;
        counter_frame=counter_frame+1;
    end
end

movie(F,10)
