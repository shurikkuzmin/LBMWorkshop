%% periodic flow my version
clear all;
clc;

% Macroscopic density and velocities
NX=16;
NY=16;
NPOP=9;
NSTEPS=2500;

rho0=1;
umax=0.001;

rho=ones(NX,NY);
ux=zeros(NX,NY);
uy=zeros(NX,NY);

uxinit=zeros(NX,NY);
uyinit=zeros(NX,NY);

feq=zeros(NPOP);
f1=zeros(NPOP,NX,NY);
f2=zeros(NPOP,NX,NY);

weights=[4/9 1/9 1/9 1/9 1/9 1/36 1/36 1/36 1/36];
cx=[0 1 0 -1 0 1 -1 -1 1];
cy=[0 0 1 0 -1 1 1 -1 -1];
omega=1.0;

for y=1:NY
    for x=1:NX
        rho(x,y)=rho0+3*0.25*umax^2*...
            (cos(4*pi*(x-1)/NX)-cos(4*pi*(y-1)/NY));
        ux(x,y)=umax*sin(2*pi*(x-1)/NX)*sin(2*pi*(y-1)/NY);
        uy(x,y)=umax*cos(2*pi*(x-1)/NX)*cos(2*pi*(y-1)/NY);
        vx=ux(x,y);
        vy=uy(x,y);
        
        %Please initialize it with the equilibrium function        
    end
end

for counter=1:NSTEPS
    for y=1:NY
        for x=1:NX
            dense=0;
            vx=0;
            vy=0;
            
            %Construct densitie and velocities by summation
            rho(x,y)=dense;
            vx = vx/dense;
            vy = vy/dense;
            ux(x,y)=vx;
            uy(x,y)=vy;
             
            for k=1:NPOP
                %Equilibrium distribution construction
                %newx,newy is for streaming
                
                newx=1+mod(x-1+cx(k)+NX,NX);
                newy=1+mod(y-1+cy(k)+NY,NY);
         
                %Collsion 
                %Streaming from f1 to f2
            end
           
        end
    end
    
    track(counter)=ux(NX/4+1,NY/4+1);
    decay(counter)=umax*exp(-1/3*(1/omega-0.5)*counter*2*(2*pi/NX)^2);
    f1=f2;
end

plot(1:NSTEPS,track,'+')
xlim([0 200])
hold on
plot(1:NSTEPS,decay,'o')
