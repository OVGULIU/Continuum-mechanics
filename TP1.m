%%%%%%%%%%%%%%%%%%%%%%
%%MMC - TP1         %%
%%Yann LE GUILLY    %%
%%%%%%%%%%%%%%%%%%%%%%

close all;clear;clc;

d=1; L=5; A=d*d; E=210e9; nu=1/2;
lam=(E*nu)/((1+nu)*(1-2*nu));
mu=E/(2*(1+nu));

%On construit les points pour tracer la barre que l on va etirer

ttM=[...
    +d/2,-d/2,0;...
    -d/2,-d/2,0;...
    -d/2,+d/2,0;...
    +d/2,+d/2,0;...
    +d/2,-d/2,0;...
    +d/2,-d/2,L;...
    +d/2,+d/2,L;...
    -d/2,+d/2,L;...
    -d/2,-d/2,L;...
    +d/2,-d/2,L;...
    ];nbM=length(ttM(:,1));

t=[-0.1:0.01:0.2];
nbt=length(t); %[s]

a0=1/2; a=a0*t;
b0=3;   b=b0*t;

figure('Position',[20 250 450 510]);figure(1);
for it=1:nbt
    ttm(:,1:2)=(1-a(it))*ttM(:,1:2);
    ttm(:,3)  =(1+b(it))*ttM(:,3);
    figure(1);
    
    hold off;
    pM=plot3(ttM(:,1),ttM(:,2),ttM(:,3),'-sk');
    hold on;
    
    pm=plot3(ttm(:,1),ttm(:,2),ttm(:,3),'-or');
    legend([pM(1),pm(1)],'reference','actuel')
    
    axis equal;view(3)
end


figure('Position',[550 550 500 200]);figure(2);hold on;
II=eye(3);
for it=1:nbt
    FF=[1-a(it),0,0;0,1-a(it),0;0,0,1+b(it)];
    %grande deformation
    CC=FF'*FF;
    EE=1/2*(CC-II);
    ENER(it)=lam/2*(trace(EE))^2+mu*trace(EE^2);
    Gener=plot(ENER,'-or');
    VTV0(it)=det(FF);
    
    %HPP
    gU=FF-II;
    ee=1/2*(gU+gU');
    ener(it)=lam/2*(trace(ee))^2+mu*trace(ee^2);
    Pener=plot(ener,'-sk');
    vtv0(it)=trace(ee)+1;
    legend([Gener(1),Pener(1)],'Grandes deformations','Petite deformations')
    title('Energies pour grandes et petites deformations');
end

figure('Position',[550 50 400 400]);figure(3);hold on;
warning('off','all') %pour eviter l alerte sur les valeurs complexe
for it=1:nbt  
    %grande deformation
    at=1-sqrt(1-nu*b(it)*(b(it)+2));    %Alpha
    Ft(it)=E*A*b(it)*(b(it)+2)/2;   %Force
    GF=plot(at, Ft(it),'--sk');
 
    %HPP
    at=nu*b(it);    %Alpha
    ft(it)=E*A*b(it);   %Force
    PF=plot(at, ft(it),'--or');
    legend([GF(1),PF(1)],'Grandes deforamtions','Petites deformations')
    title('Forces en grandes et petites defroamtions');
end
warning('on','all')

figure('Position',[1000 50 500 400]);figure(4);hold on;
for it=1:nbt  
    %grande deformation
    at=1-sqrt(1-nu*b(it)*(b(it)+2));    %Alpha
    Ft(it)=E*A*b(it)*(b(it)+2)/2;   %Force
 
    %HPP
    at=nu*b(it);    %Alpha
    ft(it)=E*A*b(it);   %Force
    Err(it)=Ft(it)-ft(it);
    PF=plot(at, Err(it), '-or');
    title('Difference entre les forces pour grandes et petites deformations');
    legend([GF(1),PF(1)],'Grandes deforamtions','Petites deformations')
end
