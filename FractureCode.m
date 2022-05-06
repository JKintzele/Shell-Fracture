% Ice Shelf/Shell Fracture Stability
% Jakob Kintzele, Princeton University Geosciences
% Last Update: 4/17/2022
% Companion Code to CoolingElasticShellStress.m (https://github.com/JKintzele/Cooling-Shell-Stress)
% NOTE: DO NOT RUN ALL AT ONCE. GO THROUGH EACH SECTION AT A TIME

strength=10^5; %ice tensile strength 0.1-1 [Pa]
Kc=150*10^3; %ice fracture toughness [Pa m^1/2]
sigmatide=10^5; %uniform tidal tension
P=-flip(sigmaW+sigmaT+sigmaV+sigmaVT); %point loads [Pa], can upload these from CoolingElasticShellStress.m
[mx,loc]=max(P); %max tension
z=flip(R-rfix(4,1:end)); %depth [m]
ds=flip(R-rfix(4,1:end)); %fracture depth [m]
H=R-rv; %shell thickness [m]
location=z(loc); %max tension depth
%%  ==== elastic layer point loads, from basal reference frame ==== %%
P_elastic=zeros(N,Nt); 
for i=1:Nt
    for j=1:N
        if Tcart(j,i)<=Te
P_elastic(j,i)=-(sigmaT(j,i)+sigmaV(j,i)+sigmaVT(j,i))+Pex(i)+PexT(i);
        end
    end 
end
%% ==== Surface Fracture Stability ==== %%

figure(576)
    grid on
    hold on
plot([0 1],[Kc Kc],'k--')
xlim([0 1])
%ylim([-Kc Kc*10])
xlabel('d_s/H')
ylabel('K_I [Pa m^{1/2}]')
title('Surface Fracture Stability')
%%
sigmatide=10^5;
num=Nt;
stabledepthS=zeros(1,num);
KIs=zeros(N);
Fn=zeros(j,1);
for i=1:1000 %# of time intervals    
    for j=2:N %loop over crack depths ds(j)<H        
zt=zeros(j,1); G=zeros(j,1);  g1=0; g2=0; g3=0; g4=0; dst=0; Fu=0;   
    if P(1,i)+sigmatide>=strength
if ds(j)<H(i)
        dst=ds(j)./H(i); %relative crack depth to current shell thickness 
            for k=1:j
        zt(k)=z(k)./ds(j); %relative depths up to the crack tip
%dimensionless weights, reset for each crack depth:        
g1=.46+3.06.*dst+.84.*(1-dst).^5+.66.*dst.^2.*(1-dst).^2;
g2=-3.52.*dst.^2;
g3=6.17-28.22.*dst+34.54.*dst.^2-14.39.*dst.^3-(1-dst).^1.5-5.88.*(1-dst).^5-2.64.*dst.^2.*(1-dst).^2;
g4=-6.63+25.16.*dst-31.04.*dst.^2+14.41.*dst.^3+2.*(1-dst).^1.5+5.04.*(1-dst).^5+1.98.*dst.^2.*(1-dst).^2;
G(k)=g1+g2.*zt(k)+g3.*zt(k).^2+g4.*zt(k).^3;

%integrate F for stress intensity factor at the tip of a crack at depth ds
Fn(k)=2.*P(k,i).*G(k)./((pi.*ds(j)).^0.5.*(1-dst)^1.5.*(1-zt(k).^2).^0.5); 
Fn(j)=0;       

%uniform tidal tension
Fu=(2/(pi*dst)*tan(0.5*pi*dst))^0.5*(0.752+2.02*dst+0.37*(1-sin(0.5*pi*dst))^3)/cos(0.5*pi*dst);

KIs(j)=trapz(z(1:j),Fn(1:j))+sigmatide*(pi*ds(j))^0.5*Fu;
            end
if j>1 && KIs(j-1)>=Kc && KIs(j)<Kc
    stabledepthS(i)=0.5*(ds(j-1)+...
    ds(j))./H(i);
elseif KIs(j)>Kc
    stabledepthS(i)=1;
end      
%plot(ds(1:j)./H(i),KIs(i,1:j))
end
    end
    end
end
%% ==== Basal Fracture Stability ==== %%
figure(577)
    grid on
    hold on
ft=plot([0 1],[Kc Kc],'k--');
xlim([0 1])
%ylim([-10^6 4*10^7])
vl=plot([1-lambdae/lambda 1-lambdae/lambda],[-10^6 4*10^7],'c--','LineWidth',2);
xlabel('d_b/H')
ylabel('K_I [Pa m^{1/2}]')
title('Basal Fracture Stability')
%%
sigmatide=10^5;
h=[linspace(0,1,10),linspace(2,H(end),100)]; %height within fracture
db=[linspace(0,1,10),linspace(2,H(end),100)]; %fracture height
num=100;
stabledepthB=zeros(1,num);
for i=1:num %time
    KIb=zeros(length(h),1);
    for j=2:length(h) %loop over crack heights        
hb=zeros(j,1); fb=0; gb=0; g1=0; g2=0; g3=0; g4=0; Fu=0; dbt=0; %reset scalars

        if db(j)<=H(i*length(t)/(num)) && db(j)>0                        
        dbt=db(j)./H(i*length(t)/(num)); %relative crack height to current shell thickness 
            for k=1:j %loop for stresses within the crack
        hb(k)=h(k)./db(j); %relative heights up to the crack tip
        %dimensionless weights, reset for each crack depth:
g1=.46+3.06.*dbt+.84.*(1-dbt).^5+.66.*dbt.^2.*(1-dbt).^2;
g2=-3.52.*dbt.^2;
g3=6.17-28.22.*dbt+34.54.*dbt.^2-14.39.*dbt.^3-(1-dbt).^1.5-5.88.*(1-dbt).^5-2.64.*dbt.^2.*(1-dbt).^2;
g4=-6.63+25.16.*dbt-31.04.*dbt.^2+14.41.*dbt.^3+2.*(1-dbt).^1.5+5.04.*(1-dbt).^5+1.98.*dbt.^2.*(1-dbt).^2;
fb=(g2 + (2*g4)/3 + ((2*g1 + g3)*pi)/4 - (dbt*(48*g1 + 32*g3 + 12*g2*pi + ... 
    9*g4*pi))/48)/((1 - dbt)^(3/2)*dbt); %ice overburden (integrated using Mathematica)
gb=(4*(12*g2 + 8*g4 + 3*(2*g1 + g3)*pi)*rho_i - dbt*(48*g1 + 32*g3 + ...%hydrostatic crack water pressure, 
    12*g2*pi + 9*g4*pi)*rho_w)/(48*(1 - dbt)^(3/2)*dbt*rho_w); %integrated using Mathematica

%uniform tidal tension
Fu=(2/(pi*dbt)*tan(0.5*pi*dbt))^0.5*(0.752+2.02*dbt+0.37*(1-sin(0.5*pi*dbt))^3)/cos(0.5*pi*dbt);
%in viscous layer:
%if db(j)/H(i*length(t)/(num))<lambda/lambdae
KIb(j)=2*g*db(j)^1.5/sqrt(pi)*(rho_w*gb-rho_i*fb)+sigmatide*(pi*db(j))^0.5*Fu; %net stress intensity factor
%in elastic layer, include ocean pressurization:
%elseif db(j)/H(i*length(t)/(num))>=lambda/lambdae 
%   KIb(j)=2*g*db(j)^1.5/sqrt(pi)*(rho_w*gb-rho_i*fb)+(sigmatide+Pex(i)+PexT(i))*(pi*db(j))^0.5*Fu;
%end

if j>1 && KIb(j-1)>=Kc && KIb(j)<Kc
    stabledepthB(i)=0.5*(db(j-1)+...
    db(j))./H(i*length(t)/(num));
end
            end
            
%plot(db(1:j)./H(i*length(t)/(num)),KIb(1:j))
        end
            if KIb(j)<-Kc %for efficiency, stop crack height loop once KI is negative
                KIb(j:end)=0;
            break
            end             
    end       
end
%%
legend([ft,vl],{'Ice Fracture Toughness','Elastic Layer Height'},'Location','NW')
%% ==== Viscous Layer Penetrating Basal Fractures ==== %%
hfix=rfix(4,:)-rv(end); %height within fracture
dbfix=rfix(4,:)-rv(end); %fracture height
%%
num=100;
stabledepthU=zeros(1,num);
for i=1:Nt%time
    KIu=zeros(length(hfix),1);
    KIe=zeros(length(hfix),1);
    for j=2:length(hfix) %loop over crack heights        
 hb=zeros(j,1); fb=0; gb=0; g1=0; g2=0; g3=0; g4=0; Fu=0; dbt=0; G=zeros(j,1); Fn=zeros(j,1); %reset 

        if dbfix(j)<=H(i)                      
        dbt=dbfix(j)./H(i); %relative crack height to current shell thickness 
            for k=1:j %loop for stresses within the crack
        hb(k)=hfix(k)./dbfix(j); %relative heights up to the crack tip
        %dimensionless weights, reset for each crack depth:
g1=.46+3.06.*dbt+.84.*(1-dbt).^5+.66.*dbt.^2.*(1-dbt).^2;
g2=-3.52.*dbt.^2;
g3=6.17-28.22.*dbt+34.54.*dbt.^2-14.39.*dbt.^3-(1-dbt).^1.5-5.88.*(1-dbt).^5-2.64.*dbt.^2.*(1-dbt).^2;
g4=-6.63+25.16.*dbt-31.04.*dbt.^2+14.41.*dbt.^3+2.*(1-dbt).^1.5+5.04.*(1-dbt).^5+1.98.*dbt.^2.*(1-dbt).^2;
fb=(g2 + (2*g4)/3 + ((2*g1 + g3)*pi)/4 - (dbt*(48*g1 + 32*g3 + 12*g2*pi + ... 
    9*g4*pi))/48)/((1 - dbt)^(3/2)*dbt); %ice overburden
gb=(4*(12*g2 + 8*g4 + 3*(2*g1 + g3)*pi)*rho_i - dbt*(48*g1 + 32*g3 + ...%hydrostatic crack water pressure
    12*g2*pi + 9*g4*pi)*rho_w)/(48*(1 - dbt)^(3/2)*dbt*rho_w);

Fu=(2/(pi*dbt)*tan(0.5*pi*dbt))^0.5*(0.752+2.02*dbt+0.37*(1-... %uniform tidal tension
    sin(0.5*pi*dbt))^3)/cos(0.5*pi*dbt);

% cooling stresses:
G(k)=g1+g2.*hb(k)+g3.*hb(k).^2+g4.*hb(k).^3; 
%integrate F for cooling component of KI
Fn(k)=2.*P_elastic(k,i).*G(k)./((pi.*db(j)).^0.5.*(1-dbt)^1.5.*(1-hb(k).^2).^0.5); 
Fn(j)=0;     
KIe(j)=trapz(hfix(1:j),Fn(1:j));

%Net KI:
%in viscous layer:
%if dbfix(j)/H(i)<(1-lambdae/lambda)
%KIu(j)=2*g*dbfix(j)^1.5/sqrt(pi)*(rho_w*gb-rho_i*fb)+...
%    sigmatide*(pi*dbfix(j))^0.5*Fu; %net stress intensity factor

%in elastic layer, include crack pressurization and tensionial stresses:
%elseif dbfix(j)/H(i)>=(1-lambdae/lambda)
   KIu(j)=2*g*dbfix(j)^1.5/sqrt(pi)*(rho_w*gb-rho_i*fb)+... % hydrostatic stress
       (sigmatide)*(pi*dbfix(j))^0.5*Fu+...% tidal load 
       KIe(j);                                               % cooling stresses +excess crack pressure
%end
% find stable crack height:
if KIu(j-1)>=Kc && KIu(j)<Kc
    stabledepthU(i)=0.5*(dbfix(j-1)+...
    dbfix(j))./H(i);
elseif KIu(j)>=Kc && dbfix(j)>=H(i)
    stabledepthU(i)=1;
end
            end
            if KIu(j)<-Kc %for efficiency, stop crack height loop once KI is negative
                KIu(j:end)=0;
            break
            end  
            
        end
                 
    end    
end


%% ==== Figures ==== %%
figure(578)
grid on
hold on
b1=plot(H(1:Nt/100:Nt)./10^3,1-stabledepthB,'b.');
b2=plot(H(1:Nt/num:Nt)./10^3,stabledepthS,'r.');
%b3=plot(H(1:Nt/num:Nt)./10^3,1-stabledepthU,'g.');
vl=plot([min(H) max(H)]./10^3,[lambdae/lambda lambdae/lambda],'c--','LineWidth',2);
%legend([b2 b1 b3],{'Surface Fracture','Basal Crevasse (unstressed shell)','Basal Crevasse (stressed shell)'},'Location', 'NE')
text(5, 0.02, 'Surface')
text(5, 0.4,'Elastic Layer')
text(5, 0.98, 'Ocean')
xlabel('Shell Thickness [km]')
ylabel('Relative Penetration Depth')
xlim([min(H), max(H)]./10^3)
set(gca, 'YDir', 'Reverse')
 iptsetpref('ImshowBorder','tight')
 
 %%
 figure(588)
 subplot(2,1,2)
 grid on
 hold on

a1=plot(H(1:Nt/1000:Nt)./10^3,H(1:Nt/1000:Nt)./10^3,'k-','LineWidth',2);
a2=plot(H(1:Nt/1000:Nt)./10^3,(lambdae/lambda).*H(1:Nt/1000:Nt)./10^3,'c--','LineWidth',2);
%a3=plot(H(1:Nt/1000:Nt)./10^3,(1-rho_i/rho_w).*H(1:Nt/1000:Nt)./10^3,'m--','LineWidth',1);
a4= plot(H(39:Nt)./10^3,(1-stabledepthU(39:Nt)).*H(39:Nt)./10^3,'b.');
a6=plot(linspace(H(2)./10^3,H(end)./10^3,100),zeros(1,100),'r.');
a7=plot(H(500:500:Nt)./10^3,location(500:500:Nt)./10^3,'m--');
legend([a1,a2,a4,a6,a7],{'Shell Depth','Elastic Layer','Basal Fracture Height','Surface Fracture Depth','Depth of Maximum Tension'},'Location','SW')
set(gca, 'YDir', 'Reverse')
xlim([5 H(end)./10^3])
xlabel('Shell Thickness [km]')
ylabel('Depth [km]')
%title('Fractures In a Cooling Ice Shell')

subplot(2,1,1)
 grid on
 hold on
plot(H(1:Nt/1000:Nt)./10^3,H(1:Nt/1000:Nt)./10^3,'k-','LineWidth',2);
plot(H(1:100:800)./10^3,location(1:100:800)./10^3,'m--');
plot(H(1:Nt/1000:Nt)./10^3,(lambdae/lambda).*H(1:Nt/1000:Nt)./10^3,'c--','LineWidth',2);
%plot(H(1:Nt/1000:Nt)./10^3,(1-rho_i/rho_w).*H(1:Nt/1000:Nt)./10^3,'m--','LineWidth',1);
plot(H(39:1000)./10^3,(1-stabledepthU(39:1000)).*H(39:1000)./10^3,'b.');
a5=plot(H(1:5:35)./10^3,zeros(size((1:5:35))),'b*');
plot([H(1)./10^3, linspace(H(2)./10^3,H(end)./10^3,100)],[(stabledepthS(1)).*H(1)./10^3, zeros(1,100)],'r.');
legend(a5,'Surface Eruption','Location','SW')
% plot(H(1:Nt/100:Nt)./10^3,(1-stabledepthB).*H(1:Nt/100:Nt)./10^3,'b.');
set(gca, 'YDir', 'Reverse')
xlim([H(1)./10^3 5])
%xlabel('Shell Thickness [km]')
ylabel('Depth [km]')
title('Fractures In a Cooling Ice Shell')