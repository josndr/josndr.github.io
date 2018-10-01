%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Supplementary Material to: Managing a Conflict %
%%%%%%%%%% By B.Balzer and J. Schneider %%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% ADR Solution Program %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This Program Calculates the optimal solution for section 2 
% of the paper as a function of the prior belief p where the only
% input is the difference between the two types: K

% The output contains a series of plots of variables of interest based 
% on the  optimal solution (b_A,b_B). If needed any of these solutions 
% can be exportert as tikz files if the Matlab2Tikz functions are included
% in the matlab path.

%NOTATION: h(igh) korresponds to type theta_i=K, and l(ow) to type theta_i=1

clear all
close all



% Turn fsolve message off
options = optimset('Display','off');



%Output Folder for Pictures
folder='matlab_pictures\';

% Grid Size
n=50; %grid p
m=500; % grid b

%%%%%%%%%%%%%%%%%%
% Input Variable %
%%%%%%%%%%%%%%%%%%
K=3;

% define bounds
p0=(K-2)/(2*(K-1));
pprime=1/(6*(K-1))*(K-8+sqrt(28-4*K+K^2));
pprimeprime=1/(2+3*K)*(2*(K-1)-sqrt(8-4*K+K^2));



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%First Step: Corner Solution
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
gll=1;
%%%%%%%%%%
%FirstPart: bB<=p and gll<=1 are binding.
%%%%%%%%%%

%Only invoke if pprime>0, otherwise non relevant and return empty matrices
if pprime>0
    pcorner1=linspace(0.001,pprime,n);
    bA=linspace(0,(1+pprimeprime)/2,m);
    bBcorner1=zeros(size(pcorner1));
    bAcorner1=bBcorner1;
    for i=1:length(pcorner1)
        p=pcorner1(i);
        bBh=p*ones(size(bA));
        %clear all bA that do not satisfy gll<=1
        const=cons(K,p,bA,bBh,gll);
        bAh=bA;
        bAh(const>1)=NaN;
        objective=obj(K,p,bAh,bBh,gll);
        [~,bAindex]=max(objective);
        bAcorner1(i)=bA(bAindex);
        bBcorner1(i)=p;
    end
    clear i bAindex p bBh
else
    pcorner1=[];
    bAcorner1=[];
    bBcorner1=[];
end

%%%%%%%%%%%%%%
%SecondPart: gll<=1 is binding but bB>p
%%%%%%%%%%%%%%%

pcorner2=linspace(max(0.001,pprime),pprimeprime,n);
bA=linspace((1-p0)/2,(1+p0)/2,m);

bAoptcorner=zeros(size(pcorner2));
bBoptcorner=bAoptcorner;


bBofbA=(1-pcorner2)/2;

for i=1:length(pcorner2)
    p=pcorner2(length(pcorner2)+1-i);
    bBofbA=fsolve(@(bB) cons(K,p,bA,bB,gll)-1, bA,options);

    objective=obj(K,p,bA,bBofbA,gll);
    
    [~,bAindex]=max(objective);
    
    bAoptcorner(length(pcorner2)+1-i)=bA(bAindex);
    bBoptcorner(length(pcorner2)+1-i)=bBofbA(bAindex);
end

if min(cons(K,p,bAoptcorner,bAoptcorner,gll))<1
    warning('Maybe symmetry optima');
    return
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Second Step: Interior Solution
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

pinterior=linspace(pprimeprime,p0,n);

bAoptinterior=(1+pinterior)./2;
bBoptinterior=(1-pinterior)./2;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Plot optimal result
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

bAopt=[bAcorner1,bAoptcorner,bAoptinterior];
bBopt=[bBcorner1,bBoptcorner,bBoptinterior];
p=[pcorner1,pcorner2,pinterior];

gll=cons(K,p,bAopt,bBopt,gll);

[U,z,g,Eg,~]=results(K,p,bAopt,bBopt,gll);
probEscalation= p.*Eg(1,:)+(1-p).*Eg(3,:);
object=obj(K,p,bAopt,bBopt,gll);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Unconstraint Benchmark
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

bAunconstraint=(1+p)./2;
bBunconstraint=(1-p)./2;
gllunconstraint=cons(K,p,bAunconstraint,bBunconstraint,gll);
[Uunconstraint,zunconstraint,gunconstraint,Egunconstraint,~]=results(K,p,bAunconstraint,bBunconstraint,gllunconstraint);
objunconstraint=obj(K,p,bAunconstraint,bBunconstraint,gllunconstraint);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%% Create Plots
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

figure
plot(p,bAunconstraint)
hold on
plot(p,bBunconstraint)
plot(p,bAopt)
plot(p,bBopt)
legend({'bAunconstraint','bBuncsontraint','bA*','bB*'},'location','southeast')
name='b' ;
tit='optimal beliefs constraint and unconstraint';
title(strcat(name,{': '},tit))
box off
print(strcat(folder,name),'-dpdf')
figure
hold on
plot(p,g)
legend({'gll','ghl','glh','ghh'},'location','northeast')
name='g';
tit='escalation probabilities';
title(strcat(name,{': '},tit))
box off
print(strcat(folder,name),'-dpdf')
figure
plot(p,Eg)
hold on
legend({'gAl','gBl','gAh','gBh'},'location','northeast')
name='Eg';
tit='optimal expected escalation probabilities';
title(strcat(name,{': '},tit))
box off
print(strcat(folder,name),'-dpdf')
figure
plot(p,z)
legend({'zAl','zBl','zAh','zBh'},'location','northeast')
tit='optimal settlement values';
name='z';
title(strcat(name,{': '},tit))
box off
print(strcat(folder,name),'-dpdf')
figure
plot(p,z./(1-Eg))
legend({'x1l','x2l','x1h','x2h'},'location','northeast')
tit='optimal settlement shares';
name='x';
title(strcat(name,{': '},tit))
box off
print(strcat(folder,name),'-dpdf')
figure
plot(p,U)
legend({'Ul','UhA','UhB'},'location','northeast')
tit='utilities';
name='U';
title(strcat(name,{': '},tit))
box off
print(strcat(folder,name),'-dpdf')
figure
plot(p,probEscalation)
legend('Pr(Gamma)')
tit='Escalation probability';
name='PrGamma';
title(strcat(name,{': '},tit))
box off
print(strcat(folder,name),'-dpdf')
figure
plot(p,z-zunconstraint)
legend({'1l','2l','1h','2h'},'location','northeast')
tit='Diferences in z (actual-unconstraint)';
name='zDiff';
title(strcat(name,{': '},tit))
box off
print(strcat(folder,name),'-dpdf')
figure
plot(p,g-gunconstraint)
legend({'ll','hl','lh','hh'},'location','northeast')
tit='Diferences in gamma  (actual-unconstraint)';
name='gammaDiff';
title(strcat(name,{': '},tit))
box off
print(strcat(folder,name),'-dpdf')
figure
plot(p,Eg-Egunconstraint)
legend({'1l','2l','1h','2h'},'location','northeast')
tit='Differences in gamma_i (actual-unconstraint)';
name='gammaiDiff';
title(strcat(name,{': '},tit))
box off
print(strcat(folder,name),'-dpdf')
figure
plot(p,(z./(1-Eg))-(zunconstraint./(1-Egunconstraint)))
legend({'1l','2l','1h','2h'},'location','northeast')
tit='Diferences in x_i (actual-unconstraint)';
name='xiDiff';
title(strcat(name,{': '},tit))
box off
print(strcat(folder,name),'-dpdf')
figure
plot(p, z + Eg.*[U(1,:);U])
legend({'1l','2l','1h','2h'},'location','northeast')
tit='Payoff from participating (Pi)';
name='Pi';
title(strcat(name,{': '},tit))
box off
print(strcat(folder,name),'-dpdf')
figure
plot(p, (z + Eg.*[U(1,:);U])-(zunconstraint + Egunconstraint.*[Uunconstraint(1,:);Uunconstraint]))
legend({'1l','2l','1h','2h'},'location','northeast')
tit='Difference in Pi (actual-unconstraint)';
name='PiDiff';
title(strcat(name,{': '},tit))
box off
print(strcat(folder,name),'-dpdf')



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%% Create TikZ plots %%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%% uncomment if needed %%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% set(gca, 'box', 'off')
% sv=4*length(p)/10;
% figure
% plot(p(sv:length(p)), z(:,sv:length(p)))
% box off
% matlab2tikz('settlement_data.tex','externalData',true)

% figure
% plot(p(sv:length(p)), Eg(:,sv:length(p)).*[U(1,sv:length(p));U(:,sv:length(p))])
% box off
% matlab2tikz('escalation_data.tex','externalData',true)

% figure
% plot(p(sv:length(p)), g(:,sv:length(p)))
% box off
% matlab2tikz('escalation_rule_data.tex','externalData',true)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%||||||||||||||||||| Below local functions only ||||||||||||||||%%%%%%%
%%%%%vvvvvvvvvvvvvvvvvvv%%%%%%%%%%%%%%%%%%%%%%%%%%%%vvvvvvvvvvvvvvvv%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%results of Interest
function [U,z,g,Eg,V]=results(K,p,bA,bB,gll)

    X=(K-1)/K;
    Ul=X.*(1-bB);
    UhA=0.*ones(size(bA));
    UhB=X.*(bA-bB);


    
    V=X.*(1-p);

   

    ghl=(p./(1-p)).*((1-bB)./bB).*gll;
    glh=(p./(1-p)).*((1-bA)./bA).*gll;

    ghh=ghl.*((1-bA)./bA).*(p./(1-p));
    
    
    gAl= p.*gll + (1-p).*glh;
    gBl= p.*gll + (1-p).*ghl;
    gAh= p.*ghl + (1-p).*ghh;
    gBh= p.*glh + (1-p).*ghh;
    

    zAl=V-gAl.*Ul;
    zBl=V-gBl.*Ul;
    zAh=zAl+gAl.*UhA-gAh.*UhA;
    zBh=zBl+gBl.*UhB-gBh.*UhB;
    
    U=[Ul;UhA;UhB];
    z=[zAl;zBl;zAh;zBh];
    Eg=[gAl;gBl;gAh;gBh];
    if length(gll)<length(ghl)
        gll=gll*ones(size(ghl));
    else
    end
    g=[gll;ghl;glh;ghh];
end

%value of the objective%
function objective=obj(K,p,bA,bB,gll)
    
    [U,~,~,~,~]=results(K,p,bA,bB,gll);
    Ul=U(1,:);
    UhA=U(2,:);
    UhB=U(3,:);
    
    r1= bB;
    r2= bA;
    
    EU= r1.*Ul + r2.*Ul + (1-r2).*UhB;
    D1=Ul-UhA;
    D2=(Ul-UhB);
    
    Epsi= r1.*D1.*(1-p)./p + r2.*D2.*(1-p)./p;
    
    objective= Epsi + EU;
    
end

%%Boundary of the set of Informaiton structures
function constraint=cons(K,p,bA,bB,gll)

    [~,z,~,~,V]=results(K,p,bA,bB,gll);

    zAl=z(1,:);
    zBl=z(2,:);
    zAh=z(3,:);
    zBh=z(4,:);

    v=2.*V-1;
    
    Q=(2.*V-(p.*zAl+(1-p).*zAh +p.*zBl+(1-p).*zBh))./gll;

    R=(p.^2)./(bA.*bB);

    constraint=v./(Q-R);
end


