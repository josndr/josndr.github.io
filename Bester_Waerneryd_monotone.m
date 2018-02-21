% ALGORITHM FOR MONOTONE LOTTERIES

clear
options = optimset('Display','off');
%{
INPUT part. Here we determine the input part. The base line model is one
similar to that discussed in Bester/Waerneryd with a fraction q<1 available
of the pie only in case of escalation.

Change to your liking.
If conditions from text are not met, a Warning message will appear.
%}

% Number of types
n=4;

% Creating the type space 
types=linspace(1,n,n);
t1= repmat(types,n,1);
t2= repmat(types',1,n);

%Match Utility

%DEFAULT GAME FORM

% Fraction of the pie available in the escalation game
q=.8;
% Type sensitivity. Extrem Values
        % If s=1 then type 1 takes the whole pie if matched with type 0. 
        % If s=0 types do not matter and the pie is always shared 50/50.
s=1;

%Utility
u= (1/2 + s*(t1-t2)./(2*(n-1))).*q;

% PRIOR DISTRIBUTION

%Uniform
%   p=pdf('Uniform',types,0,n);
%POISSON
%  p=pdf('Poisson',types,2);
%NORMAL
        p=pdf('Normal',types,n/2,1);
%Adjusting mass to truncate to entire mass on types.
%Mass open is put half half on lowest and highest type
    missingp=1-sum(p);
    p(n)=p(n)+missingp/2;
    p(1)=p(1)+missingp/2;

%{
COMPUTATION PART.
    Runs the top-down algorithm as proposed in the text.
    
    DO NOT CHANGE
%}

    % Calculate Expressions for analysis

% 1. Hazard Rate
w=zeros(n,1);
for i=1:n
    w(i)= (1-sum(p(types<=i)))/p(i);
end



% 2. Outside Option 
V=zeros(n,1);
for i=1:n
	V(i)=sum(p.*u(i,:));
end

% 3. Match ability premium and match welfare
A=ones(n,n);
W=ones(n,n);
for i=1:n-1
    for j=1:n
A(i,j)= u(i,j) - u(i+1,j);
W(i,j)= u(i,j) + u(j,i);
    end
end
W(n,:)=W(:,n);
W(n,n)=u(n,n)+u(n,n);

% 4. Virtual Valuation
VV=zeros(n,n);
for i=1:n
    for j=1:n 
    VV(i,j)= w(i)*A(i,j)+w(j)*A(j,i)+W(i,j)-W(n,n);
    end
end



% ALGORITHM (see Appendix of the paper)

% Define budget Condition
c=@(x) (1-2*u(2,2)-sum(sum(x.*VV)))*(p(1))^2+(2*V(1)-1).*(x(1,1));

% Weights for maximum of all rhos. Needed to determine how far they can be
% increased until gamma(k,t)=1. 
rhoweight=zeros(n,n);
for i=1:n
     for j=1:n
         rhoweight(i,j) = p(i)*p(j)/(p(1))^2;
     end
 end


% 1. First step (only 1-types meet)
rho=zeros(n,n);
rho(1,1)=1;

% Determine maximum given set of rhos
maxrhofun= @(x) sum(sum(rhoweight(rho>0).*x))-1;

% Naive resulat as starting point
starting_point=c(rho); 

% Step 2. 
current_condition=starting_point;
while current_condition>=0


    % Identify second highest VVs
    % Length I is the number of virtual valuations with the same value
    I=find(VV==max(max(VV(rho==0))));
    rho(I)=1;
    
    % Calculate max given set of active rhos (Step 2 a and b)
    newrho11=fsolve(maxrhofun,0,options);
    newrho=zeros(n,n);
    newrho(rho>0)=rhoweight(rho>0).*newrho11;

    % if max gives still fails condition c. increase set of active rhos
    if c(newrho)>=0
        current_condition=c(newrho);
    % Otherwise fine tune the last added ones (Step 2c)
    else
        rhoweighthelp=zeros(n,n);
        rhoweighthelp(rho>0)=rhoweight(rho>0);
        rhoweighthelp(I)=0;
        remainder=@(x) (1-sum(sum(rhoweighthelp.*x)))/length(I);
        csolver=@(r11) (1-2*u(2,2)-sum(sum(rhoweighthelp.*r11.*VV)) - length(I).*remainder(r11).*VV(I(1))).*(p(1))^2+(2*V(1)-1).*r11;
        r11opti=fsolve(csolver,1,options);
        rho=rhoweighthelp.*r11opti;
        rem=remainder(r11opti);
        rho(I)=rem;
        %If solution does not violate maximum for any pair with that VV twerminate      
        if min(rhoweight(I).*rho(I))>rem
            % In case it does redistribute among the last active ones.
        else
            Ihelp=I;
            while isempty(Ihelp)==0
                [minrhoweight,IndexinIofMin]=min(rhoweight(Ihelp));
                rhomax=rhoweight(I(IndexinIofMin))*r11opti;
                if rhomax<rho(Ihelp(IndexinIofMin))
                    rho(Ihelp(IndexinIofMin))=rhomax;
                    difference=rem-rhomax;
                    Ihelp=setdiff(Ihelp,Ihelp(IndexinIofMin));
                    rem=rem+difference/length(Ihelp);
                    rho(Ihelp)=rem;
                else
                    break
                end
            end
        end
        
        current_condition=c(rho);
        break
    end

end

%Display final distribution
disp('READING INSTRUCTIONS: Top Left Corner is all high-types')
disp('Distribution of types after Escalation')
rho
%Display distribution of escalation probabilities
disp('Esclation Probabilities')
gamma=rho./(rhoweight.*r11opti)
% Display Probability of Escalation
disp('Ex-ante Probability of Escalation')
ProbEscalation=p*gamma*p'

%plot escalation prob
% imagesc(gamma)
% colorbar


%{
Show Warnings if no monotone lottery
%}

%Test whether Distribution meets condition
if all(diff(w)<=0)==1
else
    warning(sprintf('Sufficient condition violated Distribution implies non monotonic w\n Results may not be correct.\n Verify that VV is monotone'))
end

% Test whether VV monotone otherwise a warning is displayed
if max(max(max((diff(VV))))+max(max(diff(VV))))<=0
else
    warning(sprintf('VV non-monotone. No monotone lottery.\n If w is monotone (see other warning), change u, such that\n u(t1,t2)-u(t1+1,t2) is weakly decreasing in both arguments'))
end
