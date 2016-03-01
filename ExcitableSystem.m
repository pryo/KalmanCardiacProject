%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% requena_carrion.m
%
% J Requena-Carrion, 2016
%
% 1. Simulates the stochastic model of excitable cell proposed in the manuscript: 
%
% Requena-Carrion J & Requena-Carrion VJ, "Distribution of transition times 
% in a stochastic model of excitable cell: Insights into the cell-intrinsic 
% mechanisms of randomness in neuronal inter-spike intervals"
%
%
% 2. Provides the histograms for the state probability distributions, 
% the distributions of transition times and the inter-spyke interval durations.
%
%
% 3. Plost the exact solutions of the distributions.
%
%
% 4. Can be extended to a 3-state model for reproducing subthreshold
% oscillations.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


% Simulation parameters and initialization

N = 10000;           % Population size, i.e. number of (uncoupled) systems
T = 4;               % Simulation duration
Deltat = 10^(-4);    % Time step

S = ones(1,N);          % States vector, can take the values 1 (excited), 2 (rest) and 3 (rest)
nextS = zeros(size(S)); % Identifies which systems change state

Trest = zeros(size(S)); % Last time instant a system entered S1 from S2
Texc = zeros(size(S));  % Last time instant a system entered S2
Tsub = zeros(size(S));  % Last time instant a system entered S3
Tnu = zeros(size(S));   % Last time instant a system entered S1
nextTnu = zeros(size(S));% Identifies which systems need update Tnu

RT = [];                % Collection of measured total resting times
APD = [];               % Collection of measured action potential durations
ISI = [];               % Collection of measured inter spyke intervals

P_exc  = zeros(size(S));% Probability of exiting S1
P_reset = 0;            % Probability of entering S3 after exiting S1. If 0, 
                        % the two-state variant is simulated

kqDt = 10^2*Deltat;     % Product kq
nextAPD = ones(size(S));% Next APD after excitation 

lambda1 = 0.5;          % Excitability recovery rate
lambda2 = 10;           % APD recovery rate


for t=0:Deltat:T
        
    S1 = S==1;
    S2 = S==2;
    S3 = S==3; 
    
    
    
    % Systems that will exit S1   
    P_exc(S1) = kqDt * (1-exp(-lambda1*(t-Tnu(S1))));
    nextS(S1) = P_exc(S1)>rand(size(P_exc(S1)));

    % If P_reset=0, a 2S model is implemented and Rest2Rest plays no role. 
    % Otherwise, Rest2rest is used to determine which systems exiting S1 
    % enter S2 and which enter S3
    Rest2Rest = P_reset > rand(size(S1));
    nextS12 = (nextS==1)&(S1)&(~Rest2Rest);
    nextS13 = (nextS==1)&(S1)&Rest2Rest;
    nextAPD(nextS12) = (1-exp(-lambda2*(t-Tnu(nextS12))));% 1; %
    nextAPD(nextS13) = (1-exp(-lambda2*(t-Tnu(nextS13))));% 1; %

    
        
    % Systems that will leave S2
    nextS(S2) = (t-Texc(S2))>=nextAPD(S2);
    nextS21 = (nextS==1)&(S2);
      
    
    
    % Systems that will leave S3
    nextS(S3) = (t-Tsub(S3))>=nextAPD(S3);
    nextS31 = (nextS==1)&(S3);
    
    
    
    % Measured events
    ISI = [ISI (t-Texc(nextS12 & (Texc~=0)))];
    RT = [RT (t-Trest(nextS12))];
    APD = [APD nextAPD(nextS12)];
    
    
    
    % Times update
    Texc(nextS12) = t;
    Tsub(nextS13) = t;
    
    Trest(nextS21) = t;
    Tnu(nextS21) = t;
    Tnu(nextS31) = t;

    
    % State update
    S(nextS12) = 2;
    S(nextS21) = 1;
    S(nextS13) = 3;
    S(nextS31) = 1;
    
    nextS = zeros(size(S));
    
end



%%
Nbins=60;
% Data for the steady-state distributions rho and sigma
S1 = S==1;
taus1 = t-Trest(S1);
S2 = S==2;
taus2 = t-Texc(S2);


% Analytical solutions for rho, alpha and beta (No analytical solution is
% provided for sigma)
T=2;
Delta_tau=0.001;
q1=(kqDt/(Deltat*lambda1));
q2=(kqDt/(Deltat*lambda2));

tau1end=1-Delta_tau;
tau1=0:Delta_tau:tau1end;
tau2end=tau1end;
tau2=0:Delta_tau:tau2end;
tauISI=0:Delta_tau:tau1end+tau2end;

rho=exp(-(lambda1*tau1+exp(-lambda1*tau1))*q1);
rhon=rho/sum(rho*Delta_tau);

alpha=(kqDt/Deltat)*(1-exp(-lambda1*tau1)).*exp(-(lambda1*tau1+exp(-lambda1*tau1))*q1);
alphan=alpha/(sum(alpha)*Delta_tau);

beta=(kqDt/Deltat)*((1-(1-tau2).^(lambda1/lambda2))).*((1-tau2).^q2).*exp(-q1*((1-tau2).^(lambda1/lambda2)))./(lambda2*(1-tau2));
betan=beta/(sum(beta)*Delta_tau);

ISIexact=conv(alphan,betan);
ISIexactn=ISIexact/(sum(ISIexact)*Delta_tau);



%%

figure
subplot(411)
histogram(taus1,Nbins,'Normalization','pdf')
hold on
if P_reset==0, plot(tau1,rhon,'r'), end
axis([0 1 0 inf])
title('\rho(\tau_1,\infty)')

subplot(412)
histogram(RT,Nbins,'Normalization','pdf')
hold on
if P_reset==0, plot(tau1,alphan,'r'), end
axis([0 1 0 inf])
title('f_1(\tau_1,\infty)')

subplot(413)
histogram(APD,Nbins,'Normalization','pdf')
hold on
if P_reset==0, plot(tau2,betan,'r'), end
axis([0 1 0 inf])
title('f_2(\tau_2,\infty)')

subplot(414)
histogram(ISI,Nbins,'Normalization','pdf')
hold on
if P_reset==0, plot(tauISI,ISIexactn,'r'), end
axis([0 2 0 inf])
title('f_{ISI}(\tau_{ISI},\infty)')
