%%% Neutrino in Free Space
%% Constants and Units of Neutrino Ossilation in Free Space(Three Flavor)
%Direct is for appearance of oscillation, indirect is for disappearance
Energy = 0.5;              % GeV(Using Nuclear Source for energy, use normal 
%non-relativistic formula for energy if have a momentum(p) to use
EnergyIndex = 0;
DataPoints = 500; %% Used to determine how many data points to get
PositionData = zeros(1,DataPoints); % Used to store position data
PositionIndex = 0; % Used to increment data matrices
ProbabilityMatrix = zeros(3,3); % Matrix to store oscillation probabilities
OscillationProbabilityData = zeros(DataPoints);
SurvivalProbabilityData = zeros(DataPoints); % Used to store probability data
EnergyData = zeros(1,DataPoints);
dm21 = 7.49e-5;       % eV^2
dm31 = 2.51e-3;       % eV^2
dm32 = dm31 - dm21;   % eV^2
Del_m_jk = [0 0 0;
    7.49e-5 0 0;
    2.51e-3 dm31-dm21 0];

theta12 = 33.7*pi/180; %Solar Mixing Levels
theta13 = 8.5*pi/180;  %Reactor Mixing Levels
theta23 = 48.0*pi/180; %Atmospheric Mixing Levels

%Solar neutrinos
L_coh_solar=10^8; %cm
sigma_x_solar = 10^-7; %cm
%Reactor
L_coh_reactor=10^9; %cm
%Accelarator
L_coh_accelerator=10^21; %cm
%Supernova
L_coh_supernova=100; %cm

%serphant
eta = 0.49; %m_mue^2 / m_pi^2


%Shorthanded Trig Variables
c12=cos(theta12);
s12=sin(theta12);
c13=cos(theta13);
s13=sin(theta13);
c23=cos(theta23);
s23=sin(theta23);

deltaCP = 0;  %CP-Violating Phaseck
%Flavor Representation in Rows
electron=1;
muon=2;
taon=3;
%Flavors, makes it easier to show change from a muon neutrino change to
%electron neutrino
alpha=muon;
beta=electron;
gamma=taon;
%Base Vectors for the Neutrinos
nu_e=[1;0;0];
nu_m=[0;1;0];
nu_t=[0;0;1];

%Complex Exponential for CP Violation
exp_minusi_delta=exp(-1i*deltaCP);
exp_plusi_delta=exp(1i*deltaCP);

%PMNS Mixing Matrix
Ue1=c12*c13;
Ue2=s12*c13;
Ue3=s13*exp_minusi_delta;

Umu1=-s12*c23 - c12*s23*s13*exp_plusi_delta;
Umu2=c12*c23 - s12*s23*s13*exp_plusi_delta;
Umu3= s23*c13;

Utau1=s12*s23 - c12*c23*s13*exp_plusi_delta;
Utau2=-c12*s23 - s12*c23*s13*exp_plusi_delta;
Utau3=c23*c13;

%Final Unitary Matrix
UnitaryMatrix=[Ue1,Ue2,Ue3; Umu1,Umu2,Umu3; Utau1,Utau2,Utau3];

%U*
UnitaryInverse = inv(UnitaryMatrix);

function Output = NonDirectProbability(Initial,Length,Energy,dm21,dm31,U)
%Probability of particle being in same initial state (Non-Direct)
phi1=0;
phi2=1.267*dm21*Length/Energy;
phi3=1.267*dm31*Length/Energy;
phase_vector=[exp(-1i*phi1);exp(-1i*phi2);exp(-1i*phi3)];
PhaseDiagonal=diag(phase_vector);
FlavorEvolutionMatrix=U*PhaseDiagonal*U';%This is the Full Flavor Evolution Matrix
Output=FlavorEvolutionMatrix(Initial,Initial);
end 

%Direct Survival Probability
%function Output = OscillationProbability(ini,fin,Len,Erg,U,Uinv,x)
%   Ljk_osc = 4*pi*Energy_2 / abs(del_m_jk^2);
%   Output=(symsum((abs(U(fin,j))^2*abs(U(ini,j))^2),j,1,3))...
%   + (2*real(Uinv(ini,k)*U(ini,j)*U(fin,k)*Uinv(fin,j)*exp((i*2*pi)*(x/(Ljk_osc))))
%end

sum1 = 0;
sum2 = 0;

sigma_x = sigma_x_solar;
xdih = 1:10000:400000000;
Energy2=5; %keV
initial = beta;
final = gamma;
dih = 0;
legenddih=0
function DICK(initial,final,UnitaryMatrix,UnitaryInverse,Del_m_jk,Ljk_coh,eta,sigma_x)
    sum1 = 0;
    sum2 = 0;
    Energy2=5; %keV
    xdih = 1:10000:400000000;
    dih=0;
    for i = 1:40000
        x = i*10000;
        for j = 1:3
            sum1 = sum1 + (abs((UnitaryMatrix(final,j))))^2*(abs(UnitaryMatrix(initial,j)))^2; %Baseline
            for k = 1:3
                if k>j
                    Ljk_osc = 4*pi*Energy2 / abs(Del_m_jk(j,k)^2);
                    %Ljk_coh = (4*sqrt(2)*Energy2^2*sigma_x)/(abs(Del_m_jk(j,k)^2));
                    %Ljk_coh = L_coh_solar;
                    sum2 = sum2 + (2*real(UnitaryInverse(initial,k)*UnitaryMatrix(initial,j)*UnitaryMatrix(final,k)*UnitaryInverse(final,j)*exp(1i*2*pi*x/(Ljk_osc))))*exp(-((x/Ljk_coh)^2))*exp(-2*pi^2*eta^2*(sigma_x/Ljk_osc)^2);
                end

            end
        end
        dih(1,i) = sum1+sum2;
        sum1=0;
        sum2=0;
    end
    i=1:40000;
    plot(xdih(i),dih(i),LineWidth=2)
end
Ljk_coh=L_coh_solar;
DICK(alpha,beta,UnitaryMatrix,UnitaryInverse,Del_m_jk,Ljk_coh,eta,sigma_x)
hold on
DICK(beta,gamma,UnitaryMatrix,UnitaryInverse,Del_m_jk,Ljk_coh,eta,sigma_x)
DICK(alpha,gamma,UnitaryMatrix,UnitaryInverse,Del_m_jk,Ljk_coh,eta,sigma_x)
DICK(beta,alpha,UnitaryMatrix,UnitaryInverse,Del_m_jk,Ljk_coh,eta,sigma_x)
DICK(gamma,beta,UnitaryMatrix,UnitaryInverse,Del_m_jk,Ljk_coh,eta,sigma_x)
DICK(gamma,alpha,UnitaryMatrix,UnitaryInverse,Del_m_jk,Ljk_coh,eta,sigma_x)
legend('alphabeta','betagamma','alphagamma','betaalpha','gammabeta','gammaalpha')

%figure(1);
%i=1:40000;
%plot(xdih(i),dih(i))