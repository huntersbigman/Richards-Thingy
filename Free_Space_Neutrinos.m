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
deltaCP = 0;  %CP-Violating Phaseck

%% Flavor Representation
electron=1; 
muon=2;
taon=3; %Flavor Representation in Rows

alpha=muon;
beta=electron;
gamma=taon; %Flavor Represenation in greek letters for ease of use

%% Change in neutrino masses from flavor change in ***m^2***
dm21 = 7.49e-5;       % eV^2
dm31 = 2.51e-3;       % eV^2
dm32 = dm31 - dm21;   % eV^2
dm23 = -dm32; %eV^2; this is due to the change in mass from one flavor to another will be the opposite when going the opposite way; ultimately will be squared in calculations.
dm12 = -dm21; %eV^2
dm13 = -dm31; %eV^2
Del_m_jk_squared = [0    dm12 dm13;
            dm21 0    dm23;
            dm31 dm32 0   ]; %Matrix for Del_m
% Cite "Review of Particle Physics*" *Particle Data Group* page 1334

%% Mixing Levels
theta12 = 33.7*pi/180; %Solar Mixing Levels
theta13 = 8.5*pi/180;  %Reactor Mixing Levels
theta23 = 48.0*pi/180; %Atmospheric Mixing Levels

%% More Constants
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

%% Shorthanded Trig Variables
c12=cos(theta12);
s12=sin(theta12);
c13=cos(theta13);
s13=sin(theta13);
c23=cos(theta23);
s23=sin(theta23);

%% Base Vectors for the Neutrinos
nu_e=[1;0;0];
nu_m=[0;1;0];
nu_t=[0;0;1];

%% Complex Exponential for CP Violation
exp_minusi_delta=exp(-1i*deltaCP);
exp_plusi_delta=exp(1i*deltaCP);

%% PMNS Mixing Matrix
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
UnitaryMatrix=[Ue1  ,Ue2  ,Ue3  ; 
               Umu1 ,Umu2 ,Umu3 ; 
               Utau1,Utau2,Utau3];

UnitaryMatrixInverse = inv(UnitaryMatrix); %Final Unitary Matrix Inverse

%% Unsure if this function is correct
%function Output = NonDirectProbability(Initial,Length,Energy,dm21,dm31,U)
%Probability of particle being in same initial state (Non-Direct)
%phi1=0;
%phi2=1.267*dm21*Length/Energy;
%phi3=1.267*dm31*Length/Energy;
%phase_vector=[exp(-1i*phi1);exp(-1i*phi2);exp(-1i*phi3)];
%PhaseDiagonal=diag(phase_vector);
%FlavorEvolutionMatrix=U*PhaseDiagonal*U';%This is the Full Flavor Evolution Matrix
%Output=FlavorEvolutionMatrix(Initial,Initial);
%end 

%% Direct Survival Probability function
function DICK(initial,final,UnitaryMatrix,UnitaryInverse,Del_m_jk_squared,Ljk_coh,eta,sigma_x,Energy_input)
    sum1 = 0;
    sum2 = 0;
    xdih = 1:10000:400000000;
    dih=0;
    for i = 1:40000
        x = i*10000;
        for j = 1:3
            sum1 = sum1 + (abs((UnitaryMatrix(final,j))))^2*(abs(UnitaryMatrix(initial,j)))^2; %Baseline
            for k = 1:3
                if k>j
                    Ljk_osc = 4*pi*Energy_input / abs(Del_m_jk_squared(j,k));
                    %Using L_jk_coh as a constant
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

%% For the Plotting
sigma_x = sigma_x_solar;
Energy2=5; %keV
Ljk_coh=L_coh_solar;
%Fun(input,output,U,U*,Del_m_jk,L_jk_coh,eta,sigma_x,E)
DICK(alpha,beta,UnitaryMatrix,UnitaryMatrixInverse,Del_m_jk_squared,Ljk_coh,eta,sigma_x,Energy2)
hold on
DICK(beta,gamma,UnitaryMatrix,UnitaryMatrixInverse,Del_m_jk_squared,Ljk_coh,eta,sigma_x,Energy2)
DICK(alpha,gamma,UnitaryMatrix,UnitaryMatrixInverse,Del_m_jk_squared,Ljk_coh,eta,sigma_x,Energy2)
DICK(beta,alpha,UnitaryMatrix,UnitaryMatrixInverse,Del_m_jk_squared,Ljk_coh,eta,sigma_x,Energy2)
DICK(gamma,beta,UnitaryMatrix,UnitaryMatrixInverse,Del_m_jk_squared,Ljk_coh,eta,sigma_x,Energy2)
DICK(gamma,alpha,UnitaryMatrix,UnitaryMatrixInverse,Del_m_jk_squared,Ljk_coh,eta,sigma_x,Energy2)
legend('alphabeta','betagamma','alphagamma','betaalpha','gammabeta','gammaalpha')
