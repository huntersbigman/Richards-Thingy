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
Utau3=c23*c13; %cite "Earth-Density Effects in LBL Experiments: A Comprehensive Review of Theory, Observations, and Future Directions" T. Pandit, B. S. Koranga

%Final Unitary Matrix
UnitaryMatrix=[Ue1  ,Ue2  ,Ue3  ; 
               Umu1 ,Umu2 ,Umu3 ; 
               Utau1,Utau2,Utau3];

UnitaryMatrixInverse = inv(UnitaryMatrix); %Final Unitary Matrix Inverse

%% Direct Survival Probability function
function [Prob,Dist] = DICK(initial,final,UnitaryMatrix,UnitaryInverse,Del_m_jk_squared,Ljk_coh,eta,sigma_x,Energy_input,sizing,xlim)
    sum1 = 0;
    sum2 = 0;
    xdih = 1:(xlim/sizing):xlim;
    dih=0;
    for i = 1:sizing
        x = i*(xlim/sizing);
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
    Prob = dih;
    Dist = xdih;
    i=1:sizing;
    plot(xdih(i),dih(i),LineWidth=2)
end

%cite "INTRODUCTION TO NEUTRINO PHYSICS" P. Lipari 5.5 Equation (96)

%% For the Plotting
sigma_x = sigma_x_solar;
Energy2=5; %keV

Ljk_coh=L_coh_reactor;

samplesize = 100000; %to achieve sample spacing in function
MaxDist =    10e8; %Distance to look up until

figure(2)
[P12,x12] = DICK(alpha,beta,UnitaryMatrix,UnitaryMatrixInverse,Del_m_jk_squared,Ljk_coh,eta,sigma_x,Energy2,samplesize,MaxDist);
hold on
[P23,x23] = DICK(beta,gamma,UnitaryMatrix,UnitaryMatrixInverse,Del_m_jk_squared,Ljk_coh,eta,sigma_x,Energy2,samplesize,MaxDist);
[P13,x13] = DICK(alpha,gamma,UnitaryMatrix,UnitaryMatrixInverse,Del_m_jk_squared,Ljk_coh,eta,sigma_x,Energy2,samplesize,MaxDist);
[P21,x21] = DICK(beta,alpha,UnitaryMatrix,UnitaryMatrixInverse,Del_m_jk_squared,Ljk_coh,eta,sigma_x,Energy2,samplesize,MaxDist);
[P32,x32] = DICK(gamma,beta,UnitaryMatrix,UnitaryMatrixInverse,Del_m_jk_squared,Ljk_coh,eta,sigma_x,Energy2,samplesize,MaxDist);
[P31,x31] = DICK(gamma,alpha,UnitaryMatrix,UnitaryMatrixInverse,Del_m_jk_squared,Ljk_coh,eta,sigma_x,Energy2,samplesize,MaxDist);
legend('alphabeta','betagamma','alphagamma','betaalpha','gammabeta','gammaalpha')

figure(399);

subplot(2,3,1)
plot(x12,P12,LineWidth=2,Color='cyan')
xlabel("Distance (km)")
ylabel("Probability")
title("Muon to electron")

subplot(2,3,2)
plot(x13,P13,LineWidth=2,Color='blue')
xlabel("Distance (km)")
ylabel("Probability")
title("Muon to Taon")

subplot(2,3,3)
plot(x23,P23,LineWidth=2,Color='red')
xlabel("Distance (km)")
ylabel("Probability")
title("Electron to Taon")

subplot(2,3,4)
plot(x21,P21,LineWidth=2,Color='Green')
xlabel("Distance (km)")
ylabel("Probability")
title("Electron to Muon")

subplot(2,3,5)
plot(x31,P31,LineWidth=2,Color='magenta')
xlabel("Distance (km)")
ylabel("Probability")
title("Taon to Muon")

subplot(2,3,6)
plot(x32,P32,LineWidth=2,Color='yellow')
xlabel("Distance (km)")
ylabel("Probability")
title("Taon to Electron")

%% Plotting 2: Electric Boogaloo
figure(401)
plot(x12,P12,LineWidth=2,Color='cyan')
xlabel("Distance (km)")
ylabel("Probability")
title("Muon to electron")

figure(402)
plot(x13,P13,LineWidth=2,Color='blue')
xlabel("Distance (km)")
ylabel("Probability")
title("Muon to Taon")

figure(403)
plot(x23,P23,LineWidth=2,Color='red')
xlabel("Distance (km)")
ylabel("Probability")
title("Electron to Taon")

figure(404)
plot(x21,P21,LineWidth=2,Color='Green')
xlabel("Distance (km)")
ylabel("Probability")
title("Electron to Muon")

figure(405)
plot(x31,P31,LineWidth=2,Color='magenta')
xlabel("Distance (km)")
ylabel("Probability")
title("Taon to Muon")

figure(406)
plot(x32,P32,LineWidth=2,Color='yellow')
xlabel("Distance (km)")
ylabel("Probability")
title("Taon to Electron")
