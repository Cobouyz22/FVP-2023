%% Projct 1 
%Cessena 182: Appendix B.1
    % Assuming Cruising @ Steady-State
    clear; close all; clc     
g = 32.2;

mass=2650;
b = 36; %Wing span [ft]
S = 174; % Wing Surface Area [ft^2]  and given in the tables and 2.1 Problem
AR = (b^2)/S;
M = 0.201; 
Vp1 = 220.1;
q_bar = 49.6;
  alpha_1 = 0;
alpha = .0384;
U1 = Vp1*cos(alpha_1);
W1 = Vp1*sin(alpha_1);
V1 =0;
Phi_1 = 0;

c_T = 3.85; %ft
c_R = 5.4; %ft

t_ratio = c_T/c_R;
c_bar = 4.9;

cl_i_H = 0;

Theta1 = alpha_1;
%% Delta E
Ts=.01;
t_final = 200;

% Code for Input [Copied from Test]
timeVec = [0:Ts:t_final]';
delE = zeros(length(timeVec),1);
iH = zeros(length(timeVec),1);

delE_doublet_up_Idx = find(timeVec >=1 & timeVec <= 3); %First doublet up
delE_doublet_dn_Idx = find(timeVec >= 7 & timeVec <=9); % Second doublet down

doubletMag = 2;

delE(delE_doublet_up_Idx,1) = deg2rad(doubletMag);
delE(delE_doublet_dn_Idx,1)= deg2rad(-doubletMag);

delE_vec_horzcat = horzcat(timeVec,delE);
iH_vec = horzcat(timeVec,iH);

%% Delta A

delA = zeros(length(timeVec),1);

delA_doublet_up_Idx = find(timeVec >=1 & timeVec <=3);
delA_doublet_down_Idx = find(timeVec >= 7 & timeVec <=9);

delA(delA_doublet_up_Idx,1)= deg2rad(doubletMag);
delA(delA_doublet_down_Idx,1) =deg2rad(-doubletMag);

delA_cev_horzcat = horzcat(timeVec,delA);

%% Delta R

delR =  zeros(length(timeVec),1);

delR_doublet_up_Idx = find(timeVec >=1 & timeVec <=3);
delR_doublet_down_Idx = find(timeVec >= 7 & timeVec <=9);

delR(delA_doublet_up_Idx,1)= deg2rad(doubletMag);
delR(delA_doublet_down_Idx,1) =deg2rad(-doubletMag);

delR_cev_horzcat = horzcat(timeVec,delA);


%%
% Steady State and Stability Derivatives @ Cruise 
Ixx_B = 948; %slug ft^2
Iyy_B = 1346; %slug ft^2
Izz_B = 1967; %slug ft^2
Ixz_B = 0; %slug ft^2


% Steady-State
CL_1 = .307;
CD_1 = .032;
cm1 = 0;

c_T_X1 = .032;
c_T_ZU =0;
c_T_xalpha = 0;
c_T_Z1 = 0;
c_T_Zalpha = 0;
Cm_Tu = 0;
Cm_T1 =0;
Cm_talpha =0;

% Stability Derivates
cD0 = .027;
cDu = 0;
cDalpha = 0.121;
c_TXu = -.096;
cL0 = 0.307;
cLu = 0;
cl_alpha = 4.41;
cl_alpha_dot = 1.7;
cl_q= 3.9;
cm_0= .04;
cm_u= 0;
cm_alpha= -0.613;
cm_alpha_dot = -7.27;
cm_q = -12.4;
cm_TU = 0;
cm_T1 = 0;
% Control Derivates
cD_delta_E = 0;
cl_delta_E = 0.43;
cm_delta_E = -1.122;

%Stability Derivatives (2nd Page)
Cl_Beta=-.0923;
cl_p = -.484;
cl_r = .0798;
C_YBeta = -.393;
c_Y_p = -.075;
c_Y_r = .214;
C_nBeta =.0587;
C_n_TBeta = 0;
C_n_p = -.0278;
C_n_r = -.0937;

C_yT_Beta = 0;
C_lT_Beta =0;
C_nT_Beta = 0;



cl_deltaA = .229;
cl_deltaR = .0147;
c_Y_deltaA = 0;
C_Y_deltaR = 0.187;
c_n_deltaA = -.0216;
C_n_deltaR = -.0645;

t_PID = 2;
Kp = .0010; % This adjusts response time & aggresiveness (higher = faster, but faster causes more overshoot). Adjust this first.
Ki =  0.000000; % This adjusts the final stopping value after damping. Adjust this last if needed (mine did not need it).
Kd = 0.030; %This adjusts oscillations (higher = greater frequency). Adjust this second.


%% Now Run Code
sim('Project1_2021a.slx')
load('Proj2PID2.mat')
load('PIDTime.mat')

%%  Plots
% % Graph A
%  figure; sgtitle('Elevator Doublet')
%  subplot(411); plot(timeVec,Symoutput_u);
%  hold on
%  plot(out.time,out.u)
%  title('U Vector')
%  xlabel('Time [s]')
%  ylabel('Velcoity [ft/sec]')
%  subplot(412); plot(timeVec,Symoutput_alpha); hold on 
%  plot(out.time,out.alpha)
%  title('Alpha')
%  xlabel('Time [s]')
%  ylabel('radians')
%  subplot(413); plot(timeVec,symoutput_theta); hold on
%  plot(out.time, out.theta)
%  title('Theta')
%  xlabel('Time [s]')
%  ylabel('radians')
%  subplot(414); plot(delE_vec_horzcat(:,2)); 
%  title('Input')
%  xlabel('Time [s]')
%  ylabel('Magnitude')
% 
% %Graph B
% figure; sgtitle('Aileron Doublet')
% subplot(411); plot(timeVec,Symoutput_Beta); hold on
% plot(out.time,out.Beta)
% title('Beta')
% xlabel('time [s]')
% ylabel('radians')
% subplot(412); plot(timeVec,Symoutput_Phi); hold on
% plot(out.time,out.Phi)
% title('Phi')
% xlabel('time [s]')
% ylabel('radians')
% subplot(413); plot(timeVec,Symoutput_Psi); hold on
% plot(out.time,out.Psi)
% title('Psi')
% xlabel('time [s]')
% ylabel('radians')
% subplot(414); plot(timeVec,Symoutput_delA); 
% title('Input')
% xlabel('time [s]')
% ylabel('Magnitude')
% 
% % %Graph C
% figure; sgtitle('Rudder Doublet')
% subplot(411); plot(timeVec,Symoutput_Beta); hold on
% plot(out.time,out.Beta)
% title('Beta')
% xlabel('time [s]')
% ylabel('radians')
% subplot(412); plot(timeVec,Symoutput_Phi); hold on
%  plot(out.time,out.Phi)
% title('Phi')
% xlabel('time [s]')
% ylabel('radians')
% subplot(413); plot(timeVec,Symoutput_Psi); hold on
% plot(out.time,out.Psi)
% title('Psi')
% xlabel('time [s]')
% ylabel('radians')
% subplot(414); plot(timeVec,Symoutput_delR)
% title('Input')
% xlabel('time [s]')
% ylabel('Magnitude')

% PID Controller Graph
hold on 
plot(timeVec,Symoutput_z); 
plot(out.time,out.altitude)
title('PID Controller')
xlabel('Time [s]')
ylabel('Altitude [ft]')
hold on
grid on
y_Woop = 5500;
yline(y_Woop)
legend('Project1', 'Project2','Target Line', 'Location', 'Southeast')



%% Observations & Responses [Section D]

%[Elevator Doublet only on] - Graph in Folder
% i. For the elevator doublet, when the positive doublet hits, the u vector
% spikes up and then spikes up again as the doublet ends. Once the doublet
% goes back to zero the u vectpore also starts to lead back to zero. For
% the Alpha and Theta, as the doublet increases thier amounts decrease in
% about the same fashion as the increase. They also mirror this on the
% decline as well and then taper off toward zero. I have been notified that
% the Theta is not correct and that the alpha andf Theta shouldn't be
% similar in this manner due to the meeting. Will be updating as we
% continue onto project 2 but understand that this output is incorrect. I
% am not too worried about as this as my lateral equations seem to be
% correct so I will swing back to thois before the next project and ensure
% this is correct. 


%[Aileron Doublet only on]-Graph in folder
%ii. For the aileron doubet, the Beta spiked in a positive direction for
%when the doublet would go positive. Then the it would drop down during the
% negative mangitude doublet. It would then taper out to zero once the doublet was completed. 
% For the Phi, it shaprly grow and then steady itself for a moemtn before
% coming back down to zero. It would then come down and staret to stabilize
% toward 0. For Psi, It sharpley rose over the entiretty of the doiblet
% then slowlyt stabilized toward zero. That would make sense due to the Phi
% being part of the lateral plane in which the distrubance of the angle in
% the z direction would change quickly then slowly taper back to zero in a
% Cessena. 

%[Rudder Doublet only on]-Graph in folder
%iii. For the rudder doublet, the Beta increases during the positive
%section of the doublet then decrease during the negative section of teh
%doublet while also oscilatting toward zero after the doublet was
%completed. For Phi, it almost makes a big "hump" over the eintreity of the
% doublet. However at the end it slightly osicalltes and comes to zero. For
% Psi, it goes negative during teh psoitve part of the doublet and vice
% versa. It then slowly makes it way to stability of 0. 

% TO NOTE: I noticed that for the Rudder and Aileron that they had very
% similar responses when activated. This makes sense both the Rudder and
% Aileron act in the lateral axis so with small perturbations = they should
% have similar looks and magntiudes. 

% PID Controller
% For my PID Controller, due to the nature of the problem, I went ahead and
% did a theoretical version of a climb from its regular cruise altitude. I
% did this due to my answers for derivatives were at cruise.
% I went from a 5000 cruise to a 8300 climb. I used the PID controller to try
% and attempt this in an efficent but non-sporty manner. Knowing that the
% performance of a Cessena has a theroetical rate of climb of 16.3 ft/sec from
% online sources, I wanted to ensure that I stayed below that in my climb
% to ensure its possibility. Howeevr, I understand that over 8000 is quite
% high for a Cessena but for the purposes of the PID controller I believe
% it is sufficent. Also the flight envelope for the Cesserna is at about 14000 feet. The picture given in the folder shows the output of the
% tuned PID. ALthough nto perfect in its climb  and has a taper, this was
% the best I could get with my first time to tuning. However it seems to be
% a decent climb all things considered. 
% 
% I spent time after trying to improve with some success. 
% I went ahead and edited the Kp values and was
% able to achieve a better tuned climb! Picture and given values are
% updated in the code with also the better picture in the folder as well.