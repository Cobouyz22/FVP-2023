%% Project 2
% This will be inputed to the Sim to see the reactions of the 
%% Delta E
Ts=.01;
t_final = 300;

% Code for Input [Copied from Test]
timeVec = [0:Ts:t_final]';
delE = zeros(length(timeVec),1);
iH = zeros(length(timeVec),1);

delE_doublet_up_Idx = find(timeVec >=1 & timeVec <= 3); %First doublet up
delE_doublet_dn_Idx = find(timeVec >= 7 & timeVec <=9); % Second doublet down

doubletMag = 0;

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

%% GAIN TUNING

%Kp =
%Ki =
%kd = 

%% Run Sim From Here
%run('')

%% Plotting
% % Graph A
% figure; sgtitle('Elevator Doublet')
% subplot(411); plot(timeVec,Symoutput_u);
% title('U Vector')
% xlabel('Time [s]')
% ylabel('Velcoity [ft/sec]')
% subplot(412); plot(timeVec,Symoutput_alpha);
% title('Alpha')
% xlabel('Time [s]')
% ylabel('radians')
% subplot(413); plot(timeVec,symoutput_theta);
% title('Theta')
% xlabel('Time [s]')
% ylabel('radians')
% subplot(414); plot(delE_vec_horzcat(:,2));
% title('Input')
% xlabel('Time [s]')
% ylabel('Magnitude')
% 
% %Graph B
% figure; sgtitle('Aileron Doublet')
% subplot(411); plot(timeVec,Symoutput_Beta)
% title('Beta')
% xlabel('time [s]')
% ylabel('radians')
% subplot(412); plot(timeVec,Symoutput_Phi)
% title('Phi')
% xlabel('time [s]')
% ylabel('radians')
% subplot(413); plot(timeVec,Symoutput_Psi)
% title('Psi')
% xlabel('time [s]')
% ylabel('radians')
% subplot(414); plot(timeVec,Symoutput_delA)
% title('Input')
% xlabel('time [s]')
% ylabel('Magnitude')
% 
% %Graph C
% figure; sgtitle('Rudder Doublet')
% subplot(411); plot(timeVec,Symoutput_Beta)
% title('Beta')
% xlabel('time [s]')
% ylabel('radians')
% subplot(412); plot(timeVec,Symoutput_Phi)
% title('Phi')
% xlabel('time [s]')
% ylabel('radians')
% subplot(413); plot(timeVec,Symoutput_Psi)
% title('Psi')
% xlabel('time [s]')
% ylabel('radians')
% subplot(414); plot(timeVec,Symoutput_delR)
% title('Input')
% xlabel('time [s]')
% ylabel('Magnitude')

% PID Controller Graph
%plot(timeVec,Symoutput_z)
%title('PID Controller')
%xlabel('Time [s]')
%ylabel('Altitude [ft]')
%hold on
%y_Woop = 5500;
%yline(y_Woop)

