% <<<<<<<<<<<<<<<<<<< PDE2103 Week 20 Lab Modulation and Demodulation with 2-dimensional constellation >>>>>>>>>>>>>>>>>>>
clc;
clear all;
close all;





% ******************* Digital/Binary input information ********************
x = input('Enter Digital Input Information in an array form with number of bits that can be divided by 6, e.g., [1 0 1 0], = ');   % Binary information as stream of bits (binary signal 0 or 1)
L = length(x);
while (mod(L,6))% This to make sure that the number of bits to the input can be divided by 6.
    disp('The number of bits should be even.');
    x = input('Enter Digital Input Information in an array form with number of bits, e.g., [1 0 1 0], = ');
    L = length(x); % update the length of x
end
Tb = 0.000001;   %Symbol period, data rate =1/Tb= 1Mbps
disp('Binary Input Information at Transmitter: ');
disp(x);








% ************* Represent input information with respect to time *************
nb = 100;   % Represented time duration per bit,converts bits into digital waveform
digit = [];%creates an empty array to store digit signal waveformi
for n1 = 1:1:L  %Loop from 1 to L → process each bit of x.
    if x(n1) == 1;
        sig = ones(1,nb);
    else x(n1) == 0;
        sig = zeros(1,nb);
    end
    digit = [digit sig];
end
t1=Tb/nb:Tb/nb:nb*L*(Tb/nb);   % Time period for plotting the digital input signal (not used in modulation math, only for display)
figure('Name','2-dimensional Constellation Modulation and Demodulation','NumberTitle','off');%opens a new figure window for all plots
subplot(4,1,1);%the figure has 4 rows, 1 column.
plot(t1,digit,'lineWidth',2.5);    % Plot the raw digital input bitstream over time
grid on;
axis([0 Tb*L -0.5 1.5]);
xlabel('Time(Sec)');             % X-axis = time
ylabel('Amplitude(Volts)');      % Y-axis = voltage level of bits (0 or 1)
title('Digital Input Signal at Tx');








%****************************Vector Encoder*******************************
R=20
alpha = 2 * pi / 64; % Define the angle for the 64-point constellation
X = zeros(2, 64); % Initialize the constellation matrix as 2x64

% assume R and alpha are defined, where alpha = 2*pi/64
X_0  = [R*cos(0*alpha);  R*sin(0*alpha)];
X_1  = [R*cos(1*alpha);  R*sin(1*alpha)];
X_2  = [R*cos(2*alpha);  R*sin(2*alpha)];
X_3  = [R*cos(3*alpha);  R*sin(3*alpha)];
X_4  = [R*cos(4*alpha);  R*sin(4*alpha)];
X_5  = [R*cos(5*alpha);  R*sin(5*alpha)];
X_6  = [R*cos(6*alpha);  R*sin(6*alpha)];
X_7  = [R*cos(7*alpha);  R*sin(7*alpha)];
X_8  = [R*cos(8*alpha);  R*sin(8*alpha)];
X_9  = [R*cos(9*alpha);  R*sin(9*alpha)];
X_10 = [R*cos(10*alpha); R*sin(10*alpha)];
X_11 = [R*cos(11*alpha); R*sin(11*alpha)];
X_12 = [R*cos(12*alpha); R*sin(12*alpha)];
X_13 = [R*cos(13*alpha); R*sin(13*alpha)];
X_14 = [R*cos(14*alpha); R*sin(14*alpha)];
X_15 = [R*cos(15*alpha); R*sin(15*alpha)];
X_16 = [R*cos(16*alpha); R*sin(16*alpha)];
X_17 = [R*cos(17*alpha); R*sin(17*alpha)];
X_18 = [R*cos(18*alpha); R*sin(18*alpha)];
X_19 = [R*cos(19*alpha); R*sin(19*alpha)];
X_20 = [R*cos(20*alpha); R*sin(20*alpha)];
X_21 = [R*cos(21*alpha); R*sin(21*alpha)];
X_22 = [R*cos(22*alpha); R*sin(22*alpha)];
X_23 = [R*cos(23*alpha); R*sin(23*alpha)];
X_24 = [R*cos(24*alpha); R*sin(24*alpha)];
X_25 = [R*cos(25*alpha); R*sin(25*alpha)];
X_26 = [R*cos(26*alpha); R*sin(26*alpha)];
X_27 = [R*cos(27*alpha); R*sin(27*alpha)];
X_28 = [R*cos(28*alpha); R*sin(28*alpha)];
X_29 = [R*cos(29*alpha); R*sin(29*alpha)];
X_30 = [R*cos(30*alpha); R*sin(30*alpha)];
X_31 = [R*cos(31*alpha); R*sin(31*alpha)];
X_32 = [R*cos(32*alpha); R*sin(32*alpha)];
X_33 = [R*cos(33*alpha); R*sin(33*alpha)];
X_34 = [R*cos(34*alpha); R*sin(34*alpha)];
X_35 = [R*cos(35*alpha); R*sin(35*alpha)];
X_36 = [R*cos(36*alpha); R*sin(36*alpha)];
X_37 = [R*cos(37*alpha); R*sin(37*alpha)];
X_38 = [R*cos(38*alpha); R*sin(38*alpha)];
X_39 = [R*cos(39*alpha); R*sin(39*alpha)];
X_40 = [R*cos(40*alpha); R*sin(40*alpha)];
X_41 = [R*cos(41*alpha); R*sin(41*alpha)];
X_42 = [R*cos(42*alpha); R*sin(42*alpha)];
X_43 = [R*cos(43*alpha); R*sin(43*alpha)];
X_44 = [R*cos(44*alpha); R*sin(44*alpha)];
X_45 = [R*cos(45*alpha); R*sin(45*alpha)];
X_46 = [R*cos(46*alpha); R*sin(46*alpha)];
X_47 = [R*cos(47*alpha); R*sin(47*alpha)];
X_48 = [R*cos(48*alpha); R*sin(48*alpha)];
X_49 = [R*cos(49*alpha); R*sin(49*alpha)];
X_50 = [R*cos(50*alpha); R*sin(50*alpha)];
X_51 = [R*cos(51*alpha); R*sin(51*alpha)];
X_52 = [R*cos(52*alpha); R*sin(52*alpha)];
X_53 = [R*cos(53*alpha); R*sin(53*alpha)];
X_54 = [R*cos(54*alpha); R*sin(54*alpha)];
X_55 = [R*cos(55*alpha); R*sin(55*alpha)];
X_56 = [R*cos(56*alpha); R*sin(56*alpha)];
X_57 = [R*cos(57*alpha); R*sin(57*alpha)];
X_58 = [R*cos(58*alpha); R*sin(58*alpha)];
X_59 = [R*cos(59*alpha); R*sin(59*alpha)];
X_60 = [R*cos(60*alpha); R*sin(60*alpha)];
X_61 = [R*cos(61*alpha); R*sin(61*alpha)];
X_62 = [R*cos(62*alpha); R*sin(62*alpha)];
X_63 = [R*cos(63*alpha); R*sin(63*alpha)];

% Splitting the 64 bit input into 6 Bits and mapping the grey code constellation%
C = [X_0, X_1, X_2, X_3, X_4, X_5, X_6, X_7, X_8, X_9, X_10, X_11, X_12, X_13, X_14, X_15, X_16, X_17, X_18, X_19, X_20, X_21, X_22, X_23, X_24, X_25, X_26, X_27, X_28, X_29, X_30, X_31, X_32, X_33, X_34, X_35, X_36, X_37, X_38, X_39, X_40, X_41, X_42, X_43, X_44, X_45, X_46, X_47, X_48, X_49, X_50, X_51, X_52, X_53, X_54, X_55, X_56, X_57, X_58, X_59, X_60, X_61, X_62, X_63];
% 0 - 16 %
for n2 = 1:6:L
    if x(n2:n2+5) == [0 0 0 0 0 0] % input bits are 0 0 0 0 0 0
        X(:,n2) = X_0(:,1); % Constellation mapping to X_0
    elseif x(n2:n2+5) == [0 0 0 0 0 1] % input bits 
        X(:,n2) = X_1(:,1);% Constellation mapping to X_1
    elseif x(n2:n2+5) == [0 0 0 0 1 1] % input bits 
        X(:,n2) = X_2(:,1);% Constellation mapping to X_2
    elseif x(n2:n2+5) == [0 0 0 0 1 0] % input bits 
        X(:,n2) = X_3(:,1);% Constellation mapping to X_3
    elseif x(n2:n2+5) == [0 0 0 1 1 0] % input bits 
        X(:,n2) = X_4(:,1);% Constellation mapping to X_4
    elseif x(n2:n2+5) == [0 0 0 1 1 1] % input bits 
        X(:,n2) = X_5(:,1);% Constellation mapping to X_5
    elseif x(n2:n2+5) == [0 0 0 1 0 1]% input bits 
        X(:,n2) = X_6(:,1);% Constellation mapping to X_6
    elseif x(n2:n2+5) == [0 0 0 1 0 0] % input bits ar
        X(:,n2) = X_7(:,1);% Constellation mapping to X_1
    elseif x(n2:n2+5) == [0 0 1 1 0 0] % input bits 
        X(:,n2) = X_8(:,1);% Constellation mapping to X_2
    elseif x(n2:n2+5) == [0 0 1 1 0 1] % input bits 
        X(:,n2) = X_9(:,1);% Constellation mapping to X_3
    elseif x(n2:n2+5) == [0 0 1 1 1 1] % input bits 
        X(:,n2) = X_10(:,1);% Constellation mapping to X_4
    elseif x(n2:n2+5) == [0 0 1 1 1 0] % input bits 
        X(:,n2) = X_11(:,1);% Constellation mapping to X_5
    elseif x(n2:n2+5) == [0 0 1 0 1 0]% input bits 
        X(:,n2) = X_12(:,1);% Constellation mapping to X_6
    elseif x(n2:n2+5) == [0 0 1 0 1 1] % input bits 
        X(:,n2) = X_13(:,1);% Constellation mapping to X_3
    elseif x(n2:n2+5) == [0 0 1 0 0 1] % input bits 
        X(:,n2) = X_14(:,1);% Constellation mapping to X_4
    elseif x(n2:n2+5) == [0 0 1 0 0 0] % input bits 
        X(:,n2) = X_15(:,1);% Constellation mapping to X_5
        % 16 - 31 %
    elseif x(n2:n2+5) == [0 1 1 0 0 0] % input bits 
        X(:,n2) = X_16(:,1); % Constellation mapping to X_16
    elseif x(n2:n2+5) == [0 1 1 0 0 1] % input bits 
        X(:,n2) = X_17(:,1);% Constellation mapping to X_17
    elseif x(n2:n2+5) == [0 1 1 0 1 1] % input bits 
        X(:,n2) = X_18(:,1);% Constellation mapping to X_18
    elseif x(n2:n2+5) == [0 1 1 0 1 0] % input bits 
        X(:,n2) = X_19(:,1); % Constellation mapping to X_19
    elseif x(n2:n2+5) == [0 1 1 1 1 0] % input bits 
        X(:,n2) = X_20(:,1); % Constellation mapping to X_20
    elseif x(n2:n2+5) == [0 1 1 1 1 1] % input bits 
        X(:,n2) = X_21(:,1); % Constellation mapping to X_21
    elseif x(n2:n2+5) == [0 1 1 1 0 1] % input bits 
        X(:,n2) = X_22(:,1); % Constellation mapping to X_22
    elseif x(n2:n2+5) == [0 1 1 1 0 0] % input bits
        X(:,n2) = X_23(:,1); % Constellation mapping to X_23
    elseif x(n2:n2+5) == [0 1 0 1 0 0] % input bits 
        X(:,n2) = X_24(:,1); % Constellation mapping to X_24
    elseif x(n2:n2+5) == [0 1 0 1 0 1] % input bits 
        X(:,n2) = X_25(:,1); % Constellation mapping to X_25
    elseif x(n2:n2+5) == [0 1 0 1 1 1] % input bits
        X(:,n2) = X_26(:,1); % Constellation mapping to X_26
    elseif x(n2:n2+5) == [0 1 0 1 1 0] % input bits 
        X(:,n2) = X_27(:,1); % Constellation mapping to X_27
    elseif x(n2:n2+5) == [0 1 0 0 1 0] % input bits
        X(:,n2) = X_28(:,1); % Constellation mapping to X_28
    elseif x(n2:n2+5) == [0 1 0 0 1 1] % input bits
        X(:,n2) = X_29(:,1); % Constellation mapping to X_29
    elseif x(n2:n2+5) == [0 1 0 0 0 1] % input bits 
        X(:,n2) = X_30(:,1); % Constellation mapping to X_30
    elseif x(n2:n2+5) == [0 1 0 0 0 0] % input bits 
        X(:,n2) = X_31(:,1); % Constellation mapping to X_31
        % 32 - 47 %​
    elseif x(n2:n2+5) == [1 1 0 0 0 0]% input bits
        X(:,n2) = X_32(:,1);% Constellation mapping to X_32
    elseif x(n2:n2+5) == [ 1 1 0 0 0 1]% input bits
        X(:,n2) = X_33(:,1);%Constellation mapping to X_33
    elseif x(n2:n2+5) == [ 1 1 0 0 1 1] % input bits
        X(:,n2) = X_34(:,1);%Constellation mapping to X_34
    elseif x(n2:n2+5) == [ 1 1 0 0 1 0]% input bits
        X(:,n2) = X_35(:,1);%Constellation mapping to X_35
    elseif x(n2:n2+5) == [ 1 1 0 1 1 0]% input bits
        X(:,n2) = X_36(:,1);%Constellation mapping to X_36
    elseif x(n2:n2+5) == [ 1 1 0 1 1 1]% input bits
        X(:,n2) = X_37(:,1);%Constellation mapping to X_37
    elseif x(n2:n2+5) == [ 1 1 0 1 0 1]% input bits
        X(:,n2) = X_38(:,1);%Constellation mapping to X_38
    elseif x(n2:n2+5) == [ 1 1 0 1 0 0]% input bits
        X(:,n2) = X_39(:,1);%Constellation mapping to X_39
    elseif x(n2:n2+5) == [ 1 1 1 1 0 0]% input bits
        X(:,n2) = X_40(:,1);%Constellation mapping to X_40
    elseif x(n2:n2+5) == [ 1 1 1 1 0 1]% input bits
        X(:,n2) = X_41(:,1);%Constellation mapping to X_41
    elseif x(n2:n2+5) == [ 1 1 1 1 1 1]%input bits
        X(:,n2) = X_42(:,1);%Constellation mapping to X_42
    elseif x(n2:n2+5) == [ 1 1 1 1 1 0]%input bits
        X(:,n2) = X_43(:,1);%Constellation mapping to X_43
    elseif x(n2:n2+5) == [ 1 1 1 0 1 0]%input bits
        X(:,n2) = X_44(:,1);%Constellation mapping to X_44
    elseif x(n2:n2+5) == [ 1 1 1 0 1 1]%input bits
        X(:,n2) = X_45(:,1);%Constellation mapping to X_45
    elseif x(n2:n2+5) == [ 1 1 1 0 0 1]%input bits
        X(:,n2) = X_46(:,1);%Constellation mapping to X_46
    elseif x(n2:n2+5) == [ 1 1 1 0 0 0]%input bits
        X(:,n2) = X_47(:,1);%Constellation mapping to X_47

        % 48 - 63 %
    elseif x(n2:n2+5) == [1 0 1 0 0 0] % input bits 
        X(:,n2) = X_48(:,1); % Constellation mapping to X_48
    elseif x(n2:n2+5) == [1 0 1 0 0 1] % input bits 
        X(:,n2) = X_49(:,1);% Constellation mapping to X_49
    elseif x(n2:n2+5) == [1 0 1 0 1 1] % input bits 
        X(:,n2) = X_50(:,1);% Constellation mapping to X_50
    elseif x(n2:n2+5) == [1 0 1 0 1 0] % input bits 
        X(:,n2) = X_51(:,1); % Constellation mapping to X_51
    elseif x(n2:n2+5) == [1 0 1 1 1 0] % input bits 
        X(:,n2) = X_52(:,1);% Constellation mapping to X_52
    elseif x(n2:n2+5) == [1 0 1 1 1 1] % input bits 
        X(:,n2) = X_53(:,1);% Constellation mapping to X_53
    elseif x(n2:n2+5) == [1 0 1 1 0 1] % input bits 
        X(:,n2) = X_54(:,1); % Constellation mapping to X_54
    elseif x(n2:n2+5) == [1 0 1 1 0 0] % input bits 
        X(:,n2) = X_55(:,1);% Constellation mapping to X_55
    elseif x(n2:n2+5) == [1 0 0 1 0 0] % input bits 
        X(:,n2) = X_56(:,1);% Constellation mapping to X_56
    elseif x(n2:n2+5) == [1 0 0 1 0 1] % input bits 
        X(:,n2) = X_57(:,1); % Constellation mapping to X_57
    elseif x(n2:n2+5) == [1 0 0 1 1 1] % input bits 
        X(:,n2) = X_58(:,1);% Constellation mapping to X_58
    elseif x(n2:n2+5) == [1 0 0 1 1 0] % input bits 
        X(:,n2) = X_59(:,1);% Constellation mapping to X_59
    elseif x(n2:n2+5) == [1 0 0 0 1 0] % input bits 
        X(:,n2) = X_60(:,1); % Constellation mapping to X_60
    elseif x(n2:n2+5) == [1 0 0 0 1 1] % input bits 
        X(:,n2) = X_61(:,1);% Constellation mapping to X_61
    elseif x(n2:n2+5) == [1 0 0 0 0 1] % input bits 
        X(:,n2) = X_62(:,1);% Constellation mapping to X_62
    elseif x(n2:n2+5) == [1 0 0 0 0 0] % input bits 
        X(:,n2) = X_63(:,1); % Constellation mapping to X_63
    end
end % Note that the first iteration n2=1. Since the step size is 6, the next
% iteration n2 will be 7. In the first iteration the first column of X will
% be assign to a constellation. After the first iteration, n2=7 hence the
% next value of X(:,7) will be assigned to another constellation. As the
% result, X(:,2) will be empty, The same for X(:,3),
% X(:,4).....
% The following lines of codes are to remove the empty columns in X
[r1 c1]=size(X);
for n3=c1:-1:1 % Note here we start from the last column of X and go back ward to the first column
    if (X(:,n3)==[0;0])% check if the colum n3 is empty
        X(:,n3)=[]; %removing the empty column of X
    end
end
[r2 c2]=size(X);% update the size of X after dropping empty column r1=2 and c2= N
disp('Symbol Vector at Transmitter: ');
disp(X);

X_copy = X





%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% **************************** 2-dimentional constellation Modulation ****************************
Ac = 1;      % Carrier amplitude for binary input
Fc = 10^7;   % Carrier frequency
t2 = Tb/nb:Tb/nb:Tb;          % Signal time,Time vector for a single symbol duration,Tb is the actual duration of one symbol (or bit) in real time,nb determines how many samples we use to represent 1 bit in our simulation,t2 constructs the time axis for one bit with nb samples.
phi_1 = cos(2*pi*Fc*t2);      % basis 1(I component carrier)
phi_2 = sin(2*pi*Fc*t2);      % basis 2(Q component carrier)
mod = [];
for k1 = 1:1:c2
    y = X(1,k1)*Ac*phi_1+X(2,k1)*Ac*phi_2;  % Combine I and Q components with the two basis carriers to form the transmitted signal
    mod=[mod y];      %phase of thewaveform change at every symbol interval
end
t3=Tb/nb:Tb/nb:Tb*c2;    % Time vector for the entire transmitted sequence
subplot(4,1,2);          %2nd plot in the 4-row figure (modulated signal)
plot(t3,mod);
xlabel('Time(Sec)');
ylabel('Amplitude(Volts)');
title('Modulated/Transmitted Symbol/Signal');







% ********************* Transmitted signal x ******************************
X = mod;
% ********************* Channel model h and w *****************************
h = 0.7;   % channel attenuation
[rn cn]=size(X); % matlab treats X as a 1xN row vector,eventhough x is a symbol matlab cares about its shape,and this was stored as a row
w = randn(1,cn);   % Additive White Gaussian Noise
% ********************* Received signal y *********************************
y = h.*X + w;   % AWGN channel
% ****** plot the received signal****************************
subplot(4,1,3);
plot(t3,y);
xlabel('Time(Sec)');
ylabel('Amplitude(Volts)');
title('Received Symbol/Signal');







%Tb = time duration of one symbol
% ***************************  Demodulation ***************************
s = length(t2);%number of samples in one symbol period
demod = [];
for n1 = s:s:length(y) %move through recieved signal y in steps of one symbol,each iteration process one symbol duration of y
    z1 = trapz(t2,phi_1.*y((n1-(s-1)):n1)); % intregation of the received signal with basis 1
    z2 = trapz(t2,phi_2.*y((n1-(s-1)):n1)); % intregation of the received signal with basis 1
    rz1 = ((2*z1/Tb))/Ac;% normalised by the Carrier amplitude to get the coordinate assosiated with basis 1,Divide by Tb = removes integration scaling
                         %Multiply by 2/Tb = undo basis energy
                         %Divide by Ac = undoes carrier amplitude, leaving only symbol coordinate
    rz2 = ((2*z2/Tb))/Ac;% normalised by the Carrier amplitude to get the coordinate assosiated with basis 2
    demod = [demod [rz1;rz2]]; % stack two coordinates to form 2 dimension vector of our demodulated vector symbol
end
disp('Demodulated Symbol Vector at Receiver: ');
disp(demod);







% ****************** ML detection***********************************
m = [];  %m is defined once at the top
[r4 c4]=size(demod);%  the size of demod r4=2 and c4= N

for l1=1:c4 % iterating through demod taking l1 demodulated vectors one by one
    m=[m  ML_detection(demod(:,l1),C)]; %in ML detection code computes the distance between the recieved vector and every constellation point
end

  
disp('Binary output Information at Receiver: ');
disp(m);   % m becomes the entire binary sequence, symbol after symbol.







% ********** Represent demodulated information with respect to time **********
digit_demod = [];
L=length(m);
for n1 = 1:L
    if (m(n1) == 1)
        sig_demod = ones(1,nb);% nb Represented time duration per bit,converts bits into digital waveform
    elseif (m(n1) ==0)
        sig_demod = zeros(1,nb);
    end
    digit_demod = [digit_demod sig_demod];
end
%%%%%%%%%%
t5=Tb/nb:Tb/nb:nb*L*(Tb/nb);   % Time period
%t5=Tb/nb:Tb/nb:nb*length(m)*(Tb/nb);   % Time period
subplot(4,1,4)
plot(t5,digit_demod,'LineWidth',2.5);grid on;
axis([0 Tb*L -1.5 1.5]);

%axis([0 Tb*L -0.5 1.5]);

xlabel('Time(Sec)');
ylabel('Amplitude(Volts)');
title('Digital output Signal at Rx');
%%%%%%%%%%%%%%%%










%************ Plot the transmitted vector constellation
%**********************
figure('Name','2-dimensional Constellation Modulation and Demodulation','NumberTitle','off');
Marker_Counter = 1;
Markers = {'+','o','*','x','v','d','^','s','>','<'};

subplot(2,1,1);
plot([-20 20],[0 0],'r-'); 
hold on;
plot([0 0],[-20 20],'r-'); 
hold on;
axis([-20 20 -20 20]);


for k = 1:c2
   plot(X_copy(1,k),X_copy(2,k),strcat(Markers{Marker_Counter}));
    if (Marker_Counter < length(Markers)) % assigning different mark for different points
        Marker_Counter=Marker_Counter+1;
    else
        Marker_Counter=1;
    end    
        hold on;
end

title('Transmitted constellation');
axis([-20 20 -20 20]);








%************ Plot the received vector constellation
%**********************
Marker_Counter=1;
Markers = {'+','o','*','x','v','d','^','s','>','<'};
subplot(2,1,2);
plot([-20 20],[0 0],'r-')
hold on;
plot([0 0],[-20 20],'r-')
hold on;
axis([-20 20 -20 20]);
for k=1:1:c2
    plot(demod(1,k),demod(2,k),strcat(Markers{Marker_Counter}));
    if (Marker_Counter < length(Markers)) % assigning different mark for different points
        Marker_Counter=Marker_Counter+1;
    else
        Marker_Counter=1;
    end
    hold on;
end
title('Received constellation');
axis([-20 20 -20 20]);
% ************************** End of the program ***************************