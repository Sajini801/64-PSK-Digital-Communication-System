function m=ML_detection(X,C) %m defined for for each demodulated vector
% C is the constellation which is C=[X_0,X_1,X_2...X_63]
% Mapping rule:
% bits 0 0 0 0 0 0 map to X_0; bits 0 0 0 0 0 1 map to X_1...
% This function returns bit m given the input of demodulated symbol X and
% the constellations
[ r c]=size(C);
d = zeros(1, c);           % preallocate distance array

for k1=1:c
   d(k1)=sqrt((X(1)-C(1,k1))^2+(X(2)-C(2,k1))^2);% calculate the distances between the demodulated 
   % symbol X and 64 constellations X_0, X_1, X_2...X_63. d(1) is the
   % distance between X and X_0. d(2) is the distance between X and X_1.
   % d(3) is the distance between X and X_2. d(4) is the distance between X
   % and X_3.
end

a = find(d == min(d),1); % Find the minimal distance

if a==1
    m = [0 0 0 0 0 0];
elseif a==2
    m = [0 0 0 0 0 1];
elseif a==3
    m = [0 0 0 0 1 1];
elseif a==4
    m = [0 0 0 0 1 0];
elseif a==5
    m = [0 0 0 1 1 0];
elseif a==6
    m = [0 0 0 1 1 1];
elseif a==7
    m = [0 0 0 1 0 1];
elseif a==8
    m = [0 0 0 1 0 0];
elseif a==9
    m = [0 0 1 1 0 0];
elseif a==10
    m = [0 0 1 1 0 1];
elseif a==11
    m = [0 0 1 1 1 1];
elseif a==12
    m = [0 0 1 1 1 0];
elseif a==13
    m = [0 0 1 0 1 0];
elseif a==14
    m = [0 0 1 0 1 1];
elseif a==15
    m = [0 0 1 0 0 1];
elseif a==16
    m = [0 0 1 0 0 0];
elseif a==17
    m = [0 1 1 0 0 0];
elseif a==18
    m = [0 1 1 0 0 1];
elseif a==19
    m = [0 1 1 0 1 1];
elseif a==20
    m = [0 1 1 0 1 0];
elseif a==21
    m = [0 1 1 1 1 0];
elseif a==22
    m = [0 1 1 1 1 1];
elseif a==23
    m = [0 1 1 1 0 1];
elseif a==24
    m = [0 1 1 1 0 0];
elseif a==25
    m = [0 1 0 1 0 0];
elseif a==26
    m = [0 1 0 1 0 1];
elseif a==27
    m = [0 1 0 1 1 1];
elseif a==28
    m = [0 1 0 1 1 0];
elseif a==29
    m = [0 1 0 0 1 0];
elseif a==30
    m = [0 1 0 0 1 1];
elseif a==31
    m = [0 1 0 0 0 1];
elseif a==32
    m = [0 1 0 0 0 0];
elseif a==33
    m = [1 1 0 0 0 0];
elseif a==34
    m = [1 1 0 0 0 1];
elseif a==35
    m = [1 1 0 0 1 1];
elseif a==36
    m = [1 1 0 0 1 0];
elseif a==37
    m = [1 1 0 1 1 0];
elseif a==38
    m = [1 1 0 1 1 1];
elseif a==39
    m = [1 1 0 1 0 1];
elseif a==40
    m = [1 1 0 1 0 0];
elseif a==41
    m = [1 1 1 1 0 0];
elseif a==42
    m = [1 1 1 1 0 1];
elseif a==43
    m = [1 1 1 1 1 1];
elseif a==44
    m = [1 1 1 1 1 0];
elseif a==45
    m = [1 1 1 0 1 0];
elseif a==46
    m = [1 1 1 0 1 1];
elseif a==47
    m = [1 1 1 0 0 1];
elseif a==48
    m = [1 1 1 0 0 0];
elseif a==49
    m = [1 0 1 0 0 0];
elseif a==50
    m = [1 0 1 0 0 1];
elseif a==51
    m = [1 0 1 0 1 1];
elseif a==52
    m = [1 0 1 0 1 0];
elseif a==53
    m = [1 0 1 1 1 0];
elseif a==54
    m = [1 0 1 1 1 1];
elseif a==55
    m = [1 0 1 1 0 1];
elseif a==56
    m = [1 0 1 1 0 0];
elseif a==57
    m = [1 0 0 1 0 0];
elseif a==58
    m = [1 0 0 1 0 1];
elseif a==59
    m = [1 0 0 1 1 1];
elseif a==60
    m = [1 0 0 1 1 0];
elseif a==61
    m = [1 0 0 0 1 0];
elseif a==62
    m = [1 0 0 0 1 1];
elseif a==63
    m = [1 0 0 0 0 1];
elseif a==64
    m = [1 0 0 0 0 0];
end
