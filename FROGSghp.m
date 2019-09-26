% FROGS
% ver1.2 (181011edited)
%
% main
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
clear;
close all
% 変数の宣言
global l lcg0 lcgf lcgp lcp m0 mf mp0 I0 If Ip0 n
global Cd Cnalpha Cmq Vpara1 Vpara2 Hpara lLnchr
global WindModel dt Cdv Zr thrust tThrust g
global S SIMULATION
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Choose the type of simulation(弾道or減速)
%%% 3=Ballistic fall，4=Retarding fall
SIMULATION  = 4;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

FROGSparameters;        % parameterの読み込み
FROGSthrust;            % thrustデータの読み込み
%
GHP = zeros(16,17);
for Vtemp = 1:8
    Vwaz = 0.0+Vtemp*1.0;
for k = 1:17
    Waz = 22.5*(k-1)/180*pi;

[Ve1,Ve2,Ve3,Xe1,Xe2,Xe3,omg2,omg3,q1,q2,q3,q4]  = FROGSprset;
IV  = [Ve1,Ve2,Ve3,Xe1,Xe2,Xe3,omg2,omg3,q1,q2,q3,q4];
Ve  = [IV(1); IV(2); IV(3)];            %地球座標系における機体速度ベクトル
Xe  = [IV(4); IV(5); IV(6)];            %地球座標系における機体位置ベクトル
omg = [0; IV(7); IV(8)];                %角速度
q   = [IV(9); IV(10); IV(11); IV(12)];  %クォータニオン
i = 1;                                  %ステップ数
t = 0;                                  %時間
%%%
for i = 1:n                             % 1~nまで繰り返し計算
t = i*dt;                               % time [s]

% transformation matrix(earth frame --> body frame)(座標変換行列)(地上座標系を機体座標系へ)
Aeb=[1-2*((q(2)^2)+(q(3)^2)),...
    2*(q(1)*q(2)+q(3)*q(4)),2*(q(3)*q(1)-q(2)*q(4));
    2*(q(1)*q(2)-q(3)*q(4)),1-2*((q(1)^2)+(q(3)^2)),...
    2*(q(2)*q(3)+q(1)*q(4));
    2*(q(3)*q(1)+q(2)*q(4)),2*(q(2)*q(3)-q(1)*q(4)),...
    1-2*((q(1)^2)+(q(2)^2))];

if t<tThrust                                % 燃焼中の飛翔
    T       = thrust(i);                    % thrust [N]
    mdot    = (m0-mf)/tThrust;              % weight loss [kg/sec]
    m       = m0-mdot*t;                    % weight [kg]
    I       = I0-(I0-If)/tThrust*t;         % moment of inertia [kgm^2]
    lcg     = lcg0-(lcg0-lcgf)/tThrust*t;   % center of gravity [m]
else                                        % 燃焼終了後
    T       = 0;                            % thrust [N]
    mdot    = 0;                            % weight loss [kg/sec]
    m       = mf;                           % weight [kg]
    I       = If;                           % moment of inertia [kgm^2]
    lcg     = lcgf;                         % center of gravity [m]
end

rho = 1.23-(0.11e-3)*Xe(3);             % atmospheric density
                                        % Don't use this equation over 2km

% wind
switch(WindModel)
    case 1                                                %べき乗則
        if i==1                                           %ステップが1の時は高度0mなので仕方なく計測した風速そのまま利用
            Vw = [Vwaz*cos(Waz);Vwaz*sin(Waz);0];
        else
            Vw = [Vwaz*((Xe(3)/Zr)^(1/Cdv))*cos(Waz);...  %べき法則に基づいた計算
                  Vwaz*((Xe(3)/Zr)^(1/Cdv))*sin(Waz);0];
        end
    case 2
        Vw = [Vwaz*cos(Waz);Vwaz*sin(Waz);0];             %一様風速
end

% airspeed [m/s](対気速度)
Vab = Aeb*(Ve-Vw);
Va = norm(Vab);
% Ve:地球座標系(earth frame)における機体速度ベクトル
% Vw:地球座標系(earth frame)における風速ベクトル 

% angle of attack & angle of sideslip(迎角&横滑り角)
alpha = asin(Vab(3)/Va);
bet = asin(Vab(2)/Va);
 
% drag & normal force & side force
D = 0.5*rho*Va*Va*S*Cd;
Y = 0.5*rho*Va*Va*S*Cnalpha*sin(bet);
N = 0.5*rho*Va*Va*S*Cnalpha*sin(alpha);
% D:抗力(作用方向は軸力)
% Y:機体軸(Y軸方向)に作用する力→side force
% N:法線力
% force @body frame
Fe = [T-D;-Y;-N];

% acceleration/velocity/position @earth-fixed frame
switch(SIMULATION)
    case 3                         % 弾道落下
        Ae = Aeb'*(Fe./m)-[0;0;g];
        Ve = Ae*dt+Ve;
        Xe = Ve*dt+Xe;
    case 4                         % 減速落下：下段の燃焼終了後かつ速度が負でパラ展開，終端速度に即達する
        if (Ve(3)<=0)&&(t>tThrust)&&(Xe(3)>Hpara)
            Ae = [0;0;0];
            Ve = [Vw(1); Vw(2); -Vpara1];
            Xe = Ve*dt+Xe;
        elseif (Ve(3)<=0)&&(t>tThrust)&&(Xe(3)<Hpara)   % 2段パラ
            Ae = [0;0;0];
            Ve = [Vw(1); Vw(2); -Vpara2];
            Xe = Ve*dt+Xe;          
        else                                            % 上昇中
            Ae = Aeb'*(Fe./m)-[0;0;g];
            Ve = Ae*dt+Ve;
            Xe = Ve*dt+Xe;
        end
end

% coefficient
Kj = -(Ip0/mp0+(lcg-lcgp)^2-(l-lcg)^2)*mdot;
Ka = 0.5*rho*S*(Va^2)*(l^2)/2/Va*Cmq;
% Kj:ジェットダンピング係数
% Ka:空気力による減衰モーメント係数

% angular velocity
kt1 = -((lcp-lcg)*N+(Ka+Kj)*omg(2))/I;
kt2 = -((lcp-lcg)*N+(Ka+Kj)*(omg(2)+kt1*dt/2))/I;
kt3 = -((lcp-lcg)*N+(Ka+Kj)*(omg(2)+kt2*dt/2))/I;
kt4 = -((lcp-lcg)*N+(Ka+Kj)*(omg(2)+kt3*dt))/I;

kp1 = ((lcp-lcg)*Y+(Ka+Kj)*omg(3))/I;
kp2 = ((lcp-lcg)*Y+(Ka+Kj)*(omg(3)+kp1*dt/2))/I;
kp3 = ((lcp-lcg)*Y+(Ka+Kj)*(omg(3)+kp2*dt/2))/I;
kp4 = ((lcp-lcg)*Y+(Ka+Kj)*(omg(3)+kp3*dt))/I;
% kt:y軸周りの角速度
% kp:z軸周りの角速度

switch(SIMULATION)
    case 3
        if (Xe(3)<lLnchr)&&(t<tThrust)                    %ランチャに刺さってる時は回転しない
            omg = [0;0;0];
        else                                              %飛翔中ずっとルンゲクッタ
            omg = [0;
                   1/6*(kt1+2*kt2+2*kt3+kt4)*dt+omg(2);
                   1/6*(kp1+2*kp2+2*kp3+kp4)*dt+omg(3)];
        end
        
    case 4
        if (Xe(3)<lLnchr)&&(t<tThrust)                    %ランチャに刺さってる時は回転しない
            omg = [0;0;0];
        elseif (Ve(3)<=0)&&(t>tThrust)                    %減速落下中も回転しない
            omg = [0;0;0];
        else                                              %それ以外の時(上昇中)ルンゲクッタ
            omg = [0;
                   1/6*(kt1+2*kt2+2*kt3+kt4)*dt+omg(2);
                   1/6*(kp1+2*kp2+2*kp3+kp4)*dt+omg(3)];
        end
end

% quaternion
qmat = [0 omg(3) -omg(2) 0;
        -omg(3) 0 0 omg(2);
        omg(2) 0 0 omg(3);
        0 -omg(2) -omg(3) 0];
q = (0.5*qmat*q*dt+q)/norm(q);

% attitude angle(姿勢角)
the = asin(Aeb(1,3));
psi = atan(Aeb(1,2)/Aeb(1,1));

%
if (Xe(3)<0)&&(t>tThrust)
    break
end
end
%
GHP(2*Vtemp-1,k) = real(Xe(1));
GHP(2*Vtemp,k) = real(Xe(2));
end
% plot
plot(GHP(2*Vtemp-1,:),GHP(2*Vtemp,:),'-squareb');
hold on;
end