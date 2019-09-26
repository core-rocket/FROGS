% FROGS
% ver1.3 (181117edited)
%
% main
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
clear;
close all
% �ϐ��̐錾
global l lcg0 lcgf lcgp lcp m0 mf mp0 I0 If Ip0 n
global Cd Cnalpha Cmq Vpara1 Vpara2 Hpara lLnchr
global WindModel dt Cdv Zr thrust tThrust g
global Vwaz Waz S SIMULATION
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Choose the type of simulation(�e��or����)
%%% 1=Ballistic fall 2=Retarding fall
SIMULATION  = 1;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

FROGSparameters;        % parameter�̓ǂݍ���
FROGSthrust;            % thrust�f�[�^�̓ǂݍ���
[Ve1,Ve2,Ve3,Xe1,Xe2,Xe3,omg2,omg3,q1,q2,q3,q4]  = FROGSprset;
IV  = [Ve1,Ve2,Ve3,Xe1,Xe2,Xe3,omg2,omg3,q1,q2,q3,q4];
Ve  = [IV(1); IV(2); IV(3)];            %�n�����W�n�ɂ�����@�̑��x�x�N�g��
Xe  = [IV(4); IV(5); IV(6)];            %�n�����W�n�ɂ�����@�̈ʒu�x�N�g��
omg = [0; IV(7); IV(8)];                %�p���x
q   = [IV(9); IV(10); IV(11); IV(12)];  %�N�H�[�^�j�I��
i = 1;                                  %�X�e�b�v��
t = 0;                                  %����
%%%
for i = 1:n                             % 1~n�܂ŌJ��Ԃ��v�Z
t = i*dt;                               % time [s]

% transformation matrix(earth frame --> body frame)(���W�ϊ��s��)(�n����W�n���@�̍��W�n��)
Aeb=[1-2*((q(2)^2)+(q(3)^2)),...
    2*(q(1)*q(2)+q(3)*q(4)),2*(q(3)*q(1)-q(2)*q(4));
    2*(q(1)*q(2)-q(3)*q(4)),1-2*((q(1)^2)+(q(3)^2)),...
    2*(q(2)*q(3)+q(1)*q(4));
    2*(q(3)*q(1)+q(2)*q(4)),2*(q(2)*q(3)-q(1)*q(4)),...
    1-2*((q(1)^2)+(q(2)^2))];

if t<tThrust                                % �R�Ē��̔���
    T       = thrust(i);                    % thrust [N]
    mdot    = (m0-mf)/tThrust;              % weight loss [kg/sec]
    m       = m0-mdot*t;                    % weight [kg]
    I       = I0-(I0-If)/tThrust*t;         % moment of inertia [kgm^2]
    lcg     = lcg0-(lcg0-lcgf)/tThrust*t;   % center of gravity [m]
else                                        % �R�ďI����
    T       = 0;                            % thrust [N]
    mdot    = 0;                            % weight loss [kg/sec]
    m       = mf;                           % weight [kg]
    I       = If;                           % moment of inertia [kgm^2]
    lcg     = lcgf;                         % center of gravity [m]
end

rho = 1.23-(0.11e-3)*Xe(3);                 % atmospheric density
                                            % Don't use this equation over 2km
                                            
% wind
switch(WindModel)
    case 1                                                %�ׂ��摥
        if i==1                                           %�X�e�b�v��1�̎��͍��x0m�Ȃ̂Ŏd���Ȃ��v�������������̂܂ܗ��p
            Vw = [Vwaz*cos(Waz);Vwaz*sin(Waz);0];
        else
            Vw = [Vwaz*((Xe(3)/Zr)^(1/Cdv))*cos(Waz);...  %�ׂ��@���Ɋ�Â����v�Z
                  Vwaz*((Xe(3)/Zr)^(1/Cdv))*sin(Waz);0];
        end
    case 2
        Vw = [Vwaz*cos(Waz);Vwaz*sin(Waz);0];             %��l����
end

% airspeed [m/s](�΋C���x)
Vab = Aeb*(Ve-Vw);
Va = norm(Vab);
% Ve:�n�����W�n(earth frame)�ɂ�����@�̑��x�x�N�g��
% Vw:�n�����W�n(earth frame)�ɂ����镗���x�N�g�� 

% angle of attack & angle of sideslip(�}�p&������p)
alpha = asin(Vab(3)/Va);
bet = asin(Vab(2)/Va);

% drag & normal force & side force
D = 0.5*rho*Va*Va*S*Cd;
Y = 0.5*rho*Va*Va*S*Cnalpha*sin(bet);
N = 0.5*rho*Va*Va*S*Cnalpha*sin(alpha);
% D:�R��(��p�����͎���)
% Y:�@�̎�(Y������)�ɍ�p����́�side force
% N:�@����

% force @body frame
Fe = [T-D;-Y;-N];

% acceleration/velocity/position @earth-fixed frame
switch(SIMULATION)
    case 1                         % �e������
        Ae = Aeb'*(Fe./m)-[0;0;g];
        Ve = Ae*dt+Ve;
        Xe = Ve*dt+Xe;
    case 2                         % ���������F���i�̔R�ďI���ォ���x�����Ńp���W�J�C�I�[���x�ɑ��B����
        if (Ve(3)<=0)&&(t>tThrust)&&(Xe(3)>Hpara)
            Ae = [0;0;0];
            Ve = [Vw(1); Vw(2); -Vpara1*tanh(g*(t-tmax*dt)/Vpara1)];
            Xe = Ve*dt+Xe;
        elseif (Ve(3)<=0)&&(t>tThrust)&&(Xe(3)<Hpara)   % 2�i�p��
            Ae = [0;0;0];
            Ve = [Vw(1); Vw(2); -Vpara1+(Vpara1-Vpara2)*tanh(g*(t-tpara(1,2))/(Vpara1-Vpara2))];
            Xe = Ve*dt+Xe;          
        else                                            % �㏸��
            Ae = Aeb'*(Fe./m)-[0;0;g];
            Ve = Ae*dt+Ve;
            Xe = Ve*dt+Xe;
        end
end

% coefficient
Kj = -(Ip0/mp0+(lcg-lcgp)^2-(l-lcg)^2)*mdot;
Ka = 0.5*rho*S*(Va^2)*(l^2)/2/Va*Cmq;
% Kj:�W�F�b�g�_���s���O�W��
% Ka:��C�͂ɂ�錸�����[�����g�W��

% angular velocity
kt1 = -((lcp-lcg)*N+(Ka+Kj)*omg(2))/I;
kt2 = -((lcp-lcg)*N+(Ka+Kj)*(omg(2)+kt1*dt/2))/I;
kt3 = -((lcp-lcg)*N+(Ka+Kj)*(omg(2)+kt2*dt/2))/I;
kt4 = -((lcp-lcg)*N+(Ka+Kj)*(omg(2)+kt3*dt))/I;

kp1 = ((lcp-lcg)*Y+(Ka+Kj)*omg(3))/I;
kp2 = ((lcp-lcg)*Y+(Ka+Kj)*(omg(3)+kp1*dt/2))/I;
kp3 = ((lcp-lcg)*Y+(Ka+Kj)*(omg(3)+kp2*dt/2))/I;
kp4 = ((lcp-lcg)*Y+(Ka+Kj)*(omg(3)+kp3*dt))/I;
% kt:y������̊p���x
% kp:z������̊p���x

switch(SIMULATION)
    case 1
        if (Xe(3)<lLnchr)&&(t<tThrust)                    %�����`���Ɏh�����Ă鎞�͉�]���Ȃ�
            omg = [0;0;0];
        else                                              %���Ē������ƃ����Q�N�b�^
            omg = [0;
                   1/6*(kt1+2*kt2+2*kt3+kt4)*dt+omg(2);
                   1/6*(kp1+2*kp2+2*kp3+kp4)*dt+omg(3)];
        end
        
    case 2
        if (Xe(3)<lLnchr)&&(t<tThrust)                    %�����`���Ɏh�����Ă鎞�͉�]���Ȃ�
            omg = [0;0;0];
        elseif (Ve(3)<=0)&&(t>tThrust)                    %��������������]���Ȃ�
            omg = [0;0;0];
        else                                              %����ȊO�̎�(�㏸��)�����Q�N�b�^
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

% attitude angle(�p���p)
the = asin(Aeb(1,3));
psi = atan(Aeb(1,2)/Aeb(1,1));

% log
log_t(1,i)     = t;
log_T(1,i)     = T;
log_m(1,i)     = m;
log_I(1,i)     = I;
log_lcg(1,i)   = lcg;
log_rho(1,i)   = rho;
log_Vw(:,i)    = Vw;
log_Vab(:,i)   = Vab;
log_Va(1,i)    = Va;
log_alpha(1,i) = alpha*180/pi;
log_bet(1,i)   = bet*180/pi;
log_D(1,i)     = D;
log_Y(1,i)     = Y;
log_N(1,i)     = N;
log_Fe(:,i)    = Fe;
log_Ae(:,i)    = Ae;
log_Ve(:,i)    = Ve;
log_Xe(:,i)    = Xe;
log_Kj(1,i)    = Kj;
log_Ka(1,i)    = Ka;
log_omg(:,i)   = omg*180/pi;
log_q(:,i)     = q;
log_the(1,i)   = the*180/pi;
log_psi(1,i)   = psi*180/pi;
%
[xmax,tmax] = max(log_Xe(3,:));
tmin=size(log_Xe(3,:));
MP1=log_Xe(3,tmax:tmin(1,2));
MP2=find(MP1>Hpara);
tpara = size(MP2)*dt+tmax*dt;
if (real(log_Xe(3,i))<0)&&(log_t(1,i)>tThrust)
    break
end
end

% plot
figure
plot(log_t(1,:),log_Xe(1,:),'r',log_t(1,:),log_Xe(2,:),'b',...
    log_t(1,:),log_Xe(3,:),'g')
xlabel('Time, t [sec]');
ylabel('Distance, Xe [m]');
legend('Xe(1)','Xe(2)','Xe(3)');

figure
plot(log_t(1,:),log_Ve(1,:),'r',log_t(1,:),log_Ve(2,:),'b',...
    log_t(1,:),log_Ve(3,:),'g')
xlabel('Time, t [sec]');
ylabel('Velocity, Ve [m/s]');
legend('Ve(1)','Ve(2)','Ve(3)');

figure
plot(log_t(1,:),real(log_the(1,:)),'r',log_t(1,:),real(log_psi(1,:)),'b')
xlabel('Time, t [sec]');
ylabel('Attitude Angle [deg]');
legend('the','psi');

figure
plot(log_t(1,:),log_omg(2,:),'r',log_t(1,:),log_omg(3,:),'b')
xlabel('Time, t [sec]');
ylabel('Angular Velocity, [deg/s]');
legend('pitch (omg(2))','yaw (omg(3))');

figure
plot(log_t(1,:),log_alpha(1,:),'r',log_t(1,:),log_bet(1,:),'b')
xlabel('Time, t [sec]');
ylabel('Angle [deg]');
legend('Angle of Attack (alpha)','Angle of Sideslip (bet)');

figure
plot(log_t(1,:),log_D(1,:),'r',log_t(1,:),log_N(1,:),'b',...
    log_t(1,:),log_Y(1,:),'g');
xlabel('Time, t [sec]');
ylabel('Force [N]');
legend('Drag (D)','Normal Force (N)','Side Force (Y)');

figure
plot(log_t(1,:),log_Vab(1,:),'r',log_t(1,:),log_Vab(2,:),'b',...
    log_t(1,:),log_Vab(3,:),'g')
xlabel('Time, t [sec]');
ylabel('Air Speed, Va [m/s]');
legend('Vab(1)','Vab(2)','Vab(3)');