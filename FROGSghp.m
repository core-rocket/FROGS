% FROGS
% ver1.8 (190807edited)
% for NSE15th
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
global S SIMULATION Dpara HeightH

global LeleDeg LazDeg
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Choose the type of simulation(弾道or減速)
%%% 3=Ballistic fall，4=Retarding fall 5=Delay time
SIMULATION  = 3;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

FROGSparameters;        % parameterの読み込み
FROGSthrust;            % thrustデータの読み込み
%
GHP = zeros(14,17);
Delays = zeros(7,17);
DELAY = csvread('DelayTime.csv');
%Winddata = readmatrix('Winddata.csv');

min_launchclear = 1000000.0;
max_top_vel = 0.0;
max_para_vel= 0.0;
all_xmax = 0.0;
all_xmin = 100 * 1000.0;
all_xmax_t = 0.0;
all_xmin_t = 0.0;
all_nmax = 0.0;
all_max_Va = 0.0;

fprintf("launcher:\n");
fprintf("  length   = %f\n", lLnchr);
fprintf("  azimuth  = ");
if LazDeg == 0
	fprintf("east\n");
elseif LazDeg == 90
	fprintf("north\n");
elseif LazDeg == 180
	fprintf("west\n");
elseif LazDeg == 270
	fprintf("south");
end
fprintf("  elevation= %f\n", LeleDeg);

WIND_DIR_NUM	= 16;
FLAG_RESTRICT	= 1;

tic

for Vtemp = 1:7
	Vwaz = 0+Vtemp*1.0;
	%Vwaz = 5.5;

	fprintf("  wind speed: %f m/s\n", Vwaz);

	for k = 1:(WIND_DIR_NUM+1)
		WazDeg = (360/WIND_DIR_NUM) * (k-1);
		Waz = WazDeg*pi/180;

		if FLAG_RESTRICT == 1
			if WazDeg == 360
				continue;
			end

			wdir = WazDeg + 180.0;	% 風が吹いていく方向->風向
			wdir = wdir - 90.0;		% 東0deg ->北0deg
			if wdir >= 360.0
				wdir = wdir - 360.0;
			end
			
			fprintf("    wdir = %f\n", wdir);

			if wdir == 0.0
				if Vwaz == 1 && LeleDeg > 87
					continue;
				elseif Vwaz >= 2 && Vwaz <= 3 && LeleDeg > 88
					continue;
				end
			elseif wdir == 22.5
				if LeleDeg > 87
					continue;
				elseif Vwaz == 3 && LeleDeg > 86
					continue;
				elseif Vwaz == 4 && LeleDeg > 85
					continue;
				elseif Vwaz == 5 && LeleDeg > 84
					continue;
				elseif Vwaz == 6 && LeleDeg > 83
					continue;
				elseif Vwaz == 7 && LeleDeg > 82
					continue;
				end
			elseif wdir == 45.0
				if Vwaz > 5
					continue;
				elseif Vwaz == 1 && LeleDeg > 88
					continue;
				elseif Vwaz == 2 && LeleDeg > 86
					continue;
				elseif Vwaz == 3 && LeleDeg > 84
					continue;
				elseif Vwaz == 4 && LeleDeg > 82
					continue;
				elseif Vwaz == 5 && LeleDeg > 80
					continue;
				end
			elseif wdir == 67.5
				if Vwaz > 4
					continue;
				elseif Vwaz == 1 && LeleDeg > 87
					continue;
				elseif Vwaz == 2 && LeleDeg > 85
					continue;
				elseif Vwaz == 3 && LeleDeg > 82
					continue;
				elseif Vwaz == 4 && LeleDeg > 80
					continue;
				end
			elseif wdir == 90.0
				if Vwaz > 4
					continue;
				elseif Vwaz == 1 && LeleDeg > 87
					continue;
				elseif Vwaz == 2 && LeleDeg > 85
					continue;
				elseif Vwaz == 3 && LeleDeg > 82
					continue;
				elseif Vwaz == 4 && LeleDeg > 80
					continue;
				end
			elseif wdir == 112.5
				if Vwaz > 5
					continue;
				elseif Vwaz == 1 && LeleDeg > 88
					continue;
				elseif Vwaz == 2 && LeleDeg > 86
					continue;
				elseif Vwaz == 3 && LeleDeg > 84
					continue;
				elseif Vwaz == 4 && LeleDeg > 82
					continue;
				elseif Vwaz == 5 && LeleDeg > 80
					continue;
				end
			elseif wdir == 135.0
				if LeleDeg > 88
					continue;
				elseif Vwaz == 3 && LeleDeg > 86
					continue;
				elseif Vwaz == 4 && LeleDeg > 85
					continue;
				elseif Vwaz == 5 && LeleDeg > 84
					continue;
				elseif Vwaz == 6 && LeleDeg > 83
					continue;
				elseif Vwaz == 7 && LeleDeg > 81
					continue;
				end
			elseif wdir == 157.5
				if LeleDeg > 89
					continue;
				elseif Vwaz == 2 && LeleDeg > 88
					continue;
				elseif Vwaz == 3 && LeleDeg > 87
					continue;
				elseif Vwaz == 4 && LeleDeg > 87
					continue;
				elseif Vwaz == 5 && LeleDeg > 88
					continue;
				elseif Vwaz == 6 && LeleDeg > 88
					continue;
				elseif Vwaz == 7 && LeleDeg > 89
					continue;
				end
			elseif wdir == 180
				if LeleDeg > 88
					continue;
				elseif Vwaz == 2 && LeleDeg > 88
					continue;
				elseif Vwaz == 3 && LeleDeg > 86
					continue;
				elseif Vwaz == 4 && LeleDeg > 86
					continue;
				elseif Vwaz == 5 && LeleDeg > 86
					continue;
				elseif Vwaz == 6 && LeleDeg > 87
					continue;
				elseif Vwaz == 7 && LeleDeg > 87
					continue;
				end
			elseif wdir == 202.5
				if LeleDeg > 88
					continue;
				elseif Vwaz == 2 && LeleDeg > 87
					continue;
				elseif Vwaz == 3 && LeleDeg > 86
					continue;
				elseif Vwaz == 4 && LeleDeg > 85
					continue;
				elseif Vwaz == 5 && LeleDeg > 85
					continue;
				elseif Vwaz == 6 && LeleDeg > 85
					continue;
				elseif Vwaz == 7 && LeleDeg > 85
					continue;
				end
			elseif wdir == 225.0
				if LeleDeg > 87
					continue;
				elseif Vwaz == 2 && LeleDeg > 87
					continue;
				elseif Vwaz == 3 && LeleDeg > 86
					continue;
				elseif Vwaz == 4 && LeleDeg > 85
					continue;
				elseif Vwaz == 5 && LeleDeg > 84
					continue;
				elseif Vwaz == 6 && LeleDeg > 83
					continue;
				elseif Vwaz == 7 && LeleDeg > 83
					continue;
				end
			elseif wdir == 247.5
				if LeleDeg > 87
					continue;
				elseif Vwaz == 2 && LeleDeg > 87
					continue;
				elseif Vwaz == 3 && LeleDeg > 86
					continue;
				elseif Vwaz == 4 && LeleDeg > 85
					continue;
				elseif Vwaz == 5 && LeleDeg > 85
					continue;
				elseif Vwaz == 6 && LeleDeg > 83
					continue;
				elseif Vwaz == 7 && LeleDeg > 83
					continue;
				end
			elseif wdir == 270.0
				if LeleDeg > 87
					continue;
				elseif Vwaz == 2 && LeleDeg > 87
					continue;
				elseif Vwaz == 3 && LeleDeg > 86
					continue;
				elseif Vwaz == 4 && LeleDeg > 85
					continue;
				elseif Vwaz == 5 && LeleDeg > 85
					continue;
				elseif Vwaz == 6 && LeleDeg > 84
					continue;
				elseif Vwaz == 7 && LeleDeg > 83
					continue;
				end
			elseif wdir == 292.5
				if LeleDeg > 87
					continue;
				elseif Vwaz == 2 && LeleDeg > 87
					continue;
				elseif Vwaz == 3 && LeleDeg > 86
					continue;
				elseif Vwaz == 4 && LeleDeg > 85
					continue;
				elseif Vwaz == 5 && LeleDeg > 85
					continue;
				elseif Vwaz == 6 && LeleDeg > 85
					continue;
				elseif Vwaz == 7 && LeleDeg > 84
					continue;
				end
			elseif wdir == 315.0
				if LeleDeg > 89
					continue;
				elseif Vwaz == 1 && LeleDeg > 87
					continue;
				elseif Vwaz == 2 && LeleDeg > 87
					continue;
				elseif Vwaz == 3 && LeleDeg > 86
					continue;
				elseif Vwaz == 4 && LeleDeg == 87
					continue;
				elseif Vwaz == 5 && (LeleDeg == 87 || LeleDeg == 89)
					continue;
				elseif Vwaz == 6 && (LeleDeg == 87 || LeleDeg == 89)
					continue;
				elseif Vwaz == 7 && (LeleDeg == 87 || LeleDeg == 89)
					continue;
				end
			elseif wdir == 337.5
				if LeleDeg > 89
					continue;
				elseif Vwaz == 1 && LeleDeg > 87
					continue;
				elseif Vwaz == 2 && LeleDeg > 87
					continue;
				elseif Vwaz == 3 && LeleDeg > 88
					continue;
				elseif Vwaz == 4 && LeleDeg > 89
					continue;
				elseif Vwaz == 5 && LeleDeg > 89
					continue;
				elseif Vwaz == 6 && LeleDeg > 89
					continue;
				elseif Vwaz == 7 && LeleDeg > 89
					continue;
				end
			end

			fprintf("      ");
		else
			fprintf("    wind dir: %3.1f deg ", WazDeg);
		end

		[Ve1,Ve2,Ve3,Xe1,Xe2,Xe3,omg2,omg3,q1,q2,q3,q4]  = FROGSprset;
		IV  = [Ve1,Ve2,Ve3,Xe1,Xe2,Xe3,omg2,omg3,q1,q2,q3,q4];
		Ve  = [IV(1); IV(2); IV(3)];            %地球座標系における機体速度ベクトル
		Xe  = [IV(4); IV(5); IV(6)];            %地球座標系における機体位置ベクトル
		omg = [0; IV(7); IV(8)];                %角速度
		q   = [IV(9); IV(10); IV(11); IV(12)];  %クォータニオン
		i = 1;                                  %ステップ数
		t = 0;                                  %時間
		%%%
		log_Xe	= zeros(3,n);
		log_N	= zeros(1,n);

		V_launchclear = 0.0;

		for i = 1:n                             % 1‾nまで繰り返し計算
			t = i*dt;                               % time [s]
			% transformation matrix(earth frame --> body frame)(座標変換行列)(地上座標系を機体座標系へ)
			Aeb=[1-2*((q(2)^2)+(q(3)^2)),...
				2*(q(1)*q(2)+q(3)*q(4)),2*(q(3)*q(1)-q(2)*q(4));
				2*(q(1)*q(2)-q(3)*q(4)),1-2*((q(1)^2)+(q(3)^2)),...
				2*(q(2)*q(3)+q(1)*q(4));
				2*(q(3)*q(1)+q(2)*q(4)),2*(q(2)*q(3)-q(1)*q(4)),...
				1-2*((q(1)^2)+(q(2)^2))];

			if t<tThrust								% 燃焼中の飛翔
				T       = thrust(i);					% thrust [N]
				mdot    = (m0-mf)/tThrust;				% weight loss [kg/sec]
				m       = m0-mdot*t;					% weight [kg]
				I       = I0-(I0-If)/tThrust*t;			% moment of inertia [kgm^2]
				lcg     = lcg0-(lcg0-lcgf)/tThrust*t;	% center of gravity [m]
			else										% 燃焼終了後
				T       = 0;							% thrust [N]
				mdot    = 0;							% weight loss [kg/sec]
				m       = mf;							% weight [kg]
				I       = If;							% moment of inertia [kgm^2]
				lcg     = lcgf;							% center of gravity [m]
			end

			%rho = 1.23-(0.11e-3)*Xe(3);				% atmospheric density
														% Don't use this equation over 2km
			r0=6356766;
			H = (r0*Xe(3)*0.001)/(r0+Xe(3)*0.001);		% ジオポテンシャル高度
			Temp = 15 - 6.5*H;
			P =101325 * (288.15./(Temp +273.15)).^(-5.256);
			rho = (0.0034837*P)/(Temp+273.15);			% Don't use this equation over 11km

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
					Vw = [Vwaz*cos(Waz);Vwaz*sin(Waz);0];             
				case 3                                                % 統計風
					if i==1
						Vw = [Vwaz*cos(Waz);Vwaz*sin(Waz);0];
					elseif Xe(3)>=HeightH
						Xew = ceil(Xe(3));
						VwazH = Winddata(Xew-HeightH,1);
						WazH = Winddata(Xew-HeightH,2);
						Vw = [VwazH*cos(WazH);VwazH*sin(WazH);0];
					else
						Vwazl = Vwaz+(Winddata(1,1)-Vwaz)/(HeightH-Zr)*Xe(3);
						if (Winddata(1,2)+pi) < Waz
							Wazl = Waz+(Winddata(1,2)+2*pi-Waz)/(HeightH-Zr)*Xe(3);
						elseif (Winddata(1,2)-pi) > Waz
							Wazl = Waz+(Winddata(1,2)-2*pi-Waz)/(HeightH-Zr)*Xe(3);
						else
							Wazl = Waz+(Winddata(1,2)-Waz)/(HeightH-Zr)*Xe(3);
						end
						Vw = [Vwazl*cos(Wazl); Vwazl*sin(Wazl);0];
					end
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
				case 3							% 弾道落下
					Ae = Aeb'*(Fe./m)-[0;0;g];
					Ve = Ae*dt+Ve;
					Xe = Ve*dt+Xe;
				case 4							% 減速落下：下段の燃焼終了後かつ速度が負でパラ展開，終端速度に即達する
					if (Ve(3)<=0)&&(t>tThrust)&&(Xe(3)>Hpara)
						Ae = [0;0;0];
						Ve = [Vw(1); Vw(2); -Vpara1*tanh(g*(t-tmax*dt)/Vpara1)];
						Xe = Ve*dt+Xe;
					elseif (Ve(3)<=0)&&(t>tThrust)&&(Xe(3)<Hpara)			% 2段パラ
						Ae = [0;0;0];
						Ve = [Vw(1); Vw(2); -Vpara1+(Vpara1-Vpara2)*tanh(g*(t-tpara(1,2))/(Vpara1-Vpara2))];
						Xe = Ve*dt+Xe;
					else													% 上昇中
						Ae = Aeb'*(Fe./m)-[0;0;g];
						Ve = Ae*dt+Ve;
						Xe = Ve*dt+Xe;
					end
				case 5							% 減速落下：下段の燃焼終了後かつ速度が負でパラ展開，終端速度に即達する
					if (t>=DELAY(Vtemp,k))&&(t>tThrust)&&(Xe(3)>Hpara)
						Ae = [0;0;0];
						Ve = [Vw(1); Vw(2); -Vpara1*tanh(g*(t-tmax*dt)/Vpara1)];
						Xe = Ve*dt+Xe;
					elseif (t>=DELAY(Vtemp,k))&&(t>tThrust)&&(Xe(3)<Hpara)	% 2段パラ
						Ae = [0;0;0];
						Ve = [Vw(1); Vw(2); -Vpara1+(Vpara1-Vpara2)*tanh(g*(t-tpara(1,2))/(Vpara1-Vpara2))];
						Xe = Ve*dt+Xe;
					else													% 上昇中
						Ae = Aeb'*(Fe./m)-[0;0;g];
						Ve = Ae*dt+Ve;
						Xe = Ve*dt+Xe;
					end
			end

			Veab	= norm(Ve);		% 機速
			Aeab	= norm(Ae);		% 機速

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
					if (Xe(3)<lLnchr)&&(t<tThrust)					%ランチャに刺さってる時は回転しない
						omg = [0;0;0];
					else											%飛翔中ずっとルンゲクッタ
						omg = [0;
							1/6*(kt1+2*kt2+2*kt3+kt4)*dt+omg(2);
							1/6*(kp1+2*kp2+2*kp3+kp4)*dt+omg(3)];
					end
				case 4
					if (Xe(3)<lLnchr)&&(t<tThrust)					%ランチャに刺さってる時は回転しない
						omg = [0;0;0];
					elseif (Ve(3)<=0)&&(t>tThrust)					%減速落下中も回転しない
						omg = [0;0;0];
					else											%それ以外の時(上昇中)ルンゲクッタ
						omg = [0;
							1/6*(kt1+2*kt2+2*kt3+kt4)*dt+omg(2);
							1/6*(kp1+2*kp2+2*kp3+kp4)*dt+omg(3)];
					end
				case 5
					if (Xe(3)<lLnchr)&&(t<tThrust)					%ランチャに刺さってる時は回転しない
						omg = [0;0;0];
					elseif (t>=DELAY(Vtemp,k))&&(t>tThrust)			%減速落下中も回転しない
						omg = [0;0;0];
					else											%それ以外の時(上昇中)ルンゲクッタ
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

			% ランチクリア速度
			if (norm(Xe) < lLnchr) && (t < tThrust)
				Vlc = Veab(1);
				if Vlc > V_launchclear
					V_launchclear = Vlc;
				end
			end

			log_Xe(:,i)	= Xe;
			log_Va(:,i)	= Va;
			log_N(:,i)	= N;

			[xmax,tmax] = max(log_Xe(3,:));
			tmin=size(log_Xe(3,:));
			MP1=log_Xe(3,tmax:tmin(1,2));
			MP2=find(MP1>Hpara);
			tpara = size(MP2)*dt+tmax*dt;
			tdelay = tmax*dt+Dpara;
			%
			if (Xe(3)<0)&&(t>tThrust)
				break
			end
		end

		nmax	= abs(max(log_N));
		tmp_Va	= log_Va(1,:);
		max_Va  = max(tmp_Va);
		top_vel	= tmp_Va(tmax);
		para_vel= tmp_Va(tmax+(2.5/dt));

		if V_launchclear < min_launchclear
			min_launchclear = V_launchclear;
		end
		if top_vel > max_top_vel
			max_top_vel = top_vel;
		end
		if para_vel > max_para_vel
			max_para_vel = para_vel;
		end
		if xmax > all_xmax
			all_xmax = xmax;
			all_xmax_t = tmax*dt;
		end
		if xmax < all_xmin
			all_xmin = xmax;
			all_xmin_t = tmax*dt;
		end
		if nmax > all_nmax
			all_nmax = nmax;
		end

		if max_Va > all_max_Va
			all_max_Va = max_Va;
		end

		fprintf("para vel=%f, altitude=%f, N=%f", para_vel, xmax, nmax);

		%fprintf("GHP: %f, %f", real(Xe(1)), real(Xe(2)));
		GHP(2*Vtemp-1,k) = real(Xe(1));
		GHP(2*Vtemp,k) = real(Xe(2));
		%GHP(1,k) = real(Xe(1));
		%GHP(2,k) = real(Xe(2));
		%Delays(Vtemp,k) = real(tdelay);
		fprintf("\n");
	end

	if FLAG_RESTRICT ~= 1
		fname = sprintf("ghp-%d-%d.csv", SIMULATION, LeleDeg);
        writematrix(GHP, fname);
		% plot
		plot(GHP(2*Vtemp-1,:),GHP(2*Vtemp,:),'-squareb');
		%plot(GHP(1,:),GHP(2,:),'-squareb');
		hold on;
	end
end

fprintf("simulation time: %f sec\n", toc)
fprintf("min launch clear vel: %f\n", min_launchclear);
fprintf("max air vel:  %f\n", all_max_Va);
fprintf("max top vel:  %f\n", max_top_vel);
fprintf("max para vel: %f\n", max_para_vel);
fprintf("max alitude:  %f m (%f sec)\n", all_xmax, all_xmax_t);
fprintf("min altitude: %f m (%f sec)\n", all_xmin, all_xmin_t);
fprintf("max N:        %f\n", all_nmax);
