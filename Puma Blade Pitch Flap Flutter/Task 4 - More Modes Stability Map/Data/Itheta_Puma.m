function [Ithetav,waypoints]=Itheta_Puma(xv)
%This function returns the value of Itheta in a grid given by xv, using the
%data from the Puma experiments. The second argument gives the waypoints
%that should be considered in numerical integration on the interval
%[min(xv),max(xv)]. For simplicity, I'll make these coincide with the
%datapoints in the axial positions. They must be made unique, though.
%Linear interpolation has been introduced to accurately represent the
%functional shapes given in the experimental documentation. 

% %Introduce experimental data: 
Itheta_Data = [
% Radial location	Torsional inertia
% r			(I_Theta)_u(nmodified)
% m			kg m
0.280			0.000
0.604			0.000
0.604			0.192
0.610			0.192
0.610			0.164
0.730			0.178
0.730			0.116
0.754			0.130
0.754			0.427
0.760			0.438
0.760			0.205
0.800			0.240
0.800			0.318
0.836			0.359
0.836			0.359
0.837			0.120
0.837			0.120
0.840			0.121
0.840			0.082
1.017			0.121
1.017			0.121
1.040			0.098
1.040			0.075
1.085			0.040
1.085			0.040
1.110			0.040
1.110			0.026
1.260			0.026
1.260			0.017
1.770			0.037
1.770			0.087
1.887			0.109
1.887			0.109
7.070			0.109
7.070			0.156
7.390			0.156
7.390			0.337
7.402			0.067
7.402			0.067
7.434			0.067
7.434			0.118
7.490			0.118
];

% Itheta_Data = [
% % Radial location	Torsional inertia
% % r			(I_Theta)_u(nmodified)
% % m			kg m
% 0.280			1e-6
% 0.604			1e-6
% 0.604			0.192
% 0.610			0.192
% 0.610			0.164
% 0.730			0.178
% 0.730			0.116
% 0.754			0.130
% 0.754			0.427
% 0.760			0.438
% 0.760			0.205
% 0.800			0.240
% 0.800			0.318
% 0.836			0.359
% 0.836			0.359
% 0.837			0.120
% 0.837			0.120
% 0.840			0.121
% 0.840			0.082
% 1.017			0.121
% 1.017			0.121
% 1.040			0.098
% 1.040			0.075
% 1.085			0.040
% 1.085			0.040
% 1.110			0.040
% 1.110			0.026
% 1.260			0.026
% 1.260			0.017
% 1.770			0.037
% 1.770			0.087
% 1.887			0.109
% 1.887			0.109
% 7.070			0.109
% 7.070			0.156
% 7.390			0.156
% 7.390			0.337
% 7.402			0.067
% 7.402			0.067
% 7.434			0.067
% 7.434			0.118
% 7.490			0.118
% ];

%Use the experimental data to assign values to EIv:
Ithetav=nan(size(xv));
for i=1:length(xv)
    inboard_station=find(Itheta_Data(:,1)<=xv(i),1,'last');
    %Handle the last point:
    if inboard_station==length(Itheta_Data(:,1))
        Ithetav(i)=Itheta_Data(end,2);
    else
        xl=Itheta_Data(inboard_station,1);
        xr=Itheta_Data(inboard_station+1,1);
        Ithetal=Itheta_Data(inboard_station,2);
        Ithetar=Itheta_Data(inboard_station+1,2);
        Ithetav(i)=Ithetal+(Ithetar-Ithetal)/(xr-xl)*(xv(i)-xl);
    end
end

%Waypoints: needs to create some comparison arrays due to how Matlab works.
positions=Itheta_Data(:,1);
waypoints=unique(positions(min(xv)<positions & positions<max(xv))); %unique extracts and sorts only unique elements from the positions list. Notice that without it the waypoints would be doubled. 
end