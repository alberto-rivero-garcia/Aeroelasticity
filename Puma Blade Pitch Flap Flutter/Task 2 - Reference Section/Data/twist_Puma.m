function [twistv,waypoints]=twist_Puma(xv)
%This function returns the value of twist in a grid given by xv, using the
%data from the Puma experiments. The second argument gives the waypoints
%that should be considered in numerical integration on the interval
%[min(xv),max(xv)]. For simplicity, I'll make these coincide with the
%datapoints in the axial positions. They must be made unique, though.
%Spline interpolation is used to match the functional shapes provided in
%the experimental result document. 

%Introduce experimental data: 
twist_Data = [
% Radial location	Twist
% r			Theta_u(nmodified)
% m			deg
0.280			 0.000
0.760			 0.000
1.757			 0.000
2.007			-0.100
2.257			-0.300
2.507			-0.617
2.757			-0.983
3.007			-1.283
3.257			-1.600
3.507			-1.867
6.257			-4.800
7.490			-6.117
];
%Use the experimental data to assign values to EIv:
twistv=interp1(twist_Data(:,1),twist_Data(:,2),xv,'spline');

%Waypoints: needs to create some comparison arrays due to how Matlab works.
positions=twist_Data(:,1);
waypoints=unique(positions(min(xv)<positions & positions<max(xv))); %unique extracts and sorts only unique elements from the positions list. Notice that without it the waypoints would be doubled. 
end