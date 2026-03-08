function [GJv,waypoints]=GJ_Puma(xv)
%This function returns the value of GJ in a grid given by xv, using the
%data from the Puma experiments. The second argument gives the waypoints
%that should be considered in numerical integration on the interval
%[min(xv),max(xv)]. For simplicity, I'll make these coincide with the
%datapoints in the axial positions. They must be made unique, though. 

%Introduce experimental data: 
GJ_Data = [
% Radial location	Torsional stiffness
% r			(GJ)_u(nmodified)
% m			Nm^2
0.280			 84.00e4
0.725			 84.00e4
0.725			226.00e4
0.828			226.00e4
0.828			 50.50e4
1.017			 50.50e4
1.017			  8.50e4
7.241			  8.50e4
7.241			  8.70e4
7.490			  8.70e4
];
%Use the experimental data to assign values to EIv:
GJv=nan(size(xv));

for i=1:length(xv) %We assume that xv is larger than the first point on the grid, since it is before the flap hinge. 
    %We will assign the value corresponding to the inboard station. In
    %cases where discontinuities happen, we want to assign the value on the
    %correct side of the discontinuity. That's why the max of the indices
    %of points smaller than or equal to xv is used. 
    inboard_station=find(GJ_Data(:,1)<=xv(i),1,'last');
    GJv(i)=GJ_Data(inboard_station,2);
end

%Waypoints: needs to create some comparison arrays due to how Matlab works.
positions=GJ_Data(:,1);
waypoints=unique(positions(min(xv)<positions & positions<max(xv))); %unique extracts and sorts only unique elements from the positions list. Notice that without it the waypoints would be doubled. 
end