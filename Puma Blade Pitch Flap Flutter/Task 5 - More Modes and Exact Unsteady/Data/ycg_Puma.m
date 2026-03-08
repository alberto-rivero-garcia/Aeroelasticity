function [ycgv,waypoints]=ycg_Puma(xv)
%This function returns the value of ycg in a grid given by xv, using the
%data from the Puma experiments. The second argument gives the waypoints
%that should be considered in numerical integration on the interval
%[min(xv),max(xv)]. For simplicity, I'll make these coincide with the axial
%position of the datapoints. They must be made unique, though. 

%Introduce experimental data: 
ycg_Data = [
% Radial location	C.G. offset (wrt. quarter chord, positive forward)
% r			(x/c)_u(nmodified)
% m			adim.
0.280			 0.000
1.887			 0.000
1.887			-0.010
7.070			-0.010
7.070			 0.030
7.390			 0.030
7.390			 0.147
7.402			 0.147
7.402			 0.000
7.490			 0.000
];
%Use the experimental data to assign values to EIv:
ycgv=nan(size(xv));

for i=1:length(xv) %We assume that xv is larger than the first point on the grid, since it is before the flap hinge. 
    %We will assign the value corresponding to the inboard station. In
    %cases where discontinuities happen, we want to assign the value on the
    %correct side of the discontinuity. That's why the max of the indices
    %of points smaller than or equal to xv is used. 
    inboard_station=find(ycg_Data(:,1)<=xv(i),1,'last');
    ycgv(i)=ycg_Data(inboard_station,2);
end

%Waypoints: needs to create some comparison arrays due to how Matlab works.
positions=ycg_Data(:,1);
waypoints=unique(positions(min(xv)<positions & positions<max(xv))); %unique extracts and sorts only unique elements from the positions list. Notice that without it the waypoints would be doubled. 
end