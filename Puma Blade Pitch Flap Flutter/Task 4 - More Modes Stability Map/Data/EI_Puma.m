function [EIv,waypoints]=EI_Puma(xv)
%This function returns the value of EI in a grid given by xv, using the
%data from the Puma experiments. The second argument gives the waypoints
%that should be considered in numerical integration on the interval
%[min(xv),max(xv)]. For simplicity, I'll make these coincide with the
%datapoints in the axial positions. They must be made unique, though. 

%Introduce experimental data: 
EI_Data = [
% Radial location	Flap stiffness
% r			(EI_f)_u(nmodified)	
% m			Nm^2			
0.280			178.00e4		
0.600			178.00e4		
0.600			178.00e4		
0.610			137.00e4		
0.610			137.00e4		
0.800			137.00e4		
0.800			137.00e4
0.810			 41.20e4	
0.810			 41.20e4		
1.240			 41.00e4		
1.240			 41.00e4		
1.250			  8.10e4		
1.250			  8.10e4		
7.300			  8.10e4		
7.300			  8.10e4		
7.310			  8.20e4		 
7.310			  8.20e4		 
7.490			  8.20e4
];

%Use the experimental data to assign values to EIv:
EIv=nan(size(xv));

for i=1:length(xv) %We assume that xv is larger than the first point on the grid, since it is before the flap hinge. 
    %We will assign the value corresponding to the inboard station. In
    %cases where discontinuities happen, we want to assign the value on the
    %correct side of the discontinuity. That's why the max of the indices
    %of points smaller than or equal to xv is used. 
    inboard_station=find(EI_Data(:,1)<=xv(i),1,'last');
    EIv(i)=EI_Data(inboard_station,2);
end

%Waypoints: needs to create some comparison arrays due to how Matlab works.
positions=EI_Data(:,1);
waypoints=unique(positions(min(xv)<positions & positions<max(xv))); %unique extracts and sorts only unique elements from the positions list. Notice that without it the waypoints would be doubled. 
end