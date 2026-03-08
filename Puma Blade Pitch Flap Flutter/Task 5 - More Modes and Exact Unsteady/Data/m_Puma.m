function [mv,waypoints]=m_Puma(xv)
%This function returns the value of m in a grid given by xv, using the
%data from the Puma experiments. The second argument gives the waypoints
%that should be considered in numerical integration on the interval
%[min(xv),max(xv)]. For simplicity, I'll make these coincide with the datapoints.
%They must be made unique, though. Linear interpolation is used to match
%the graph provided in the Puma experimental results. 

%Introduce experimental data: 
m_Data = [
% Radial location	Running mass
% r			m_u(nmodified)
% m			kg/m
0.280			58.400
0.610			58.400
0.610			50.000
0.730			50.000
0.730			16.175
0.754			16.175
0.754			53.333
0.760			53.333
0.760			24.949
0.800			24.949
0.800			33.100
0.840			33.100
0.840			22.400
1.040			22.400
1.040			17.110
1.110			17.110
1.110			11.225
1.260			11.225
1.260			 7.150
1.770			 7.150
1.770			 7.150
1.887			 8.929
1.887			 8.929
7.070			 8.929
7.070			12.754
7.390			12.754
7.390			34.600
7.402			34.600
7.402			 6.930
7.490			 6.930
];
%Use the experimental data to assign values to mv: linear interpolation
mv=nan(size(xv));
for i=1:length(xv)
    inboard_station=find(m_Data(:,1)<=xv(i),1,'last');
    %Handle the last point:
    if inboard_station==length(m_Data(:,1))
        mv(i)=m_Data(end,2);
    else
        xl=m_Data(inboard_station,1);
        xr=m_Data(inboard_station+1,1);
        ml=m_Data(inboard_station,2);
        mr=m_Data(inboard_station+1,2);
        mv(i)=ml+(mr-ml)/(xr-xl)*(xv(i)-xl);
    end
end

%Waypoints: needs to create some comparison arrays due to how Matlab works.
positions=m_Data(:,1);
waypoints=unique(positions(min(xv)<positions & positions<max(xv))); %unique extracts and sorts only unique elements from the positions list. Notice that without it the waypoints would be doubled. 
end