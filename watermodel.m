Function FYDP
% Input data
%% input water quality
%% input constraints
filename = 'water in and out.xlsx';
wq_o=xlsread(filename,'B:B');      %required outlet water quality (column B)
wq_i =xlsread(filename,'C:C');     % input water quality (column C)
wq_i =xlsread(filename,'D:D');     % input constraints (column D)

%identify which modules are being used
filename1 = 'moduleinfo.xlsx';
dataRO=xlsread(filename1,RO);%import entire sheet
dataIEC=xlsread(filename1,IEC);%import entire sheet
dataUltra=xlsread(filename1,ULTRA);%import entire sheet
dataFilter=xlsread(filename1,FILTER)%import entire sheet

%define parameter for single unit operation
%% Filtration 1-3
FILTER1=dataFilter(2,:); %second row
FILTER2=dataFilter(3,:); %third row
FILTER3=dataFilter(4,:); %fourth row
%% RO 1-3
RO1=dataRO(2,:); %second row
RO2=dataRO(3,:); %third row
RO3=dataRO(4,:); %fourth row
%% IEC 1-3
IEC1=dataIEC(2,:); %second row
IEC2=dataIEC(3,:); %third row
IEC3=dataIEC(4,:); %fourth row
%% Ultra-purification 1-3
ULTRA1=dataUltra(2,:); %second row
ULTRA2=dataUltra(3,:); %third row
ULTRA3=dataUltra(4,:); %fourth row

%%conductivity equation
%valid for NaCl only
%conductivity=2.4472ppm; %from CHE480

%define operating condition
Global T F Di Pi conductivity
T=273+20; %temperature in K
F=300/(1000*2*3600); %flow rate in m^3/s, 300L/min or 80 gpm
Di= 0.1016; %inner diameter of pipe (schedule 40), m
Pi= 45*6894.757; %domestic water pressure 45 psi to abs
conductivity=316; %microS/cm
electricity=0; %keeps track of electricity use

%https://www.researchgate.net/post/How_is_it_possible_to_convert_conductivity_of_NaCl_solution_in_uS_to_its_salinity_NaCl_concentration_in_ppm
end

% Define sub-functions for each module
%%Accumulator
function a=accumulator ()
Global T F Di Pi conductivity
%%increase inlet water pressure
%negligible height, pipe length
Po=Pi+1000*9.8*(0.55*1864.25/9.8)%in Pa, Head as 62 ft, no friction


end

%%Filtration
function f= filtration ()
	%similar to RO equation except different diffusivity
	
end

% Ultrafiltraion
function u= ultrafiltration ()

end

% Ion exchange
Function c= IEC ()

fnd

%Reverse osmosis
Function r=RO ()

end

%Mixing
Function m=mixing ()

end

%Cost calculation
Function cost=cost()
Global T F Di Pi conductivity
%calculate annual op cost based on water and power consumption, salt, filter replacement, maintenance, and other contingencies
AnnualCost=Electricity*0.089+1.5*F*1000*2*3600+
%calculate capital cost based on piping, and equipment installation

end
