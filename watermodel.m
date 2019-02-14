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
Global T F Di Pi conductivity C_Ca C_Cl C_Mg C_Na C_S
T=273+20; %temperature in K
F=300/(1000*2*3600); %flow rate in m^3/s, 300L/min or 80 gpm
Di= 0.1016; %inner diameter of pipe (schedule 40), m
Pi= 45*6894.757; %domestic water pressure 45 psi to abs
conductivity=316; %microS/cm
electricity=0; %keeps track of electricity use
C_Ca=0.01; %concentration of species Ca in mol/m^3
C_Cl=0.01; %concentration of species Cl in mol/m^3
C_Mg=0.01; %concentration of species Mg in mol/m^3
C_Na=0.01; %concentration of species Na in mol/m^3
C_S=0.01; %concentration of species S in mol/m^3
C_Organism=0.03; %concentration in mol/m^3
%https://www.researchgate.net/post/How_is_it_possible_to_convert_conductivity_of_NaCl_solution_in_uS_to_its_salinity_NaCl_concentration_in_ppm
End

% Define sub-functions for each module
%%Accumulator
function a=accumulator (Po)
Global T F Di Pi conductivity Po
%%increase inlet water pressure
%negligible height, pipe length
Po=Pi+1000*9.8*(0.55*1864.25/9.8)%in Pa, Head as 62 ft, no friction
end

% Ultrafiltration
Function u= ultrafiltration (C_treated_organism)
Global T F Di Pi conductivity C_Organism Po
%removes viruses, bacteria, etc.
%requires air (air pump)


%good upto 6.25 bar, operating TMP is 2.1bar, filtrate flux is 40-120 L/m^2/hr
%Laminar fluid through the capillary (Hagen-Poiseuille Model for Laminar Flow through Channel)

kB=1.38064852e-23; %in J/K
Dp=0.03*10^-6; %mean pore diameter in meter
Do= (0.2e-6+5e-9)/2; %mean diameter of solute (ORGANISMS) (bacteria from 0.2 to 10 microns)(virus from 5-300 nm)
N=3*10^9*100^2; %number of pores per m^2
eps=N*(pi/4)*Dp^2; *surface porosity unitless
Dh=(4*(pi*0.7^2/4)/(pi*0.7))/1000; %Dh=4A/P generic formula for hydraulic diameter, 0.7mm inner diameter
mu=8.9e-4; %dynamic viscosity of permeate, that of water
L_pore=0.3/1000; %thickness of the membrane in m; 0.3mm
Po=101325;%??what would this be? DUMMY VARIABLE
D21=kB*T/(6*pi*Do*mu); %diffusivity of organism in liquid (Stokes-Einstein, only valid in Low Reynolds)
J1=((eps*Dh^2)/(32*mu*L_pore))*(Pi-Po); %initial guess of permeate flux (or can we just use this???)

%iterative process
error=1;
While error>0.01
	v=J1*4/((0.7/1000)^2); %flow velocity (of which one??)

	%Re
	Re=1000*v*Dh/(8.90*10^-4);

	%Sc
	%phi=Do/Dp; %solute diameter/pore diameter
	%F1=(1-phi)^2; F2=1-2.104*phi+2.09*phi^3-0.95*phi^5; %from CHE312
	%D21=D21*F1*F2; %solute diffusing through pore
	Sc=(8.90*10^-4)/(1000*D21);

	%Sh (is this the right one?? this one is for solvent flow 
	L=0.05; %in meters, channel length; (from CHE 490)
	Lstar=0.029*Dh;
	If Re<2100 && L>Lstar
	Sh=0.664*(Re^0.5)*(Sc^0.33)*(Dh/L)^0.33;
		If Re<2100 && L<Lstar
		Sh=1.86*(Re^0.33)*(Sc^0.33)*(Dh/L)^0.33;
			Else
			%turbulent regime
			Sh=0.023*(Re^0.83)*(Sc^0.33);
			End
		End
	End

	%(C_2s-C_2p)=(C_2f-C_2p)exp(J*/k)
	%C_2s: surface concentration of solute (unknown)
	%C_2p: permeate concentration of solute (unknown)
	%C_2f: bulk concentration of solute
	%THIS IS ONLY FOR BACTERIA, VIRUS, etc.
	C_2f=C_organism;
	k=Sh*D21/Dh;
	alpha=2.084555/0.001204 %P2/P1 from CHE 490 
	MolarVolumeWater=(18.02/1000)/1000; %in m^3/mol
	osmoticP=(8.314*T*(C_Ca+C_Cl+C_Mg+C_Na+C_S+C_Organism)); %lump sum for osmotic pressure (is this okay??)
	deltaX=((Pi-Po)-osmoticP)/(8.314*T);

	%need to solve for both C_2p and C_2s (iterative process)
	fun-@ROfunction
	C0=[0.01 0.01]; %C_2p and C_2s initial guess
	Answer=@ROfunction(F,C0);

	%new osmotic pressure
	Answer(2)=C_2s;
	osmoticP=(8.314*T*(C_Ca+C_Cl+C_Mg+C_Na+C_S+C_2s));

	%calculate new flux
	deltaX=((Pi-Po)-osmoticP)/(8.314*T);
	C_1s=(1000/18.02)/1000; %in mol/m^3
	J1=0.001204*C_1s*deltaX; %solvent flux
End
C_treated_organism=C_2p;
J1=F; %permeate flow is the new flowrate to the next module
End

function F=ROfunction(C_2f,C_2p)
F=[(alpha/(alpha+deltaX))*C_2s-C_2p;
	(C_2f-C_2p)exp(J*/k)-(C_2s-C_2p)];	
End

% Ion exchange
Function c= IEC ()

fnd

%Reverse osmosis
Function r=RO ()
%refer to CHE 490 RO Lab

end

%Mixing
Function m=mixing ()
%mixing assembly power consumption
%diffusion of salt in water

end

%Cost calculation
Function cost=cost()
Global T F Di Pi conductivity
%calculate annual op cost based on water and power consumption, salt, filter replacement, maintenance, and other contingencies
AnnualCost=Electricity*0.089+1.5*F*1000*2*3600+
%calculate capital cost based on piping, and equipment installation

end
