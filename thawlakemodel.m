% Thaw Lake Model-1D

% This model a 1-D numerical model of permafrost and subsidence processes. 
% It aims to investigate the subsurface thermal impact of thaw lakes of various depths, 
% and to evaluate how this impact might change in a warming climate. 

% Key paper: Matell, N., Anderson, R.S., Overeem, I., Wobus, C., Urban, F.,
% Clow, G., in review 2011. Modeling the subsurface thermal impact of Arctic thaw lakes in a warming climate. Computers and Geosciences. 
 
% Copyright (C) <2011> <Nora Matell, Irina Overeem, Robert Anderson, Cameron Wobus>

% Developer can be contacted by irina.overeem@colorado.edu

% Dr. Irina Overeem
% CSDMS Community Surface Dynamics Modeling System
% INSTAAR, University of Colorado at Boulder
% PO Box 450, 80309-0450
% Boulder, CO, USA


% This program is free software; you can redistribute it and/or modify it under the terms of the GNU General Public License as published by the Free Software Foundation; either version 2 of the License, or (at your option) any later version.
% This program is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details.
% You should have received a copy of the GNU General Public License along with this program; if not, write to the Free Software Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301 USA. 

% This model a 1-D numerical model of permafrost and subsidence processes. 
% It aims to investigate the subsurface thermal impact of thaw lakes of various depths, 
% and to evaluate how this impact might change in a warming climate.  

% The model was designed for the Alaskan Arctic Coastal Plain
% The model uses average and amplitude of temperature from observed permafrost temperatures at 5cm in the subsurface at Drew Point data (USGS)
% Clow, G.D., 2008a. Continued permafrost warming in northern Alaska, 2008 update. NOAA/ESRL Global Monitoring Annual Conference, Boulder.  
% Clow, G. D., 2008b. USGS Polar Temperature Logging System, Description and Measurement Uncertainties. U.S. Geological Survey Techniques and Methods 2–E3, 24 p.  

% The model uses observed heat flow for thermal gradient from:
% Lachenbruch, A.H., Sass, J.H., Marshall, B.V., Moses, T.H., Jr., 1982. Temperatures, heat flow, and the geothermal regime at Prudhoe By, Alaska. Journal of Geophysical Research 87, 9301-9316.  


% in 1-D, subsurface temperature changes through time and space of a semi-infinite
% half-space, then plots result.  The model code incorporates phase change by using
% the apparent heat capacity scheme for temperatures within the
% "phase-change envelope".  

% lake-permafrost model - lake freezes and thaws, permafrost changes temperature,
% when permafrost thaw it subsides, assuming that if excess ice melts all 
% the water leaves the permafrost and thus the subsurface volume decreases.  

% boundary condition of ice bottom and lake water top = Tf
% if Ts<Tf, top boundary condition for calculation is Ts.  if Ts>Tf, top
% boundary condition --> Tf and excess energy is used to melt ice.  When
% ice melted from the top, all other boundary just moved up by equivelant
% amount, if moved up a full control volume then bottom control volume
% added that is same temp as (end-1) control volume.  Assumes that when ice
% free, lake is completely mixed.  (ie Liston & Hall)



clear all
figure(1)
clf
figure(2)
clf
figure(3)
clf
figure(4)
clf
figure(5)
clf

% set up 1-d grid
dz0 = 0.05;                             % cell size in top 10 m
dz2 = 1;								% cell size below 10 m
deptht = 500;                           % model depth [m]
z = [0:dz0:10,11:dz2:deptht];           % vertical grid (distance below surface) [m]
dz = diff(z);
zbetween = z(2:end)+diff(z)/2;
dzbetween = diff(zbetween);
numcells = length(z);                   % number of grid cells
iceinit = 0.0001;                       % initial ice depth [m]
celltype = int16(1*ones(1,numcells));   % substrate type - 1=regular soil with excess ice, 2=compacted soil without excess ice, 8=water, 9=ice
celltype(find(z<6.7)) = 2;               % depth of subsided cells, this is a function of the initial lake depth (seed lake depth/excess ice content)
celltype(find(z<2.0)) = 8;              % seed lake depth
celltype(find(z<=iceinit)) = 9;         % ice thickness
celltype(find(z>=50)) = 2;               % underlying cells with no excess ice
thawed = zeros(1,numcells);
thawedspec = zeros(1,numcells);

deptht = max(z);
depthsubside = 2.0;			% seed lake depth, has to match with ln 79, because the lake contains water
depthtalik = 0;				% usually unkown, so the model will trend toward equilibrium in the first 100-200 years of a simulation
depthi = iceinit; 

% set time steps and simulation duration
periodyear = 3600*24*365;       % period (1 year) [s]
dt = 3600*24;                   % length of timestep [s]
years = 1000;					% simulation duration
nt = ((periodyear/dt)*years)+1; % number of timesteps
t = 0:dt:dt*nt;                 % time at each timestep [s]
tday = t/(3600*24); 
tyear = t/(3600*24*365);
tyearsonly = 1:years;
yearnum = 0;

% declare arrays

depthin = zeros(1,nt);
depthwn = zeros(1,nt);
depthsubsiden = zeros(1,nt+1);
depthtalikn = zeros(1,nt+1);
icethickness = zeros(1,nt+1);
icegrowth = zeros(1,nt+1);
soilthickness = zeros(1,nt+1);
waterthickness = zeros(1,nt+1);
Tbot = zeros(1,nt+1);
T3m = zeros(1,nt+1);
T5m = zeros(1,nt+1);
T10m = zeros(1,nt+1);
T25m = zeros(1,nt+1);
T50m = zeros(1,nt+1);
T100m = zeros(1,nt+1);
Tbotavgann = zeros(1,years);
T3mavgann = zeros(1,years);
T5mavgann = zeros(1,years);
T10mavgann = zeros(1,years);
T25mavgann = zeros(1,years);
T50mavgann = zeros(1,years);
T100mavgann = zeros(1,years);
Trecord = zeros(365,numcells);
watts = zeros(size(t));


% define constants
% Model input and boundary conditions (compiled from field and modeling studies by West & Plug, 2008; Ling and Zhang, 2003;.....
% Ling & Zhang, 2004; Zhou and Huang, 2004; Romanovsky and Osterkamp, 2000; Hinzman, 1998; and French, 2007).

kelvin = 273.15;

rhoi = 917;         % ice density [kg/m3]
ci = 2108;          % ice specific heat capacity [J/kg/K]
ki = 2.18;          % ice thermal conductivity [J/s/m/K]
kappai = ki/(ci*rhoi);
rhow = 1000;        % water density [kg/m3]
cw = 4210;          % water specific heat capacity [J/kg/K]
kw = 0.58;          % water thermal conductivity [J/s/m/K]
kappaw = kw/(cw*rhow);
rhorock = 1200;     % mineral soil density [kg/m3]
crock = 1000;       % mineral soil specific heat capacity [J/kg/degree C]
krock = 1.5;        % rock thermal conductivity [J/s/m/degree C]

totalice = 0.65;     % ice content of original substrate, from bulk samples in top 5m at Drew Point
excessice = 0.3;	% excess ice content - gone once melted
porewater = totalice-excessice;     % porewater ice - still there if refreezes
W = totalice;       % total water content percent of material by mass
Wu = totalice*0.05;  % unfrozen percent water content of material by mass at 
                    % temperature T = Tpc-pcenv
Ws = porewater;
Wus = porewater*0.05;
totalice2 = totalice-Wu;
porewater2 = porewater-Wus;

cf = (ci*totalice2)+(cw*Wu)+(crock*(1-totalice));        
cu = (cw*totalice)+(crock*(1-totalice)); 
cfs = (ci*porewater2)+(cw*Wu)+(crock*(1-porewater));        
cus = (cw*porewater)+(crock*(1-porewater));
rhof = (rhoi*totalice2)+(rhow*Wu)+(rhorock*(1-totalice));
rhou = (rhow*totalice)+(rhorock*(1-totalice));
rhofs = (rhoi*porewater2)+(rhow*Wus)+(rhorock*(1-porewater));
rhous = (rhow*porewater)+(rhorock*(1-porewater));
Cf = cf*rhof;       % frozen volumetric heat capacity [J/m3/degree C]
Cu = cu*rhou;       % thawed volumetric heat capacity [J/m3/degree C]
Cfs = cfs*rhofs;
Cus = cus*rhous;
kf = (ki^totalice2)*(kw^Wu)*(krock^(1-totalice));
ku = (kw^totalice)*(krock^(1-totalice));
kfs = (ki^porewater2)*(kw^Wus)*(krock^(1-porewater));
kus = (kw^porewater)*(krock^(1-porewater));

L = 334*1000;       % latent heat of fusion for water [J/kg]
hw = 0.56;          % convective transfer coefficient =[J/s/m2/K]
Kh = 0.0;           % turbulent diffusivity for heat
Tpc = kelvin+0;     % freezing temperature [K]
pcenv = 1;          % width of phase change envelope [m]


% set up surface and initial temperatures, including geothermal gradient
% theoretical warming scenario is set here
Tbar = kelvin-11;    % MAAT [K]
warming = 5;
Tbarplus = Tbar+(warming/36500)*(1:36500);
Tamp = 17.5;          % amplitude of air temperature fluctuations [K]
%Tsurface = Tbar*ones(1,nt);
Tsurface = Tbar-Tamp*sin(2*pi*(t)/periodyear);          % surface temperatures at each timestep [K]
Tsurface(182500:218999) = Tbarplus(1:end)-Tamp*sin(2*pi*(t(182500:218999))/periodyear);
Tsurface(219000:end) = Tbar+warming-Tamp*sin(2*pi*(t(219000:end))/periodyear);
q = 0.056;          % mantle heat flow [J/m2/s]
dTdzbase = q/kf;
dTdzbase2 = q/kfs;
Tgrad = [dTdzbase.*z(find(z<=10)),dTdzbase2.*z(find(z>10))];
Twi = Tsurface(1);  % start up water temperature [K]
Tsi = Tbar;         % start up permafrost temperature [K]
Tinit = ones(size(z));     % initial temperature grid     
Tinitlake = [Tpc*Tinit(find(celltype==9)),Twi*Tinit(find(celltype==8))];
water = length(Tinitlake);
for n = 1:numcells
    if celltype == 9
        Tinit(n) = Tpc;
    elseif celltype == 8
        Tinit(n) = Twi;
    else
        Tinit(n) = Tsi+Tgrad(n);
    end
end
% load daily average radiation, function from Drew Point meteorological data
load('radin_dailyavg'); 
        
%Tinitperm = [Tsi*Tinit(find(celltype==1))+Tgrad(find(celltype==1)),Tsi*Tinit(find(celltype==2))+Tgrad(find(celltype==2))];
%Tinit = [Tinitlake, Tinitperm];
T = Tinit;
zstarthaw = sqrt((ku/Cu)*periodyear/(pi));
zstarfrozen = sqrt((kfs/Cf)*periodyear/(pi));
Tbase = Tbar+(dTdzbase2*deptht)+Tamp.*exp(-deptht/zstarthaw).*sin((2*pi*t./periodyear)...
        -(deptht/zstarthaw));
Tright = -kelvin+Tbar+(dTdzbase2.*z)+Tamp.*exp(-z/zstarthaw); %outer edges of the funnel
Trightf = -kelvin+Tbar+(dTdzbase2.*z)+Tamp.*exp(-z/zstarfrozen); %outer edges of the funnel
Tleftf = -kelvin+Tbar+(dTdzbase2.*z)-Tamp.*exp(-z/zstarfrozen);
Tleft = -kelvin+Tbar+(dTdzbase2.*z)-Tamp.*exp(-z/zstarthaw);
%Ts0 = Tsurface(1);
          


count = 1;
daynum = 0;
numplots = 50;
countprint = dt*nt/numplots;
nplot=0;
figure(1)
zero = zeros(size(z));
plot(zero,-z,'k','linewidth',2)
hold on
plot(Tright,-z,'g',Tleft,-z,'g',Trightf,-z,'g',Tleftf,-z,'g')
hold on
plot(Tinit-kelvin,-z)
axis([Tbar-Tamp-kelvin Tbar+Tamp-kelvin -deptht 0]);
xlabel('Temperature (C)','fontsize',18)
ylabel('Depth (m)','fontsize',18)
%hold on

% model equations
for n = 1:nt;
    thawedspec = zeros(1,numcells);
    thawed = zeros(1,numcells);
    clear 'Tcalc' 'Tcalci' 'Tcalcw' 'Tcalcperm' 'Tit' 'Titi' 'Titw' 'Ti' 'Tw' 'Twater' 'Tice' 'T2'
    
    count = count+1;
    time = n*dt;    % time into run [s]
    error = 1;
    Told = T;       % remember temperatures from last timestep... 
    Ts0 = Tsurface(n);
    T(1)=Ts0; 
    if depthi>0 && Tsurface(n)>Tpc
        Ts0 = Tpc;
    end

    water = 0;
    ice = 0;
    subsided = 0;
    regsoil = 0;
    for i=1:numcells
        if celltype(i)==8
            water = water+1;
        elseif celltype(i)==9
            ice = ice+1;
        elseif celltype(i)==2
                subsided = subsided+1;
        elseif celltype(i)==1
            regsoil = regsoil+1;  
        end
    end

    day = rem(n,365)+1;
    solarrad = radin_dailyavg(day);
    
    % compute temperatures in water layer if no ice exists
    if celltype(1)==8
       Twater = T(1:water);%T(find(celltype==8));
       %Tmix = mean(Twater);  % thoroughly mix lake
       
       %Tmix-kelvin
       %Twater = Tmix*ones(size(Twater));
       %Twater(1) = Ts0;
            %%% explicit
       %for i = 2:water-1
       %    acoeff = cw*rhow*dzbetween(i-1)/dt;
       %    bcoeff = kw/dz(i);
       %    ccoeff = kw/dz(i-1);
       %    dcoeff = acoeff-bcoeff-ccoeff;
       %    Tcalc(i) = (bcoeff*Told(i+1)+ccoeff*Told(i-1)+dcoeff*Told(i))/acoeff;
       %end
       %Twater = [Ts0,Tcalc(2:end),Twater(end)];
       Tmix = mean(Twater); % mix again
%       Tmix-kelvin
       %T(find(celltype==8)) = Tmix*ones(size(Twater));
            %%%implicit
%       fprintf('no ice on lake,should get radiation')
       %solarrad = 100;      % incoming solar radiation;
       sextinc = 0.6;       % solar extinction constant
       albedo = 0.06;
       qrad = (1-albedo)*solarrad*exp(-sextinc*zbetween(1:water));
       watts(count) = solarrad;
       
       error = 1;
       while error>0.0001
           Tit = Twater;
           for i = 2:water-1
               coeff1 = kw/dz(i);
               coeff2 = kw/dz(i-1);
               coeff3 = cw*rhow*dzbetween(i-1)/dt;
               coeff4 = coeff1+coeff2+coeff3;
               Tcalc(i) = (coeff1*T(i+1)+coeff2*T(i-1)+coeff3*Told(i)+qrad(i-1)-qrad(i))/coeff4;
           end
           Twater=[Ts0,Tcalc(2:end),Tcalc(end)];
           error = max(abs(Twater-Tit));
       end

       Twater2 = Twater(1:end-1);
       Tmix = mean(Twater2);  % thoroughly mix lake
%       Tmix-kelvin
       %Tmix = mean(Twater) % mix again
       T(2:water) = Tmix*ones(size(Twater2));
       
          
       
    % compute temperatures in ice layer if completely frozen to bottom
    elseif water==0
            %%%explicit
       %for i = 2:ice-1
       %    acoeff = ci*rhoi*dzbetween(i-1)/dt;
       %    bcoeff = ki/dz(i);
       %    ccoeff = ki/dz(i-1);
       %    dcoeff = acoeff-bcoeff-ccoeff;
       %    Tcalc(i) = (bcoeff*Told(i+1)+ccoeff*Told(i-1)+dcoeff*Told(i))/acoeff;
       %end
       %Tice=[Ts0,Tcalc(2:end),Tice(end)];     
       %T(find(celltype==9)) = Tice;
       Tice = T(find(celltype==9));
            %%%implicit
       while error>0.0001
           Tit = Tice;
           for i = 2:ice
           %%%
           %for i = 2:ice-1
               coeff1 = ki/dz(i);
               coeff2 = ki/dz(i-1);
               coeff3 = ci*rhoi*dzbetween(i-1)/dt;
               coeff4 = coeff1+coeff2+coeff3;
               Tcalc(i) = (coeff1*T(i+1)+coeff2*T(i-1)+coeff3*Told(i))/coeff4;
           end
           Tice = [Ts0,Tcalc(2:end)];
           %%%%%
           %Tice=[Ts0,Tcalc(2:end),Tcalc(end)];
           error = max(abs(Tice-Tit));
       end    
       T(find(celltype==9)) = Tice; 
   
    % compute temperatures in ice and water layers if both exist
    else
       % ice layer with top boundary Ts0, bottom boundary Tpc if freezing,
       % isothermal at Tpc if melting
       if Tsurface(n)>=Tpc
           T(1:ice) = Tpc;
       elseif ice<3
           T(1) = Tsurface(n);
           T(ice) = Tpc;   
       else
                %%% explicit
           %for i = 2:ice-1
           %    acoeff = ci*rhoi*dzbetween(i-1)/dt;
           %    bcoeff = ki/dz(i);
           %    ccoeff = ki/dz(i-1);
           %    dcoeff = acoeff-bcoeff-ccoeff;
           %    Tcalci(i) = (bcoeff*Told(i+1)+ccoeff*Told(i-1)+dcoeff*Told(i))/acoeff;
           %end
           %Ti=[Ts0,Tcalci(2:end),Tpc];
           %T(1:ice) = Ti;
           Ti = T(1:ice);
                %%% implicit
           error = 1;
           while error>0.0001
               Titi = Ti;
               for i = 2:ice-1
                   coeff1 = ki/dz(i);
                   coeff2 = ki/dz(i-1);
                   coeff3 = ci*rhoi*dzbetween(i-1)/dt;
                   coeff4 = coeff1+coeff2+coeff3;
                   Tcalci(i) = (coeff1*T(i+1)+coeff2*T(i-1)+coeff3*Told(i))/coeff4;
               end
               Ti=[Ts0,Tcalci(2:end),Tpc];
               error = max(abs(Ti-Titi));
           end
           T(1:ice) = Ti;
        end
       
       % water layer with top boundary Tf
       error = 1;
       if water==1
        T(find(celltype==8)) = Tpc;
       %elseif maxcelli>(numcells-2)
       %    T(end) = Tpc;
       else
           %%% explicit
          % for i = ice+1:ice+water-1
          %     acoeff = cw*rhow*dzbetween(i-1)/dt;
          %    bcoeff = kw/dz(i);
          %    ccoeff = kw/dz(i-1);
          %    dcoeff = acoeff-bcoeff-ccoeff;
          %    Tcalc(i) = (bcoeff*Told(i+1)+ccoeff*Told(i-1)+dcoeff*Told(i))/acoeff;
          %end
            %%% implicit
          error = 1;
          while error>0.0001
              Tit = T;
              for i = ice+1:ice+water
                coeff1 = kw/dz(i);
                coeff2 = kw/dz(i-1);
                coeff3 = cw*rhow*dzbetween(i-1)/dt;
                coeff4 = coeff1+coeff2+coeff3;
                Tcalcw(i) = (coeff1*T(i+1)+coeff2*T(i-1)+coeff3*Told(i))/coeff4;
              end
              Tw=[Tpc,Tcalcw(ice+1:end)];
              T(ice:ice+water) = Tw;
              error = max(abs(T-Tit));
          end
       end
    end
    
    % Calculate temperatures in underlying permafrost
    clear 'z2' 'zbetween2' 'Told2' 'T2' 'celltype2'
    %numcells2 = subsided+regsoil+1;
    z2 = z(ice+water:end);       
    zbetween2 = zbetween(ice+water:end);
    Told2 = Told(ice+water:end);
    T2 = T(ice+water:end);
    numcells2 = length(T2);
    celltype2 = celltype(ice+water:end);    
    C = zeros(1,numcells2);
    k = zeros(1,numcells2);
    Tcalc = 0;
    for i=1:numcells2;
        if celltype2(i)==8
            C(i) = cw*rhow;
            k(i) = kw;
        elseif celltype2(i)==9
            C(i) = ci*rhoi;
            k(i) = ki;
        elseif celltype2(i)==1 && T2(i)<(Tpc-pcenv)
            C(i) = Cf;
            k(i) = kf;
            %fprintf('real permafrost, frozen')
        elseif celltype2(i)==1 && T2(i)>=(Tpc-pcenv) && T2(i)<=Tpc
            C(i) = Cf+L*rhof*((W-Wu)/pcenv);
            k(i) = kf+((ku-kf)/pcenv)*(T2(i)-(Tpc-pcenv));
            %fprintf('real permafrost, thawing')
        elseif celltype2(i)==1 && T2(i)>Tpc
            C(i) = Cu;
            k(i) = ku;
            %fprintf('real permafrost, thawed')
        elseif celltype2(i)==2 && T2(i)<(Tpc-pcenv)
            C(i) = Cfs;
            k(i) = kfs;
        elseif celltype2(i)==2 && T2(i)>=(Tpc-pcenv) && T2(i)<=Tpc
            C(i) = Cfs+L*rhofs*((Ws-Wus)/pcenv);
            k(i) = kfs+((kus-kfs)/pcenv)*(T2(i)-(Tpc-pcenv));
        elseif celltype2(i)==2 && T2(i)>Tpc
            C(i) = Cus;
            k(i) = kus;
        end
    end

    kbetween = k(2:end)-(diff(k)/2);
    dz2 = diff(z2);
    dzbetween2 = diff(zbetween2);
        %%% explicit central difference
    %for i = 2:numcells2-1
    %    acoeff = C(i)*dzbetween2(i-1)/dt;
    %    bcoeff = kbetween(i)/dz2(i);
    %    ccoeff = kbetween(i-1)/dz2(i-1);
    %    dcoeff = acoeff-bcoeff-ccoeff;
    %    Tcalcperm(i) = (bcoeff*Told2(i+1)+ccoeff*Told2(i-1)+dcoeff*Told2(i))/acoeff;
    %end
    %T(ice+water:end) = [T(ice+water),Tcalcperm(2:end),Tbase(n)];
        %%% implicit
    error = 1;
    while error>0.0001
        Tit = T2;
        for i = 2:numcells2-1
            coeff1 = kbetween(i)/dz2(i);
            coeff2 = kbetween(i-1)/dz2(i-1);
            coeff3 = C(i-1)*dzbetween2(i-1)/dt;
            coeff4 = coeff1+coeff2+coeff3;
            Tcalcperm(i) = (coeff1*T2(i+1)+coeff2*T2(i-1)+coeff3*Told2(i))/coeff4;
        end
        T2 = [T(ice+water),Tcalcperm(2:end),Tbase(n)];
        error = max(abs(T2-Tit));
    end
    T(ice+water:end) = T2;
    
    % subside permafrost
    for i = water+ice+1:numcells;
        if T(i)<Tpc
            break 
        elseif T(i)>Tpc && celltype(i) == 1
            thawedspec(i) = 1;
            thawed(i) = 1;
            celltype(i) = 2;
        elseif T(i)>Tpc && celltype(i) == 2
            thawed(i) = 1;
            %celltype(i) = 2;
        end
    end
    numthawedspec =sum(thawedspec)+1;
    depthsubsidenew = z(numthawedspec)*excessice;
    depthsubside = depthsubside+depthsubsidenew;
    depthsubsiden(n) = depthsubside;
    cellsubside = int16(floor((depthsubside)/dz0));
    maxcellsubside = cellsubside;
    if depthsubside>0%depthsubsidenew>0
        for i = ice+1:maxcellsubside
            celltype(i) = 8;
        end
    end
    
    %depthtalik = sum(thawed)*dz0;
    %depthtalikn(n) = depthtalik;
    
    
    % find ice and water depths    
    if Tsurface(n)<Tpc;
        icecheckdepth = int16(0.33/dz0)+ice;    % z location of ~0.33 m into water column
        %icecheckdepthREM(n) = icecheckdepth;
        if icecheckdepth>ice+water
            icecheckdepth = ice+water-1;
        end
        Twater = T(icecheckdepth);
        if water==1
            Twater = T(find(celltype==8));
            Twater = Twater(1);
        end
        if Twater<Tpc && ice==0
            depthi = 0.01;
            %%%%depthw = deptht-depthi;
        end

        depthmix = z(icecheckdepth)-depthi;
        ddepthidt = (((Tpc-Ts0)*((depthi/ki))^(-1))-hw*(Twater-Tpc))./(rhoi*L);
        ddepthi = ddepthidt*dt;
        %%%%% ADJUSTMENT HERE
        if ddepthi>0.05
            ddepthi = 0.05;
            fprintf('had to adjust')
        end
        icegrowth(n) = ddepthi;
        depthi = depthi+ddepthi;
        if depthi<0.00001
            depthi = 0;
        end
    
        depthw = depthsubside-depthi;
        if depthi>depthsubside
            depthi = depthsubside;
            depthw = 0;
        end
    
    % check if Ts>Tpc (if there is ice there already)
    elseif Tsurface(n)>Tpc && depthi>0
        Qm = (Tsurface(n)-Tpc)*(((dz0/ki))^(-1));
        dMdt = Qm/(L*rhoi);
        dM = dMdt*dt;
        icegrowth(n) = -dM;
        depthi = depthi-dM;
        if depthi<0
            depthi = 0;
        end
        %depthw = deptht-depthi;
    end

    %maxcelli = int16(ceil(depthi/dz0));
    maxcelli = int16(round(depthi/dz0));
    celltype(1:ice+water) = 8;
    celltype(1:maxcelli) = 9;
    
    icethickness(n) = depthi;
    soilthickness(n) = deptht-depthsubside;
    waterthickness(n) = deptht-icethickness(n)-soilthickness(n);
    
    
    Tbot(n) = T(ice+water+1)-kelvin;
    T3m(n) = T(60)-kelvin;
    T5m(n) = T(100)-kelvin;
    T10m(n) = T(200)-kelvin;
    T25m(n) = T(215)-kelvin;
    T50m(n) = T(240)-kelvin;
    T100m(n) = T(290)-kelvin;
    if rem(tday(n),365)==0 && tday(n)>364
        yearnum = yearnum+1;
        Tbotavgann(yearnum) = mean(Tbot(n-364:n));
        T3mavgann(yearnum) = mean(T3m(n-364:n));
        T5mavgann(yearnum) = mean(T5m(n-364:n));
        T10mavgann(yearnum) = mean(T10m(n-364:n));
        T25mavgann(yearnum) = mean(T25m(n-364:n));
        T50mavgann(yearnum) = mean(T50m(n-364:n));
        T100mavgann(yearnum) = mean(T100m(n-364:n));
    end

    if n>=nt-364
    %if rem(n,10)==0
        daynum = daynum+1;
        Trecord(daynum,1:end) = T;
    end
    
        
    %print temperature profile with depth at intermediate time steps

    if(rem (n,(365*50)) ==0)   
    %    nplot = nplot+1
        figure(1)
    %    ice
        icethickness(n) 
    %    water
        waterthickness(n)
    %    %soilthickness(n) 
    %    %thickness = icethickness(n)+waterthickness(n)+soilthickness(n)
        lakedepthx = [-20:1:25];
        lakedepth = icethickness(n)+waterthickness(n)*ones(size(lakedepthx));
        Tplot = T-kelvin;
        plot(Tplot,-z,'m')
        %hold on
        timeplot = tday(n)/365 %will stamp time in years on your screen
    %    %pause(0.5)
        xlabel('temperature','fontname','arial','fontsize',18)
        ylabel('depth','fontname','arial','fontsize',18)
    %    axis([-30 25 -5 0])
    %
        figure(2)
        clf(2)
        icedepth = icethickness(n)*ones(size(lakedepthx));
        plot(zero,-z,'k','linewidth',2)
        hold on
        plot(lakedepthx,-lakedepth,'b',lakedepthx,-icedepth,'g',Tplot,-z,'m')
        xlabel('temperature','fontname','arial','fontsize',18)
        ylabel('depth','fontname','arial','fontsize',18)
        axis([-30 25 -5 0])
        %celltype(1:200)
    end
    
    
end
 

% Plot output parameters after run is completed

% plot surface temperature fluctuations
figure(2)
TsK = Tsurface-kelvin;
plot(tday,TsK)
xlabel('time(days)','fontname','arial','fontsize',18)
ylabel('surface temperature(degrees C)','fontname','arial','fontsize',18)    
    
    
% plot thickness of lake ice
%figure(3)
%plot(tday,icethickness)
%xlabel('time(days)','fontname','arial','fontsize',18)
%ylabel('lake ice thickness (m)','fontname','arial','fontsize',18)

% plot thickness
figure(3)
plot(tyear(1:end-1),icethickness(1:end-1),tyear(1:end-1),waterthickness(1:end-1),...
    tyear(1:end-1),depthsubsiden(1:end-1))
xlabel('time(years)','fontname','arial','fontsize',18)
ylabel('thickness (m)','fontname','arial','fontsize',18)
legend('ice thickness','water thickness','depth subside')

% plot temps at various depths
figure(4)
T0C = zeros(size(tyear));
plot(tyearsonly,Tbotavgann,tyearsonly,T3mavgann,tyearsonly,T5mavgann,tyearsonly,T10mavgann,tyearsonly,T25mavgann,...
    tyearsonly,T50mavgann,tyearsonly,T100mavgann,tyear(1:end-1),T0C(1:end-1))
xlabel('time(years)','fontname','arial','fontsize',18)
ylabel('temperature (C)','fontname','arial','fontsize',18)
legend('at lake bottom','at 3m depth','at 5m depth','at 10m depth','at 25m depth','at 50m depth','at 100m depth','0C')

% plot temps at lake bottom
figure(5)
T0C = zeros(size(tyear));
plot(tyear(1:end-1),Tbot(1:end-1),tyearsonly,Tbotavgann,tyear(1:end-1),T0C(1:end-1))
xlabel('time(years)','fontname','arial','fontsize',18)
ylabel('temperature (C)','fontname','arial','fontsize',18)
legend('at lake bottom','lake bottom average annual temp','0C')




