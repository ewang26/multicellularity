%% Diffusion Model
%% Cellular Automata model for growth of yeast communities

% This code is free to use for any purposes, provided any publications resulting from the use of this code reference the original code/author.
% Author:  Babak Momeni (bmomeni@gmail.com)
% Date:    01/2013
% The code is not guaranteed to be free of errors. Please notify the author of any bugs, and contribute any modifications or bug fixes back to the original author.

% Diffusion in 3D; diffusion inside the community slower than agarose
% 3D diffusion for transfer of nutrients
% nonhomogeneous media
% space occupied by dead cells
% FCU: fast uptake
% CI: confinement independent
% ES3: Expansion to sides to fill the confinement neighborhood; once confined, bud to top with 70% chance or otherwise to the sides, and push up all the cells in the column above
% EB: Boundary condition (continuous flow) explicitly applied at the agar-community interface
% EBP: Periodic boundary condition on the sides (in addition to EB)
% LRA: Release of adenine by live R cells
% GD: Glucose-dependent growth

clear
rseed = 2417;
rand('twister',rseed)
savename = 'CA3DCOOPCHEAT_GD_FCU_CI_ES3_EBP_LRA_c5g15d60_395N3300_35N3300_415N3300_t1000_dt1_5s_rg5Nc144_Ndz20Naz200_D0_360_D1_20v1';

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Parameters
% 1: cooperator (R),	2: partner (G),		3: cheater (C)
r1 = 0.395; % reproduction rate, per hour
r2 = 0.35; % reproduction rate, per hour
r3 = 0.415; % reproduction rate, per hour
d1 = 0.054; % death rate, per hour
d2 = 0.018; % death rate, per hour
d3 = 0.054; % death rate, per hour
alphaA = 1; % reproduction nutrient consumption factor; fmole per cell
alphaL = 2; % reproduction nutrient consumption factor; fmole per cell
alphaG = 1; % reproduction nutrient consumption factor; fmole per cell
gammaA = 0.08; % ade release rate; fmole per cell per hour
betaL = 15; % death nutrient release factor; fmole per cell
KA = 1e5; % Michaelis-Menten constant; fmole/ml
KL = 1e6; % Michaelis-Menten constant; fmole/ml
KG = 1e6; % Michaelis-Menten constant; fmole/ml

T1 = log(2)/r1; % minimum division time
T2 = log(2)/r2; % minimum division time
T3 = log(2)/r3; % minimum division time
vmL1 = alphaL*r1/log(2); % maximum uptake rate, fmole/hour
vmA = alphaA*r2/log(2); % maximum uptake rate, fmole/hour
vmL3 = alphaL*r3/log(2); % maximum uptake rate, fmole/hour
vmG = alphaG/T1; % maximum uptake rate, fmole/hour

D0 = 360; % diffusion constant in agarose/agar, microns2/s
D1 = 20; % diffusion constant in cells, microns2/s
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Definition of solution domain
% Length scale
c = 5; % grid size for cells, microns
SC = 3; % ratio of diffusion to cell grids
SD = 4; % ratio of diffusion grid (agar to cells)
g = SC*c; % grid size for diffusion in cells, microns
f = SD*g; % grid size for diffusion in agar, microns
dV = c^3 * 1e-12; % volume of each grid cube for cells, ml
dVc = g^3 * 1e-12; % volume of each grid cube for diffusion in cells, ml
dVa = f^3 * 1e-12; % volume of each grid cube for diffusion in agar, ml
h = 0.5*(g+f); % grid size at the interface between the two grids

rg0 = 5; % neghborhood radius for expansion to sides

%% Diffusion simulation domain
Nax = 12; % agar domain size along x
Nay = 12; % agar domain size along y
Naz = 200; % agar domain height
Ndx = Nax*SD; % cell domain size for diffusion
Ndy = Nay*SD; % cell domain size for diffusion
Ndz = 20; % cell domain height for diffusion

%% Cell simulation domain
Ncx = Ndx*SC; % cell domain size for cells
Ncy = Ndy*SC; % cell domain size for cells
Ncz = Ndz*SC; % cell domain height for cells
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Initialization of nutrient concentration functions
SLa = zeros(Nax,Nay,Naz); % concentration of lys in agarose
SLd = zeros(Ndx,Ndy,Ndz); % concentration of lys in community, diffusion-grid resolution
SLde = zeros(Ncx,Ncy,Ncz); % concentration of lys in community, cell-grid resolution
SAa = zeros(Nax,Nay,Naz); % concentration of ade in agarose
SAd = zeros(Ndx,Ndy,Ndz); % concentration of ade in community, diffusion-grid resolution
SAde = zeros(Ncx,Ncy,Ncz); % concentration of ade in community, cell-grid resolution
SG0 = 2e8; % fmole/ml
SGa = SG0*ones(Nax,Nay,Naz); % concentration of glucose in agarose
SGd = zeros(Ndx,Ndy,Ndz);
SGde = zeros(Ncx,Ncy,Ncz);
SLb = zeros(Ndx,Ndy); % concentration of lys at the boundary of agar and community; 0.5 or 1 or 0 makes no difference
SAb = zeros(Ndx,Ndy); % concentration of ade at the boundary of agar and community; 0.5 or 1 or 0 makes no difference
SGb = zeros(Ndx,Ndy); % concentration of glucose at the boundary of agar and community; 0.5 or 1 or 0 makes no difference

Dd = zeros(Ndx,Ndy,Ndz); % Diffusion coefficient throughout the community region

[xxd,yyd] = meshgrid((SD+1)/2:Ndx+(SD-1)/2,(SD+1)/2:Ndy+(SD-1)/2);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Initial cell distrbution
CSD = 3300; % initial cell surface density per type, 1/mm2
N1 = round(CSD*c*Ncx*c*Ncy/1e6); % initial number of type 1 cells, (R)
N2 = round(CSD*c*Ncx*c*Ncy/1e6); % initial number of type 2 cells, (Y)
N3 = round(CSD*c*Ncx*c*Ncy/1e6); % initial number of type 3 cells, (C)
X = zeros(Ncx,Ncy,Ncz); % Cells
U = zeros(Ncx,Ncy,Ncz); % Nutrients taken up by cells
UG = zeros(Ncx,Ncy,Ncz); % Glucose taken up by cells
T = zeros(Ncx,Ncy,Ncz); % Last division time for each cell
NB = zeros(Ncx,Ncy,Ncz); % Number of budding events at each location
nn = 0;
while nn < N1,
    i = floor(Ncx*rand(1))+1;
    j = floor(Ncy*rand(1))+1;
    k = 1;
    if X(i,j,k) == 0;
        nn = nn+1;
        X(i,j,k) = 1;
        T(i,j,k) = -rand(1)*T1;
        U(i,j,k) = rand(1)*alphaL;
        UG(i,j,k) = U(i,j,k)*alphaG/alphaL; %to account for the asynchronous nature of cells - how much nutrient accumulated up to t=0
    end
end
nn = 0;
while nn < N2,
    i = floor(Ncx*rand(1))+1;
    j = floor(Ncy*rand(1))+1;
    k = 1;
    if X(i,j,k) == 0;
        nn = nn+1;
        X(i,j,k) = 2;
        T(i,j,k) = -rand(1)*T2;
        U(i,j,k) = rand(1)*alphaA;
        UG(i,j,k) = U(i,j,k)*alphaG/alphaA; %to account for the asynchronous nature of cells - how much nutrient accumulated up to t=0
    end
end
nn = 0;
while nn < N3,
    i = floor(Ncx*rand(1))+1;
    j = floor(Ncy*rand(1))+1;
    k = 1;
    if X(i,j,k) == 0;
        nn = nn+1;
        X(i,j,k) = 3;
        T(i,j,k) = -rand(1)*T3;
        U(i,j,k) = rand(1)*alphaL;
        UG(i,j,k) = U(i,j,k)*alphaG/alphaL; %to account for the asynchronous nature of cells - how much nutrient accumulated up to t=0
    end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Cell-growth time-course
tau0 = 0; % in hours
tauf = 1000; % in hours
dtau = 0.1; % in hours, cell growth update and uptake timescale
ts = 2; % in hours, sampling time for snapshots of sections

r1e = r1*dtau; % probability of reproduction in a time step
r2e = r2*dtau; % probability of reproduction in a time step
r3e = r3*dtau; % probability of reproduction in a time step

d1e = d1*dtau; % probability of death in a time step
d2e = d2*dtau; % probability of death in a time step
d3e = d3*dtau; % probability of death in a time step

taurng = tau0:dtau:tauf;

% Initialization
X1m = zeros(size(taurng));
X2m = zeros(size(taurng));
X3m = zeros(size(taurng));
X1lm = zeros(size(taurng));
X2lm = zeros(size(taurng));
X3lm = zeros(size(taurng));
SPLa = zeros(Naz,length(taurng));
SPLd = zeros(Ndz,length(taurng));
SPAa = zeros(Naz,length(taurng));
SPAd = zeros(Ndz,length(taurng));
SPGa = zeros(Naz,length(taurng));
SPGd = zeros(Ndz,length(taurng));
UAccL = zeros(size(taurng)); % Accumulated nutrients in live cells
UAccA = zeros(size(taurng)); % Accumulated nutrients in live cells
SLAccA = zeros(size(taurng)); % Total lys nutrients in the agar region
SLAccC = zeros(size(taurng)); % Total lys nutrients in the cell region
SAAccA = zeros(size(taurng)); % Total ade nutrients in the agar region
SAAccC = zeros(size(taurng)); % Total ade nutrients in the cell region

ct = 0;
cS = 0;

Q = [1 0; -1 0; 0 1; 0 -1; 1 -1; -1 1; 1 1; -1 -1]; % locations of eight neighbor grids
Qc = [1 0; -1 0; 0 1; 0 -1]; % locations of bud for confined cells

for tau = taurng,
    ct = ct+1;
    tic

    % Update cell activity
    [zm,zi] = sort(rand(1,Ncz-1)); % random order for cells at different heights
    for z = zi,
        zd = floor((z-1)/SC)+1;
        [I,J] = find((X(:,:,z)==1)+(X(:,:,z)==2)+(X(:,:,z)==3)); % find all live cells at height z
        Ncc = length(I);
        [SS,OS] = sort(rand(1,Ncc));
        I = I(OS);
        J = J(OS);
        for cc = 1:Ncc;

            Id = floor((I(cc)-1)/SC)+1; % diffusion index along I
            Jd = floor((J(cc)-1)/SC)+1; % diffusion index along J
            xd = I(cc); % location of cell along x
            yd = J(cc); % location of cell along y

            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
			% cell live/dead decision
			live = 1;
            cd = rand(1);
            if (X(xd,yd,z) == 1)&&(cd < d1e),
                X(xd,yd,z) = 0.5;
                live = 0;
            end
            if (X(xd,yd,z) == 2)&&(cd < d2e),
                X(xd,yd,z) = 1.5;
                SLd(Id,Jd,zd) = SLd(Id,Jd,zd)+betaL/dVc;
                live = 0;
            end
            if (X(xd,yd,z) == 3)&&(cd < d3e),
                X(xd,yd,z) = 2.5;
                live = 0;
            end
			%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % Cell division and rearrangement
            if (X(xd,yd,z) == 1),
                Ti = T1;
                alpha = alphaL;
                celltype = 1;
            end
            if (X(xd,yd,z) == 2),
                Ti = T2;
                alpha = alphaA;
                celltype = 2;
            end
            if (X(xd,yd,z) == 3),
                Ti = T3;
                alpha = alphaL;
                celltype = 3;
            end
            if (live==1)&&(U(xd,yd,z) >= alpha)&&(UG(xd,yd,z) >= alphaG),
                Budded = 0;
                Natt = 0;
                U(xd,yd,z) = U(xd,yd,z)-alpha;
                UG(xd,yd,z) = UG(xd,yd,z)-alphaG;
                T(xd,yd,z) = tau;

                [qm1,qi1] = sort(rand(1,4)); % random index for neighboring grids
                [qm2,qi2] = sort(rand(1,4)); % random index for neighboring grids
                qi = [qi1 4+qi2];
                while (Budded==0)&&(Natt<8), % try immediate neighboring grids
                    Natt = Natt+1;
                    xb = I(cc) + Q(qi(Natt),1);
                    if xb > Ncx, xb = xb-Ncx; end
                    if xb < 1, xb = xb+Ncx; end
                    yb = J(cc) + Q(qi(Natt),2);
                    if yb > Ncy, yb = yb-Ncy; end
                    if yb < 1, yb = yb+Ncy; end

                    if X(xb,yb,z)==0, % if available space, bud into it
                        zb = z;
                        Budded = 1;
                        X(xb,yb,zb) = X(xd,yd,z);
                        T(xb,yb,zb) = tau;
                        U(xb,yb,zb) = 0;
                        UG(xb,yb,zb) = 0;
                    end
                end
                if (Budded==0), % try extended neighborhood in the same plane
                    rg = rg0;
                    xdp = xd;
                    ydp = yd;
                    if xd+rg > Ncx,
                        Xp = [X(:,:,z); X(1:xd+rg-Ncx,:,z)];
                    else
                        if xd-rg < 1,
                            Xp = [X(Ncx+xd-rg:Ncx,:,z); X(:,:,z)];
                            xdp = rg+1;
                        else
                            Xp = X(:,:,z);
                        end
                    end
                    if yd+rg > Ncy,
                        Xp = [Xp, Xp(:,1:yd+rg-Ncy)];
                    else
                        if yd-rg < 1,
                            Xp = [Xp(:,Ncy+yd-rg:Ncy), Xp];
                            ydp = rg+1;
                        end
                    end
                    Ng = Xp(xdp-rg:xdp+rg,ydp-rg:ydp+rg);
                    if sum(sum(abs(Ng)>0.25))<(2*rg+1)^2,
                        [Ie,Je] = find(Ng==0);
                        [SSe,OSe] = min(sqrt((Ie-rg-1).^2+(Je-rg-1).^2));
                        if SSe <= rg,
                            xe = Ie(OSe);
                            ye = Je(OSe);
                            TP = TracePath(rg+1,rg+1,xe,ye);
                            NT = size(TP,1);
                            TP(:,1) = TP(:,1) + xd-rg-1;
                            TP(:,1) = TP(:,1) + Ncx*((TP(:,1)<1)-(TP(:,1)>Ncx));
                            TP(:,2) = TP(:,2) + yd-rg-1;
                            TP(:,2) = TP(:,2) + Ncy*((TP(:,2)<1)-(TP(:,2)>Ncy));
                            for ctr = NT-1:-1:2,
                                X(TP(ctr,1),TP(ctr,2),z) = X(TP(ctr-1,1),TP(ctr-1,2),z);
                                U(TP(ctr,1),TP(ctr,2),z) = U(TP(ctr-1,1),TP(ctr-1,2),z);
                                UG(TP(ctr,1),TP(ctr,2),z) = UG(TP(ctr-1,1),TP(ctr-1,2),z);
                                T(TP(ctr,1),TP(ctr,2),z) = T(TP(ctr-1,1),TP(ctr-1,2),z);
                            end
                            zb = z;
                            X(TP(NT,1),TP(NT,2),zb) = X(TP(NT-1,1),TP(NT-1,2),z);
                            U(TP(NT,1),TP(NT,2),zb) = U(TP(NT-1,1),TP(NT-1,2),z);
                            UG(TP(NT,1),TP(NT,2),zb) = UG(TP(NT-1,1),TP(NT-1,2),z);
                            T(TP(NT,1),TP(NT,2),zb) = T(TP(NT-1,1),TP(NT-1,2),z);
                            X(TP(1,1),TP(1,2),z) = X(xd,yd,z);
                            U(TP(1,1),TP(1,2),z) = 0;
                            UG(TP(1,1),TP(1,2),z) = 0;
                            T(TP(1,1),TP(1,2),z) = tau;
                            Budded = 1;
                        end
                    end
                end
                if (Budded == 0), % bud to top or sides and push cells above
                    cq = rand(1);
                    if cq < 0.7, % bud straight up with 70% chance
                        X(xd,yd,z+1:Ncz) = X(xd,yd,z:Ncz-1);
                        U(xd,yd,z+1:Ncz) = U(xd,yd,z:Ncz-1);
                        UG(xd,yd,z+1:Ncz) = UG(xd,yd,z:Ncz-1);
                    else % bud to one of the sides and push cells above the new bud upward
                        [qm3,qi3] = sort(rand(1,4)); % random index for neighboring grids
                        xb = I(cc) + Qc(qi3(1),1);
                        if xb > Ncx, xb = xb-Ncx; end
                        if xb < 1, xb = xb+Ncx; end
                        yb = J(cc) + Qc(qi3(1),2);
                        if yb > Ncy, yb = yb-Ncy; end
                        if yb < 1, yb = yb+Ncy; end
                        X(xb,yb,z+1:Ncz) = X(xb,yb,z:Ncz-1);
                        U(xb,yb,z+1:Ncz) = U(xb,yb,z:Ncz-1);
                        UG(xb,yb,z+1:Ncz) = UG(xb,yb,z:Ncz-1);
                        X(xb,yb,z) = celltype;
                        U(xb,yb,z) = 0;
                        UG(xb,yb,z) = 0;
                        xd = xb;
                        yd = yb;
                    end
                    zt = Ncz;
                    while (zt>z)&&(X(xd,yd,zt)==0),
                        zt = zt-1;
                    end
                    Xe = [X(Ncx,:,zt-1); X(:,:,zt-1); X(1,:,zt-1)];
                    Xe = [Xe(:,Ncy), Xe, Xe(:,1)];
                    if sum(sum(Xe(xd:xd+2,yd:yd+2)>0.1))<=6,
                        [xe,ye] = find(Xe(xd:xd+2,yd:yd+2)==0);
                        et = floor(rand(1)*length(xe))+1;
                        xet = xd + xe(et) - 2;
                        xet = xet + Ncx*((xet<1)-(xet>Ncx));
                        yet = yd + ye(et) - 2;
                        yet = yet + Ncy*((yet<1)-(yet>Ncy));
                        zn = zt-1;
                        while (zn>=2)&&(X(xet,yet,zn-1)==0),
                            zn = zn-1;
                        end
                        X(xet,yet,zn) = X(xd,yd,zt);
                        U(xet,yet,zn) = U(xd,yd,zt);
                        UG(xet,yet,zn) = UG(xd,yd,zt);
                        X(xd,yd,zt) = 0;
                        U(xd,yd,zt) = 0;
                        UG(xd,yd,zt) = 0;
                    end
                end
            end
        end
        X1lm(z,ct) = sum(sum(X(:,:,z)==1));
        X2lm(z,ct) = sum(sum(X(:,:,z)==2));
        X3lm(z,ct) = sum(sum(X(:,:,z)==3));
        X1m(z,ct) = X1lm(z,ct) + sum(sum(X(:,:,z)==0.5));
        X2m(z,ct) = X2lm(z,ct) + sum(sum(X(:,:,z)==1.5));
        X3m(z,ct) = X3lm(z,ct) + sum(sum(X(:,:,z)==2.5));
    end
    UAccL(ct) = sum(sum(sum(U.*((X==1)+(X==3)))));
    UAccA(ct) = sum(sum(sum(U.*(X==2))));
    SLAccC(ct) = sum(sum(sum(SLd)));
    SLAccA(ct) = sum(sum(sum(SLa)));
    SAAccC(ct) = sum(sum(sum(SAd)));
    SAAccA(ct) = sum(sum(sum(SAa)));
    disp([tau  sum(X1lm(:,ct)) sum(X2lm(:,ct)) sum(X3lm(:,ct)) sum(X1m(:,ct)+X2m(:,ct)+X3m(:,ct)) sum(X2m(:,ct)-X2lm(:,ct))*betaL/alphaL sum(X1m(:,ct)+X3m(:,ct))-N1-N3])

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Update diffusion constant
    Dd = zeros(Ndx,Ndy,Ndz);
    for ii = 1:SC,
        for jj = 1:SC,
            for kk = 1:SC,
                Dd = Dd + D1/SC^3*(X(ii:SC:Ncx,jj:SC:Ncy,kk:SC:Ncz)>0.1);
            end
        end
    end
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	
	dt0 = min(0.15*f^2/D0,0.15*g^2/D1); % time step for diffusion simulation, in seconds
    trng = linspace(0,3600*dtau-dt0,round(3600*dtau/dt0)); % in seconds
    dt = trng(2)-trng(1);

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% Diffusion loop, simulate the diffusion of nutrients within each time step for cell activities
    for t = trng,
        % Update concentration at the boundary, SLb
        SLam = [SLa(Nax,:,Naz); SLa(:,:,Naz); SLa(1,:,Naz)];
        SLam = [SLam(:,Nay), SLam, SLam(:,1)];
        SLae = interp2(SD*(0:Nax+1),SD*(0:Nay+1),SLam',xxd,yyd,'linear')';
        SLb = SLb - dt/h * (D0/f*(SLb - SLae) + 1/g*Dd(:,:,1).*(SLb - SLd(:,:,1)));

        % Diffusion in the agar region
        SLbm = zeros(Nax,Nay);
        for ii = 1:SD,
            for jj = 1:SD,
                SLbm = SLbm + 1/SD^2*SLb(ii:SD:Ndx,jj:SD:Ndy);
            end
        end
        SLab = SLa(:,:,3:Naz);
        SLab(:,:,Naz-1) = SLbm;
        dSLa = dt/f^2 * D0*([SLa(Nax,:,2:Naz); SLa(1:Nax-1,:,2:Naz)] + [SLa(2:Nax,:,2:Naz); SLa(1,:,2:Naz)] + ...
            [SLa(:,Nay,2:Naz), SLa(:,1:Nay-1,2:Naz)] + [SLa(:,2:Nay,2:Naz), SLa(:,1,2:Naz)] + ...
            SLa(:,:,1:Naz-1) + SLab - 6*SLa(:,:,2:Naz));
        SLa(:,:,2:Naz) = SLa(:,:,2:Naz) + dSLa;
        SLa(:,:,1) = SLa(:,:,2);

        % Diffusion in the community region
        SLdb = zeros(Ndx,Ndy,Ndz-1);
        SLdb(:,:,1) = SLb;
        SLdb(:,:,2:Ndz-1) = SLd(:,:,1:Ndz-2);
        dSLd = dt/g^2 * Dd(:,:,1:Ndz-1).*([SLd(Ndx,:,1:Ndz-1); SLd(1:Ndx-1,:,1:Ndz-1)] + ...
            [SLd(2:Ndx,:,1:Ndz-1); SLd(1,:,1:Ndz-1)] + [SLd(:,Ndy,1:Ndz-1), SLd(:,1:Ndy-1,1:Ndz-1)] + ...
            [SLd(:,2:Ndy,1:Ndz-1), SLd(:,1,1:Ndz-1)] + SLdb + SLd(:,:,2:Ndz) - 6*SLd(:,:,1:Ndz-1)) + dt/g^2 * (...
            ([Dd(2:Ndx,:,1:Ndz-1); Dd(1,:,1:Ndz-1)]-Dd(:,:,1:Ndz-1)).*([SLd(2:Ndx,:,1:Ndz-1); SLd(1,:,1:Ndz-1)]-SLd(:,:,1:Ndz-1))+...
            ([Dd(:,2:Ndy,1:Ndz-1), Dd(:,1,1:Ndz-1)]-Dd(:,:,1:Ndz-1)).*([SLd(:,2:Ndy,1:Ndz-1), SLd(:,1,1:Ndz-1)]-SLd(:,:,1:Ndz-1))+...
            (Dd(:,:,2:Ndz)-Dd(:,:,1:Ndz-1)).*(SLd(:,:,2:Ndz)-SLd(:,:,1:Ndz-1)));
        SLd(:,:,1:Ndz-1) = SLd(:,:,1:Ndz-1) + dSLd;
        SLd(:,:,Ndz) = SLd(:,:,Ndz-1);

        % Update concentration at the boundary, SAb
        SAam = [SAa(Nax,:,Naz); SAa(:,:,Naz); SAa(1,:,Naz)];
        SAam = [SAam(:,Nay), SAam, SAam(:,1)];
        SAae = interp2(SD*(0:Nax+1),SD*(0:Nay+1),SAam',xxd,yyd,'linear')';
        SAb = SAb - dt/h * (D0/f*(SAb - SAae) + 1/g*Dd(:,:,1).*(SAb - SAd(:,:,1)));

        % Diffusion in the agar region
        SAbm = zeros(Nax,Nay);
        for ii = 1:SD,
            for jj = 1:SD,
                SAbm = SAbm + 1/SD^2*SAb(ii:SD:Ndx,jj:SD:Ndy);
            end
        end
        SAab = SAa(:,:,3:Naz);
        SAab(:,:,Naz-1) = SAbm;
        dSAa = dt/f^2 * D0*([SAa(Nax,:,2:Naz); SAa(1:Nax-1,:,2:Naz)] + [SAa(2:Nax,:,2:Naz); SAa(1,:,2:Naz)] + ...
            [SAa(:,Nay,2:Naz), SAa(:,1:Nay-1,2:Naz)] + [SAa(:,2:Nay,2:Naz), SAa(:,1,2:Naz)] + ...
            SAa(:,:,1:Naz-1) + SAab - 6*SAa(:,:,2:Naz));
        SAa(:,:,2:Naz) = SAa(:,:,2:Naz) + dSAa;
        SAa(:,:,1) = SAa(:,:,2);

        % Diffusion in the community region
        SAdb = zeros(Ndx,Ndy,Ndz-1);
        SAdb(:,:,1) = SAb;
        SAdb(:,:,2:Ndz-1) = SAd(:,:,1:Ndz-2);
        dSAd = dt/g^2 * Dd(:,:,1:Ndz-1).*([SAd(Ndx,:,1:Ndz-1); SAd(1:Ndx-1,:,1:Ndz-1)] + ...
            [SAd(2:Ndx,:,1:Ndz-1); SAd(1,:,1:Ndz-1)] + [SAd(:,Ndy,1:Ndz-1), SAd(:,1:Ndy-1,1:Ndz-1)] + ...
            [SAd(:,2:Ndy,1:Ndz-1), SAd(:,1,1:Ndz-1)] + SAdb + SAd(:,:,2:Ndz) - 6*SAd(:,:,1:Ndz-1)) + dt/g^2 * (...
            ([Dd(2:Ndx,:,1:Ndz-1); Dd(1,:,1:Ndz-1)]-Dd(:,:,1:Ndz-1)).*([SAd(2:Ndx,:,1:Ndz-1); SAd(1,:,1:Ndz-1)]-SAd(:,:,1:Ndz-1))+...
            ([Dd(:,2:Ndy,1:Ndz-1), Dd(:,1,1:Ndz-1)]-Dd(:,:,1:Ndz-1)).*([SAd(:,2:Ndy,1:Ndz-1), SAd(:,1,1:Ndz-1)]-SAd(:,:,1:Ndz-1))+...
            (Dd(:,:,2:Ndz)-Dd(:,:,1:Ndz-1)).*(SAd(:,:,2:Ndz)-SAd(:,:,1:Ndz-1)));
        SAd(:,:,1:Ndz-1) = SAd(:,:,1:Ndz-1) + dSAd;
        SAd(:,:,Ndz) = SAd(:,:,Ndz-1);

        % Update concentration at the boundary, SGb
        SGam = [SGa(Nax,:,Naz); SGa(:,:,Naz); SGa(1,:,Naz)];
        SGam = [SGam(:,Nay), SGam, SGam(:,1)];
        SGae = interp2(SD*(0:Nax+1),SD*(0:Nay+1),SGam',xxd,yyd,'linear')';
        SGb = SGb - dt/h * (D0/f*(SGb - SGae) + 1/g*Dd(:,:,1).*(SGb - SGd(:,:,1)));

        % Diffusion in the agar region
        SGbm = zeros(Nax,Nay);
        for ii = 1:SD,
            for jj = 1:SD,
                SGbm = SGbm + 1/SD^2*SGb(ii:SD:Ndx,jj:SD:Ndy);
            end
        end
        SGab = SGa(:,:,3:Naz);
        SGab(:,:,Naz-1) = SGbm;
        dSGa = dt/f^2 * D0*([SGa(Nax,:,2:Naz); SGa(1:Nax-1,:,2:Naz)] + [SGa(2:Nax,:,2:Naz); SGa(1,:,2:Naz)] + ...
            [SGa(:,Nay,2:Naz), SGa(:,1:Nay-1,2:Naz)] + [SGa(:,2:Nay,2:Naz), SGa(:,1,2:Naz)] + ...
            SGa(:,:,1:Naz-1) + SGab - 6*SGa(:,:,2:Naz));
        SGa(:,:,2:Naz) = SGa(:,:,2:Naz) + dSGa;
        SGa(:,:,1) = SGa(:,:,2);

        % Diffusion in the community region
        SGdb = zeros(Ndx,Ndy,Ndz-1);
        SGdb(:,:,1) = SGb;
        SGdb(:,:,2:Ndz-1) = SGd(:,:,1:Ndz-2);
        dSGd = dt/g^2 * Dd(:,:,1:Ndz-1).*([SGd(Ndx,:,1:Ndz-1); SGd(1:Ndx-1,:,1:Ndz-1)] + ...
            [SGd(2:Ndx,:,1:Ndz-1); SGd(1,:,1:Ndz-1)] + [SGd(:,Ndy,1:Ndz-1), SGd(:,1:Ndy-1,1:Ndz-1)] + ...
            [SGd(:,2:Ndy,1:Ndz-1), SGd(:,1,1:Ndz-1)] + SGdb + SGd(:,:,2:Ndz) - 6*SGd(:,:,1:Ndz-1)) + dt/g^2 * (...
            ([Dd(2:Ndx,:,1:Ndz-1); Dd(1,:,1:Ndz-1)]-Dd(:,:,1:Ndz-1)).*([SGd(2:Ndx,:,1:Ndz-1); SGd(1,:,1:Ndz-1)]-SGd(:,:,1:Ndz-1))+...
            ([Dd(:,2:Ndy,1:Ndz-1), Dd(:,1,1:Ndz-1)]-Dd(:,:,1:Ndz-1)).*([SGd(:,2:Ndy,1:Ndz-1), SGd(:,1,1:Ndz-1)]-SGd(:,:,1:Ndz-1))+...
            (Dd(:,:,2:Ndz)-Dd(:,:,1:Ndz-1)).*(SGd(:,:,2:Ndz)-SGd(:,:,1:Ndz-1)));
        SGd(:,:,1:Ndz-1) = SGd(:,:,1:Ndz-1) + dSGd;
        SGd(:,:,Ndz) = SGd(:,:,Ndz-1);

        %% Uptake
        dtu = dt/3600; % in hours
        % Nutrient uptake
        for i1 = 1:SC,
            for i2 = 1:SC,
                for i3 = 1:SC,
                    SLde(i1:SC:Ncx,i2:SC:Ncy,i3:SC:Ncz) = SLd;
                    SAde(i1:SC:Ncx,i2:SC:Ncy,i3:SC:Ncz) = SAd + 1/dV*gammaA*dtu*(X(i1:SC:Ncx,i2:SC:Ncy,i3:SC:Ncz)==1);
                    SGde(i1:SC:Ncx,i2:SC:Ncy,i3:SC:Ncz) = SGd;
                end
            end
        end

        UcL = min(dV*SLde,(dtu*(vmL1*(X==1) + vmL3*(X==3)).*SLde./(SLde+KL))).*(SLde>0).*(U<=1.2*alphaL); % Lysine uptake
        UcA = (X==2).* min(dV*SAde,(dtu*vmA*SAde./(SAde+KA))).*(SAde>0).*(U<=1.2*alphaA); % Adenine uptake
        UcG = ((X==1)+(X==2)+(X==3)).* min(dV*SGde,(dtu*vmG*SGde./(SGde+KG))).*(SGde>0).*(UG<=1.2*alphaG); % Glucose uptake
        U = U+UcL+UcA;
        UG = UG+UcG;
        SLde = SLde-1/dV*UcL;
        SAde = SAde-1/dV*UcA;
        SGde = SGde-1/dV*UcG;

        SLd = zeros(Ndx,Ndy,Ndz);
        SAd = zeros(Ndx,Ndy,Ndz);
        SGd = zeros(Ndx,Ndy,Ndz);
        for i1 = 1:SC,
            for i2 = 1:SC,
                for i3 = 1:SC,
                    SLd = SLd + 1/SC^3*SLde(i1:SC:Ncx,i2:SC:Ncy,i3:SC:Ncz);
                    SAd = SAd + 1/SC^3*SAde(i1:SC:Ncx,i2:SC:Ncy,i3:SC:Ncz);
                    SGd = SGd + 1/SC^3*SGde(i1:SC:Ncx,i2:SC:Ncy,i3:SC:Ncz);
                end
            end
        end

    end
    SPLa(:,ct) = shiftdim(mean(mean(SLa,2),1));
    SPLd(:,ct) = shiftdim(mean(mean(SLd,2),1));
    SPAa(:,ct) = shiftdim(mean(mean(SAa,2),1));
    SPAd(:,ct) = shiftdim(mean(mean(SAd,2),1));
    SPGa(:,ct) = shiftdim(mean(mean(SGa,2),1));
    SPGd(:,ct) = shiftdim(mean(mean(SGd,2),1));

    if mod(tau+0.03,ts)<0.1,
        cS = cS+1;
        XS1(:,:,cS) = sum((X>0.25).*(X<1.25),3);
        XS2(:,:,cS) = sum((X>1.25).*(X<2.25),3);
        XS3(:,:,cS) = sum((X>2.25).*(X<3.25),3);
        ss = round(Ncx/2);
        for ii = 1:Ncy,
            for jj = 1:Ncz,
                Sec(ii,jj,cS) = X(ss,ii,jj);
            end
        end
        save(strcat(savename,'_ts2_cs',num2str(cS)),'X','tau')
    end
    toc
    
end

figure
semilogy(linspace(0,tau-dtau,ct-1),sum(X1m(:,1:ct-1)),'r')
hold on
semilogy(linspace(0,tau-dtau,ct-1),sum(X2m(:,1:ct-1)),'g')
semilogy(linspace(0,tau-dtau,ct-1),sum(X3m(:,1:ct-1)),'b')
semilogy(linspace(0,tau-dtau,ct-1),sum(X1lm(:,1:ct-1)),'r:')
semilogy(linspace(0,tau-dtau,ct-1),sum(X2lm(:,1:ct-1)),'g:')
semilogy(linspace(0,tau-dtau,ct-1),sum(X3lm(:,1:ct-1)),'b:')
xlabel('Time (hours)')
ylabel('Population (# of cells)')
xlim([0 tau-dtau])


save(savename)
