function [Output,paras] = MultiConMatrix_L(x,y,a,paras,option,type,gap,group)
%   Impedance Matrix Calculation for Round Conductors with Boundaries
%
%   Input variables
%   x           X coordinates of centre points m
%   y           Y coordinates of centre points m
%   a           Radius m
%
%   para.sigma  Conductivity S/m
%   para.freq   Frequency list
%   para.bx     Coordinates of vertical boundary m
%   para.by     Coordinates of horizental boundary m
%   para.b      Outer radius m
%   para.V      Voltage V
%   para.I      Current A
%   para.Hext   External magnetic field A/m
%
%   option.r0   Reference point m
%   option.Nord Truncated order, default N = 3
%   option.k    Reflection coefficient, default k = 1
%   option.Nref Reflection time, defalt Nref = 2
%   option.a0   Boundary potential
%
%   type.result Type of calculation
%   type.correct Flag for boundary compensation
%   type.ExtH   Flag for external magnetic field
%   type.gap    Flag for air gap
%
%   gap.xg      Gap's centre X coordinate
%   gap.yg      Gap's centre Y coordinate
%   gap.len     Length of gap
%   gap.dire    Direction of gap
%   gap.kgap    Coefficient to decide anti source strength
%   gap.gtype   Type of anti source
%
%   Output variable
%   Output:     Output based on input
%   paras       all parameters in the calculation
%
%   Reference:
%   ...
%
%   Authors:
%   Tianmingluo
%   Delft University of Technology
%
%   Email:
%   T.Luo-1@tudelft.nl
%
%   Version 1.0.238

arguments
    x  (:,1) double  % x coordinates of centre points
    y  (:,1) double  % y coordinates of centre points
    a  (:,1) double  % radius of conductors
    
    paras.bx        (1,:) double =[];
    paras.by        (1,:) double =[];
    paras.sigma     (:,1) double {mustBeNumeric} = 5.96e7; % copper conductivity
    paras.freq      (1,:) double {mustBeNumeric} = 1; % frequency list
    paras.b         (:,1) double = [];
    paras.mur       (:,1) double = 1; % permeability of conductor
    paras.V         (:,1) double =[];
    paras.I         (:,1) double =[];
    paras.Hext      (1,2) double = [0,0];
    
    option.r0       double {mustBeNumeric} = 1;
    option.Nord     double {mustBeNumeric} = 3;
    option.k        double {mustBeNumeric} = 1;
    option.Nref     double {mustBeNumeric} = 2;
    option.a0       double {mustBeNumeric} = 0;
    
    type.result     string {mustBeMember(type.result,["case","matrix"])} = "case";
    type.correct    logical = false;
    type.ExtH       logical = false;
    type.gap        logical = false;
    type.group      logical = false;
    
    gap.xg          (:,1) double = [];
    gap.yg          (:,1) double = [];
    gap.len         (:,1) double = [];
    gap.dire        (:,1) char   = [];
    gap.kgap        (1,1) double = -1;
    gap.gtype        string {mustBeMember(gap.gtype,["points","lines"])} = "points";
    
    group.gplist    (:,2) single = repmat((1:length(x))',1,2);
end
%% Constant define
mu0 = 4e-7*pi; % unit H/m

%% Input check
flag = 0;

warnid = 'MATLAB:nearlySingularMatrix';
warning('off',warnid);

num = length(x);

if ~isequal(length(y),num)
    flag = 1;
    errorinfo = "Coordinates must have the same length as x";
end
if ~isequal(length(a),num)
    flag = 1;
    errorinfo = "Radius must have the same length as x";
end

if isempty(paras.freq)
    flag = 1;
    errorinfo = "Found empty freqeuncy array";
end
if length(paras.mur) == 1
    paras.mur = ones(num,1)*paras.mur;
end
if length(paras.sigma) ~= 1
    if length(paras.sigma) ~= num
        flag = 1;
        errorinfo = "Please assign conductivity to each conductors or the same conductivity";
    end
end
if type.result == "case"
    if isempty(paras.I)&&isempty(paras.V)
        flag = 1;
        errorinfo = "Found empty current or voltage array";
    elseif ~isequal(length(paras.I),num)&&~isequal(length(paras.V),num)
        flag = 1;
        errorinfo = "Current or voltage must have the same length as x";
    elseif ~isempty(paras.I)&&~isempty(paras.V)
        flag = 1;
        errorinfo = "Found double input array";
    end
end


if type.gap
    if isempty(gap.xg) || isempty(gap.yg)
        type.gap = 0;
    else
        if length(gap.xg) ~= length(gap.yg)
            flag = 1;
            errorinfo = "yg must have the same length as xg";
        end
        if gap.gtype == "lines"
            if length(gap.len) ~= length(gap.xg)
                flag = 1;
                errorinfo = "len must have the same length as xg";
            end
            if length(gap.dire) ~= length(gap.xg)
                flag = 1;
                errorinfo = "dire must have the same length as xg";
            end
            if nnz(gap.dire~='y' & gap.dire~='x')
                flag = 1;
                errorinfo = "dire must be 'x' or 'y'";
            end
        end
    end
end

if type.correct
    if isempty(paras.bx) && isempty(paras.by)
        type.correct = false;
    end
    if option.k == -1 && numel(paras.bx)<=1 && numel(paras.by)<=1
        type.correct = false;
    end
end

if type.group
    group.gplist = array2table(group.gplist);
    group.gplist.Properties.VariableNames={'NoC','NoG'};
    if  max(group.gplist.NoC)>num
        flag = 1;
        errorinfo = "the No. of conductor exceed the existing conductor";
    end
    group.gplist.NoG= categorical(group.gplist.NoG);
    gpkind = categories(group.gplist.NoG);
    gpcount = countcats(group.gplist.NoG);
    idcount = gpcount ~= 1;
    if ~nnz(idcount)
        type.group = 0;
    end
    iddel = ismember(group.gplist.NoG,gpkind(~idcount));
    group.gplist(iddel,:) = [];
    group.gplist.NoG = removecats(group.gplist.NoG,gpkind(~idcount));
    gpkind = gpkind(idcount);
end

%%
if flag == 0
    
    lnr0 = log(option.r0);
    N = option.Nord;
    k = option.k;
    Nref = option.Nref;
    a0 = option.a0;
    bx = paras.bx;
    by = paras.by;
    sigma = paras.sigma;
    w = 2*pi*paras.freq;
    units = 2*N+1;
    n1 = 1:N;
    totlen = units*num;
    
    [matA0,matB0,countD] = mat_basic_v2(x,y,bx,by,N,Nref,k);
    
    if type.gap
        numg = numel(gap.xg);
        [matB2C,countD2C,matA2g,matB2g] = mat_gap_v2(gap.xg,gap.yg,x,y,bx,by,N,Nref,k,gap.len,gap.dire,gap.gtype);
        matR2C = zeros(totlen,numg);
        matR2C(1:units:totlen,:)= countD2C;
        matR2g = countD2C';
    end
    
    matR0 = zeros(totlen,num);
    matR0(1:units:totlen,:)= countD;
    
    if type.group % disable correction when group is used
        type.correct = false;
    end
    
    if type.correct
        [mark,bound,ind] = point_comp(x,y,bx,by,type.gap,gap.xg,gap.yg);
    else
        matB0 = matB0-matR0*lnr0;
        if type.gap
            matB2C = matB2C-matR2C*lnr0;
            matB2g = matB2g-matR2g*lnr0;
        end
    end
    
    
    %%
    
    mur = paras.mur;
    num2 = length(w);
    
    k2 = sqrt(-1j*sigma*mu0*w);
    id0 = k2 == 0;
    k2(id0) = NaN;
    
    [matf,matbessel] = matf_L(a,mur,N,k2);
    munit = eye(totlen);
    matA0 = matA0.*matf+munit;
    
    [matVc,matVd,matVr,aveAp2] = matV_L(mur,a,k2,w,N,sigma);
    
    if type.correct
        [Acomp,Bcomp] = matComp(mark,ind,bound,x,y,matf(1,(ind-1)*units+1:ind*units,:),totlen);
        %fprintf('Compensation for Boundary Potential Applied. \n');
    end
    
    idC = false(totlen,1);
    idC(1:units:totlen)=true;
    coef = repmat(n1,1,2);
    if type.gap
        if gap.gtype == "lines"
            input_len = sum(gap.len);
        else
            input_len = numg;
        end
    end
    
    if type.result == "case"
        paras.FactorAB = zeros(totlen,1,num2);
        paras.aveA = zeros(num,1,num2);
        paras.Poyn_in = zeros(num,1,num2);
        paras.Poyn_out = zeros(num,1,num2);
        
        if ~isempty(paras.I)
            vorq = false;
            paras.V = zeros(num,1,num2);
        elseif ~isempty(paras.V)
            vorq = true;
            paras.I = zeros(num,1,num2);
            if type.correct
                type.correct = false;
                %fprintf('Potential Correction only Work with Current. \n');
            end
        end
        nloop = 1;
    else
        paras.FactorAB = zeros(totlen,num,num2);
        paras.aveA = zeros(num,num,num2);
        paras.Poyn_in = zeros(num,num,num2);
        paras.Poyn_out = zeros(num,num,num2);
        paras.V = zeros(num,num,num2);
        paras.I = eye(num);
        nloop = num;
        vorq = false;
    end
    
    if type.gap
        if ~vorq
            matG = matG_L(matA2g,matf(1,:,:),vorq);
        end
    end
    
    for idx = 1:nloop
        if ~vorq
            I = paras.I(:,idx);
            D = -mu0*I/2/pi;
            if type.group
                [matA,matB] = matGroup(matA0,matB0,D,matVc,matVd-matVr*lnr0,group.gplist);
                paras.I = zeros(num,nloop,num2);
            else
                matB = matB0*D;
                if type.gap % average allocated the courter source
                    Dgap = sum(D)*gap.kgap*ones(numel(gap.xg),1)/input_len;
                    matBg = matB2C*Dgap;
                    matB = matB + matBg;
                end
                if type.correct
                    matB = [matB;a0+D(ind)*Bcomp];
                    if type.gap
                        matR = matR0*D+matR2C*Dgap;
                    else
                        matR = matR0*D;
                    end
                    matR = repmat(matR,1,1,num2);
                    matA = [matA0,matR;Acomp,repmat(-D(ind),1,1,num2)];
                else
                    matA = matA0;
                end
            end
        else
            matA = [matA0,-repmat(matB0,1,1,num2);matVc,matVd-matVr*lnr0];
            matB = [zeros(totlen,1);paras.V];
            
            if type.gap
                [matG,rowD,columC] = matG_L(matA2g,matf(1,:,:),vorq,gap.kgap,input_len,num);
                matA2C = [-matB2C;zeros(num,numg);-eye(numg);zeros(numg)];
                matA = [[matA;rowD;[matG,repmat(matB2g,1,1,num2)]],repmat(matA2C,1,1,num2),columC];
                matB = [matB;zeros(2*numg,1)];
            end
        end
        
        %--- Adding External field strenght ---
        if type.ExtH
            matH = mat_ext(mu0,paras.Hext(1),paras.Hext(2),N,num,x,y);
            if numel(matB) > numel(matH)
                matH(numel(matB))=0;
            end
            matB = matB+matH;
            %fprintf('Uniform External Magnetic Field Applied. \n');
        end
        %--- Part end ---
        
        for idx3 = 1:num2
            temp_result = matA(:,:,idx3)\matB;
            if vorq
                if ~type.gap
                    D(:,idx3) = temp_result(totlen+1:end);
                else
                    D(:,idx3) = temp_result(totlen+1:totlen+num);
                    Cg0(:,idx3) = temp_result(totlen+num+numg+1:end);
                end
            else
                if type.correct
                    lnr0(1,idx3) = real(temp_result(end));
                    if ~isreal(temp_result(end))
                        temp_result(1:totlen) = temp_result(1:totlen)+1i*imag(temp_result(end))*matR(:,:,idx3);
                    end
                end
                if type.group
                    D(group.gplist.NoC,idx3) = temp_result(totlen+1:end);
                end
            end
            paras.FactorAB(:,idx,idx3) = temp_result(1:totlen);
        end
        
        C = permute(paras.FactorAB(idC,idx,:),[1,3,2]);
        AB = squeeze(paras.FactorAB(~idC,idx,:));
        
        AB = permute(reshape(AB,2*N,[],num2),[2,1,3]);
        
        if vorq || type.group
            I = -2*pi/mu0*D;
        end
        
        if type.gap
            if ~vorq
                Cg0 = sum(matG.*permute(paras.FactorAB(:,idx,:),[2,1,3]),2)+matB2g*D;
                Cg0 = permute(Cg0,[1,3,2]);
            end
            gap.Cg = sum(Cg0,1)/input_len;
            normI = vecnorm(I,2,1);
            Izero = normI == 0;
            normI(Izero) = 1;
            
            % loss part allocation is based on the current value
            % and phase
            % ! allocation need to recheck %
            C = C + gap.kgap * I.*(real(I)*1i.*imag(gap.Cg)-1i*imag(I).*real(gap.Cg))./normI.^2;
        end
        
        paras.aveA(:,idx,:) = C+D.*(log(a)-lnr0 - aveAp2);
        
        if ~vorq || type.gap
            paras.V(:,idx,:) = I./sigma/pi./a.^2+1j*w.*permute(paras.aveA(:,idx,:),[1,3,2]);
            if type.gap && vorq
                fprintf('Please note the input voltage is changed. \n');
            end
        end
        
        paras.Poyn_in(:,idx,:) = I.*conj(I).*(1./sigma/pi./a.^2+1j*w*mu0/2/pi.*aveAp2);
        
        Pp2 = coef.*a.^(2*coef).*(1+matbessel).*(1-conj(matbessel)).*AB.*conj(AB);
        paras.Poyn_out(:,idx,:) =  1j*permute(w,[1,3,2])*pi/mu0.*sum(Pp2 ,2);
        
        if nnz(id0)
            if vorq
                paras.Poyn_in(:,idx,id0) =  I(:,id0).*conj(I(:,id0)).*(1./sigma/pi./a.^2+1j*mu0/2/pi*aveAp2(:,id0));
            else
                paras.Poyn_in(:,idx,id0) =  I.*conj(I).*(1./sigma/pi./a.^2+1j*mu0/2/pi*aveAp2(:,id0));
            end
            Pp3 =  coef.*a.^(2*coef).*(1+matbessel(:,:,id0)).^2.*AB(:,:,id0).^2;
            paras.Poyn_out(:,idx,id0) =  1j*pi/mu0./mur.*permute(sum(Pp3 ,2),[1,3,2]);
        end
        
        paras.Poyn = paras.Poyn_in+paras.Poyn_out;
        
        paras.I(:,idx,:) = permute(I,[1,3,2]);
    end
    
    if type.ExtH
        paras.P = real(paras.Poyn);
        paras.Q = [];
        paras.Wm = imag(paras.Poyn)./permute(w,[1,3,2])/2; % Energy stored in the conductors
        if nnz(id0)
            paras.Wm(:,id0) = imag(paras.Poyn(:,id0))/2;
        end
        %             fprintf('Reactive power is not applied under external field. \n');
    else
        if ~nnz(paras.I==0)
            Imp = paras.V./paras.I;
        else
            Imp = paras.V;
        end
        S = paras.V.*conj(paras.I);
        paras.Q = imag(S);
        paras.P = real(S);
        paras.Wm = imag(paras.Poyn)./permute(w,[1,3,2])/2;
        if nnz(id0)
            paras.Wm(:,id0) = imag(paras.Poyn(:,id0))/2;
        end
    end
    
    %%
    if type.result == "matrix"
        Output = Imp;
    else
        if exist('Imp','var')
            Output = Imp;
        else
            Output = paras.P;
        end
    end
    
    paras.x = x;
    paras.y = y;
    paras.a = a;
    paras.option = option;
    paras.type = type;
    paras.gap = gap;
else
    error(errorinfo);
end

end

function matH = mat_ext(mu,Hx,Hy,N,num,x,y)
units = 2*N+1;
Bext = zeros(units,1);
Bext(2) = -mu*Hy;
Bext(N+2) = mu*Hx;

matH = cell(num,1);
matH(:) = {Bext};
matH = cell2mat(matH);
for k = 1:num
    matH(1+(k-1)*units) = -mu*Hy*x(k)+mu*Hx*y(k);
end
end

function [matf,matbessel] = matf_L(a,mur,N,k2)
%--- Matrix converting A'' and B'' coefficients to A' and B' in L situation---

num = length(a);
[num_s,num2] = size(k2);
units = 2*N+1;
id0 = isnan(k2);
flag = 1;
matf = cell(num,num);
um = ones(units,1,num2);
matbessel = zeros(num,N,num2);
[matn,~,matk] = meshgrid((1:N)',(0:N)',k2(1,:));
if num_s ==1
    flag = 0;
end
for id = 1:num
    if flag
        matk = permute(repmat(k2(id,:),N+1,1,N),[1,3,2]);
    end
    if ~nnz(mur ~= 1)
        matka = matk*a(id);
        pbessel = besselj(matn+1,matka)./besselj(matn-1,matka);
        if nnz(id0)
            idnan = isnan(pbessel);
            pbessel(idnan) = 0;
        end
    else
        matka = matk*sqrt(mur(id))*a(id);
        part1 = besselj(matn+1,matka);
        part2 = besselj(matn-1,matka);
        pbessel = ((mur(id)+1).*part1+(mur(id)-1).*part2)./((mur(id)+1).*part2+(mur(id)-1).*part1);
        if nnz(id0)
            idnan = isnan(pbessel);
            pbessel(idnan) = (mur(id)-1)./(mur(id)+1);
        end
    end
    
    pf = a(id).^(2*matn).*pbessel;
    pf2 = pf(2:end,:,:);
    
    matf(:,id)= {[um,[repmat(pf,1,2);repmat(pf2,1,2)]]};
    matbessel(id,:,:) = pbessel(1,:,:);
    
end

matbessel = repmat(matbessel,1,2);
matf = cell2mat(matf);
end

function [matVc,matVd,matVr,aveAp2] = matV_L(mur,a,k2,w,N,sigma)
mu0 = 4e-7*pi;
[num_s,num2] = size(k2);
num = length(a);
units = 2*N+1;
flag = 1;
if num_s ==1
    flag = 0;
end
matVc = cell(num,num,num2);
matVc(:) = {zeros(1,units)};
for k = 1:num
    matVc(k,k,:) = {eye(1,units)};
end
matVc = cell2mat(matVc);
matVc = matVc .* permute(1j*w,[1,3,2]);

chi_dc = diag(-2./a.^2/mu0./sigma);
matVd = repmat(chi_dc,1,1,num2);

if ~nnz(mur ~= 1)
    if flag
        k2a = a.*k2;
    else
        k2a = a*k2;
    end
    chi = log(a)-besselj(2,k2a)./(k2a)./besselj(1,k2a);
else
    if flag
        k2a = sqrt(mur).*a.*k2;
    else
        k2a = sqrt(mur).*a*k2;
    end
    chi = log(a)-mur.*besselj(2,k2a)./(k2a)./besselj(1,k2a);
end

aveAp2 = log(a)-chi;
chi = chi.*(1j*w);

if nnz(isnan(k2))
    idnan = isnan(chi);
    chi(idnan) = 0;
    if ~nnz(mur ~= 1)
        aveAp2(idnan) = 1/4;
    else
        aveAp2(idnan) = mur/4;
    end
end
chi = permute(chi,[1,3,2]);

matVr = zeros(num,num,num2);
matr0 = permute(ones(num,1)*(1j*w),[1,3,2]);
for k = 1:num2
    matVd(:,:,k) = matVd(:,:,k)+ diag(chi(:,:,k));
    matVr(:,:,k) = diag(matr0(:,1,k));
end
end

function [matG,rowD,columC] = matG_L(matA2g,matf,vorq,k,len,num)
arguments
    matA2g
    matf
    vorq    = 0;
    k       = -1;
    len     = 1;
    num   	= 0;
end
[~,~,num2] = size(matf);
[numg,totlen] = size(matA2g);
if ~vorq
    matG = matA2g.*matf;
else
    columC = zeros(totlen+num+numg,numg);
    columC = [columC;-eye(numg)];
    columC = repmat(columC,1,1,num2);
    rowD = zeros(numg,totlen);
    rowD = [rowD, ones(numg,num)*k/len];
    rowD = repmat(rowD,1,1,num2);
    matG = matA2g.*matf;
end
end

function [mark,bound,ind] = point_comp(x,y,bx,by,flag,xg,yg)
%--- Select cell for compensation ---
arguments
    x
    y
    bx      (1,:) = [];
    by      (1,:) = [];
    flag    = false;
    xg      (1,:) = [];
    yg      (1,:) = [];
end
if ~isempty(bx)
    dx = abs(x-bx);
    valx = min(dx,[],'all');
    midx = mean(bx);
end
if ~isempty(by)
    dy = abs(y-by);
    valy = min(dy,[],'all');
    midy = mean(by);
end
if flag
    dxgap = min(abs(x-xg),[],2);
    dygap = min(abs(y-yg),[],2);
end

if ~isempty(bx)&&~isempty(by)
    if valx<valy
        mark = 'x';
        if ~flag
            [~,ind,bound_id] = process(dx,valx,abs(y-midy),1);
        else
            [~,ind,bound_id] = process(dx,valx,dygap,0);
        end
        bound = bx(bound_id);
    elseif valx>valy
        mark = 'y';
        if ~flag
            [~,ind,bound_id] = process(dy,valy,abs(x-midx),1);
        else
            [~,ind,bound_id] = process(dy,valy,dxgap,0);
        end
        bound = by(bound_id);
    else
        note = 0;
        if ~flag
            [valx2,indx,bound_idx] = process(dx,valx,abs(y-midy),1);
            [valy2,indy,bound_idy] = process(dy,valy,abs(x-midx),1);
            if valx2 <= valy2
                note = 1;
            end
        else
            [valx2,indx,bound_idx] = process(dx,valx,dygap,0);
            [valy2,indy,bound_idy] = process(dy,valy,dxgap,0);
            if valx2 >= valy2
                note = 1;
            end
        end
        if note
            mark = 'x';
            ind = indx;
            bound = bx(bound_idx);
        else
            mark = 'y';
            ind = indy;
            bound = by(bound_idy);
        end
    end
elseif ~isempty(bx)
    mark = 'x';
    if ~flag
        [~,ind,bound_id] = process(dx,valx,abs(y-mean(y)),1);
    else
        [~,ind,bound_id] = process(dx,valx,dygap,0);
    end
    bound = bx(bound_id);
elseif ~isempty(by)
    mark = 'y';
    if ~flag
        [~,ind,bound_id] = process(dy,valy,abs(x-mean(x)),1);
    else
        [~,ind,bound_id] = process(dy,valy,dxgap,0);
    end
    bound = by(bound_id);
end
    function [val,C_id,bound_id] = process(a,valin,clue,type)
        ind_1 = find(a==valin);
        [row,col] = ind2sub(size(a),ind_1);
        if type
            [val,ind_2] = min(clue(row));
        else
            [val,ind_2] = max(clue(row));
        end
        C_id = row(ind_2);
        bound_id = col(ind_2);
    end
end

function [matcomp,Bcomp] = matComp(mark,ind,bound,x,y,matf,totlen)
%--- parameters in compensation equation ---
[~,units,num2] = size(matf);
N =(units-1)/2;
n1 = 1:N;
matcomp = zeros(1,totlen,num2);

switch mark
    case 'x'
        xm = bound-x(ind);
        A21p1 = xm.^(n1);
        A21p2 = xm.^(-n1);
        Bcomp = -log(abs(xm));
        
    case 'y'
        ym = bound-y(ind);
        A21p1 = (1i*ym).^(n1);
        A21p2 = (-1i*ym).^(-n1);
        Bcomp = -log(abs(ym));
end

A21p1 = [0,real(A21p1),imag(A21p1)];
A21p2 = [1,real(A21p2),imag(A21p2)];

matcomp(1,(ind-1)*units+1:ind*units,:) = A21p1+A21p2.*matf;
end

function [matA, matB] = matGroup(matA0,matB0,D,matVc,matVd,gplist)
% num = length(D);
[totlen,~,num2] = size(matA0);

numgp = numel(gplist.NoC);

gpkind = categories(gplist.NoG);
nog = numel(gpkind);

matA = [matA0,-repmat(matB0(:,gplist.NoC),1,1,num2)];
matA = [matA;zeros(numgp,totlen+numgp,num2)];

matB0(:,gplist.NoC) = [];

tempD = D;
tempD(gplist.NoC) = [];

matB = matB0*tempD;
matB = [matB;zeros(numgp,1)];

id = totlen+1;

for k = 1:nog
    kind = gpkind(k);
    list = ismember(gplist.NoG,kind);
    matB(id) = sum(D(list));
    rowD = zeros(1,numgp);
    rowD(list) = 1;
    matA(id,totlen+1:end,:) = repmat(rowD,1,1,num2);
    
    idx = find(list == 1);
    rowVc = matVc(idx(2:end),:,:) - matVc(idx(1),:,:);
    rowVd = matVd(idx(2:end),gplist.NoC,:) - matVd(idx(1),gplist.NoC,:);
    
    matA(id+1:id+nnz(list)-1,:,:) = [rowVc,rowVd];
    id = id + nnz(list);
end

end