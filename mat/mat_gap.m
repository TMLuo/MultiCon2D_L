function [matB,countD,matA,matBg] =  mat_gap(xg,yg,x,y,bx,by,N,Nref,k,gap,flag,type)
%   Matrix Geometric Supplementary Parts for Anti-potential Source
%   for 2D Round Conductors Approach
%   Cgap = matA*xc +matBg*Dc
%
%   Input variables
%   N           Truncated order, default N = 3
%   xg          X coordinates of centre points m
%   yg          Y coordinates of centre points m
%   x           X coordinates of centre points m
%   y           Y coordinates of centre points m
%   bx          Coordinates of vertical boundary m
%   by          Coordinates of horizental boundary m
%   Nref        Number of reflection
%   k           Reflection coefficient, default k = 1
%   gap         Gap length m
%   flag        Direction  of linear anti source 'x'/'y'
%   type        Type of anti source
%
%   Output variable
%   matA:       Geometric part on the left side
%   matB:       Geometric part on the right side
%   countD:     Contribution times of each round cell
%   matBg       Geometric part contributing to anti source
%
%   Reference:
%   ...
%
%   Authors:
%   Tianmingluo
%   Delft University of Technology
%
%   Email:
%   T.Luo-1@tudelft.n

arguments
    xg      (:,1) double  % x coordinates of centre points
    yg      (:,1) double  % y coordinates of centre points
    x       (:,1) double  % x coordinates of centre points
    y       (:,1) double  % y coordinates of centre points
    bx      (1,:) double
    by      (1,:) double
    N       (1,1) double {mustBeNumeric};
    Nref    (1,1) double {mustBeNumeric}
    k       (1,1) double {mustBeNumeric}
    gap     (:,1) double =[];
    flag    (:,1) char =[];
    type    string {mustBeMember(type,["points","lines"])} = "points";
end

if nnz(isnan(bx))
    id = isnan(bx);
    bx(id) = [];
end

if nnz(isnan(by))
    id = isnan(by);
    by(id) = [];
end

units = 2*N+1;

num = length(x);

numg = length(xg);

matB = cell(num,numg);
matB(:) = {zeros(units,1)};

countD = zeros(num,numg);

matAij = cell(numg,num);
matAij(:) = {zeros(1,units)};

matBg = zeros(numg,num);

n1 = (1:N)';

nx = length(bx);
ny = length(by);

ind = 0;

lenin = 0;
flagin = 'y';

if type=="lines"
    ind = 1;
end

if nx == 2
    w = bx(2)-bx(1);
    w(1,1,2) = -w(1,1,1);
end

if ny == 2
    h = by(2)-by(1);
    h(1,1,2) = -h(1,1,1);
end

for idx = 1:numg
    
    dx = xg(idx)-x;
    dy = yg(idx)-y;
    
    if ind
        lenin = gap(idx);
        flagin = flag(idx);
    else
        lenin = 0;
        flagin = 'y';
    end
    
    iterate_zero(dx,dy,idx,k^0);
    
    if (nx+ny)~=0
        
        if nx~= 0
            dx1 = 2*bx(1)-x-xg(idx);
            if nx == 2
                dx2 = 2*bx(2)-x-xg(idx);
            end
        end
        if ny~=0
            dy1 = 2*by(1)-y-yg(idx);
            if ny ==2
                dy2 = 2*by(2)-y-yg(idx);
            end
        end
        
        if nx~=0
            iterate_zero(dx1,dy,idx,k^1);
            if nx == 2
                iterate_zero(dx2,dy,idx,k^1);
                iterate_two(dx,dx2,dy,w,"xnone",idx,k);
            end
        end
        
        if ny~=0
            iterate_zero(dx,dy1,idx,k^1);
            if ny == 2
                iterate_zero(dx,dy2,idx,k^1);
                iterate_two(dy,dy2,dx,h,"ynone",idx,k);
            end
        end
        
        if ny==1 && nx ==1
            iterate_zero(dx1,dy1,idx,k^2);
        end
        
        if ny ==1 && nx ==2
            iterate_two(dx,dx2,dy1,w,"yxy",idx,k);
        end
        
        if nx == 1 && ny ==2
            iterate_two(dy,dy2,dx1,h,"xxy",idx,k);
        end
        
        if nx+ny == 4
            iterate_four(dx,dx2,dy,dy2,w,h,idx,k);
        end
    end
end

matA = cell2mat(matAij);

matB = cell2mat(matB);

    function iterate_zero(dx,dy,idx,k0)
        for idx2 = 1:num
            xin = dx(idx2);
            yin = dy(idx2);
            
            [Bplus,Dplus,Bgplus,Aplus] = computblock(xin,yin,k0,n1,ind,lenin,flagin);
            
            adding(idx,idx2,Bplus,Dplus,Bgplus,Aplus);
        end
    end

    function iterate_two(dx,dx2,dy,w,type,idx,k0)
        for idx2 = 1:num
            yin = dy(idx2);
            for time = 2:Nref
                if ~mod(time,2)
                    if contains(type,"none")
                        xin = dx(idx2)+time*w;
                    else
                        vmove = reshape([time-2,time],[1,1,2]);
                        xin = dx2(idx2)+vmove.*w;
                    end
                else
                    if contains(type,"none")
                        vmove = reshape([time-1,time+1],[1,1,2]);
                        xin = dx2(idx2)+vmove.*w;
                    else
                        xin = dx(idx2)+(time-1)*w;
                    end
                end
                
                switch type
                    case {"xnone","yxy"}
                        [Bplus,Dplus,Bgplus,Aplus] = computblock(xin,yin,k0^time,n1,ind,lenin,flagin);
                    case {"ynone","xxy"}
                        [Bplus,Dplus,Bgplus,Aplus] = computblock(yin,xin,k0^time,n1,ind,lenin,flagin);
                end
                
                adding(idx,idx2,Bplus,Dplus,Bgplus,Aplus);
            end
        end
    end

    function iterate_four(dx,dx2,dy,dy2,w,h,idx,k0)
        for idx2 = 1:num
            for time = 2:Nref
                
                mx = 1:time-1;
                my = time-mx;
                
                for idx3 = 1:time-1
                    
                    if mod(mx(idx3),2) == 1 && mod(my(idx3),2) == 1
                        xmove = reshape([mx(idx3)-1,mx(idx3)+1],[1,1,2]);
                        ymove = reshape([my(idx3)-1,my(idx3)+1],[1,1,2]);
                        xin = dx2(idx2)+repmat(xmove,1,1,2).*repmat(w,1,1,2);
                        yin = dy2(idx2)+reshape(repmat(ymove,1,2).*repmat(h,1,2),1,1,4);
                        
                    elseif mod(mx(idx3),2) == 0 && mod(my(idx3),2) == 1
                        ymove = reshape([my(idx3)-1,my(idx3)+1],[1,1,2]);
                        xin = dx(idx2)+ mx(idx3)*repmat(w,1,1,2);
                        yin = dy2(idx2)+reshape(repmat(ymove,1,2).*repmat(h,1,2),1,1,4);
                        
                    elseif mod(mx(idx3),2) == 1 && mod(my(idx3),2) == 0
                        xmove = reshape([mx(idx3)-1,mx(idx3)+1],[1,1,2]);
                        xin = dx2(idx2)+repmat(xmove,1,1,2).*repmat(w,1,1,2);
                        yin = dy(idx2)+reshape(my(idx3)*repmat(h,1,2),1,1,4);
                        
                    else
                        xin = dx(idx2)+ mx(idx3)*repmat(w,1,1,2);
                        yin = dy(idx2)+reshape(my(idx3)*repmat(h,1,2),1,1,4);
                        
                    end
                    [Bplus,Dplus,Bgplus,Aplus] = computblock(xin,yin,k0^time,n1,ind,lenin,flagin);
                    adding(idx,idx2,Bplus,Dplus,Bgplus,Aplus);
                end
            end
        end
    end

    function adding(i,j,A,B,C,D)
        matB{j,i} = matB{j,i}+ A;
        countD(j,i) = countD(j,i)+B;
        matBg(i,j) = matBg(i,j)+C;
        matAij{i,j} = matAij{i,j}+D;
    end
end

function [Bplus,Dplus,Bgplus,Aplus] = computblock(xin,yin,k,n1,ind,len,flag)

arguments
    xin
    yin
    k
    n1
    ind
    len
    flag char {mustBeMember(flag,['x','y'])} = 'x';
end

num = max(length(xin),length(yin));

if ind
    [B,rfgeo,ifgeo] = intFgeo(xin,yin,len,n1,flag);
    
    Dplus = k*num*abs(len);
    
    Aplus = intcontri_ij(-xin,-yin,len,k,n1,flag);
else
    [B,rfgeo,ifgeo] = Fgeo(xin,yin,n1);
    
    Dplus = k*num;
    
    Aplus = contri_ij(-xin,-yin,k,n1);
end

Bplus = -k*[B;rfgeo;ifgeo];

Bgplus = - k*B;

end

function Mij = contri_ij(dx,dy,k,n)
% contribution from Aj to Agap
arguments
    dx  (:,1) double
    dy  (:,1) double
    k   (:,1) double
    n   (1,:) double
end
pgeo = sum(k.*(-1).^n./(dx-1i*dy).^n,1);
mr = real(pgeo);
mi = imag(pgeo);
Mij = [0,mr,mi];
end

function Mij = intcontri_ij(x,y,len,k,n,flag)
% contribution from Aj to Agap
arguments
    x   (:,1) double
    y   (:,1) double
    len (1,1) double
    k
    n   (1,:) double
    flag
end

if flag == 'x'
    x1 = min([x+len/2,x-len/2],[],2);
    x2 = max([x+len/2,x-len/2],[],2);
end
if flag == 'y'
    y1 = min([y+len/2,y-len/2],[],2);
    y2 = max([y+len/2,y-len/2],[],2);
end

n1 = n(2:end);

if flag == 'x'
    Fac1 = sum(log(abs(x2-1i*y))+1i*atan(-y./x2)-log(abs(x1-1i*y))-1i*atan(-y./x1),1);
    Fac2 = sum(-((x2-1i*y).^(1-n1)-(x1-1i*y).^(1-n1))./(n1-1),1);
    Fac = k.*[Fac1,Fac2].*(-1).^n;
    
elseif flag =='y'
    Fac1 = sum(1i*(log(abs(x-1i*y2))+1i*atan(-y2./x)-log(abs(x-1i*y1))-1i*atan(-y1./x)),1);
    Fac2 = sum(-1i*((x-1i*y2).^(1-n1)-(x-1i*y1).^(1-n1))./(n1-1),1);
    Fac = k.*[Fac1,Fac2].*(-1).^n;
end

mr = real(Fac);
mi = imag(Fac);
Mij = [0,mr,mi];
end

function [B0,fgeor,fgeoi] = Fgeo(x,y,n)

arguments
    x  (1,:) double
    y  (1,:) double
    n  (:,1) double
end

Fac = sum((x-1i*y).^-n./n,2);

fgeor = real(Fac);
fgeoi = imag(Fac);

B0 = sum(-1/2*log(x.^2+y.^2));

end

function [B0,fgeor,fgeoi] = intFgeo(x,y,len,n,flag)

arguments
    x  (1,:) double
    y  (1,:) double
    len(1,1) double
    n  (:,1) double
    flag
end

if flag == 'x'
    x1 = min([x+len/2;x-len/2],[],1);
    x2 = max([x+len/2;x-len/2],[],1);
end
if flag == 'y'
    y1 = min([y+len/2;y-len/2],[],1);
    y2 = max([y+len/2;y-len/2],[],1);
end

n1 = n(2:end);

if flag == 'x'
    Fac1 = sum(log(abs(x2-1i*y))+1i*atan(-y./x2)-log(abs(x1-1i*y))-1i*atan(-y./x1),2);
    Fac2 = sum(-((x2-1i*y).^(1-n1)-(x1-1i*y).^(1-n1))./(n1.*(n1-1)),2);
    if isempty(n1)
        Fac = Fac1;
    else
        Fac = [Fac1;Fac2];
    end
elseif flag =='y'
    Fac1 = sum(1i*(log(abs(x-1i*y2))+1i*atan(-y2./x)-log(abs(x-1i*y1))-1i*atan(-y1./x)),2);
    Fac2 = sum(-1i*((x-1i*y2).^(1-n1)-(x-1i*y1).^(1-n1))./(n1.*(n1-1)),2);
    if isempty(n1)
        Fac = Fac1;
    else
        Fac = [Fac1;Fac2];
    end
end

fgeor = real(Fac);
fgeoi = imag(Fac);

if flag =='y'
    B0 = sum(-1/2*(y2.*log(y2.^2+x.^2)-2*y2+2*x.*atan(y2./x)...
        -y1.*log(y1.^2+x.^2)+2*y1-2*x.*atan(y1./x)));
elseif flag =='x'
    B0 = sum(-1/2*(x2.*log(x2.^2+y.^2)-2*x2+2*y.*atan(x2./y)...
        -x1.*log(x1.^2+y.^2)+2*x1-2*y.*atan(x1./y)));
end

end