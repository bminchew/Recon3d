
%  3D vector displacement estimation using weighted least squares


clear all
%close all
%clc
format long g


weighted = 1;    %  = 1 (true) for weighted least squares (must provide correlation files)

% LOS file names
losnames = ['./Data/LOS/Icelnd_02003_09039_09040.los.grd';
    './Data/LOS/Icelnd_20003_09039_09040.los.grd';
    './Data/LOS/Icelnd_32011_09039_09040.los.grd';
    './Data/LOS/Icelnd_14012_09039_09040.los.grd'];

% Unwrapped interferogram  or defomation file names (in radians/unit_time or deformation/unit_time...must specify which)
typunw = 'phase';  % 'phase' or 'deformation'  (phase is in units radians/time deformation is units distance/time)
unwnames = ['Data/HH/Icelnd_02003_09039_09040.vel_day.flat.grd';
    'Data/HH/Icelnd_20003_09039_09040.vel_day.flat.grd';
    'Data/HH/Icelnd_32011_09039_09040.vel_day.flat.grd';
    'Data/HH/Icelnd_14012_09039_09040.vel_day.flat.grd'];

% Correlatin file names (can be left empty to do unweighted least squares)
cornames = ['Data/HH/Icelnd_02003_09039_09040.cor.grd';
    'Data/HH/Icelnd_20003_09039_09040.cor.grd';
    'Data/HH/Icelnd_32011_09039_09040.cor.grd';
    'Data/HH/Icelnd_14012_09039_09040.cor.grd'];
if ~weighted
    cornames = [];
end

if isempty(cornames) || ~weighted
    disp(' Estimating velocity field using')
    disp('Unweighted least squares');disp(' ');disp('Line = ')
else
    disp(' Estimating velocity field using')
    disp('Weighted least squares based on given InSAR correlation files');disp(' ');disp('Line = ')
end

% open output files
outerror = 1;     % 1 = (true) output model error estimates 
outfold = 'Output';
if ~weighted
    outpref = 'Icelnd_09039_09040_unweighted';%  Icelnd_09039_09040
else
    outpref = 'Icelnd_09039_09040_weighted';
end

cols = [8106; 8210; 11181; 11030];  % samples in each scene

master = 1;                     % master scene
nlooks = 36*ones(size(cols));   % number of looks in each scene (must be same size as cols)
nulu   = 0;                     % null value for displacement files
nulc   = 0;                     % null value for correlation files
nullos = -100;                  % null value for LOS files
nulout = -100;                  % null value for output

corn = [65.12298720    -19.37436182;   % upper left corner lat/lon for each scene
     65.16343488    -19.53124912;
     65.09248476    -19.74435812;
     65.19327060    -19.45969428];
spc  = [-0.000055560     0.000111120;  % spacing deg/pixel for each scene
    -0.000055560     0.000111100; 
    -0.000055560     0.000111120;
    -0.000055560     0.000111100];

timeconv = 365;     % converts time from orig time units to new time units
dispconv = 1;       % converts from one distance unit to another (not necessary if unw is in radians...just convert lambda in that case)
lambda   = 0.24;    % radar wavelength in preferred distance units (lambda [m] -> displacement in [m/time_unit])

%%%%    End of input parameters
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

fideast = fopen([outfold,'/',outpref,'.east'],'w','ieee-le');
fidnorth= fopen([outfold,'/',outpref,'.north'],'w','ieee-le');
fidup   = fopen([outfold,'/',outpref,'.up'],'w','ieee-le');
if outerror == 1
    fiderre = fopen([outfold,'/',outpref,'.east.error'],'w','ieee-le');
    fiderrn = fopen([outfold,'/',outpref,'.north.error'],'w','ieee-le');
    fiderru = fopen([outfold,'/',outpref,'.up.error'],'w','ieee-le');
end

numscenes = length(cols);   % number of scenes
if numscenes ~= size(corn,1) || numscenes ~= size(spc,1) 
    disp('You must enter corner lat/lon and pixel spacing for each scene')
    stop
end

% initialize
fidlos = zeros(1,numscenes);
fidunw = zeros(1,numscenes);
fidcor = zeros(1,numscenes);
lines  = zeros(size(cols));
offset = zeros(size(corn));
eline  = zeros(size(cols));
ccol   = zeros(size(cols));
acol   = zeros(size(cols));
scol   = zeros(size(cols));
bcol   = zeros(size(cols));

% open input files and setup conditions for solver
for ii = 1:numscenes
    fidlos(ii) = fopen(losnames(ii,:),'r','ieee-le');
    fidunw(ii) = fopen(unwnames(ii,:),'r','ieee-le');
    if ~isempty(cornames)
        fidcor(ii) = fopen(cornames(ii,:),'r','ieee-le');
    end
    dirinfo = dir(unwnames(ii,:));
    lines(ii) = dirinfo.bytes/(4*cols(ii));
    
    offset(ii,1) = -round((corn(master,1) - corn(ii,1))/spc(ii,1));   %  positive southward
    offset(ii,2) = -round((corn(master,2) - corn(ii,2))/spc(ii,2));   %  positive eastward 
    if offset(ii,1) < 0
        fseek(fidunw(ii),-offset(ii,1)*4*cols(ii),'bof');            %  place markers for fread 
        fseek(fidlos(ii),-offset(ii,1)*12*cols(ii),'bof');
        if ~isempty(cornames)
            fseek(fidcor(ii),-offset(ii,1)*4*cols(ii),'bof');
        end
    end
    
    eline(ii) = lines(ii) + offset(ii,1);  % end line for each data set 
    
    % set column offset conditions:  4 cases (ccol)
    if offset(ii,2) <= 0
        bcol(ii)  = 0;                         % left side buffer
        scol(ii)  = -offset(ii,2) + 1;         % starting column
        dcol      = cols(ii) - scol(ii) - 1;
        if dcol < cols(master)
            acol(ii) = cols(master) - dcol;    % right side buffer
            ccol(ii) = 1;    % left shifted column that needs buffering on the right
        else
            acol(ii) = 0;
            ccol(ii) = 2;    % left shifted column that doesn't need any buffering
        end
    else
        bcol(ii)  = offset(ii,2);
        scol(ii)  = 1;
        dcol      = cols(ii);
        if dcol < cols(master)
            acol(ii) = cols(master) - dcol;
            ccol(ii) = 3;    % right shifted column that needs buffering on the right (acol)
        else
            acol(ii) = 0;
            ccol(ii) = 4;    % right shifted column that doesn't need any buffering
        end
    end
end
                       

%  Parameters
if strcmp(typunw,'phase') || strcmp(typunw,'p')
    r2dsp = lambda/(4*pi)*timeconv*dispconv;     % convert from radians/orig_time_unit to displacement/new_time_unit
elseif strcmp(typunw,'deformation') || strcmp(typunw,'d') 
    r2dsp = timeconv*dispconv;
end
ccoef = 1./sqrt(2*nlooks);           % coefficient for sigma of correlation 
elinesrt = sort(eline,'descend');

% loop over rows

for ii=1:min([lines(master) elinesrt(3)])
    
    if mod(ii,100) == 0
        disp(ii)
    end
    
    coloff = zeros(1,numscenes);
    
    unw = nulu*ones(numscenes,cols(master));
    cor = nulc*ones(numscenes,cols(master));
    los = nullos*ones(numscenes,3*cols(master));
    
    count = 0;
    for mm = 1:numscenes
        if ii > offset(mm,1) && ii <= eline(mm) % read in files
            dumunw = (fread(fidunw(mm),cols(mm),'float'))';
            if ~isempty(cornames)
                dumcor = (fread(fidcor(mm),cols(mm),'float'))';
            end
            dumlos = (fread(fidlos(mm),3*cols(mm),'float'))';
            % Reference all vectors to the master before loading matrices:
            % 
            % start by buffering on the left and right, then trim 
            dumunw = [nulu*ones(1,bcol(mm)) dumunw nulu*ones(1,acol(mm))];
            dumunw = dumunw(scol(mm):cols(master)+scol(mm)-1);
            if ~isempty(cornames)
                dumcor = [nulc*ones(1,bcol(mm)) dumcor nulc*ones(1,acol(mm))];
                dumcor = dumcor(scol(mm):cols(master)+scol(mm)-1);
            end
            dumlos = [nullos*ones(1,3*bcol(mm)) dumlos nullos*ones(1,3*acol(mm))];
            dumlos = dumlos(3*(scol(mm)-1)+1:3*(scol(mm)-1)+3*cols(master));
            % now load the matrices
            unw(mm,:) = dumunw;
            if ~isempty(cornames)
                cor(mm,:) = dumcor;
            end
            los(mm,:) = dumlos;
            count = count + 1;
        end
    end
    
    east  = nulout*ones(1,cols(master));
    north = nulout*ones(1,cols(master));
    up    = nulout*ones(1,cols(master));
    
    erre = nulout*ones(1,cols(master));
    errn = nulout*ones(1,cols(master));
    erru = nulout*ones(1,cols(master));
    
    if count >= 3   % loop over cols
        
        for jj = 1:cols(master)
            nonul = 0;
            logic = zeros(1,numscenes);
            
            for mm = 1:numscenes
                if unw(mm,jj) ~= nulu && any(los(mm,(1+(jj-1)*3):(3+(jj-1)*3))==nullos) == 0
                    nonul = nonul + 1;
                    logic(mm) = 1;
                end
            end
            if nonul >= 3
                G = zeros(nonul,3);
                W = zeros(nonul,nonul);
                d = zeros(nonul,1);
                kk = 1;
                for mm = 1:numscenes
                    if logic(mm) == 1
                        if ~isempty(cornames)
                            W(kk,kk) = (ccoef(mm)*sqrt(1-cor(mm,jj)^2)/cor(mm,jj))^2;
                        else
                            W(kk,kk) = 1;
                        end
                        G(kk,:) = los(mm,(1+(jj-1)*3):(3+(jj-1)*3));
                        d(kk)   = unw(mm,jj)*r2dsp;
                        kk = kk + 1;
                    end
                end

                if ~isempty(cornames)
                    [enuv,stdx,mse,S] = lscov(G,d,W);
                else
                    [enuv,stdx] = lscov(G,d,W);
                end
                    
                east(jj)  = enuv(1);
                north(jj) = enuv(2);
                up(jj)    = enuv(3);
                
                erre(jj)  = stdx(1);
                errn(jj)  = stdx(2);
                erru(jj)  = stdx(3);
            end
        end
    end
    fwrite(fideast,east,'float');
    fwrite(fidnorth,north,'float');
    fwrite(fidup,up,'float');
    
    if outerror == 1
        fwrite(fiderre,erre,'float');
        fwrite(fiderrn,errn,'float');
        fwrite(fiderru,erru,'float');
    end        
end

fclose('all');

disp('Done with inversion')
    


%   looking
disp(' ')
disp('looking')
disp(' ')



   
fideast = fopen([outfold,'/',outpref,'.east'],'r','ieee-le');
fidnorth= fopen([outfold,'/',outpref,'.north'],'r','ieee-le');
fidup   = fopen([outfold,'/',outpref,'.up'],'r','ieee-le');
fide    = fopen([outfold,'/',outpref,'.east.looks'],'w','ieee-le');
fidn    = fopen([outfold,'/',outpref,'.north.looks'],'w','ieee-le');
fidu    = fopen([outfold,'/',outpref,'.up.looks'],'w','ieee-le');
dirinf  = dir([outfold,'/',outpref,'.east']);

col = cols(1); 
line = dirinf.bytes/(4*col);

winsz = [5 5];                  % window size in [rows columns]
brkline = 2360;                 % stop looking at this line (set brkline = huge# to do all lines)
otline  = brkline;              % print null values until you reach this line (set to 0 to turn this option off)
nulin  = -100;                  % null values
nulout = NaN;

side  = floor(winsz/2);
colo  = floor(col/prod(winsz));

me = zeros(1,colo);
mn = me;
mu = me;

lncnt = 1;

east  = zeros(winsz(1),col);
north = east;
up    = east;

for ii = side(1)+1:winsz(1):line-side(1)
    
    if mod(lncnt,100) == 0
        disp(ii)
    end
    easti  = (fread(fideast,col*winsz(1),'float'))';
    northi = (fread(fidnorth,col*winsz(1),'float'))';
    upi    = (fread(fidup,col*winsz(1),'float'))';
    for mm = 1:winsz(1)
        east(mm,:)  = easti(1+(mm-1)*col:mm*col);
        north(mm,:) = northi(1+(mm-1)*col:mm*col);
        up(mm,:)    = upi(1+(mm-1)*col:mm*col);
    end
    
    cnt = 0;
    for jj = side(2)+1:winsz(2):col-side(2)
        cnt = cnt + 1;
        es = east(:,jj-side(2):jj+side(2));
        es = es(:);
        es = es(es~=nulin);
        if isempty(es)
            me(cnt) = nulout;
        else
            me(cnt) = mean(es);
        end
        no = north(:,jj-side(2):jj+side(2));
        no = no(:);
        no = no(no~=nulin);
        if isempty(no)
            mn(cnt) = nulout;
        else
            mn(cnt) = mean(no);
        end
        u  = up(:,jj-side(2):jj+side(2));
        u  = u(:);
        u  = u(u~=nulin);
        if isempty(u)
            mu(cnt) = nulout;
        else
            mu(cnt) = mean(u);
        end
    end
    fwrite(fide,me,'float');
    fwrite(fidn,mn,'float');
    fwrite(fidu,mu,'float');
    if lncnt == 2360
        break
    else
        lncnt = lncnt + 1;
    end
end
if lncnt < otline
    me(:) = nulout;
    mn(:) = nulout;
    mu(:) = nulout;
    for ii=lncnt:otline
        fwrite(fide,me,'float');
        fwrite(fidn,mn,'float');
        fwrite(fidu,mu,'float');
    end
end
fclose('all');


