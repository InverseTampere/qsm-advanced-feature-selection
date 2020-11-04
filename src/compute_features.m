function [FeatureValues,FeatureNames] = compute_features(QSMs)

% -------------------------------------------------------------------------
% COMPUTE_FEATURES.M      Computes thousands of features from QSMs. The
%                           features are different kind of tree attributes,
%                           ratios of these attributes and comparison of
%                           tree distributions with reference
%                           distributions. This version assumes the
%                           structure of QSMs from the TREEQSM version
%                           2.4.0.
%
% Version 1.0.0
% Latest update     3 Now 2020
%
% Copyright (C) 2020 Pasi Raumonen
% ---------------------------------------------------------------------

% Inputs:
% QSMs    Structure array (the output of the qsmtree) that needs the
%           following fields 
%   cylinder  Must have the following fields:
%             radius, length, start, axis, parent, extension, branch, 
%             BranchOrder, PositionInBranch
%   branch    Must have the following fields:
%             order, parent, volume, length, angle, height, azimuth, diameter
%   treedata  Must contain the fields computed with treedata.m version
%             3.0.0 or later
%           
% Output:
% FeatureValues   The feature values, (n_features x n_qsms)-matrix
% FeatureNames    The names of the features, (n_features x 1)-cell
% ---------------------------------------------------------------------

% Computes the features used in Åkerblom et al. 2017, Remote Sensing of
% Environment. In addition to these the following types of features are
% computed: 
% 1) Treedata (= QSM.treedata) and treedata divided by other treedata
% 2) Tree attributes (cylinder volumes, areas and lengths between certain 
%     diameter and height classes and branch orders) divided by treedata 
% 3) Percentiles (5,10,15,..,95%) of the cylinder and branch distributions 
% 4) Comparison of the cylinder and branch distributions to reference
%     distributions (uniform, triangle, normal):
%     - Cylinder volume, area, length as function of their height, diameter, 
%       azimuth and zenith. 
%     - Branch (all branches and 1st-order) volume, area, length, and number 
%       as function of their height, diameter, azimuth, zenith, angle, and
%       order
% 5) Basic statistics (medians, averages, minimums, maximums and their 
%     ratios) of branch data (order, volume, area, length, angle, azimuth, 
%     zenith, diameter)
% 6) Branch tip and base path lengths based features (lengths along the 
%     branching structure from the branch's tip and base to the tree's base): 
%     - Basic statistics (medians, averages, minimums, maximums and their 
%       ratios) of:
%       - Tip/base path lengths relative to the tip/base height
%       - Tip/base path lengths relative to the tree height
%       - Tip/base path lengths relative to the maximum tip path length
%     - Percentile values of:
%       - absolute tip/base path length
%       - tip/base path length relative tip/base height
%       - tip/base path length relative tree height
%       - relative tip/base path length
%     - distribution comparisons to reference distributions (uniform, 
%       triangle, normal)
%     - Relations between different branching orders
% Reference distributions are uniform (step), triangle and normal 
% distributions with the following parameters:
% - 15 uniform (step) distributions:
%   The values divided into 5 equal length section and the following relative
%   starting and ending section:
%   Su = [0 0 0 0 0 1 2 3 4 1 2 3 1 1 2]; relative start (0/5, 1/5, 2/5, etc)
%   Eu = [1 2 3 4 5 5 5 5 5 4 4 4 3 2 3]; relative ends (1/5, 2/5, 5/5, etc)
% - 29 triangle distributions:
%   The values divided into 4 equal length section and the following relative
%   starting, mid/high point, ending section:
%   St = [0 0 0 0 0 0 0 0 0 0 0 0 0 0 1 1 1 1 1 1 1 1 1 2 2 2 2 2 3];
%   Mt = [0 0 0 0 1 1 1 1 2 2 2 3 3 4 1 1 1 2 2 2 3 3 4 2 2 3 3 4 4];
%   Et = [1 2 3 4 1 2 3 4 2 3 4 3 4 4 2 3 4 2 3 4 3 4 4 3 4 3 4 4 4];
% - 6 normal distributions:
%   The values divided into 4 equal length section and the following relative
%   mean sections and standard deviations (relative to 5):
%   Mn = [1 1 2 2 3 3]; % means (relative, quadrants, i.e. 1/4 2/4, etc)
%   Sn = [1 2 1 2 1 2]; % standard deviations (relative, 1/5 and 2/5)


% COMPUTE_FEATURES is free software: you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation, either version 3 of the License, or
% (at your option) any later version.
% 
% COMPUTE_FEATURES is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
% GNU General Public License for more details, 
% see <http://www.gnu.org/licenses/>.


%% The field names in QSM.treedata computed with treedata.m (version 3.0.0)
% that are needed to compute all the features:

% TotalVolume
% TrunkVolume
% BranchVolume
% TreeHeight
% TrunkLength
% BranchLength
% TotalLength
% NumberBranches    Total number of branches
% MaxBranchOrder 
% TrunkArea 
% BranchArea 
% TotalArea 
% DBHqsm        From the cylinder of the QSM at the right heigth
% DBHcyl        From the cylinder fitted to the section 1.1-1.5m
% CrownDiamAve
% CrownDiamMax
% CrownAreaConv
% CrownAreaAlpha
% CrownBaseHeight
% CrownLength
% CrownRatio
% CrownVolumeConv
% CrownVolumeAlpha
% location      (x,y,z)-coordinates of the base of the tree
% StemTaper     Stem taper function/curve from the QSM
% VerticalProfile
% spreads
% VolCylDia     Distribution of the total volume in diameter classes
% AreCylDia     Distribution of the total area in diameter classes
% LenCylDia     Distribution of the total length in diameter classes
% VolCylHei     Distribution of the total volume in height classes
% AreCylHei     Distribution of the total area in height classes
% LenCylHei     Distribution of the total length in height classes
% VolCylAzi     Distribution of the total volume in azimuth angle classes
% AreCylAzi     Distribution of the total area in azimuth angle classes
% LenCylAzi     Distribution of the total length in azimuth angle classes
% VolCylZen     Distribution of the total volume in zenith angle classes
% AreCylZen     Distribution of the total area in zenith angle classes
% LenCylZen     Distribution of the total length in zenith angle classes
% VolBranchOrd     Branch volume per branching order
% AreBranchOrd     Branch area per branching order
% LenBranchOrd     Branch length per branching order
% NumBranchOrd     Number of branches per branching order
% VolBranchDia
% VolBranch1Dia
% AreBranchDia
% AreBranch1Dia
% LenBranchDia
% LenBranch1Dia
% NumBranchDia
% NumBranch1Dia
% VolBranchAng
% VolBranch1Ang
% AreBranchAng
% AreBranch1Ang
% LenBranchAng
% LenBranch1Ang
% NumBranchAng
% NumBranch1Ang
% VolBranchAzi 
% VolBranch1Azi
% AreBranchAzi
% AreBranch1Azi
% LenBranchAzi
% LenBranch1Azi
% NumBranchAzi
% NumBranch1Azi
% VolBranchHei
% VolBranch1Hei
% AreBranchHei
% AreBranch1Hei
% LenBranchHei
% LenBranch1Hei
% NumBranchHei
% NumBranch1Hei
% VolBranchZen 
% VolBranch1Zen 
% AreBranchZen 
% AreBranch1Zen 
% LenBranchZen 
% LenBranch1Zen 
% NumBranchZen 
% NumBranch1Zen 


%% Define the needed objects
Nmodels = max(size(QSMs)); % number of QSMs
F = zeros(20000,Nmodels); % feature values
FN = cell(20000,1); % feature names

NameB = fieldnames(QSMs(1).branch); % field names in "branch"
mb = size(NameB,1); % number of fields in "branch"

NameT = fieldnames(QSMs(1).treedata); % field names in "treedata"
mt = size(NameT,1); % number of fields in "treedata"
% Determine the indexes of certain fields of the treedata
for i = 1:mt
    if strcmp(NameT{i},'CrownVolumeAlpha')
        nf = i;
    end
    if strcmp(NameT{i},'VolCylDia')
        ncs = i;
    end
    if strcmp(NameT{i},'LenCylZen')
        nce = i;
    end
    if strcmp(NameT{i},'VolBranchOrd')
        nbs = i;
    end
    if strcmp(NameT{i},'NumBranch1Zen')
        nbe = i;
    end
end

% Shedding ratios
SR = [0 1 2 3 4 0 0 0 0 1 1 1 2 2 3]/5;
ER = [1 2 3 4 5 2 3 4 5 3 4 5 4 5 5]/5;

% relative height/diameter/zenith of bottom 5, 10, 15, etc.
% volume/area/length
RHN = (5:5:95); % for names
RH = RHN/100; % for computations

% Cylinder height/diameter classes (Sc (m/cm) - Ec (m/cm)) and branch
% orders (Sb - Eb)
Sc = [1 1 1 1 1  3 3 3 3  5 5 5  7 7  9]/10;
Ec = [2 4 6 8 10 4 6 8 10 6 8 10 8 10 10]/10;
Sb = [1 2 3 4 5 1 1 1 1 2 2 2 3 3 4];
Eb = [1 2 3 4 5 2 3 4 5 3 4 5 4 5 5];

% Uniform and step distributions
Su = [0 0 0 0 0 1 2 3 4 1 2 3 1 1 2]; % relative start (0/5, 1/5, 2/5, etc)
Eu = [1 2 3 4 5 5 5 5 5 4 4 4 3 2 3]; % relative ends (1/5, 2/5, 5/5, etc)

% normal distributions (mean,std)
Mn = [1 1 2 2 3 3]; % means (relative, quadrants, i.e. 1/4 2/4, etc)
Sn = [1 2 1 2 1 2]; % standard deviations (relative, 1/5 and 2/5)
        
% Triangle distributions (start, mid/high point, end points):
% start points:
St = [0 0 0 0 0 0 0 0 0 0 0 0 0 0 1 1 1 1 1 1 1 1 1 2 2 2 2 2 3];
% mid/high points
Mt = [0 0 0 0 1 1 1 1 2 2 2 3 3 4 1 1 1 2 2 2 3 3 4 2 2 3 3 4 4];
% end points
Et = [1 2 3 4 1 2 3 4 2 3 4 3 4 4 2 3 4 2 3 4 3 4 4 3 4 3 4 4 4];


%% Compute the features and create the feature name strings
for i = 1:Nmodels
  f = 0; % feature index
  t = QSMs(i).treedata; % treedata of the model
  if ~isempty(t)
    
    %% Features from the 2017 RSE-paper
    % Stem branch angle (median branching angle of 1st-ord. branches)
    b = QSMs(i).branch;
    Ord1 = b.order == 1;
    nb = length(b.order);
    indb = (1:1:nb)';
    Ord1 = indb(Ord1);
    m1 = length(Ord1);
    f = f+1;
    F(f,i) = median(b.angle(Ord1));
    if i == 1
      FN{f} = 'StemBranchAngle';
    end
    % Stem branch cluster size (Mean number of 1branches inside 40cm height
    % interval. Each branch can only belong to one interval).
    if m1 > 0
      h1 = b.height(Ord1);
      Int = zeros(m1,1);
      k = 10;
      j = 1;
      while j <= m1
        h0 = h1(j);
        jk = j;
        while jk <= m1 && h1(jk) < h0+0.4
          jk = jk+1;
        end
        jk = jk-1;
        k = k+1;
        Int(k) = jk-j+1;
        j = jk+1;
      end
      CS = mean(Int(1:k));
    else
      CS = 0;
    end
    f = f+1;
    F(f,i) = CS;
    if i == 1
      FN{f} = 'StemBranchClusterSize';
    end
    % Stem branch radius (Mean ratio between the 10 largest 1branches and
    % stem radius at respective height).
    r1 = b.diameter(Ord1)/2;
    [r,I] = sort(r1,'descend');
    if m1 >= 10
      r = r(I(1:10));
      Ind = Ord1(I(1:10));
    else
      Ind = Ord1;
    end
    m = length(r);
    c = QSMs(i).cylinder;
    R = r;
    for j = 1:m
      a = c.branch == Ind(j) & c.PositionInBranch == 1;
      p = c.parent(a);
      R(j) = c.radius(p);
    end
    f = f+1;
    F(f,i) = mean(r./R);
    if i == 1
      FN{f} = 'StemBranchRadius';
    end
    % Stem branch length (Mean length of 1branches normalized by DBH)
    f = f+1;
    F(f,i) = mean(b.length(Ord1))/t.DBHcyl;
    if i == 1
      FN{f} = 'StemBranchLength';
    end
    % Stem branch distance (Mean distance between 1branches from moving
    % average with 1m window. If empty window, set value to half the
    % window width. Normilize with dbh).
    if m1 > 0
      D = zeros(m1^2,1);
      k = 1;
      for j = 1:m1
        d = abs(h1-h1(j));
        I = d <= 0.5;
        I(j) = false;
        d = d(I);
        m = length(d);
        if m > 0
          D(k:k+m-1) = d;
        else
          D(k) = 0.5;
          m = 1;
        end
        k = k+m;
      end
      k = k-m;
      D = mean(D(1:k));
    else
      D = 0;
    end
    f = f+1;
    F(f,i) = D/t.DBHcyl;
    if i == 1
      FN{f} = 'StemBranchDistance';
    end
    % Crown start height (Height of first stem branch in crown relative to
    % tree height)
    f = f+1;
    F(f,i) = t.CrownBaseHeight/t.TreeHeight;
    if i == 1
      FN{f} = 'CrownStartHeight';
    end
    % Crown height (Crown length divided by the tree height)
    f = f+1;
    F(f,i) = t.CrownRatio;
    if i == 1
      FN{f} = 'CrownHeight';
    end
    % Crown evenness (Crown cylinders divided into 8 angular bins. Ratio
    % between extreme minimum heights in bins)
    nc = length(c.radius);
    indc = (1:1:nc)';
    bot = min(c.start(:,3));
    hcs = c.start(:,3)-bot;
    E = c.start+[c.length.*c.axis(:,1) c.length.*c.axis(:,2) ...
      c.length.*c.axis(:,3)];
    hce = E(:,3)-bot;
    I = (hcs >= t.CrownBaseHeight | hce >= t.CrownBaseHeight) & ...
      c.BranchOrder > 0;
    crown = indc(I);
    if ~isempty(crown)
      hs = hcs(crown);
      he = hce(crown);
      Cen = mean(E(crown,:),1);
      V = mat_vec_subtraction(E(crown,:),Cen);
      ang = atan2(V(:,2),V(:,1))+pi;
      H = zeros(8,1);
      for j = 1:8
        I = ang >= (j-1)*pi/4 & ang < j*pi/4;
        if any(I)
          H(j) = min(min(hs(I)),min(he(I)));
        end
      end
    end
    f = f+1;
    F(f,i) = min(H(H > 0))/max(H);
    if i == 1
      FN{f} = 'CrownEvenness';
    end
    % Crown diameter/height (Ratio between crown diameter and height)
    f = f+1;
    F(f,i) = t.CrownDiamAve/t.CrownLength;
    if i == 1
      FN{f} = 'CrownDiameter/Height';
    end
    % DBH/height ratio (Ratio between DBH and tree height)
    f = f+1;
    F(f,i) = t.DBHcyl/t.TreeHeight;
    if i == 1
      FN{f} = 'DBH/TreeHeight';
    end
    % DBH/tree volume ratio (Ratio between DBH and total tree volume)
    f = f+1;
    F(f,i) = t.DBHcyl/t.TotalVolume;
    if i == 1
      FN{f} = 'DBH/TotalVolume';
    end
    % DBH/minimum tree radius (Ratio between DBH and the minimum of the
    % vertical bins radius estimate)
    height= max(hce);
    CylVol = 1000*pi*c.radius.^2.*c.length;
    R = zeros(3,1);
    for j = 1:3
      I = hce >= (j-1)*height/3 & hce < j*height/3;
      L = indc(I);
      V = CylVol(L);
      stem = c.branch(L) == 1;
      stem = L(stem);
      if ~isempty(stem)
        Cen = mean(c.start(stem,:),1);
      end
      d = distances_to_line(E(L,:),[0 0 1],Cen);
      [d,I] = sort(d);
      v = cumsum(V(I))/sum(V);
      k = 1;
      while v(k) < 0.9
        k = k+1;
      end
      R(j) = 2*d(k);
    end
    f = f+1;
    F(f,i) = t.DBHcyl/min(R);
    if i == 1
      FN{f} = 'DBH/MinimumTreeRadius';
    end
    % Volume below 55% of height (Relative volume below 55%)
    I = hce/height < 0.55;
    f = f+1;
    F(f,i) = sum(CylVol(I))/t.TotalVolume;
    if i == 1
      FN{f} = 'VolumeBelow55%';
    end
    % Cylinder length/tree volume (Ratio between total length and total volume)
    f = f+1;
    F(f,i) = t.TotalLength/t.TotalVolume;
    if i == 1
      FN{f} = 'TotalLength/TotalVolume';
    end
    % Shedding ratio (Number of branches without children divided by the number
    % of branches in the bottom third)
    B = indb(b.height < t.TreeHeight/3 & b.order == 1);
    m = length(B);
    k = 0;
    for j = 1:m
      if ~any(b.parent == B(j))
        k = k+1;
      end
    end
    f = f+1;
    if m > 0
      F(f,i) = k/m;
    end
    if i == 1
      FN{f} = 'SheddingRatio';
    end
    
    % Shedding ratios (Number of branches without children divided by the number
    % of branches in different height layers)
    h = b.height;
    o = b.order;
    th = t.TreeHeight;
    for l = 1:15
      B = indb(h < ER(l)*th & h >= SR(l)*th & o == 1);
      m = length(B);
      k = 0;
      for j = 1:m
        if ~any(b.parent == B(j))
          k = k+1;
        end
      end
      f = f+1;
      if m > 0
        F(f,i) = k/m;
      end
      if i == 1
        FN{f} = ['SheddingRatio_',num2str(SR(l)),'_',num2str(ER(l))];
      end
    end
    
    
    %% Treedata and treedata divided by other treedata
    for j = 1:nf
      f = f+1;
      F(f,i) = t.(NameT{j});
      if i == 1
        FN{f} = NameT{j};
      end
      for k = 1:nf
        f = f+1;
        if t.(NameT{k}) ~= 0
          F(f,i) = t.(NameT{j})/t.(NameT{k});
        end
        if i == 1
          FN{f} = [NameT{j},'/',NameT{k}];
        end
      end
    end
    
    
    %% Tree attributes divided by treedata
    % Cylinder volumes, areas and lengths between certain diameter and
    % height classes and branch orders divided by treedata
    for j = 1:nf
      a = t.(NameT{j});
      for l = ncs:ncs+5
        D = t.(NameT{l});
        if length(D) < 10
          D(10) = 0;
        end
        b = length(D);
        for k = 1:length(Sc)
          f = f+1;
          F(f,i) = sum(D(ceil(Sc(k)*b):floor(Ec(k)*b)))/a;%Vol(S-E Diam)/a
          if i == 1
            FN{f} = [NameT{l},'_',num2str(Sc(k)-0.1),'_',...
              num2str(Ec(k)),'/',NameT{j}];
          end
        end
      end
      for l = nbs:nbs+3
        D = t.(NameT{l});
        if length(D) < 5
          D(5) = 0;
        end
        for k = 1:length(Sb)
          f = f+1;
          F(f,i) = sum(D(Sb(k):Eb(k)))/a; % #branch_k / a
          if i == 1
            FN{f} = [NameT{l},'_',num2str(Sb(k)),'_',num2str(Eb(k)),...
              '/',NameT{j}];
          end
        end
      end
    end
    
    
    %% Tree segments (cylinder) distributions
    % Volume, area, length of cylinder as functions of dia, hei, azi, zen
    for j = ncs:nce
      d = t.(NameT{j}); % distribution
      if strcmp(NameT{j}(1:3),'Vol')
        a = t.TotalVolume;
      elseif strcmp(NameT{j}(1:3),'Are')
        a = t.TotalArea;
      elseif strcmp(NameT{j}(1:3),'Len')
        a = t.TotalLength;
      end
      dr = d/a; % relative distribution
      dc = cumsum(dr); % cumulative relative distribution
      
      % relative height/diameter/zenith of 5, 10, 15 etc % of bottom
      % volume/area/length
      N = relative_height(dc,RH);
      for k = 1:length(RH)
        f = f+1;
        F(f,i) = N(k);
        if i == 1
          FN{f} = ['Rel_cyl_',NameT{j}(end-2:end),'_bottom_',...
            num2str(RHN(k)),'%_',NameT{j}(1:3)];
        end
      end
      
      % distribution comparisons:
      m = length(dr);
      % differences to the triangle distributions (start,mid,end)
      % start points:
      Sti = ceil(St/4*m);
      Mti = ceil(Mt/4*m);
      Eti = ceil(Et/4*m);
      for k = 1:length(St)
        di = abs(dr-triad(Sti(k),Mti(k),Eti(k),m));
        % mean and max difference to the comparison distribution:
        f = f+1;
        if ~isempty(di)
          F(f,i) = mean(di);
          F(f+1,i) = max(di);
        end
        if i == 1
          FN{f} = [NameT{j},'_triad_',num2str(St(k)),'_',...
            num2str(Mt(k)),'_',num2str(Et(k)),'_mean'];
          FN{f+1} = [NameT{j},'_triad_',num2str(St(k)),'_',...
            num2str(Mt(k)),'_',num2str(Et(k)),'_max'];
        end
        f = f+1;
      end
      
      % differences to the normal distributions (mean,std)
      Mni = ceil(Mn/4*m);
      Sni = ceil(Sn/5*m);
      for k = 1:length(Mn)
        di = abs(dr-normd(Mni(k),Sni(k),m));
        % mean and max difference to the comparison distribution:
        f = f+1;
        if ~isempty(di)
          F(f,i) = mean(di);
          F(f+1,i) = max(di);
        end
        if i == 1
          FN{f} = [NameT{j},'_normd_',num2str(Mn(k)),'_',...
            num2str(Sn(k)),'_mean'];
          FN{f+1} = [NameT{j},'_normd_',num2str(Mn(k)),'_',...
            num2str(Sn(k)),'_max'];
        end
        f = f+1;
      end
      
      % difference to the uniform and step distributions
      Sui = ceil(Su/5*m);
      Eui = ceil(Eu/5*m);
      for k = 1:length(Su)
        di = abs(dr-unid(Sui(k),Eui(k),m));
        % mean and max difference to the comparison distribution:
        f = f+1;
        if ~isempty(di)
          F(f,i) = mean(di);
          F(f+1,i) = max(di);
        end
        if i == 1
          FN{f} = [NameT{j},'_unid_',num2str(Su(k)),'_',...
            num2str(Eu(k)),'_mean'];
          FN{f+1} = [NameT{j},'_unid_',num2str(Su(k)),'_',...
            num2str(Eu(k)),'_max'];
        end
        f = f+1;
      end
    end
    
    
    %% Branch distributions
    % Volume, area, length, number (of 1st-order) of branches as functions
    % of height, diameter, angle
    for j = nbs:nbe
      d = t.(NameT{j}); % distribution
      if strcmp(NameT{j}(1:3),'Vol')
        a = t.BranchVolume;
        if strcmp(NameT{j}(end-3),'1')
          a = t.VolBranchOrd(1);
        end
      elseif strcmp(NameT{j}(1:3),'Are')
        a = t.BranchArea;
        if strcmp(NameT{j}(end-3),'1')
          a = t.AreBranchOrd(1);
        end
      elseif strcmp(NameT{j}(1:3),'Len')
        a = t.BranchLength;
        if strcmp(NameT{j}(end-3),'1')
          a = t.LenBranchOrd(1);
        end
      elseif strcmp(NameT{j}(1:3),'Num')
        a = t.NumberBranches;
        if strcmp(NameT{j}(end-3),'1')
          a = t.NumBranchOrd(1);
        end
      end
      dr = d/a; % relative distribution
      dc = cumsum(dr); % cumulative relative distribution
      
      % relative height/diameter/zenith of 5, 10, 15 etc % of bottom
      % volume/area/length
      N = relative_height(dc,RH);
      for k = 1:19
        f = f+1;
        F(f,i) = N(k);
        if i == 1
          if strcmp(NameT{j}(end-3),'h')
            FN{f} = ['Rel_branch_',NameT{j}(end-2:end),...
              '_bottom_',num2str(RHN(k)),'%_',NameT{j}(1:3)];
          else
            FN{f} = ['Rel_branch1_',NameT{j}(end-2:end),...
              '_bottom_',num2str(RHN(k)),'%_',NameT{j}(1:3)];
          end
        end
      end
      
      % distribution comparisons:
      m = length(dr);
      % differences to the triangle distributions (start,mid,end)
      Sti = ceil(St/4*m);
      Mti = ceil(Mt/4*m);
      Eti = ceil(Et/4*m);
      for k = 1:length(St)
        di = abs(dr-triad(Sti(k),Mti(k),Eti(k),m));
        % mean and max difference to the comparison distribution:
        f = f+1;
        if ~isempty(di)
          F(f,i) = mean(di);
          F(f+1,i) = max(di);
          
        end
        if i == 1
          FN{f} = [NameT{j},'_triad_',num2str(St(k)),...
            '_',num2str(Mt(k)),'_',num2str(Et(k)),'_mean'];
          FN{f+1} = [NameT{j},'_triad_',num2str(St(k)),'_',...
            num2str(Mt(k)),'_',num2str(Et(k)),'_max'];
        end
        f = f+1;
      end
      
      % differences to the normal distributions (mean,std)
      Mni = ceil(Mn/4*m);
      Sni = ceil(Sn/5*m);
      for k = 1:length(Mn)
        di = abs(dr-normd(Mni(k),Sni(k),m));
        % mean and max difference to the comparison distribution:
        f = f+1;
        if ~isempty(di)
          F(f,i) = mean(di);
          F(f+1,i) = max(di);
        end
        if i == 1
          FN{f} = [NameT{j},'_normd_',num2str(Mn(k)),'_',num2str(Sn(k)),'_mean'];
          FN{f+1} = [NameT{j},'_normd_',num2str(Mn(k)),'_',num2str(Sn(k)),'_max'];
        end
        f = f+1;
      end
      
      % difference to the uniform and step distributions
      Sui = ceil(Su/5*m);
      Eui = ceil(Eu/5*m);
      for k = 1:length(Su)
        di = abs(dr-unid(Sui(k),Eui(k),m));
        % mean and max difference to the comparison distribution:
        f = f+1;
        if ~isempty(di)
          F(f,i) = mean(di);
          F(f+1,i) = max(di);
        end
        if i == 1
          FN{f} = [NameT{j},'_unid_',num2str(Su(k)),'_',...
            num2str(Eu(k)),'_mean'];
          FN{f+1} = [NameT{j},'_unid_',num2str(Su(k)),'_',...
            num2str(Eu(k)),'_max'];
        end
        f = f+1;
      end
    end
    
    
    %% Branch-data
    % Medians, averages, minimums, maximums and their ratios of branch-data
    b = QSMs(i).branch;
    I1 = b.order == 1;
    I2 = b.order == 2;
    I3 = b.order == 3;
    for j = [1 3:mb]
      d = double(b.(NameB{j}));
      for o = 0:3
        if o == 0
          a = d; % all branches
        elseif o == 1
          a = d(I1); % 1st-order branches
        elseif o == 2
          a = d(I2); % 2nd-order branches
        elseif o == 3
          a = d(I3); % 3rd-order branches
        end
        f = f+1;
        if ~isempty(a)
          F(f:f+6,i) = [median(a); mean(a); min(a); max(a); ...
            mean(a)/max(a); min(a)/max(a); min(a)/mean(a)];
        end
        if i == 1
          FN{f} = ['branch_or',num2str(o),'_',NameB{j},'_median'];
          FN{f+1} = ['branch_or',num2str(o),'_',NameB{j},'_mean'];
          FN{f+2} = ['branch_or',num2str(o),'_',NameB{j},'_min'];
          FN{f+3} = ['branch_or',num2str(o),'_',NameB{j},'_max'];
          FN{f+4} = ['branch_or',num2str(o),'_',NameB{j},'_mean/max'];
          FN{f+5} = ['branch_or',num2str(o),'_',NameB{j},'_min/max'];
          FN{f+6} = ['branch_or',num2str(o),'_',NameB{j},'_min/mean'];
        end
        f = f+6;
        for k = [1 3:mb]
          if k ~= j
            c = double(b.(NameB{k}));
            for or = 0:3
              if or == 0
                a = c;
              elseif or == 1
                a = c(I1);
              elseif or == 2
                a = c(I2);
              elseif or == 3
                a = c(I3);
              end
              f = f+1;
              if ~isempty(a) && ~isempty(c)
                F(f,i) = mean(c)/mean(a);
              end
              if i == 1
                FN{f} = ['branch_or',num2str(or),'_',NameB{k},...
                  '_mean/','or',num2str(o),'_',NameB{j},'_mean'];
              end
            end
          end
        end
      end
    end
    
    
    %% Branch azimuth
    % distribution comparisons
    m = t.NumberBranches;
    a = t.NumBranchAzi/m; % relative branch-angle distribution
    di = abs(a-unid(0,36,36)); % difference to the 0-360 uniform distribution
    f = f+1;
    if ~isempty(di)
      F(f,i) = mean(di);   % mean difference to the uniform angle distribution
      F(f+1,i) = max(di);   % maximum difference to the uniform angle distribution
    end
    if i == 1
      FN{f} = 'NumBranchAzi_unid_0_360_mean';
      FN{f+1} = 'NumBranchAzi_unid_0_360_max';
    end
    f = f+1;
    
    m1 = t.NumBranchOrd(1);
    a = t.NumBranch1Azi/m1; % relative branch-angle distribution (1st-branches)
    di = abs(a-unid(0,36,36));   % difference to the 0-360 uniform distribution
    f = f+1;
    if ~isempty(di)
      F(f,i) = mean(di); % mean difference to the uniform angle distribution
      F(f+1,i) = max(di); % maximum difference to the uniform angle distribution
    end
    if i == 1
      FN{f} = 'NumBranch1Azi_unid_0_360_mean';
      FN{f+1} = 'NumBranch1Azi_unid_0_360_max';
    end
    f = f+1;
    
    
    %% Path-length based features
    c = QSMs(i).cylinder;
    nc = length(c.radius);
    % Determine the child-cylinders for each cylinder:
    Chi = cell(nc,1);
    for j = 2:nc
      Chi{c.parent(j)} = [Chi{c.parent(j)}; j];
    end
    % Compute the path lengths for each cylinder
    PL = zeros(nc,1); % path lengths for every cylinder
    C = Chi{1};
    PL(1) = c.length(1);
    while ~isempty(C)
      PL(C) = PL(c.parent(C))+c.length(C);
      C = vertcat(Chi{C});
    end
    % Compute the heights of the cylinder tips:
    E = c.start+[c.length.*c.axis(:,1) c.length.*c.axis(:,2) ...
      c.length.*c.axis(:,3)];
    H = E(:,3)-min(c.start(:,3));
    
    % Branch tip path lengths to the tree's base
    for j = 1:4
      % Path lengths for tip cylinders:
      if j == 1
        TPL = PL(c.extension == 0);
        TH = H(c.extension == 0);
        str = '_all';
      else
        TPL = PL(c.extension == 0 & c.BranchOrder == j-1);
        TH = H(c.extension == 0 & c.BranchOrder == j-1);
        str = ['_ord',num2str(j-1)];
      end
      
      % Basic statistics:
      f = f+1;
      if ~isempty(TPL)
        F(f:f+6,i) = [mean(TPL); min(TPL); max(TPL); min(TPL)/max(TPL);...
          mean(TPL)/max(TPL); std(TPL); std(TPL)/mean(TPL)];
      end
      if i == 1
        FN{f} = ['mean_tip_path_length',str];
        FN{f+1} = ['min_tip_path_length',str];
        FN{f+2} = ['max_tip_path_length',str];
        FN{f+3} = ['min/max_tip_path_length',str];
        FN{f+4} = ['mean/max_tip_path_length',str];
        FN{f+5} = ['std_tip_path_length',str];
        FN{f+6} = ['std/max_tip_path_length',str];
      end
      f = f+6;
      
      % Tip path lengths relative to the tip height:
      f = f+1;
      if ~isempty(TPL)
        LH = TPL./TH;
        F(f:f+6,i) = [mean(LH); min(LH); max(LH); min(LH)/max(LH);...
          mean(LH)/max(LH); std(LH); std(LH)/mean(LH)];
      end
      if i == 1
        FN{f} = ['mean(tip_path_length/tip_height)',str];
        FN{f+1} = ['min(tip_path_length/tip_height)',str];
        FN{f+2} = ['max(tip_path_length/tip_height)',str];
        FN{f+3} = ['min/max(tip_path_length/tip_height)',str];
        FN{f+4} = ['mean/max(tip_path_length/tip_height)',str];
        FN{f+5} = ['std(tip_path_length/tip_height)',str];
        FN{f+6} = ['std/max(tip_path_length/tip_height)',str];
      end
      f = f+6;
      
      % Tip path lengths relative to the tree height:
      f = f+1;
      if ~isempty(TPL)
        LTH = TPL/t.TreeHeight;
        F(f:f+6,i) = [mean(LTH); min(LTH); max(LTH); min(LTH)/max(LTH);...
          mean(LTH)/max(LTH); std(LTH); std(LTH)/mean(LTH)];
      end
      if i == 1
        FN{f} = ['mean(tip_path_length/tree_height)',str];
        FN{f+1} = ['min(tip_path_length/tree_height)',str];
        FN{f+2} = ['max(tip_path_length/tree_height)',str];
        FN{f+3} = ['min/max(tip_path_length/tree_height)',str];
        FN{f+4} = ['mean/max(tip_path_length/tree_height)',str];
        FN{f+5} = ['std(tip_path_length/tree_height)',str];
        FN{f+6} = ['std/max(tip_path_length/tree_height)',str];
      end
      f = f+6;
      
      % Tip path lengths relative to the maximum tip path length:
      f = f+1;
      if ~isempty(TPL)
        LML = TPL/max(TPL);
        F(f:f+6,i) = [mean(LML); min(LML); max(LML); min(LML)/max(LML);...
          mean(LML)/max(LML); std(LML); std(LML)/mean(LML)];
      end
      if i == 1
        FN{f} = ['mean(relative_tip_path_length)',str];
        FN{f+1} = ['min(relative_tip_path_length)',str];
        FN{f+2} = ['max(relative_tip_path_length)',str];
        FN{f+3} = ['min/max(relative_tip_path_length)',str];
        FN{f+4} = ['mean/max(relative_tip_path_length)',str];
        FN{f+5} = ['std(relative_tip_path_length)',str];
        FN{f+6} = ['std/max(relative_tip_path_length)',str];
      end
      f = f+6;
      
      % Percentile values:
      TPL1 = sort(TPL); % absolute tip path length
      n = length(TPL);
      I = floor(RH*n)+1;
      for k = 1:19
        f = f+1;
        if ~isempty(TPL1)
          F(f,i) = TPL1(I(k));
        end
        if i == 1
          FN{f} = [num2str(RHN(k)),'%_tip_path_length',str];
        end
      end
      
      LH1 = sort(LH); % tip path length relative tip height
      for k = 1:19
        f = f+1;
        if ~isempty(LH1)
          F(f,i) = LH1(I(k));
        end
        if i == 1
          FN{f} = [num2str(RHN(k)),'%_tip_path_length/tip_height',str];
        end
      end
      
      LTH1 = sort(LTH); % tip path length relative tree height
      for k = 1:19
        f = f+1;
        if ~isempty(LTH1)
          F(f,i) = LTH1(I(k));
        end
        if i == 1
          FN{f} = [num2str(RHN(k)),'%_tip_path_length/tree_height',str];
        end
      end
      
      LML1 = sort(LML); % relative tip path length
      for k = 1:19
        f = f+1;
        if ~isempty(LML1)
          F(f,i) = LML1(I(k));
        end
        if i == 1
          FN{f} = [num2str(RHN(k)),'%_relative_tip_path_length',str];
        end
      end
      
      
      % distribution comparisons:
      dr = histcounts(LML,0:0.05:1,'Normalization','probability');
      m = length(dr);
      % differences to the triangle distributions (start,mid,end)
      Sti = ceil(St/4*m);
      Mti = ceil(Mt/4*m);
      Eti = ceil(Et/4*m);
      for k = 1:length(St)
        di = abs(dr-triad(Sti(k),Mti(k),Eti(k),m));
        % mean and max difference to the comparison distribution:
        f = f+1;
        if ~isempty(di)
          F(f,i) = mean(di);
          F(f+1,i) = max(di);
        end
        if i == 1
          FN{f} = ['Tip_path_length_triad_',num2str(St(k)),...
            '_',num2str(Mt(k)),'_',num2str(Et(k)),'_mean'];
          FN{f+1} = ['Tip_path_length_triad_',num2str(St(k)),'_',...
            num2str(Mt(k)),'_',num2str(Et(k)),'_max'];
        end
        f = f+1;
      end
      
      % differences to the normal distributions (mean,std)
      Mni = ceil(Mn/4*m);
      Sni = ceil(Sn/5*m);
      for k = 1:length(Mn)
        di = abs(dr-normd(Mni(k),Sni(k),m));
        % mean and max difference to the comparison distribution:
        f = f+1;
        if ~isempty(di)
          F(f,i) = mean(di);
          F(f+1,i) = max(di);
        end
        if i == 1
          FN{f} = ['Tip_path_length_normd_',num2str(Mn(k)),'_',...
            num2str(Sn(k)),'_mean'];
          FN{f+1} = ['Tip_path_length_normd_',num2str(Mn(k)),'_',...
            num2str(Sn(k)),'_max'];
        end
        f = f+1;
      end
      
      % difference to the uniform and step distributions
      Sui = ceil(Su/5*m);
      Eui = ceil(Eu/5*m);
      for k = 1:length(Su)
        di = abs(dr-unid(Sui(k),Eui(k),m));
        % mean and max difference to the comparison distribution:
        f = f+1;
        if ~isempty(di)
          F(f,i) = mean(di);
          F(f+1,i) = max(di);
        end
        if i == 1
          FN{f} = ['Tip_path_length_unid_',num2str(Su(k)),'_',...
            num2str(Eu(k)),'_mean'];
          FN{f+1} = ['Tip_path_length_unid_',num2str(Su(k)),'_',...
            num2str(Eu(k)),'_max'];
        end
        f = f+1;
      end
      
      
      % Relations between different branching orders:
      for k = j:3
        TPLI = PL(c.extension == 0 & c.BranchOrder == k+1);
        stri = ['ord',num2str(k)];
        f = f+1;
        if ~isempty(TPLI)
          F(f:f+6,i) = [mean(TPL)/mean(TPLI); max(TPL)/mean(TPLI); ...
            mean(TPL)/max(TPLI); min(TPL)/mean(TPLI); mean(TPL)/min(TPLI);...
            min(TPL)/max(TPLI); max(TPL)/min(TPLI)];
        end
        if i == 1
          FN{f} = ['mean',str,'_tip_path_length/mean_',stri,...
            '_tip_path_length'];
          FN{f+1} = ['max',str,'_tip_path_length/mean_',stri,...
            '_tip_path_length'];
          FN{f+2} = ['mean',str,'_tip_path_length/max_',stri,...
            '_tip_path_length'];
          FN{f+3} = ['min',str,'_tip_path_length/mean_',stri,...
            '_tip_path_length'];
          FN{f+4} = ['mean',str,'_tip_path_length/min_',stri,...
            '_tip_path_length'];
          FN{f+5} = ['min',str,'_tip_path_length/max_',stri,...
            '_tip_path_length'];
          FN{f+6} = ['max',str,'_tip_path_length/min_',stri,...
            '_tip_path_length'];
        end
        f = f+6;
      end
    end
    
    % Branch base path lengths to the tree's base
    PLB = PL-c.length;
    % Compute the heights of the cylinder bases:
    H = c.start(:,3)-min(c.start(:,3));
    for j = 1:4
      % Path lengths for base cylinders:
      if j == 1
        BPL = PLB(c.PositionInBranch == 1);
        BH = H(c.PositionInBranch == 1);
        str = '_all';
      else
        BPL = PLB(c.PositionInBranch == 1 & c.BranchOrder == j-1);
        BH = H(c.PositionInBranch == 1 & c.BranchOrder == j-1);
        str = ['_ord',num2str(j-1)];
      end
      
      % Basic statistics:
      f = f+1;
      if ~isempty(BPL)
        F(f:f+6,i) = [mean(BPL); min(BPL); max(BPL); min(BPL)/max(BPL);...
          mean(TPL)/max(BPL); std(BPL); std(BPL)/mean(BPL)];
      end
      if i == 1
        FN{f} = ['mean_base_path_length',str];
        FN{f+1} = ['min_base_path_length',str];
        FN{f+2} = ['max_base_path_length',str];
        FN{f+3} = ['min/max_base_path_length',str];
        FN{f+4} = ['mean/max_base_path_length',str];
        FN{f+5} = ['std_base_path_length',str];
        FN{f+6} = ['std/max_base_path_length',str];
      end
      f = f+6;
      
      % Base path lengths relative to the base height:
      f = f+1;
      if ~isempty(BPL)
        LH = BPL./BH;
        LH(1) = mean(LH(2:end));
        F(f:f+6,i) = [mean(LH); min(LH); max(LH); min(LH)/max(LH);...
          mean(LH)/max(LH); std(LH); std(LH)/mean(LH)];
      end
      if i == 1
        FN{f} = ['mean(base_path_length/base_height)',str];
        FN{f+1} = ['min(base_path_length/base_height)',str];
        FN{f+2} = ['max(base_path_length/base_height)',str];
        FN{f+3} = ['min/max(base_path_length/base_height)',str];
        FN{f+4} = ['mean/max(base_path_length/base_height)',str];
        FN{f+5} = ['std(base_path_length/base_height)',str];
        FN{f+6} = ['std/max(base_path_length/base_height)',str];
      end
      f = f+6;
      
      % Base path lengths relative to the tree height:
      f = f+1;
      if ~isempty(BPL)
        LTH = BPL/t.TreeHeight;
        F(f:f+6,i) = [mean(LTH); min(LTH); max(LTH); min(LTH)/max(LTH);...
          mean(LTH)/max(LTH); std(LTH); std(LTH)/mean(LTH)];
      end
      if i == 1
        FN{f} = ['mean(base_path_length/tree_height)',str];
        FN{f+1} = ['min(base_path_length/tree_height)',str];
        FN{f+2} = ['max(base_path_length/tree_height)',str];
        FN{f+3} = ['min/max(base_path_length/tree_height)',str];
        FN{f+4} = ['mean/max(base_path_length/tree_height)',str];
        FN{f+5} = ['std(base_path_length/tree_height)',str];
        FN{f+6} = ['std/max(base_path_length/tree_height)',str];
      end
      f = f+6;
      
      % Base path lengths relative to the maximum base path length:
      f = f+1;
      if ~isempty(BPL)
        LML = BPL/max(BPL);
        F(f:f+6,i) = [mean(LML); min(LML); max(LML); min(LML)/max(LML);...
          mean(LML)/max(LML); std(LML); std(LML)/mean(LML)];
      end
      if i == 1
        FN{f} = ['mean(relative_base_path_length)',str];
        FN{f+1} = ['min(relative_base_path_length)',str];
        FN{f+2} = ['max(relative_base_path_length)',str];
        FN{f+3} = ['min/max(relative_base_path_length)',str];
        FN{f+4} = ['mean/max(relative_base_path_length)',str];
        FN{f+5} = ['std(relative_base_path_length)',str];
        FN{f+6} = ['std/max(relative_base_path_length)',str];
      end
      f = f+6;
      
      % Percentile values:
      BPL1 = sort(BPL); % Absolute base path length
      n = length(BPL);
      I = floor(RH*n)+1;
      for k = 1:19
        f = f+1;
        if ~isempty(BPL1)
          F(f,i) = BPL1(I(k));
        end
        if i == 1
          FN{f} = [num2str(RHN(k)),'%_base_path_length',str];
        end
      end
      
      LH1 = sort(LH); % base path length relative base height
      for k = 1:19
        f = f+1;
        if ~isempty(LH1)
          F(f,i) = LH1(I(k));
        end
        if i == 1
          FN{f} = [num2str(RHN(k)),'%_base_path_length/base_heigth',str];
        end
      end
      
      LTH1 = sort(LTH); % base path length relative tree height
      for k = 1:19
        f = f+1;
        if ~isempty(LTH1)
          F(f,i) = LTH1(I(k));
        end
        if i == 1
          FN{f} = [num2str(RHN(k)),'%_base_path_length/tree_height',str];
        end
      end
      
      LML1 = sort(LML); % relative base path length
      for k = 1:19
        f = f+1;
        if ~isempty(LML1)
          F(f,i) = LML1(I(k));
        end
        if i == 1
          FN{f} = [num2str(RHN(k)),'%_relative_base_path_length',str];
        end
      end
      
      % distribution comparisons:
      dr = histcounts(LML,0:0.05:1,'Normalization','probability');
      m = length(dr);
      % differences to the triangle distributions (start,mid,end)
      Sti = ceil(St/4*m);
      Mti = ceil(Mt/4*m);
      Eti = ceil(Et/4*m);
      for k = 1:length(St)
        di = abs(dr-triad(Sti(k),Mti(k),Eti(k),m));
        % mean and max difference to the comparison distribution:
        f = f+1;
        if ~isempty(di)
          F(f,i) = mean(di);
          F(f+1,i) = max(di);
        end
        if i == 1
          FN{f} = ['Base_path_length_triad_',num2str(St(k)),...
            '_',num2str(Mt(k)),'_',num2str(Et(k)),'_mean'];
          FN{f+1} = ['Base_path_length_triad_',num2str(St(k)),'_',...
            num2str(Mt(k)),'_',num2str(Et(k)),'_max'];
        end
        f = f+1;
      end
      
      % differences to the normal distributions (mean,std)
      Mni = ceil(Mn/4*m);
      Sni = ceil(Sn/5*m);
      for k = 1:length(Mn)
        di = abs(dr-normd(Mni(k),Sni(k),m));
        % mean and max difference to the comparison distribution:
        f = f+1;
        if ~isempty(di)
          F(f,i) = mean(di);
          F(f+1,i) = max(di);
        end
        if i == 1
          FN{f} = ['Base_path_length_normd_',num2str(Mn(k)),'_',...
            num2str(Sn(k)),'_mean'];
          FN{f+1} = ['Base_path_length_normd_',num2str(Mn(k)),'_',...
            num2str(Sn(k)),'_max'];
        end
        f = f+1;
      end
      
      % difference to the uniform and step distributions
      Sui = ceil(Su/5*m);
      Eui = ceil(Eu/5*m);
      for k = 1:length(Su)
        di = abs(dr-unid(Sui(k),Eui(k),m));
        % mean and max difference to the comparison distribution:
        f = f+1;
        if ~isempty(di)
          F(f,i) = mean(di);
          F(f+1,i) = max(di);
        end
        if i == 1
          FN{f} = ['Base_path_length_unid_',num2str(Su(k)),'_',...
            num2str(Eu(k)),'_mean'];
          FN{f+1} = ['Base_path_length_unid_',num2str(Su(k)),'_',...
            num2str(Eu(k)),'_max'];
        end
        f = f+1;
      end
      
      % Relations between different branching orders:
      for k = j:3
        BPLI = PLB(c.PositionInBranch == 1 & c.BranchOrder == k+1);
        stri = ['ord',num2str(k)];
        f = f+1;
        if ~isempty(BPLI)
          F(f:f+6,i) = [mean(BPL)/mean(BPLI); max(BPL)/mean(BPLI); ...
            mean(BPL)/max(BPLI); min(BPL)/mean(BPLI); mean(BPL)/min(BPLI);...
            min(BPL)/max(BPLI); max(BPL)/min(BPLI)];
        end
        if i == 1
          FN{f} = ['mean',str,'_base_path_length/mean_',stri,...
            '_base_path_length'];
          FN{f+1} = ['max',str,'_base_path_length/mean_',stri,...
            '_base_path_length'];
          FN{f+2} = ['mean',str,'_base_path_length/max_',stri,...
            '_base_path_length'];
          FN{f+3} = ['min',str,'_base_path_length/mean_',stri,...
            '_base_path_length'];
          FN{f+4} = ['mean',str,'_base_path_length/min_',stri,...
            '_base_path_length'];
          FN{f+5} = ['min',str,'_base_path_length/max_',stri,...
            '_base_path_length'];
          FN{f+6} = ['max',str,'_base_path_length/min_',stri,...
            '_base_path_length'];
        end
        f = f+6;
      end
    end
    
    % Tip path length divided by base path lengths:
    for j = 1:4
      
      % Basic statistics:
      f = f+1;
      if ~isempty(TPL)
        TB = TPL./BPL;
        TB(1) = mean(TB(2:end));
        F(f:f+6,i) = [mean(TB); min(TB); max(TB); min(TB)/max(TB);...
          mean(TB)/max(TB); std(TB); std(TB)/mean(TB)];
      end
      if i == 1
        FN{f} = ['mean(tip_path_length/base_path_length)',str];
        FN{f+1} = ['min(tip_path_length/base_path_length)',str];
        FN{f+2} = ['max(tip_path_length/base_path_length)',str];
        FN{f+3} = ['min/max(tip_path_length/base_path_length)',str];
        FN{f+4} = ['mean/max(tip_path_length/base_path_length)',str];
        FN{f+5} = ['std(tip_path_length/base_path_length)',str];
        FN{f+6} = ['std/max(tip_path_length/base_path_length)',str];
      end
      f = f+6;
      
      % Distribution values:
      TB1 = sort(TB);
      n = length(TB);
      I = floor(RH*n)+1;
      for k = 1:19
        f = f+1;
        if ~isempty(TB1)
          F(f,i) = TB1(I(k));
        end
        if i == 1
          FN{f} = [num2str(RHN(k)),...
            '%_(tip_path_length/base_path_length)',str];
        end
      end
    end
  end
end
f
% Define the outputs:
FeatureValues = F(1:f,:);
FeatureNames = FN(1:f);

end % End of main function



function d = normd(m,v,n)
x = (1/n:1/n:1);
G = normpdf(x,m/n,v/n);
d = G/sum(G);
end

function d = triad(a,c,b,n)
x = (0.5:1:n-0.5);
T = zeros(1,n);
T(a+1:c) = 2*(x(a+1:c)-a)/(c-a)/(b-a);
T(c+1:b) = 2*(b-x(c+1:b))/(b-a)/(b-c);
d = T/sum(T);
end

function d = unid(a,b,n)
c = b-a;
d = [zeros(1,a) 1/c*ones(1,c) zeros(1,n-a-c)];
end

function N = relative_height(d,N)
h = length(d);
k = 1;
n = length(N);
for j = 1:h
  for i = k:n
    if k <= n && d(j) > N(k)
      N(k) = j;
      k = k+1;
    end
  end
end
N = N/h;
end
