function raincloudPlotForPaper(ax,data,center,width,col,kde)

% kde is a struct containing:
%%% kde.bw - bandwidth
%%% kde.xi - values of x at which kde should be assessed
%%% kde.support - if kde unbounded or not

% Extract the lines for the boxplot
data = sort(data);
q = quantile(data,[0.25 0.50 0.75]);
iqr = q(3)-q(1);
lower_fence = q(1)-(1.5*iqr);
upper_fence = q(3)+(1.5*iqr);
lower_tail = data(find(data>lower_fence,1,'first'));
upper_tail = data(find(data<upper_fence,1,'last'));
clear iqr

quartileWidth = 0.5;
tailWidth = 0.25;

% Run kde on data
if ~exist('kde','var')
    kde = struct;
end
if ~isfield(kde,'support')
    kde.support = 'unbounded';
end
if strcmpi(kde.support,'positive')  % In this case, we need to check that there aren't any negative values, and we should (as a "hack") set any 0 values to eps
    if any(data<0)
        error('Support was specified as "positive" but data contained negative values...')
    end
    if any(data==0)
        disp('Warning: some zero-values inside data.  These are being set to eps, in order to be compatible with "positive" support specified.');
        data(data==0) = data(data==0)+eps;
    end
end
if ~isfield(kde,'xi')
    [~,xi,~] = ksdensity(data,'support',kde.support);
    kde.xi = xi; clear xi;
end
if ~isfield(kde,'bw')
    [~,~,u] = ksdensity(data,kde.xi,'support',kde.support);
    kde.bw = u; clear u;
end

[f,~,~] = ksdensity(data,kde.xi,'bandwidth',kde.bw,'support',kde.support);


    
% Plot the lines of the boxplot
if isempty(ax)
    figure;
    ax = gca;
else 
end
hold on;


% Run the scatter
jitter = (rand(size(data))*width*0.9)+(center-(width*0.5*1.8));
ps = scatter(jitter,data);
ps.CData = col;
ps.Marker = '.';
ps.SizeData = 300;

% Plot KDE
fplot = f./(max(f)*(1/width)); % Scale f so that max(f) = width;
ppatch = patch([center+fplot(:);center;center;center+fplot(1)],[kde.xi(:);kde.xi(end);kde.xi(1);kde.xi(1)],col); % Make patch object
ppatch.FaceAlpha = 0.6;
pkde = plot(fplot+center,kde.xi);
pkde.Color = col; pkde.LineWidth = 2;
pcenter = plot([center,center],[kde.xi(1),kde.xi(end)]);
pcenter.Color = col; pcenter.LineWidth = 2;
pline = plot([center;center+fplot(1)],[kde.xi(1);kde.xi(1)]);
pline.Color = col; pline.LineWidth = 2;

%p1 = plot(ax,[center,center-(width*quartileWidth)],[q(1),q(1)]);
%p2 = plot(ax,[center,center-(width*quartileWidth)],[q(2),q(2)]);
%p3 = plot(ax,[center,center-(width*quartileWidth)],[q(3),q(3)]);
%p4 = plot(ax,[center,center-(width*tailWidth)],[lower_tail,lower_tail]);
%p5 = plot(ax,[center,center-(width*tailWidth)],[upper_tail,upper_tail]);
%p6 = plot(ax,[center,center],[q(1),q(3)],'k','LineWidth',3);
%p7 = plot(ax,[center,center],[lower_tail,upper_tail]);

% Change the properties of the lines
%p1.Color = col; p1.LineWidth = 2;
%p2.Color = col; p2.LineWidth = 4;
%p3.Color = col; p3.LineWidth = 2;
%p4.Color = col; p4.LineWidth = 2;
%p5.Color = col; p5.LineWidth = 2;
%p6.Color = 'k'; p6.LineWidth = 3; p6.LineJoin = 'round';
%p7.Color = col; p7.LineWidth = 2;

p6 = rectangle(ax,'Position',[center,q(1),0,q(3)-q(1)],'LineWidth',3,'Curvature',0.5);
median = scatter(center,q(2),50,'k','filled');
