% Name: Correlation Function Backend; corr_bend.m
% Usage: backend for corr_fe.m ; corr_bend(filepath, xlnumbers, temporary_set_index, corr_vars(saved in running dir), counter1, counter2)
% Outputs: Savefile jpg of scatter plot and line plot for param A and B
% Author: Dhruv
% Last Edit: 030515
% Calculates correlation between two sets of points.

function varargout = corr_bend(filepath, xlnumbers, tidx, corr_vars, ct1_alph, ct2) %tidx is the current running set of indicies. ct1 and ct2 are for associative naming.
%--1.1--Hard Coded Values (to be added to GUI in the future):
%corr_vars(3)=2;             %No Normalization
corr_vars(4)=0;              %Don't print/save graphs

%--1.2--Error Handling
if isnan(corr_vars(1))||isnan(corr_vars(2))
    errordlg('Incorrect Inputs For Correlation Function!', 'Correlation Funtion')
end

%--1.3--Raw Values
ColA = xlnumbers(tidx, corr_vars(1));
ColB = xlnumbers(tidx, corr_vars(2));

%--2.1--Normalized Data
switch corr_vars(3)
    case 1
        normA = min(nonzeros(ColA));
        normB = min(nonzeros(ColB));
    case 2
        normA = max(ColA);
        normB = max(ColB);
    case 3 
        normA = 1;
        normB = 1;
end

if length(normA)<1
    normA = 0;
end
if length(normB)<1
    normB = 0;
end

ColAn = ColA./normA;
ColBn = ColB./normB;

%--2.2--Getting rid of zeros to prevent false correlation rates:
ColBn(ColAn==0)=0;
ColAn(ColBn==0)=0;
ColAn = nonzeros(ColAn);
ColBn = nonzeros(ColBn);

%--3.1--Correlation Function:
%Legacy Code
corren = @(x,y) sum((x-mean(x(:))).*(y-mean(y(:))));
corr_v1 = corren(ColAn, ColBn)/(std(ColAn)*std(ColBn));
corr_val = corr_v1/(length(ColAn-1));

%--4.1--Plotting both curves:
%This may not even be necessary!! Graphing code below (2 blocks) may contain errors,
%has not been tested, and possibly will be discarded.
%Legacy Code

varargout{1} = corr_val;



%Legacy Code:
%--3.1--Correlation Function:
%Legacy Code
% corrn = @(x) (x-mean(x))./std(x);
% corr_val = (corrn(ColAn)'*corrn(ColBn))/sum(corrn(ColBn).^2);

%--4.1--Plotting both curves:
%if corr_vars(4)==1
%    Col(:,1) = ColAn;
%    Col(:,2) = ColBn;
%    col2 = sortrows(Col,1);
%   figure, plot(col2), set(gcf, 'Visible', 'off');
%   title([ct1_alph num2str(ct2)])
%    hgexport(gcf, [filepath 'CorrPlot_' ct1_alph{1} num2str(ct2) '.jpg'],hgexport('factorystyle'), 'Format', 'jpeg');
%   close gcf

    %Plotting normalized scatter, straight line fit and output slope:
%    ColAs = col2(:,1);
 %   ColBs = col2(:,2);
 %   figure, scatter(ColAs, ColBs), set(gcf, 'Visible', 'off');
 %   h = lsline;
 %   p = polyfit(get(h,'xdata'),get(h,'ydata'),1);
 %   title([ct1_alph num2str(ct2)])
 %   hgexport(gcf, [filepath 'CorrScat_' ct1_alph{1} num2str(ct2) '.jpg'],hgexport('factorystyle'), 'Format', 'jpeg');
 %   close gcf
%end