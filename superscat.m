%Author: Dhruv Raina; Description: Scatter Plot function
%Last Edit: 100515
%Function: Plots a scatter between multiple treatment IDs. Accepts several treatments simultaneously, then
%quantifies the difference between the scatters.(Called from excelconv_bend.m)

%InputVariables: SuperScatVars: [SelectionFlag, colX, colY, treatmentID1, treatmentID2...treatmentIDn]

function superscat(xlnumbers, xlstrings, SuperScatVars, pmap, filepath, filename)
%% 1. Inputs and Init:
%keyboard
%1.1 Disp
disp('Plotting Scatter...')

%1.2 Inputs:
outlierperc = 0.98;                                                        %percentile values to include in the graph, eliminates highest 1% of data values.
gridblock = 20;                                                            %size of the binning for graph gridblock point ratio calculations.

%1.3 96 wp wellnames:
alph = 'A':'H';
for ct1 = 1:8
    for ct2 = 1:12
        platealph{ct1,ct2} = [alph(ct1) num2str(ct2)];
    end
end

%1.4 Accounting for Header Column:
minusvec = zeros(size(SuperScatVars));
[minusvec(2), minusvec(3)] = deal(1);
SuperScatVars = SuperScatVars-minusvec;

%1.5 Reformatting Data: Pooling Treatment Groups:
for ct3 = 4:length(SuperScatVars)                                          %Superscatvars changes size depending on the number of treatment groups. 1:4 is reserved for other vals.
    ct31 = ct3-3;
    tidx = find(pmap==SuperScatVars(ct3));
    treatmentname{ct31} = ['Treatment Number: ' num2str(SuperScatVars(ct3))]; %delete? 
    welllist = platealph(tidx);
    concat_idx2 = 0;
    for ct4 = 1:length(welllist)
        [concat_idx1 ~] = find(strcmp(welllist{ct4}, xlstrings));          %10x speed improvement over legacy code 1.5!
        concat_idx2 = [concat_idx2; concat_idx1];
    end
    scatX{ct31} = xlnumbers(nonzeros(concat_idx2)-1, SuperScatVars(2));
    scatY{ct31} = xlnumbers(nonzeros(concat_idx2)-1, SuperScatVars(3));
end
try % this is the original for pre-filtered contour
    labelX = {xlstrings{1, SuperScatVars(2)+1}};
    labelY = {xlstrings{1, SuperScatVars(3)+1}};
catch %This is tester for post-filtered contour
    labelX = xlstrings{1,1} {SuperScatVars(2)+1};
    labelY = xlstrings{1,1} {SuperScatVars(3)+1};
end

%1.6 Catching erroneous inputs:
MTcatcher = cell2mat(cellfun(@length, scatX, 'UniformOutput', 0));         %Can also be done with scatY
if min(MTcatcher)==0
    errordlg('Scatter Plot Error: Treatment ID not found. Please check the treatment ID inputs or Platemap you provided!', 'Scatter Plot Error')
end


%% 2. Plotting Graphs:

%2.1 Merging all treatment groups to determine the best scale:
xlim3 = cellfun(@(x) nonzeros(sort(x)), scatX, 'UniformOutput', 0);        %For Upper X bound - Calculation
ylim3 = cellfun(@(x) nonzeros(sort(x)), scatY, 'UniformOutput', 0);        %For Upper Y bound - Calculation
%xlim1 = min(cellfun(@(x) min(x)-ceil(0.25*min(x)), xlim3));                %Lower X bound - ONLY FOR HIST3
%ylim1 = min(cellfun(@(x) min(x)-ceil(0.25*min(x)), ylim3));                %Lower Y bound - ONLY FOR HIST3

%2.1.1 Calculating upper bounds:
xlim3 = cellfun(@(x) x(nonzeros(floor(outlierperc*length(x)))), xlim3, 'UniformOutput', 0);  %For Upper X bound - Calculation
xlim4 = xlim3(cellfun(@(x) length(x), xlim3)==1);                                            %For Upper X bound - Calculation
xlim2 = max(cellfun(@max, xlim4));                                                           %Higher X bound
ylim3 = cellfun(@(x) x(nonzeros(floor(outlierperc*length(x)))), ylim3, 'UniformOutput', 0);  %For Upper Y bound - Calculation
ylim4 = ylim3(cellfun(@(x) length(x), ylim3)==1);                                            %For Upper Y bound - Calculation
ylim2 = max(cellfun(@max, ylim4));                                                           %Higher Y bound

% xlim1 = 0;                                                               %Debug
% ylim1 = 0;

%2.2 Finding Maximum Distance Gridblock:
[gridsizex, gridsizey] = deal(gridblock);                                  %Gridblock is arbitrary
gridxstep = floor((xlim2/gridsizex)*10)/10;
gridystep = floor((ylim2/gridsizey)*10)/10;
xfiltlow = 0; %xlim1;
yfiltlow = 0; %ylim1;

for ct8 = 1:ct31                                                           %length of superscatvars from 4:end
    
    for ct6 = 1:gridsizex
        xfilthigh = xfiltlow + gridxstep;
        filtidx1 =(scatX{ct8}>xfiltlow&scatX{ct8}<xfilthigh);
        
        for ct7 = 1:gridsizey
            yfilthigh = yfiltlow + gridystep;
            %--2.21-Filtering Data--%
            filtidx2 = scatY{ct8}>yfiltlow&scatY{ct8}<yfilthigh;
            respoints{ct8}(ct7, ct6) = sum(filtidx1.*filtidx2);
            reslimitmap(ct7, ct6) = xfiltlow;                              %Used to input xBins and yBins into xls sheet
            reslimitmap2(ct7, ct6) = yfiltlow;
            yfiltlow = yfilthigh;
            %xdat = scatX{ct8}(filtidx);                                   %%Recovers Actual Data Points if needed
            %ydat = scatY{ct8}(filtidx);
        end %Y loop
        
        yfiltlow = 0;
        xfiltlow = xfilthigh;
    end %X loop
    
    %--2.22-Reformating Data: arranging respoints to intuitively look like the graphs-%
    respoints{ct8} = flipud((respoints{ct8}./sum(respoints{ct8}(:))*100));
    save([filepath 'superscatpoints.mat'], 'respoints', 'reslimitmap')
    xfiltlow = 0;
end %treatment loop

%Bins for the Hist3 graphs:
% bin{1} = xlim1:gridxstep:xlim2;
% bin{2} = ylim1:gridystep:ylim2;

%2.3 Setting Graph Export Props:
%Under Legacy Code

%2.4 Drawing Graphs:
for ct5 = 1:length(scatX)
    %figure, hist3([scatX{ct5}, scatY{ct5}], bin) For Hist3 Usage
    figure('Visible', 'off'); 
    contour(reslimitmap, reslimitmap2, flipud(respoints{ct5}),20) %flipud because its been flipped once before for excel printing!
    xlabel(labelX)
    ylabel(labelY)
    title(['Treatment ID: ' num2str(SuperScatVars(ct5+3))])
    grid on
    print(gcf,[filepath filename '_Contour_tID' num2str(SuperScatVars(ct5+3))], '-dpng','-r200') %Print .png @ 200dpi
    close gcf
end

%2.5 Unpackaging data for xlswrite:
resmat = zeros(1,gridblock);
paddington = ones(2, gridblock);
resname = ones(1, gridblock);
resmat = [resmat; resmat; resmat; reslimitmap(1,:);resmat; reslimitmap2(:,1)']; %this format is for xls writing
for ct9 = 1:length(respoints)
    resmat = [resmat; paddington; resname; respoints{ct9}];
end

%2.6 Conversion to Cell and Adding Labels:
resmat = num2cell(resmat);
resmat{3,1} = 'X-Axis Bin Sizes';                                          %Note: not {'X-Axis...sizes'};
resmat{5,1} = 'Y-Axis Bin Sizes';
ct11 = 9;
for ct10 = 1:ct31
    resmat{ct11,1} = ['Treatment ID: ' num2str(SuperScatVars(3+ct10))];
    ct11 = ct11+23;
end

disp('...100%')

%2.7 Write to Excel
try
    xlswrite([filepath 'scatterpoints' filename '.xlsx'], resmat);
catch
    uiwait(msgbox('The Scatter Plot output file ''scatterpoints.xlsx'' is locked or being written by another process. Please close all instances and click OK', 'Scatter Plot Error'))
    xlswrite([filepath 'scatterpoints.xlsx'], resmat);
end

%Legacy Code:
%1.5
% concat_idx1 = strmatch( welllist{ct4}, xlstrings, 'Exact');

%2.1.05 Make sure all indicies in 2.1.1 are nonzero: (May not be required after
%nonzeros in 2.1.1
% xlim4 = find(cellfun(@(x) max(x)<=0, xlim3));                           %return indicies of length<=1
% ylim4 = find(cellfun(@(x) length(x)<=1, ylim3));                           
% if max(xlim4(:))>0 %is all length>1?
% [xlim3{xlim4}] = deal([max(xlim3);max(xlim3{xlim4})]);                     %Doubling values within length<=1
% [ylim3{ylim4}] = deal([max(ylim3{ylim4});max(ylim3{ylim4})]);
% end

%2.3 Setting Graph Export Props:
% fStyle = hgexport('factorystyle');
% fStyle.Format = 'png';
% fStyle.Width = 5;
% fStyle.Height = 4;
% fStyle.Resolution = 300;
% fStyle.Units = 'inch';
% hgexport(gcf, [filepath filename '_Contour_tID' num2str(SuperScatVars(ct5+3)) '.png'], fStyle, 'Format', 'png')
% s.fig.PaperUnits = 'inches';
% s.fig.PaperPosition = [0 0 6 3];
% s.fig.PaperPositionMode = 'manual'
