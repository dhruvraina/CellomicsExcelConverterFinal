%Author: Dhruv Raina; Back-End for Excel Converter.
%Last Edit: 051215
%Function: backend for excelconv.m GUI. Reads a Cellomics Excel file and
%          outputs averages of each column. Also prints heatmaps, histograms
%          and '% cells > n spots'

%Variables:
%calc_inputs     = (avg/var flag(1), Print Histogram(2) , Heatmap flag(3), Print in
%                    Excel(4), Run Correlation Fn(5), cmdMode(in loop)(6), cmd Mode
%                    (out loop)(7), <100 cells(8), SpotNum predictor(9), Custom Histogram Column(10), SparseFieldFilter(11),
%                     p-ValueFlag(12), pValueColumnA (13), pValueColumnB (14))
%reader_inputs   = (SpotCount Col(1), Param1(2), Param2(3), Param3(4), SpotCountNumber(5))
%pf_inputs       = (pf Flag(1), pf Col Num(2), pf Cutoff Value(3))
%SuperScatInputs = (SelectionFlag(1), colX(2), colY(3), treatmentID1(4),
%                    treatmentID2...treatmentIDn)

%To Do List:  3. Sparse Fields Filter; 4. Batch Mode

function excelconv_bend(calc_inputs, reader_inputs, pf_inputs, SuperScatInputs , fullfilepath, matflag, PmapMemFlag, bat_pmapname)
%% 1. Inputs and Init:
clc;
calc_inputs(11) = 0; %Temporary; sparse fields filter is set to 'off'
histtype = 2; %(1) - fourth root transformed hist, (2) - untransformed hist; NOTE:- No GUI inputs for this.

%1.1 Sparse Fields Cutoff:
sparsenum = 20;

%1.2 Figure Properties for saving images (only for hgexport, not with
%'print'):
%Look in Legacy Code
             
%1.3 'Running' Msgbox and Spot Number Predictor:
if reader_inputs(1)==0
    handle_msgbox1  = msgbox('''% cells > n spots'' will not be recorded. Running...   ', 'ExcelConv');
    spotpredict_col = 0;
else
    handle_msgbox1  = msgbox('Running.. Please Wait', 'ExcelConv');
    spotpredict_col = reader_inputs(1)-1;                                  % remove header
end

%1.4 Getting excel file paths:
disp('Loading Excel File..')
[xlnumbers, xlstrings] = xlsread(fullfilepath);
pos1                   = strfind(fullfilepath, '/');
pos1                   = pos1(end);
filename               = fullfilepath(pos1+1:end-5);
filepath               = fullfilepath(1:pos1);
fpsave                 = filepath;                                         %to save filepath in the pmap savefile. Saving the var. directly causes loaded filepath to overwrite runtime filepath
disp('...')

%1.5 Getting the platemap using 'pmap_ed.m'
switch PmapMemFlag(1)
    case(1)%Running from GUI
        if exist([filepath 'Platemap_save.mat'], 'file')==2
            load([filepath 'Platemap_save.mat']);
            pmap              = platemap_savefile;
        else
            pmap              = pmap_ed;
            msgbox('Platemap does not exist!', 'ExcelConv')
            save([filepath 'Platemap_Save.mat'], 'platemap_savefile', 'fpsave')
        end
        
    case(2)%BatchMode - No GUI feed
        if exist([filepath bat_pmapname], 'file')==2
            load([filepath bat_pmapname]);
            pmap              = platemap_savefile;
        else
            pmap              = pmap_ed;
            platemap_savefile = pmap;
            msgbox('Platemap does not exist!', 'ExcelConv')
            save([filepath 'Platemap_Save.mat'], 'platemap_savefile', 'fpsave')
        end
        
    otherwise
        pmap                 = pmap_ed;
        platemap_savefile    = pmap;
        save([filepath 'Platemap_Save.mat'], 'platemap_savefile', 'fpsave')
end

%1.6 Initializing Input Variables:
[platemap, platemap1, platemap2, platemap3, platemap4, platemap5] = deal(zeros(8, 12));
debug_raw_recorder  = cell(8, 12);
histmap1{8,12}      = 0;
platealph           = {'A', 'B', 'C', 'D', 'E', 'F', 'G', 'H'};
zerocol             = size(xlnumbers,2)+1;                                 % Pad xlnumbers to avoid 'out of bound' errors when reader_inputs=0;
xlnumbers(:,zerocol)= zeros;
reader_inputs1      = reader_inputs;
reader_inputs       = reader_inputs-1;                                     % 1st column is a string (well ID), so it gets removed by xlsread.
reader_inputs(5)    = reader_inputs1(5);                                   % reader_inputs(5) is the spotNum position; you don't want to -1 from that.
pf_colnum           = pf_inputs(2)-1;
spotNum             = reader_inputs(5);
hist_flag           = calc_inputs(2);
CorrFlag            = calc_inputs(5);
excludeFlag         = calc_inputs(8);
debug_filt_idx      = [];

%1.7 Enter the loop flag:
runplate_flag = sum(reader_inputs1(1:4))+ hist_flag;

%1.8 Error Handling and Warnings for inputs:
if reader_inputs(1)<=0 || isnan(reader_inputs(1))                          % 'Inactive' buttons get assigned zero values to prevent errors
    reader_inputs(1) = zerocol;
    reader_inputs1(1)= zerocol;
end
if reader_inputs(2)<=0 || isnan(reader_inputs(2))
    reader_inputs(2) = zerocol;
    reader_inputs1(2)= zerocol;
end
if reader_inputs(3)<=0 || isnan(reader_inputs(3))
    reader_inputs(3) = zerocol;
    reader_inputs1(3)= zerocol;
end
if reader_inputs(4)<=0 || isnan(reader_inputs(4))
    reader_inputs(4) = zerocol;
    reader_inputs1(4)= zerocol;
end

if max(isnan(SuperScatInputs(1:3)))
    errordlg('Scatter Plot column numbers not entered!', 'Excel Converter')
end

%1.9 Finding the 'paramter-agnostic' column to use in defining cell number
%(wells with 0 cells are not recorded to avoid NaNs)
CellNoCol = reader_inputs(find(reader_inputs(1:4)<zerocol, 1));

disp(['Plate: ' filename ' loaded... \n Filtering Data...'])

%1.10 Don't display excess information in batchMode (PmapMemFlag==2)
if PmapMemFlag~=2 
msgbox('Please note: algorithm is changed from  ''%CellsWithSpots>spotNumber''   to now   ''%CellsWithSpots>=spotNumber''   . You will have to change the standardized No. of Foci cutoffs by +1 :::Message Dated: 8.11.15 - Dhruv', 'WARNING!')
end

%% 2. Main
%2.1 Scatter Plot Function:
%Removed; put in *after* recording data for pre-filter application

%2.2 Reading Excel File:
if runplate_flag~=0                                                            %Enter this loop either to print histograms or record information
    for ct1 = 1:length(platealph)
        for ct2 = 1:size(platemap,2)
            tID = [platealph{ct1}, num2str(ct2)];
            pID{ct1,ct2} = tID;
            tidx = strmatch(tID, xlstrings, 'Exact');                          %SLOW STEP: Match each alphanumeric : 'A1, B1'..etc in the excel sheet
            tidx = tidx-1;
            disp(tID)
            try                                                                %1st row in xl sheet is Header, gets removed during xlsread
                cellNo = xlnumbers(tidx, CellNoCol);                           %Estimate No. of Cells - parameter agnostic, looks for any 'ON' parameter
                pf_param_raw = xlnumbers(tidx, :);
            catch                                                              %#ok<*CTCH>
                errordlg('Please make sure your excel file has a header! Empty header row detected!')
            end
%TBA: Workaround to eliminate cellNo variable and reduce memory usage

         if ~isempty(cellNo)   
             
              %---2.21-PreFilter--%
                if pf_inputs(1)==1
                    pfvals = xlnumbers(tidx, pf_colnum);
                    pf1 = find(pfvals>pf_inputs(3));     %>                    %Used 'find' to return indicies
                    Filt_tidx = tidx(pf1);
                else
                    pf1 = 1;
                    pf_param_raw(pf1, :) = 1;
                    Filt_tidx = tidx;
                end
                
                %--2.22-Correlation Function--%                                %Inside loop because it calculates correlation on a 'per-well' basis
                ct1_alph = platealph(ct1);
                if CorrFlag==1
                    if exist('Corr_Vars.mat', 'file')==2
                        load('Corr_Vars.mat');
                    else
                        msgbox('Corr Vars do not exist!', 'Correlation Function')
                    end
                    corrslope(ct1, ct2) = corr_bend(filepath, xlnumbers, Filt_tidx, corr_vars, ct1_alph, ct2);
                end
                                
                %--2.23-Sparse Fields filter--% Not yet implemented
                if calc_inputs(11)==1
                    fieldvec = [xlnumbers(Filt_tidx)';xlnumbers(Filt_tidx, 2)'];
                    fieldvec = fieldvec';                                  %TBA: Capture the areas instead; express sparsity as a fraction of cellular area vs. image area
                    fieldveclist = unique(fieldvec(:,1));
                    for ct3 = fieldveclist'
                        fieldfiltMaxCell(ct3) = max(fieldvec((fieldvec(:,1)==ct3),2));
                    end
                    fieldfilt2 = find(fieldfiltMaxCell<sparsenum);
                    fieldvec = fieldvec(:,1);
                    fieldfilt3 = fieldvec.*(~ismember(fieldvec, fieldfilt2));
                    xlnumbers(~xlnumbers(Filt_tidx,1)==fieldfilt3,1)=NaN;               %Using NaN to avoid mean calculations getting affected
                    Filt_tidx = xlnumbers(~isnan(xlnumbers(Filt_tidx,1)));              %Updating 'Filt_tidx' - temp. running index to exclude NaNs
                end
                
                %--2.24-Saving data in platemaps--%
                    tRawspot = xlnumbers(Filt_tidx, reader_inputs(1));
                    tFiltspot = find(tRawspot>=spotNum);                                 %NOTE:- GREATER OR EQUAL TO spotnum!
                    platemap(ct1, ct2) = (length(tFiltspot)/length(tRawspot))*100;       %Param1: Perc. Cells; Note: Fraction WITHIN the prefiltered set (doesn't incl. filtered out data)
                 switch calc_inputs(1)
                  case(1) %Average Calculations
                    platemap1(ct1, ct2) = mean(xlnumbers(Filt_tidx, reader_inputs(2)));            %Param2 - Filtered
                    platemap2(ct1, ct2) = mean(xlnumbers(Filt_tidx, reader_inputs(3)));            %Param3 - Filtered
                    platemap3(ct1, ct2) = mean(nonzeros((xlnumbers(Filt_tidx, reader_inputs(4)))));%Param4 - Filtered + NONZEROS Filter%
                    platemap5(ct1, ct2) = mean((xlnumbers(Filt_tidx, reader_inputs(4))));          %Param4-  Filtered
                    %platemap5(ct1, ct2) = mean((xlnumbers(tidx, reader_inputs(4))));                        %Param4 - Unfiltered
                  case(2) %Variance Calculations
                    platemap1(ct1, ct2) = var(xlnumbers(Filt_tidx, reader_inputs(2)));             %Param2 - Filtered
                    platemap2(ct1, ct2) = var(xlnumbers(Filt_tidx, reader_inputs(3)));             %Param3 - Filtered
                    platemap3(ct1, ct2) = var(nonzeros((xlnumbers(Filt_tidx, reader_inputs(4))))); %Param4-  Filtered + NONZEROS Filter%
                    platemap5(ct1, ct2) = var((xlnumbers(Filt_tidx, reader_inputs(4))));           %Param4-  Filtered
                    %platemap5(ct1, ct2) = var((xlnumbers(tidx,reader_inputs(4))));                          %Param4 - Unfiltered
                 end
                    %platemap5(ct1, ct2) = (mean((xlnumbers(Filt_tidx, reader_inputs(4)))))/(var((xlnumbers(Filt_tidx, reader_inputs(4)))));
                    platemap4(ct1, ct2) = size(Filt_tidx,1);                                       %Number of Cells after prefilter
                    debug_raw_recorder{ct1, ct2} = xlnumbers(tidx,:);                              %Records ALL DATA in each well
                    debug_filt_idx = [debug_filt_idx; Filt_tidx];                                  %Records all indicies of all filtered data
                
                %--2.25-cmd mode: for testing new code in loop--%
                if calc_inputs(6)==1                                               %If cmd Mode 'in loop' is ON and cell number >1
                    keyboard
                end
                
                %--2.26-Assembling All Histogram Data--%
                if ~isdir([filepath 'histimg/'])                                   %Cannot create directory at hgexport
                    mkdir(filepath, 'histimg');
                end
                
                switch hist_flag
                  case (0) %Histogram is off
                        
                  case (1) %'Per-Well' histogram of 'Spot Count' Column
                              if max(tRawspot)>0
                                  for ct3 = 1:max(tRawspot)                         %DOESN'T RECORD CELLS WITH '0' SPOTS
                                      tempvec = find(tRawspot==ct3);
                                      t_histvec(ct3) = length(tempvec);
                                  end
                                  hist.figure = figure('Visible', 'off');
                                  
                                  switch histtype %No GUI input, defined at the start of this file
                                      case(1)%Fourth Root Transformed Hist:
                                          powerraw = tRawspot.^0.25;
                                          hist(powerraw)
 
                                      case(2)%Non-Transformed Hist:
                                          bar(t_histvec)
                                          set(gca, 'YLim', [1 40])
                                          set(gca, 'XLim', [1 30])
                                  end
                                  
                                  title([platealph(ct1) num2str(ct2)])
                                  print(hist.figure,[filepath 'histimg/' filename 'Hist_' char(platealph(ct1)) num2str(ct2)], '-dpng','-r200')
                                  close gcf
                              end

                              
                   case (2) %Record data for 'Consolidated' Histogram
                       try
                        histcolval = xlnumbers(Filt_tidx, calc_inputs(10)-1);
                       catch
                           errordlg('No ''Column Number'' entered for Consolidated Histogram! Please check the ''Print Histogram'' inputs.')
                       end
                        debug_hist_raw{ct1, ct2} = histcolval;                   %Storage variable
                        histcolval = 0;
                       
                 end %Histogram
                
         end%Empty Well Checker
        end%platemap COLUMN (loop)
    end%platemap ROW (loop)
    
    keyboard
    
    %% 3. Working with all recorded information, out of loop:
    %3.1 Excluding <100 cells from analysis:
    if excludeFlag==1
        pmap(platemap4<100)=0;                                             %Note that this changes the original runtime pmap variable! (the saved one obv. doesn't change)
    end
    
    %3.2 Calculating p-Values:
    switch calc_inputs(12) %PVal flag
        case 0
            pvalA = 0;
            pvalB = 0;
        case 1
            pvalA = platemap(pmap==calc_inputs(13));
            pvalB = platemap(pmap==calc_inputs(14));
        case 2
            pvalA = platemap1(pmap==calc_inputs(13));
            pvalB = platemap1(pmap==calc_inputs(14));
        case 3
            pvalA = platemap2(pmap==calc_inputs(13));
            pvalB = platemap2(pmap==calc_inputs(14));
        case 4
            pvalA = platemap3(pmap==calc_inputs(13));
            pvalB = platemap3(pmap==calc_inputs(14));
    end
    
   [~, pValRes, pval_ci, pvalstats] = ttest2(pvalA, pvalB, [], [], 'unequal');
   xlpval_d = {calc_inputs(13), mean(pvalA), pvalstats.sd(1), calc_inputs(14), mean(pvalB), ...
                pvalstats.sd(2), pval_ci(1), pval_ci(2), pValRes, pvalstats.tstat, pvalstats.df,0}; %vector length = 12
   xlpval_s = {'Treatment ID - Dataset1' 'Mean - Dataset1' 'St. Dev - Dataset1' ... 
               'Treatment ID - Dataset2' 'Mean - Dataset2' 'St. Dev - Dataset2' ...
               '95% Confidence Interval Between Difference of Means - Lower' ...
               '95% Confidence Interval Between Difference of Means - Upper' ...
               'p-Value' 't-Test Statistic' 'Degrees of Freedom' ''}; %vector length = 12; 
                      
   %3.3  Graphs:
    %--3.31 Applying Prefilter to Data for Graphs--%:
    filt_xlnumbers = xlnumbers(debug_filt_idx, :);
    filt_xlstrings = {xlstrings{debug_filt_idx,1}}';
    filt_xlstrings{1,:} = {xlstrings{1,:}};
    
    %--3.32 Scatter Plot Function--%
    switch SuperScatInputs(1)
        case (2) %All treatmentIDs
            tID_vec_ss = nonzeros(unique(pmap(:)));
            SuperScatInputs = [SuperScatInputs'; tID_vec_ss]';
            superscat(filt_xlnumbers, filt_xlstrings, SuperScatInputs, pmap, filepath, filename);
        case (3) %User-Defined Treatment IDs
            superscat(filt_xlnumbers, filt_xlstrings, SuperScatInputs, pmap, filepath, filename);
    end
     %Conf. File
    load('F:\Dhruv\CCBT2\CellomicsExcelConverter\conf_g.txt');
    if conf_g~=145967
        errordlg('Improper Installation', 'Excel Converter')
        exit
    end

    %--%3.33 Drawing Consolidated Histogram--%:                                  %'Out of loop' to scale graph between all treatments
    if hist_flag==2
        
        %----3.33.1-Scaling Histogram and Bin Widths--%
        binno = 25;
        hist_vals_mat = debug_hist_raw(pmap>0);
        [binwidth binmax] = gscale(hist_vals_mat, binno, 'cell');          %gscale.m scales values appropriately for Consolidated Hist
        bins =0:binwidth:binmax;
        ct_nxt = 1;
        ct_histvec = unique(nonzeros(pmap))';                              %Stored in var to later retrieve treatment IDs
        
        %----3.33.2-Saving Hist. values for each treatment--%
        for ct_hist = ct_histvec
            hist_pmap = debug_hist_raw{pmap==ct_hist};                     %Note: Debug_hist_raw is after 'pre-filter' application
            counts(ct_nxt, :) = histc(hist_pmap,bins);
            counts_Norm(ct_nxt, :) = (counts(ct_nxt,:)./sum(counts(ct_nxt,:)).*100); %Normalized values if you so choose
            ct_nxt = ct_nxt+1;
        end
              
       %----3.33 - Individual Histograms after filtering--% 
       %LEGACY CODE
        
        %----3.33.3-Setting figure paramters, and saving--%
        chist.figure = figure('Visible', 'off')
        bar3(counts_Norm', 0.5)
        colormap(flipud(colormap('Hot')))
        ylabel(xlstrings(1, calc_inputs(10)))
        title([filename ' Consolidated Histogram for ' xlstrings{1, calc_inputs(10)}])
        set(gca, 'YTick', 0:5:26, 'XTickLabel',num2str(ct_histvec'), 'YTickLabel', bins(1:5:end)', 'FontSize', 6)
        childHandle = get(gca,'Children');                                 %Get GCA child handles, 'FaceAlpha' is the child that controls transparency
        set(childHandle,'FaceAlpha',0.5);                                  %0 = transparent, 1 = opaque.
        view([-53.5 30])                                                   %setting Azimuth and Elevation
        if exist([filepath filename '_Consolidated_Hist.png'], 'file')==2  %'Print' is not overwriting existing files, and doesn't throw an error to indicate this behaviour!
            delete([filepath filename '_Consolidated_Hist.png']);
        end
        print(chist.figure,[filepath filename '_Consolidated_Hist'], '-dpng','-r200')
        close(gcf)
    end
    
    %3.4 Associating values with treatment IDs from the platemap
    if max(pmap(:))>0
        resmat_con = {'Treatment Number' 'P1:Mean Perc' 'P1: Std Perc' 'Mean P2_prefiltered' 'Std P2' 'Mean P3_prefiltered' 'Std P3' 'Mean P4_NonZeros' 'Std P4' 'Mean P4_Unfiltered' 'Std P4_Unfiltered' '-'}; %extra spaces for vertcat
        for ct4 = 1:max(pmap(:))                                           %Max number of conditions
            ct3 = ct4+1;                                                   %Header
            pmap_idx = find(pmap==ct4);
            resmat_con{ct3,1} = ct4; 
            resmat_con{ct3,2} = mean(nonzeros(platemap(pmap_idx)));        % Spot Perc
            resmat_con{ct3,3} = std(nonzeros(platemap(pmap_idx)));         
            resmat_con{ct3,4} = mean(nonzeros(platemap1(pmap_idx)));       % Param2
            resmat_con{ct3,5} = std(nonzeros(platemap1(pmap_idx)));
            resmat_con{ct3,6} = mean(nonzeros(platemap2(pmap_idx)));       % Param3
            resmat_con{ct3,7} = std(nonzeros(platemap2(pmap_idx)));
            resmat_con{ct3,8} = mean(nonzeros(platemap3(pmap_idx)));       % Param4 - Filtered
            resmat_con{ct3,9} = std(nonzeros(platemap3(pmap_idx)));
            resmat_con{ct3,10} = mean(nonzeros(platemap5(pmap_idx)));      % Param4 - Unfiltered
            resmat_con{ct3,11} = std(nonzeros(platemap5(pmap_idx)));
            debug_raw_reshuffle{ct4} = vertcat(debug_raw_recorder{pmap_idx});  %Putting all treatment group values together
        end
    else
        resmat_con = num2cell(zeros(1,12));                                    %For concat @ ExcelWriting
        debug_raw_reshuffle = {(zeros(1,12))};                                 %In case nothing is recorded (e.g running Scatter Plot only) - otherwise errors out at spotnum predictor
    end
    
    %% 4. Outputs:
    disp('Saving Outputs')
    
    %4.1 Printing Heatmaps
    if calc_inputs(3)==1                                                       %Heatmap flag
       
        %--4.1.1-Number of Cells--%
        h.figure = figure('Visible', 'off');
        imagesc(platemap4)
        set(gca, 'YTickLabel', platealph, 'Xlim', [0.5 12.5])
        set(gcf, 'Position', [0 0 600 300])
        colorbar
        title('Number of Cells per Well')
        print(h.figure,[filepath filename 'reformat_CellsPerWell'], '-dpng','-r200')
        close gcf
        
        %--4.1.2-Perc. Spots.--%
        if reader_inputs(1)~=zerocol
            p.figure = figure('Visible', 'off');
            imagesc(platemap)
            set(gca, 'YTickLabel', platealph, 'Xlim', [0.5 12.5])
            set(gcf, 'Position', [0 0 600 300])
            colorbar
            title(['% cells >' num2str(spotNum) ' spots'])
            if exist([filepath filename 'reformat_spotperc'], 'file')==2    %Print is not overwriting existing files, doesn't throw error to indicate this behaviour
                delete([filepath filename 'reformat_spotperc.png']);
            end
            print(p.figure,[filepath filename 'reformat_spotperc'], '-dpng','-r200')
            close gcf
        end
        
        %--4.1.3-Param1--%
        if reader_inputs(2)~=zerocol
            p1.figure = figure('Visible', 'off');
            imagesc(platemap1)
            set(gca, 'YTickLabel', platealph, 'Xlim', [0.5 12.5])
            set(gcf, 'Position', [0 0 600 300])
            colorbar
            title(xlstrings{1,reader_inputs1(2)})
            if exist([filepath filename 'reformat_param1'], 'file')==2      %Print is not overwriting existing files
                delete([filepath filename 'reformat_param1.png']);
            end
            print(p1.figure,[filepath filename 'reformat_param1'], '-dpng','-r200')
            close gcf
        end
        
        %--4.1.4-Param2--%
        if reader_inputs(3)~=zerocol
            p2.figure = figure('Visible', 'off');
            imagesc(platemap2)
            set(gca, 'YTickLabel', platealph, 'Xlim', [0.5 12.5])
            set(gcf, 'Position', [0 0 600 300])
            colorbar
            title(xlstrings{1,reader_inputs1(3)})
            if exist([filepath filename 'reformat_param2'], 'file')==2      %Print is not overwriting existing files
                delete([filepath filename 'reformat_param2.png']);
            end
            print(p2.figure,[filepath filename 'reformat_param2'], '-dpng','-r200')
            close gcf
        end
        
        %--4.1.5-Param3--%
        if reader_inputs(4)~=zerocol
            p3.figure = figure('Visible', 'off');
            imagesc(platemap3)
            set(gca, 'YTickLabel', platealph, 'Xlim', [0.5 12.5])
            set(gcf, 'Position', [0 0 600 300])
            colorbar
            title(xlstrings{1,reader_inputs1(4)})
            if exist([filepath filename 'reformat_param3'], 'file')==2      %Print is not overwriting existing files
                delete([filepath filename 'reformat_param3.png']);
            end
            print(p3.figure,[filepath filename 'reformat_param3'], '-dpng','-r200')
            close gcf
        end
    end
    
    %4.2 Saving m-files
    if matflag==1
        save([filepath filename 'excelconv_vars.mat'], 'platemap', 'platemap1', 'platemap2', 'platemap3', 'platemap4', 'resmat_con', 'pmap' )
    end
    
    %4.3 Writing Files to Excel
    if calc_inputs(4)==1                                                   %xlswrite flag
        xlspacer = ones(4, ct2);
        xldata = [xlspacer; platemap; xlspacer; platemap1; xlspacer; platemap2; xlspacer; platemap3; xlspacer; platemap4; xlspacer; platemap5; xlspacer; pmap];
        tnan = isnan(xldata);
        xldata(tnan)=0;
        xldata = num2cell(xldata);
        xldata{4,1}  = ['P1: Perc. Cells >' num2str(spotNum) 'spots'];
        xldata{52,1} = ('Number Of Cells Recorded');
        xldata{64,1} = ('P4 - Unfiltered');
        xldata{76,1} = ('Platemap');
        try
            xldata{16,1} = ['P2: ' xlstrings{1,reader_inputs1(2)}];
            xldata{28,1} = ['P3: ' xlstrings{1,reader_inputs1(3)}];
            xldata{40,1} = ['P4: ' xlstrings{1,reader_inputs1(4)}];
        catch
            errordlg('Wrong Column Numbers entered. Please check the excel sheet again and verify the numbering!', 'Excel Converter')
        end
        keyboard
        xldata = [xldata; num2cell(xlspacer); xlpval_s; xlpval_d; num2cell(xlspacer); resmat_con];
        try
            xlswrite([filepath 'reformat_' filename '.xlsx'], xldata);
        catch
            uiwait(errordlg('The .xls file is either write protected or currently in use. Please close all instances and try again.', 'Excel Converter Error'))
            xlswrite([filepath 'reformat_' filename '.xlsx'], xldata);
        end
    end
    
    %% 5. End + Cleanup:
    
    %5.1 Spot Number Predictor:
    %--5.1.1-Finding the 'Control' conditions - typically either xx0 or xx1,(last digit ==0 or ==1)
    try                                                                    %Catch spotPredictor errors
        if calc_inputs(9)>0                                                %Spot Predictor Flag
            ctrl_finder = floor(length(debug_raw_reshuffle)/10);
            if ctrl_finder==0
                ctrl_vec   = 1;
            else
                ctrl_vec(1)= 1;
                
                for ct5 = 1:ctrl_finder
                    idxctrl  = ct5*10;
                    tempctrl = debug_raw_reshuffle(idxctrl);
                    
                    try
                        while sum(tempctrl{1,1}(:))<1 && idxctrl+1<=size(debug_raw_reshuffle,2)
                            idxctrl  = idxctrl+1;
                            tempctrl = debug_raw_reshuffle(idxctrl);
                        end
                    catch
                        disp('Unidentified Error in Spot Predictor. Type ''dbquit + Enter'' or  ''Ctrl+C'' to exit')
                        keyboard
                    end
                    
                    ctrl_vec(ct5+1) = idxctrl;                             % ctrl_vec contains the indicies for all 'control' treatment groups
                end
            end
            
            %--5.1.2-Showing you the closest value to the top 1 percentile:
            for ct6 = ctrl_vec
                normcontrol = sort(debug_raw_reshuffle{1,ct6}(:,reader_inputs(1))); %Getting raw spot data for tID = 1
                cval_1perc = ceil(0.01*length(normcontrol));
                cval_spot1 = length(normcontrol)-cval_1perc;                        %Can also be done as floor(0.99) directly
                ControlNormVal = normcontrol(cval_spot1);
                cval_threshperc = max(find(normcontrol==ControlNormVal));
                controlperc = (1-cval_threshperc/length(normcontrol))*100;
                TreatmentID = ct6                                                   %For Display
                Prediction = ControlNormVal                                         %For Display
                PredictedPercentage = controlperc                                   %For Display
            end
        end
    catch
        errordlg('Unidentified error in Spot Predictor. Please turn the prediction option OFF and try again', 'Excel Converter')
    end
end

%5.2 msgbox
if ishandle(handle_msgbox1)                                                % Check if msgbox is still open or user has closed it
    delete(handle_msgbox1);
    clear('handle_msgbox1');
end
handle_msgbox = msgbox(['File: ' filename '  Done.'],'Excel Converter');
disp(['File: ' filename '  Done.'])

%5.3 Analyze xCorrelation values
if CorrFlag==1 || calc_inputs(7)==1
figure, imagesc(corrslope)
colorbar;
keyboard                                                               % Analyzing Correlation inputs done manually
end

if PmapMemFlag==2 %Close extra windows on batchMode
close all
end

end %End


%LEGACY CODE:
% %1.2 Figure Properties for saving images (only for hgexport, not with 'print'):
% fStyle = hgexport('factorystyle');
% fStyle.Format = 'png';
% fStyle.Width = 5;
% fStyle.Height = 4;
% fStyle.Resolution = 300;
% fStyle.Units = 'inch';
% hgexport(h.figure, [filepath filename 'reformat_CellsPerWell'], fStyle);


%  %--3.21-Scaling Histogram and Bin Widths--%
%         binno = 25;
%         hist_vals_mat = cell2mat(debug_hist_raw(pmap>0));                 %Converts it to double, 10x slower
%         [binwidth binmax] = gscale(hist_vals_mat, binno, 'double');

%  %--3.23-Print Individual Histograms after filtering--%
%         for tet = 1:size(counts_Norm,1)
%             figure, bar(counts_Norm(tet,:))
%             set(gca, 'XTick', 0:5:26,'XTickLabel', bins(1:5:end)')
%             print(gcf,[filepath filename '_Ind_Hist' num2str(ct_histvec(tet))], '-dpng','-r200')
%             close(gcf)
%         end
        