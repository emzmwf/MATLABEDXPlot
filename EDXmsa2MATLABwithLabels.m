%   EDX read and replot script
%   This version uses the labelpoints script of Adam Danz
%   http://www.mathworks.com/matlabcentral/fileexchange/46891-labelpoints

%Data handling - clear everything else
clc
clear all
close all

%   Input - either .txt or .msa file, dialog box will ask for file

%-----
%Ask for file name of msa format file, 
[filename,dirname] = uigetfile('*.*','Select msa format file');
fullpath = fullfile(dirname,filename);
%fid = fopen(fullpath,'rt');
%fclose(fid);
%clear fullpath fid datacell info

%   read the whole file in
data = fileread(fullpath);

%   Confirm first line matches expected specification
%   #FORMAT      : EMSA/MAS Spectral Data File
LnOne = '#FORMAT      : EMSA/MAS Spectral Data File';

fid = fopen(fullpath);
%   read in the next line, which will be first in this case
line_ex = fgetl(fid);

if strcmpi(LnOne, line_ex)
    Check = true
else
    Check = false
    f = msgbox('Incompatible format', 'Error','error');
    return
    %Currently ends the program, rather than ask them to reselect
end

%   Ask for the x range to be plotted
answer = inputdlg('keV range to plot to:',...
             'X axis range', [1])
kevlim = str2double(answer{:});
clear answer;
         



%   Read in the text file, identify data, and put in appropriate
%   arrays

%   Labels info all in format ##OXINSTLABEL: 13, 1.487, Al
%   Identify line from start ##OXINSTLABEL: 
lab = '##OXINSTLABEL: ';

%   remove first part 
%   put other three parts into an array
%   Currently putting into a cell, and running
%   str2double() when we need it - may not be optimal?

%   Get the number of OXINSTLABEL lines
Lablines = count(data, lab);

%   Set up the labels array, 3 by lablines
ElLab = cell(Lablines, 3);

%   Format is atomic no, keV, element

%   i = 1
%   ElLab{i, 1} = 13
%   ElLab{i, 2} = 1.487
%   ElLab{i, 3} = 'Al'

%   Set up Spectrum Array SpecArr
%   Could look up NPOINTS value from data
NPOINTS = 1024;
%   DAT = zeros(NPOINTS, 2);
%   This will have all the array values read in from the file
DAT = [];

%   Now parse through the file
%   If the line starts with ##OXINSTLABEL: send to the ElLab array
%   If there is no #, send to SpecArr

%   reiterate to the end of the text
done = 0;   %Stop when we get the end of spectrum line
labno = 1;  %Number of labels, starting at one
dno = 1;    %Number of data lines read in, starting at one
SpecTitle = "EDX Spectrum";

while done == 0
    line_next = fgetl(fid);
    if numel(line_next) >15
        LabChck = line_next(1:15);
        if line_next(1:6) == "#TITLE"
            SpecTitle = extractAfter(line_next,'#TITLE       : ');
        end
    else
        LabChck = '#';
    end
    if LabChck == lab
        %put the rest of the data into ElLab
        [a, b, c] = LabIn(line_next);
        ElLab{labno, 1} = a;
        ElLab{labno, 2} = b;
        ElLab{labno, 3} = c;
        labno = labno+1;
    elseif line_next(1:1) ~= '#'
        [d, e] = DIn(line_next);
        DAT{dno,1} = d;
        DAT{dno,2} = e;
        %put the data into the DAT array
        dno = dno+1;
    elseif strcmpi(line_next, '#ENDOFDATA   : ')
        %all is done, end the loop
        done = 1;
    end
end     % end of while loop

%   ----------
%   Plot the data

%       Format for  labelpoints.m is 
%         x = [0 .1 .2]; y = [.8 .7 .6];
%         xlabs = [0 .1 .2]; ylabs = [.8 .7 .6];
%         labels = {'label 1' 'label 2' 'label 3'};
%         plot(x, y, 'o'); axis([-1 1 0 1]);
%         labelpoints(xlabs, ylabs, labels);


%   All the data in the DAT array should be convertible to doubles
DAT = cell2mat(DAT);

x =  (DAT(:, 1));
y =  (DAT(:, 2));

%   Set up label array
%   labels = string(ElLab(:,3)); would do all of them, but we want to
%   select, so define blank arrays and fill selectively
xlabs = [];
ylabs = [];
labels = [];

%   Set up the figure and plot data
figure
EDXplot = plot(x,y)
title(SpecTitle)
xlabel('keV')
ylabel('counts')
xlim([-0.2 kevlim])
%ylim - max should be highest value after x = 0.1keV
xa = find(DAT==0.1,1);
L = length(y);
%   Find the max y value past the zero peak for plotting reasons
ymax = max(DAT(xa:L, 2:2));
ylim([0 ymax*1.1])

%   Find the min y value past the zero peak for labelling reasons
ymin = min(DAT(xa:L/2, 2:2));

%   Apply the Labels
alab = [];
for i = 1:Lablines
    t1 = ElLab{i,1};            % atomic no
    t2 = ElLab{i,2};            % keV
    t3 = string(ElLab(i,3));    %label
    yesno = 1;

    %ind1 = DAT(:,1) == t2;
    ind1 = DAT(:,2) == t2;
    A0 = DAT(ind1);
    
    row = find(DAT(:,1) >= t2);
    rowa = row(1);
    yval = DAT(rowa,2);
    %   if the value isn;t a lot bigger than the lowest noise value,
    %   there's no point labelling it
    yvmin = 20*ymin;
    if yval < yvmin
        yesno = 0;
    end
       
    if t2<0.1
        yesno = 0;
    end
      
    if yesno ==1 
        % append x location - t2 - to xlabs
        xlabs = [xlabs, t2];
        % append y location - yval - to ylabs
        ylabs = [ylabs, yval];
        % append label name - t3 - to labels
        labels = [labels, t3];
        blab = {t2, yval, t3};
        alab = [alab;blab];
%        tlab = text(t2, yval+yvmin, t3);
%        tlab.FontSize = 12;
        end
end

%Future - do some malarky to check if the label is the same and the keV is
%within ~150eV there's no point labelling both so remove the lower count
%line from the array
% Or group the data by half eV values and use the 'stacked', 'down' option
% on each group so they don't overlap
% would look very odd when zoomed in though


%   Sort the label array by ascending  keV
%alabsort = sortrows(alab, 1);
%TESTx = str2double(alabsort(:,1));
%TESTy = str2double(alabsort(:,2));
%TESTlab = string(alabsort(:,3));
%labelpoints(TESTx, TESTy, TESTlab, "N", 0.2);

%   labelpoints - N is directly above point, 0.2 is value offset in
%   direction
labelpoints(xlabs, ylabs, labels, "N", 0.2);

%fclose(fid);

%   ----------
%   Set up functions - 

%   Label reading function
function [ElNo, ElKv, ElName] = LabIn(linLab)
    %remove the label from the start
    newStr = extractAfter(linLab,'##OXINSTLABEL: ');
    % now make a 3x1 string
    outStr = split(newStr,[","]);
    ElNo = str2double(outStr(1));
    ElKv = str2double(outStr(2));
    ElName = outStr(3);
end


%   Data reading function
%   From #SPECTRUM    : Spectral Data Starts Here is two column data
%   If it has a # at the start, it isn't data
function [DkV, DCounts] = DIn(linDat)
doutStr = split(linDat,[","]);
    %   now make it numbers
dVals = str2double(doutStr);
DkV = dVals(1);
DCounts = dVals(2);
end
