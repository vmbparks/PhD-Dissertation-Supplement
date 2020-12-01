close;clc;clear;
StartFolder=dir('/Users/victoriabuford/Downloads/McQ02N/'); % First 3 blah

% Set default font parameters for figures & text labels
% doc AxesProperties
set(0,'DefaultAxesFontSize',12)
set(0,'DefaultTextFontSize',11)
set(0,'DefaultLegendFontSize',11)

set(0,'DefaultAxesFontName','Arial')
set(0,'DefaultTextFontName','Arial')
set(0,'DefaultLegendFontName','Arial')



method      = 'movmean';
window_size = 1;
XMax        = 1050;
XMin        = 330;


TchronDataFiles ={'/Users/victoriabuford/Box Sync/Research/MoveStuff/McQ02/McQ02N2_Tchron.txt'
    '/Users/victoriabuford/Box Sync/Research/MoveStuff/McQ02/McQ02N9_Tchron.txt'};

% fig = figure;
% %u = fig.Units;
% fig.Units = 'inches';
% set(fig,'visible','off')
% fig.Position= [0 -5 7.7778*2 5.8333*10];
% saveas(fig,'/Users/victoriabuford/Documents/MATLAB/Figure1.png')
% close

fig=figure(1);
set(gcf,'Position',[360   347   482   351])   


for kk=4:length(StartFolder)
    Folder_Ages=strcat(StartFolder(kk).folder,'/',StartFolder(kk).name);
    
    Model   = StartFolder(kk).name(1:7);
    ID      = StartFolder(kk).name(9:end);
    if strcmp(Model,'McQ02N1')==1
        TchronDataFile =TchronDataFiles{1};
        Model='McQ02N10';
        ID      = StartFolder(kk).name(10:end);
    elseif strcmp(Model,'McQ02N9')==1
        TchronDataFile =TchronDataFiles{2};
        
    end
    %%
    % %% Goodness of fit (prev analysis, from paper)
    GOF=table();
    GOF.velocity={'i'
        'ii'
        'iii'
        'iv'
        'v'
        'vi'
        'vii'
        'viii'
        'ix'
        'x'};
    GOF.model={'const'
        'VV7'
        'VV8'
        'VV17'
        'VV10'
        'VV5'
        'VV3'
        'VV4'
        'VV9'
        'VV6'};
    GOF.percent=[52.38	47.62
        47.62	52.38
        57.14	71.43
        52.38	61.90
        42.86	52.38
        47.62	71.43
        57.14	71.43
        52.38	71.43
        57.14	71.43
        57.14	61.90];
    
    
    if strcmp(Model,'McQ02N9')==1
        ModNum=4;
        ModCol=2;
    elseif strcmp(Model,'McQ02N10')==1
        ModNum=1;
        ModCol=1;
    else
        disp('model does not exist')
    end
    % assign velocity # & GOF%
    for IDD=1:length(GOF.model)
        if strcmp(GOF.model{IDD},ID)==1
            n=IDD;
            MODID=GOF.velocity{n};
            PERCFIT=GOF.percent(n,ModCol);
            break
        end
    end
    %set(gcf,'Size'
    %SubplotLoc=[1 2; 3 4; 5 6; 7 8; 9 10; 11 12; 13 14; 15 16; 17 18; 19 20];
    %subplot(10,2,SubplotLoc(n,ModCol))% instead of 10,length(StartFolder)-4
    %% %%%%%
    i=39;
    File_Age = importdata([Folder_Ages '/Ages_tec00' num2str(i) '.dat'], ' ', 4);
    File_Age.data((File_Age.data(:,4) > XMax),:) = [];
    File_Age.data((File_Age.data(:,4) < XMin),:) = [];
    File_Age.data((File_Age.data(:,5) ~= 2),:) = [];
    A=size(File_Age.data);
    for j=1:A(2)
        File_Age.data(isnan(File_Age.data(:,j)),:) = [];
    end
    
    % Sort Ages for line plot
    File_Age.data = sortrows(File_Age.data,4); % x.
    
    % Compute moving mean
    Age_Avg=[zeros(size(File_Age.data))];
    % id ; x ; y ; real x ; real y ; real z ; AHe ; AFT ; ZHe ; ZFT ; MAr;
    Age_Avg(:,4) = File_Age.data(:,4);
    %window_size=7;
    for k=7:A(2)
        Age_Avg(:,k)=smoothdata(File_Age.data(:,k),method,window_size);
    end
    %Age_Avg=File
    TchronData=readtable(TchronDataFile,'Delimiter','\t');
    
    % Define which data is which in the Tchron file
    MArind = strcmp(TchronData.Method,'MAr');
    ZFTind = strcmp(TchronData.Method,'ZFT');
    ZHeind = strcmp(TchronData.Method,'ZHe');
    AFTind = strcmp(TchronData.Method,'AFT');
    AHeind = strcmp(TchronData.Method,'AHe');
    % Matrix of thermochronometers sorted by Method
    index = [MArind ZFTind ZHeind AFTind AHeind];
    %                 if islogical(index)==1
    %                     index=double(index);
    %                 end
    [d,e]=size(index);
    
    %%% Goodness Of Fit
    %GOF_version=1.0;
    g=5;
    samplelocll=zeros(d,g);
    ModelVal=zeros(d,g);
    %Matches=zeros(d,g);
    ModelError=1; % in [Ma]
    %for mm=1:2
    XsecError=2; % in [km]
    CrossPlot=nan([size(ModelVal) XsecError*2/0.5+1]);
    
    for g=1:5 % AFT then ZHe
        if g==5 %AHe
            ii=7;
        elseif  g==4 %AFT
            ii=8; % which index value tchronometer I want
        elseif g==3 % ZHe
            ii=9;
        elseif g==2 %ZFT
            ii=10;
        elseif g==1 %MAr
            ii=11;
        end
        
        % for each sample in table
        for jj=1:d
            % if that sample exists in this thermochronometer
            if index(jj,g)==1
                %Matches(jj,g)=3;
                %%%%%% Variables
                %%%% Measured Tchron Value
                %%% TchronData.Age(jj)
                
                %%%% Measured Tchron Error
                %%% TchronData.Error(jj)
                
                %%%% Sample Location (Adjusted to Pecube Grid)
                %%% TchronData.xdist(jj)
                sampleloc(jj,g)=round(TchronData.xdist(jj)*2)/2;
                
                % Creating Variable for Cross Plot
                %  (jj,g,1) is Measured Age
                %  (jj,g,2) is Measured Error
                %  (jj,g,3:end) is Modelled ages
                CrossPlot(jj,g,1)=(TchronData.Age(jj));
                CrossPlot(jj,g,2)=(TchronData.Error(jj));
                
                % Loop over X-error in cross-section
                % Without looping over 1km wide window, N10const had 28.57%
                % fit. With 1km error, 42.86% fit.
                ll=[sampleloc(jj,g)-(XsecError): 0.5: sampleloc(jj,g)+(XsecError)];
                for lll= 1:length(ll)
                    samplelocll(jj,g)=ll(lll);
                    % Find Sample Location in MODELLED
                    IndexSampleLocMod=find(Age_Avg(:,4)==samplelocll(jj,g));
                    % Get Modelled Value at that location
                    ModelVal(jj,g)=Age_Avg(IndexSampleLocMod,ii);
                    CrossPlot(jj,g,2+lll)=Age_Avg(IndexSampleLocMod,ii);
                end
                %
            end
        end
    end
    
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%Cloud Distribution of Cross Plot %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    % !make sure max and min are included from modelled data
    % /change text/word colors to match with cloud/best fit line
    % !add legend
    % make loop to auto-go through all models
    % !add GOF val? (make sure coords with table from paper!)
    
    markcol=[0.6350 0.0780 0.1840
        0.4940, 0.1840, 0.5560
        0.8500, 0.3250, 0.0980
        0, 0.4470, 0.7410
        0.4660, 0.6740, 0.1880];
    
    % compute min/max/std for modelled data, create normal distribution based
    % on it
    %
    count=0;
    XX=NaN((length(find(index(:,3:4)==1))-1)*1002,2);
    % -1 in above statement because one of the values is >100 and I don't want
    % to include it!
    r1=NaN(1002,2);
    
    for g=3:4
        for jj=1:d
            if index(jj,g)==1 && CrossPlot(jj,g,1)<101
                S=std(CrossPlot(jj,g,3:end));
                MEAN=mean(CrossPlot(jj,g,3:end));
                MAX=max(CrossPlot(jj,g,3:end));
                MIN=min(CrossPlot(jj,g,3:end));
                mu1=[CrossPlot(jj,g,1), MEAN];
                sigma1=[2*CrossPlot(jj,g,2), 2*S];
                r1(1:1000,:)= mvnrnd(mu1,sigma1,1000);
                %r1(:,find(r1(:,2)>=MAX|r1(:,2)<=MIN))=[]; % remove values
                % outside of range... CURRENTLY BROKEN/doesn't work on whole
                % loop
                r1(1001:1002,:)=[CrossPlot(jj,g,1),MIN;CrossPlot(jj,g,1),MAX];
                count=count+1;
                XX(((count-1)*1002+1):((count)*1002),:)=r1;
            end
        end
    end
    %
    % clean data
    [row,col]=find(XX<=0 | XX(:,1)>=100 | XX(:,2)<=0 | XX(:,2)>=100);
    XX(unique(row),:)=[];
    % fit data (linearly)
    xxx=[0:100];
    [pfit,Sfit]=polyfit(XX(:,1),XX(:,2),1);
    %ypfit=polyval(pfit,xxx);
    ypfit=polyval(pfit,[0,xxx])
    X=xxx.';                         % design "matrix" as column
    b0=X./ypfit.'                       % solve for slope, zero intercept
    yhat0=b0*[0 xxx];                % evaluate over same range from origin
    plot([0 xxx],yhat0,'r-')
    pause
    % plot cloud data
    sc=scatter(XX(:,1),XX(:,2),10,'.y'); % Scatter plot with points of size 10
    hold on
    % plot fitted data
    plot(xxx,ypfit)
    % r2 of the fit
    Rsqfit=1-(Sfit.normr/norm(XX(:,2)-mean(XX(:,2))))^2
    
    % plot 1:1 line
    plot([0 100],[0 100],'-k')
    
    % plot data w/ errorbars on top of clod data
    for g=4:-1:3
        for jj=3:d
            if index(jj,g)==1
                S=std(CrossPlot(jj,g,3:end));
                MEAN=mean(CrossPlot(jj,g,3:end));
                MAX=max(CrossPlot(jj,g,3:end));
                MIN=min(CrossPlot(jj,g,3:end));
                dd=errorbar(CrossPlot(jj,g,1),MEAN,...
                    S,S,...
                    CrossPlot(jj,g,2),CrossPlot(jj,g,2),...
                    'ok','MarkerFaceColor',markcol(g,:),...
                    'MarkerEdgeColor',markcol(g,:));
                dd.Bar.LineStyle='dotted';
                text(CrossPlot(jj,g,1),MEAN, TchronData.sample_num(jj), 'VerticalAlignment','top', ...
                    'HorizontalAlignment','left')
            end
        end
    end
    
    axis([0 100 0 100])
    xlabel('Measured');ylabel('Modelled')
    
    
    %['Temperature is ',num2str(c),' C']
    title(['Model ' num2str(ModNum) ': velocity \it' num2str(MODID)]);% %s',ModNum,MODID))
    Fitlabel={sprintf('Goodness of fit=%2.2f%%',PERCFIT);...
        sprintf('R^2=%0.2f',Rsqfit);...
        sprintf('slope=%0.2f',pfit(1))};
    text(67,30,Fitlabel)
    box on
    hold off
    legend('Gaussian PDF','Best fit line','1:1 line','data','Location','southeast')
    pause 
    % Save Figure
    SubplotLoc=[1 2; 3 4; 5 6; 7 8; 9 10; 11 12; 13 14; 15 16; 17 18; 19 20];
    %subplot(10,2,SubplotLoc(n,ModCol))% instead of 10,length(StartFolder)-4
    saveas(gcf,strcat('/Users/victoriabuford/Box Sync/Research/Grants-writeups-presentations/MyPaper2018/FigsTables/Supplementary/S5-CrossPlots/', num2str(SubplotLoc(n,ModCol)), '_Gauss.eps'),'epsc')
    % remove gaussian cloud
    delete(sc)
    saveas(gcf,strcat('/Users/victoriabuford/Box Sync/Research/Grants-writeups-presentations/MyPaper2018/FigsTables/Supplementary/S5-CrossPlots/', num2str(SubplotLoc(n,ModCol)), '.eps'),'epsc')
    
end

%%%%% Then, add all of the plots to a multi-page PDF in illustrator (eg,
%%%%% place), and reduce to 67%

% If I want to fit the cloud data?
%gmPDF = @(x,y)reshape(pdf(gm,[x(:) y(:)]),size(x));
%fcontour(gmPDF)
%%%%%%%%
