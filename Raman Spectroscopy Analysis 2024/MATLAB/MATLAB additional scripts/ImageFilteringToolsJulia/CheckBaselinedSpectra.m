%[CorrMeanSpectra,Baseline]= removeBaseline(PCARun.xj,mean(PCARun.dataforPCA(PCARun.goodspectra==1,:)),PCARun.PolynomialOrder,0,50,1,'none');
clear CorrMeanSpectra Baseline
[CorrMeanSpectra,Baseline]=f_baseline_corr(PCARun.xj,...
        mean(PCARun.dataforPCA(PCARun.goodspectra==1,:)),5,600,5,'plot2D');
figure;
plot(PCARun.xj,CorrMeanSpectra);
hold on;
plot(PCARun.xj,mean(PCARun.dataforPCA(PCARun.goodspectra==1,:)));
plot(PCARun.xj,Baseline);
legend('Corrected Spectra','Uncorrected Spectra','Fitted Baseline')

 ButtonName = questdlg('Does the corrected spectra look ok?', ...
                         'Look good Question', ...
                         'Yes', 'No', 'Yes');
   switch ButtonName,
     case 'Yes',
      disp('Looks good');
      looksgood=1;
     case 'No',
      disp('Does not look good. Repeat')
      close(gcf)
      looksgood=0;
      InputforPCA
      CheckBaselinedSpectra
   end % switch
