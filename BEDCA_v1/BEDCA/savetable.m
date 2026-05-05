function savetable(tab,fname)
  % savetable(tab,fname)
%
  % sabes a csv file from tab structure (Variable stats)
  % NJB Oct 2020


  % Summary stats Rows mu, taux, taux, sdx, sdz

  summarystats=[tab.means; tab.medians; tab.sd];
  
writetable(array2table(summarystats','VariableNames',{'Mean','Median','StandardDeviation'},'RowNames',tab.params),fname,'WriteRowNames',true);
