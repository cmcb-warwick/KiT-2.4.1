function algparams=preRunchecks(exptnum,runtreatment,specnam,filterpars,options)
  % algparams=preRunchecks(exptnum,runtreatment,specnam,filterpars,options)
  % This ouputs all info on data and set up of a run, but doesnt start the run. 
  %
  % NJB March 2020

  if ~isempty(options)
  
  options.checkonly=1;

[algparams convergencediagnos]=Run_BEDCA(exptnum,runtreatment,specnam,filterpars,options);

  else

    options.checkonly=1;

[algparams convergencediagnos]=Run_BEDCA(exptnum,runtreatment,specnam,filterpars,options);  

    end
  
    disp(' ');
  disp('*************** RUN IS NOT STARTED. Use Run_BEDCA to run *****************');
    disp(' ');
