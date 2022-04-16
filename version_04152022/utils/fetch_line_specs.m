function [linespecs] = fetch_line_specs

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Now, do some analysis and make some plots. Fun!
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

linespecs.lh.style = '--';
linespecs.rh.style = '-.';
linespecs.all.style = '-';
linespecs.lh.width = 2.5;
linespecs.rh.width = 2.5;
linespecs.all.width = 2.5;


%tissues
linespecs.wm.color = 'r';
linespecs.wm.style = ':';

linespecs.ctx.color = 'b';
linespecs.ctx.style = ':';

linespecs.csf.color = 'g';
linespecs.csf.style = ':';

%bundles
linespecs.fminor.color = 'm';
linespecs.fminor.style = '-';

linespecs.fmajor.color = [0.4940, 0.1840, 0.5560];
linespecs.fmajor.style = '-';

linespecs.cst.color = 'c';
linespecs.cst.style = '-';

linespecs.slfp.color = [0 0 1];
linespecs.slfp.style = '-';

linespecs.slft.color = [0.3010, 0.7450, 0.9330];
linespecs.slft.style = '-';

linespecs.ilf.color = [0 0.4 0.8];
linespecs.ilf.style = '-';

linespecs.unc.color = [0.9290, 0.6940, 0.1250];
linespecs.unc.style = '-';

linespecs.cab.color = [0.8500, 0.3250, 0.0980];
linespecs.cab.style = '-';

linespecs.ccg.color = [0.6350, 0.0780, 0.1840];
linespecs.ccg.style = '-';


%other fs labels
linespecs.cerebellum.color = [0.75, 0.75, 0];
linespecs.cerebellum.style = '-.';

linespecs.brainstem.color = [0.25, 0.25, 0.25];
linespecs.brainstem.style = '-.';

linespecs.cc.color = [1 0 0];
linespecs.cc.style = '-.';


%subcortical
linespecs.putamen.color = [255;193;193]/255;
linespecs.putamen.style = '--';

linespecs.thalamus.color = [184;184;184]/255;
linespecs.thalamus.style = '--';

linespecs.caudate.color = 'b';
linespecs.caudate.style = '--';

linespecs.pallidum.color = 'm';
linespecs.pallidum.style = '--';

linespecs.amygdala.color = [0.9290, 0.6940, 0.1250];
linespecs.amygdala.style = '--';

linespecs.hippocampus.color = [0.75, 0, 0.75];
linespecs.hippocampus.style = '--';

linespecs.accumbens.color = [0, 0.75, 0.75];
linespecs.accumbens.style = '--';

linespecs.ventraldc.color = [0.25, 0.25, 0.25];
linespecs.ventraldc.style = '--';


end

