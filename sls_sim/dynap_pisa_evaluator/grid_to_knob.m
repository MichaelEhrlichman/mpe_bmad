clear

grid = importdata('knob_grid',' ',1);

ql1fit  = fit([grid.data(:,1),grid.data(:,2)],grid.data(:,3),'poly33');
ql2fit  = fit([grid.data(:,1),grid.data(:,2)],grid.data(:,4),'poly33');
ql3fit  = fit([grid.data(:,1),grid.data(:,2)],grid.data(:,5),'poly33');
ql4fit  = fit([grid.data(:,1),grid.data(:,2)],grid.data(:,6),'poly33');
ql5hfit = fit([grid.data(:,1),grid.data(:,2)],grid.data(:,7),'poly33');

%plot(ql1fit,[grid.data(:,1),grid.data(:,2)],grid.data(:,3))

quad1=grid.colheaders{3};
quad2=grid.colheaders{4};
quad3=grid.colheaders{5};
quad4=grid.colheaders{6};
quad5=grid.colheaders{7};

scaling = 1e-4/3.0  %

ql1str = MatlabFit2BmadKnob(ql1fit,scaling);
ql2str = MatlabFit2BmadKnob(ql2fit,scaling);
ql3str = MatlabFit2BmadKnob(ql3fit,scaling);
ql4str = MatlabFit2BmadKnob(ql4fit,scaling);
ql5hstr = MatlabFit2BmadKnob(ql5hfit,scaling);

ql1knob = [quad1 '[k1]:' ql1str];
ql2knob = [quad2 '[k1]:' ql2str];
ql3knob = [quad3 '[k1]:' ql3str];
ql4knob = [quad4 '[k1]:' ql4str];
ql5hknob = [quad5 '[k1]:' ql5hstr];

knob = sprintf('grl1: group = {%s, &\n%s, &\n%s, &\n%s, &\n%s}, var = {a,b}',ql1knob,ql2knob,ql3knob,ql4knob,ql5hknob);

knob