clear

grid = importdata('knob_grid',' ',1);

q1fit  = fit(grid.data(:,1),grid.data(:,3),'poly2')
q2fit  = fit(grid.data(:,1),grid.data(:,4),'poly2')
q3fit  = fit(grid.data(:,1),grid.data(:,5),'poly2')

%plot(ql1fit,[grid.data(:,1),grid.data(:,2)],grid.data(:,3))

quad1=grid.colheaders{3};
quad2=grid.colheaders{4};
quad3=grid.colheaders{5};

scaling = 5e-7/6.0  % Scaling divided by number of occurances cell.

q1str = MatlabFit2BmadKnob(q1fit,scaling,2);
q2str = MatlabFit2BmadKnob(q2fit,scaling,2);
q3str = MatlabFit2BmadKnob(q3fit,scaling,2);

q1knob = [quad1 '[k1]:' q1str];
q2knob = [quad2 '[k1]:' q2str];
q3knob = [quad3 '[k1]:' q3str];

knob = sprintf('grs1: group = {%s, &\n%s, &\n%s}, var = {a}',q1knob,q2knob,q3knob);

knob