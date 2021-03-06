
function bmadKnob = MatlabFit2BmadKnob_1D(qfit,sc,order)
    if order==3
        bmadKnob=sprintf('a*%0.3f*%0.1e+a^2*%0.3f*%0.1e+a^3*%0.3f*%0.1e', ...
                  qfit.p3,sc, ...
                  qfit.p2,sc^2, ...
                  qfit.p1,sc^3);
    elseif order==2              
        bmadKnob=sprintf('a*%0.3f*%0.1e+a^2*%0.3f*%0.1e+a^3*%0.3f*%0.1e', ...
                  qfit.p3,sc, ...
                  qfit.p2,sc^2, ...
                  qfit.p1,sc^3);
    end
end