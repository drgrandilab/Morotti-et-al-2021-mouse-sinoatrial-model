%% Calculate timecourse for currents and other intermediates

function currents = calcCurrents(t,y,p,output)
% After running a simulation, feed the time vector and state variables into
% this function to compute ionic currents, etc.
% currents: [I_Na I_Catot];

currents=[];
for i = 1:size(t)
    if ceil(i/1000) == i/1000
        disp(['t = ',num2str(ceil(t(i)))]);
    end
    currents = [currents;mouse_SAM_eccODEfile(t(i),y(i,:),p,output)];
end
% end calcCurrents