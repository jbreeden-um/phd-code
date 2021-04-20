function stop = PlotOutputs(t, x, extra)
persistent Constraints ... % the constraints loaded from file once at beginning of sim
    n ... % the number of constraints
    ConPlots ... % the handles to the plots of the constraint h values
    ConPhysPlots ... % the handles to the plots of the constraint values in physical units
    ComparePlots ... % the handles to the plots of what the constraint values should be
    ConData ... % the actual constraint h values
    ConPhysData ... % the constraint physical interpretations
    CompareData ... % the constraint requirements (0 if active, NaN if inactive)
    tData rData ... % time and position data
    posPlot ... % the handle to the 3D position plot
    uData ... % control input data, generated here instead of from ode45
    uPlots ... % the handles to the control input plots (two plots)
    TargetPlots ... % the handles to the plots locating the targets in the position plot
    nTarget ... % the number of targets total
    TargetPlotNums ... % the constraint indices corresponding to the TargetPlots handles
    posPlotIndividual ... % the position plots individually
    viewVectorPlot ... % quiver plot of where the spacecraft is pointing
    vecData % time history of pointing vector

if isequal(extra, 'init')
    file = load('InData/SimData.mat');
    Constraints = file.Constraints;
    n = length(Constraints);
    
    ConPlots = [];
    ConData = [];
    ConPhysPlots = [];
    ConPhysData = [];
    ComparePlots = [];
    CompareData = [];
    TargetPlots = [];
    rData = [];
    uData = [];
    vecData = [];
    TargetPlotNums = [];
    posPlotIndividual = [];
    viewVectorPlot = [];
    
    tData = t(1);
    figure(1); clf;
    rData(1,:) = x(1:3)';
    posPlot = plot3(rData(:,1), rData(:,2), rData(:,3));
    xlabel 'x'; ylabel 'y'; zlabel 'z'; hold on;
    axis equal;
    
    nTarget = 0;
    for index=1:n
        C = Constraints(index);
        if isequal(C.mode, 'proximity')
            nTarget = nTarget + 1;
            h = plot3(C.vector(1), C.vector(2), C.vector(3), 'ko');
            if nTarget==1
                TargetPlots = h;
            else
                TargetPlots(nTarget) = h;
            end
            TargetPlotNums(nTarget) = index;
        elseif isequal(C.mode, 'avoidance')
            [s1, s2, s3] = sphere(20);
            surf(s1*C.scalar, s2*C.scalar, s3*C.scalar, 'FaceColor', [.7; .7; .7], ...
                'FaceAlpha', 0.4, 'EdgeAlpha', 0.2, 'EdgeColor', [.8; .8; .8]);
        elseif isequal(C.mode, 'asteroid')
            if isequal(C.vector, 'Eros')
                shape = load('InData/Eros_Shape.mat');
                vertices = shape.vertices;
                plates = shape.plates + 1;
                trisurf(plates, vertices(:,1), vertices(:,2), vertices(:,3), ...
                    'FaceColor', [.7; .7; .7], 'FaceAlpha', 0.6, 'EdgeAlpha', 0.8);
            end
        end
    end
    
    mrp = x(7:9);
    q = PtoQ(mrp);
    e3 = [0; 0; 1];
    uhat = RotQ(e3, q);
    vecData(1,:) = uhat';
    viewVectorPlot = quiver3(rData(end,1), rData(end,2), rData(end,3), uhat(1), uhat(2), uhat(3), ...
        'r', 'MaxHeadSize', .5, 'LineWidth', 2);
    
    figure(5); clf;
    subplot(2,3,1);
    posPlotIndividual = plot(tData, rData(:,1));
    xlabel 'Time (s)'; ylabel 'Position (km)'; title 'X State';
    subplot(2,3,2);
    posPlotIndividual(2) = plot(tData, rData(:,2));
    xlabel 'Time (s)'; ylabel 'Position (km)'; title 'Y State';
    subplot(2,3,3);
    posPlotIndividual(3) = plot(tData, rData(:,3));
    xlabel 'Time (s)'; ylabel 'Position (km)'; title 'Z State';
    subplot(2,3,4);
    posPlotIndividual(4) = plot(tData, rData(:,1), 'r'); hold on;
    posPlotIndividual(5) = plot(tData, rData(:,2), 'g');
    posPlotIndividual(6) = plot(tData, rData(:,3), 'b');
    xlabel 'Time (s)'; ylabel 'Position (km)'; title 'State'; legend x y z;
    subplot(2,3,5);
    posPlotIndividual(7) = plot(tData, vecData(:,1), 'r'); hold on;
    posPlotIndividual(8) = plot(tData, vecData(:,2), 'g');
    posPlotIndividual(9) = plot(tData, vecData(:,3), 'b');
    xlabel 'Time (s)'; ylabel 'Component'; title 'View State'; legend x y z;

    w = ceil(sqrt(n)); % number of columns
    h = ceil(n/w); % number of rows
    
%     w = 5; h= 2;
    figure(2); clf;
    for index=n:-1:1
        [ConData(1,index), ConPhysData(1, index), CompareData(1, index)] = ConstraintValue(Constraints(index), t(1), x);
        subplot(h, w, index);
        ConPlots{index} = plot(tData, ConData(:, index), 'b'); hold on;
        ComparePlots{index} = plot(tData, CompareData(:, index), 'r');
        xlabel 'Time'; ylabel 'h(x, t)'; 
        title(['Constraint ' num2str(index)]);
        if isequal(Constraints(index).mode, 'asteroid')
            if ConPhysData(1,index) > 0
                error('Constraint is violated');
            end
        end
    end
    figure(4);
    for index=n:-1:1
        subplot(h, w, index);
        ConPhysPlots{index} = plot(tData, ConPhysData(:, index), 'b');
        xlabel 'Time';
        if isequal(Constraints(index).mode, 'proximity')
            ylabel 'Distance (km)';
            title 'Proximity Constraint';
        elseif isequal(Constraints(index).mode, 'boresight')
            ylabel 'Angle (deg)';
            title 'Boresight Constraint';
        elseif isequal(Constraints(index).mode, 'avoidance')
            ylabel 'Distance (km)';
            title 'Avoidance Constraint';
        elseif isequal(Constraints(index).mode, 'asteroid')
            ylabel 'Distance (km)';
            title 'Asteroid Constraint';
        end
    end
    
    figure(3); clf;
    uData(1,:) = CalculateU(t(1), x, Constraints);
    subplot(1, 3, 1);
    uPlots(1:3) = plot(tData, uData(1,1:3));
    xlabel 'Time'; ylabel 'u'; title 'Thrust';
    subplot(1, 3, 2);
    uPlots(4:6) = plot(tData, uData(1,4:6));
    xlabel 'Time'; ylabel 'u'; title 'Torque';
    subplot(1, 3, 3);
    uPlots(7:14) = plot(tData, uData(1,7:14));
    xlabel 'Time'; ylabel 'u'; title 'Slack';
    
    % Add a rule for plotting the control input if possible
elseif isequal(extra, 'done')
    % Figure out if you need any cleanup here
else
    for tstep=1:length(t)
        tn = t(tstep);
        xn = x(:,tstep);
        s = length(tData)+1;
        tData(s) = tn;
        rData(s,:) = xn(1:3)';

        for index=n:-1:1
            [ConData(s,index), ConPhysData(s,index), CompareData(s,index)] = ConstraintValue(Constraints(index),tn,xn);
        end
        
        [uData(s,:), flag] = CalculateU(tn, xn, Constraints);
        if flag < 1
            disp(['Input Error at time ' num2str(tn) '. Code was ' num2str(flag)]);
        end
        
        mrp = xn(7:9);
        q = PtoQ(mrp);
        e3 = [0; 0; 1];
        uhat = RotQ(e3, q);
        vecData(s,:) = uhat';
    end
    
    for index=n:-1:1
        set(ConPlots{index}, 'XData', tData, 'YData', ConData(:, index));
        set(ComparePlots{index}, 'XData', tData, 'YData', CompareData(:, index));
        set(ConPhysPlots{index}, 'XData', tData, 'YData', ConPhysData(:, index));
    end
    
    for index=1:nTarget
        if isnan(CompareData(s, TargetPlotNums(index)))
             set(TargetPlots(index), 'MarkerFaceColor', [.2; .2; .2]);
        else
             set(TargetPlots(index), 'MarkerFaceColor', [1; 0; 0]);
        end
    end
    
    set(posPlot, 'XData', rData(:,1), 'YData', rData(:,2), 'ZData', rData(:,3));
    for index=n:-1:1
        set(ConPlots{index}, 'XData', tData, 'YData', ConData(:, index));
        set(ComparePlots{index}, 'XData', tData, 'YData', CompareData(:, index));
    end
    for index=1:3
        set(posPlotIndividual(index), 'XData', tData, 'YData', rData(:,index));
    end
    for index=1:3
        set(posPlotIndividual(index+3), 'XData', tData, 'YData', rData(:,index));
    end
    for index=1:3
        set(posPlotIndividual(index+6), 'XData', tData, 'YData', vecData(:,index));
    end
    for index = 1:14
        set(uPlots(index), 'XData', tData, 'YData', uData(:,index));
    end
    
    set(viewVectorPlot, 'XData', xn(1), 'YData', xn(2), 'ZData', xn(3), ...
        'UData', uhat(1), 'VData', uhat(2), 'WData', uhat(3));
    
    drawnow;
end

stop = 0;
end

function [h, h_phys, req] = ConstraintValue(con, t, x)
persistent Eros
rvec = x(1:3);
pvec = con.vector;
if isequal(con.mode, 'proximity')
    rho = con.scalar;
    h = (rvec-pvec)'*(rvec-pvec) - rho^2;
    h_phys = norm(rvec-pvec) - rho;
elseif isequal(con.mode, 'boresight')
    phi = con.scalar;
    yvec = pvec - rvec;
    yhat = yvec/norm(yvec);
    mrp = x(7:9);
    q = PtoQ(mrp);
    e3 = [0; 0; 1];
    uhat = RotQ(e3, q);
    h = cos(phi) - uhat'*yhat;
    h_phys = acosd(uhat'*yhat) - rad2deg(phi);
    if imag(h_phys) ~= 0
        disp 'Imaginary number in plots. Simulation crash eminent.';
    end
elseif isequal(con.mode, 'avoidance')
    rhobar = con.scalar;
    h = rhobar^2 - (rvec-pvec)'*(rvec-pvec);
    h_phys = rhobar - norm(rvec-pvec);
elseif isequal(con.mode, 'asteroid')
    if isequal(pvec, 'Eros')
        if isempty(Eros)
            file = load('InData/Eros_Shape.mat');
            Eros = (file.vertices)'; % Note that Eros's shape model file is in units of km
        end
        vecs = rvec - Eros;
        dists = vecnorm(vecs);
        [~,index] = min(dists);
        pvec = Eros(:,index);
    else
        error('Asteroid not defined');
    end
    rhobar = con.scalar;
    h = rhobar^2 - (rvec-pvec)'*(rvec-pvec);
    h_phys = rhobar - norm(rvec-pvec);
end

% Sets the required value at the current time instance
if t >= con.t1 && t <= con.t2
    req = 0; % required value is zero
else
    req = NaN; % no requirements
end
end