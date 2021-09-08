function [t_arr, y_arr, out_arr] = ode01(span, init, opts)
global outdata
dt = opts.dt;
len = ceil(diff(span)/dt)+1;
t_arr = zeros(len, 1);
y_arr = zeros(len, length(init));
out_arr(len) = SavedData;
y = init(:);
t = span(1);
t_arr(1) = t;
y_arr(1,:) = y';
i = 1;
outdata = SavedData;

%% Parse Options
try
    opts.OutputFcn = opts.OutputFcn;
    output = 1;
catch
    disp 'No Output Function';
    output = 0;
end
if output 
    feval(opts.OutputFcn, span, y', 'init');
end
try
    debug = opts.debug;
catch
    debug = 0;
end

%% First Function Call
CalculateU(t, y);
out_arr(1) = outdata;

%% Waitbar
SimDur = span(end);
h = waitbar(0, 'Simulating Scenario');

%% Other Function Calls
while t < SimDur
    u = CalculateU(t, y);
    y = UpdateX(t, y, u, dt);
    t = t+dt;
    
    % Assign outputs to arrays once loop ends with no errors
    i = i+1;
    t_arr(i) = t;
    y_arr(i,:) = y(:)';    
    out_arr(i) = outdata;
    
    % If waitbar was closed, end sim early
    try
        waitbar(t/SimDur, h);
    catch
        warning('Exited early due to waitbar being closed.');
        t_arr = t_arr(1:i);
        y_arr = y_arr(1:i,:);
        out_arr = out_arr(1:i);
        return
    end
    if output
        if debug
            feval(opts.OutputFcn, t_arr(1:i), y_arr(1:i,:), []);
        else
            try
                feval(opts.OutputFcn, t_arr(1:i), y_arr(1:i,:), []);
            catch
                disp 'Output Function Failed';
            end
        end
    end
    if i==length(t_arr) & t < SimDur % Preallocate more steps if necessary
        t_arr = [t_arr; zeros(100, 1)];
        y_arr = [y_arr; zeros(100, length(init))];
        out_arr(end+100) = SavedData;
    end
end
t_arr = t_arr(1:i);
y_arr = y_arr(1:i,:);
out_arr = out_arr(1:i);
end