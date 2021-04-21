function [t_arr, y_arr, out_arr] = ode01(span, init, opts)
global is_error outdata
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
is_error = 0;
u = CalculateU(t, y, dt);
if is_error
    error('Initial conditions result in impossible constraints');
end
out_arr(1) = outdata;
% From this point forward, is_error is just a debugging tool and does not affect the
% operation of the simulation.

%% Waitbar
SimDur = span(end);
h = waitbar(0, 'Simulating Scenario');

%% Other Function Calls
while t < SimDur
    is_error = 0; % ensure the while loop is entered
    
    u = CalculateU(t, y, dt);
    y = UpdateX(t, y, u, dt);
    t = t+dt;
    
    if is_error
        disp(['Warning: Something is wrong at time t = ' num2str(t) '. Error code ' num2str(is_error)]);
    end
    
    % Assign outputs to arrays once loop ends with no errors
    i = i+1;
    t_arr(i) = t;
    y_arr(i,:) = y(:)';    
    out_arr(i) = outdata;
    
    % If waitbar was closed, end sim early
    try
        waitbar(t/SimDur, h);
    catch e
        warning('Exited Early. Error Occurred');
        disp(getReport(e, 'extended', 'hyperlinks', 'on'));
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
    if i==length(t_arr) % Preallocate more steps if necessary
        t_arr = [t_arr; zeros(100, 1)];
        y_arr = [y_arr; zeros(100, length(init))];
        out_arr(end+100) = SavedData;
    end
end
t_arr = t_arr(1:i);
y_arr = y_arr(1:i,:);
out_arr = out_arr(1:i);
end