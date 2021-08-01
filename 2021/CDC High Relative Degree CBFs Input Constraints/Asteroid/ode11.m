function [t_arr, y_arr, out_arr] = ode11(func, span, init, opts)
global is_error outdata enforce_errors
try
    dt0 = opts.dt;
catch
    dt0 = diff(span)/1000;
end
len = ceil(diff(span)/dt0)+1;
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
    RelTol = opts.RelTol;
    if isempty(RelTol), RelTol = 1e-4; end
catch
    RelTol = 1e-4;
end
try
    enforce_errors = opts.enforce_errors;
catch
    enforce_errors = 0;
end
try
    max_error_count = opts.max_error_count;
catch
    max_error_count = 20;
end
try
    debug = opts.debug;
catch
    debug = 0;
end

%% First Function Call
is_error = 0;
ydot = feval(func, t, y);
if is_error
    error('Initial conditions result in impossible constraints');
end
out_arr(1) = outdata;

%% Other Function Calls
while t < span(end)
    error_count = 0;
    is_error = 1; % ensure the while loop is entered
    dt = min(dt0, RelTol/norm(ydot));
    t0 = t;
    y0 = y;
    ydot0 = ydot;
    while is_error && error_count < max_error_count
        is_error = 0;
        t = t0 + dt;
        y = y0 + ydot0*dt;
        try
            ydot = feval(func, t, y);
        catch e
            warning('Exited Early. Error Occurred');
            disp(getReport(e, 'extended', 'hyperlinks', 'on'));
            t_arr = t_arr(1:i);
            y_arr = y_arr(1:i,:);
            out_arr = out_arr(1:i);
            return
        end
        error_count = error_count+1;
        dt = dt/4;
        if ~enforce_errors
            is_error = 0;
        end
    end
    if is_error
        disp(['Warning: Constraint violated at time t = ' num2str(t) '. Minimum step size used']);
            % Treat violations as errors
%             t_arr = t_arr(1:i);
%             y_arr = y_arr(1:i,:);
%             out_arr = out_arr(1:i);
%             return
    end
    i = i+1;
    t_arr(i) = t;
    y_arr(i,:) = y(:)';    
    out_arr(i) = outdata;
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