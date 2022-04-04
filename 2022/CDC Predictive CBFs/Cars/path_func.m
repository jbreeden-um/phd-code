function [p, dtau, dt, dx] = path_func(tau,t,x)
k = 1;
vdes = 12;

p = [z_future(tau-t,x(1),x(2));
    zdot_future(tau-t,x(2));
    z_future(tau-t,x(3),x(4));
    zdot_future(tau-t,x(4))];

if nargout > 1
    dtau = [z_future_dtau(tau-t,x(2));
            zdot_future_dtau(tau-t,x(2));
            z_future_dtau(tau-t,x(4));
            zdot_future_dtau(tau-t,x(4))];

    dt = -dtau;

    dx = [z_future_dx(tau-t),    0, 0;
          zdot_future_dx(tau-t), 0, 0;
          0, 0,                  z_future_dx(tau-t);
          0, 0,                  zdot_future_dx(tau-t)];
end

    function out = z_future(delta_t,z,zdot)
    out = z + vdes*delta_t + (zdot-vdes)/k*(1 - exp(-k*delta_t));
    end

    function out = zdot_future(delta_t,zdot)
    out = vdes + (zdot-vdes)*exp(-k*delta_t);
    end

    function out = z_future_dtau(delta_t,zdot)
    out = vdes + (zdot-vdes)*exp(-k*delta_t);
    end

    function out = zdot_future_dtau(delta_t,zdot)
    out = -k*zdot*exp(-k*delta_t);
    end

    function out = z_future_dx(delta_t)
    out = [1, (1 - exp(-k*delta_t))/k];
    end

    function out = zdot_future_dx(delta_t)
    out = [0, exp(-k*delta_t)];
    end

end


