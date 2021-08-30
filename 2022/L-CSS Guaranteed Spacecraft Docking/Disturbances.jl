

module Disturbances

using LinearAlgebra

wx_data = nothing;
wu_data = nothing;

w_xmax = nothing;
w_umax = nothing;

d_u = 3;
d_x = 3;

len = 100001;
t_inc = 10;

function write_data()
    f = open("WX.dat", "w");
    write(f, rand(len*d_x).*sign.(randn(len*d_x)));
    close(f);

    f = open("WU.dat", "w");
    write(f, rand(len*d_u).*sign.(randn(len*d_u)));
    close(f);
end

function read_data()
    f = open("WX.dat", "r");
    global wx_data = zeros(len, d_x);
    for i=1:len
        for j=1:d_x
            global wx_data[i,j] = read(f, Float64);
        end
        if norm(wx_data[i,:]) > 1
            wx_data[i,:] = wx_data[i,:] / norm(wx_data[i,:]);
        end
    end
    close(f);
    
    f = open("WU.dat", "r");
    global wu_data = zeros(len, d_u);
    for i=1:len
        for j=1:d_u
            global wu_data[i,j] = read(f, Float64);
        end
        if norm(wu_data[i,:]) > 1
            wu_data[i,:] = wu_data[i,:] / norm(wu_data[i,:]);
        end
    end
    close(f);
end

function Matched(t)
    index1 = Int(floor(t/t_inc));
    index2 = Int(ceil(t/t_inc));
    if index2 > len
        println("WARNING: Requested noise outside time domain.");
        index2 = len;
    end
    delta = (t - index1*t_inc)/t_inc;
    return w_umax*(wu_data[index1+1,:]*(1-delta) + wu_data[index2+1,:]*delta);
end

function Unmatched(t)
    index1 = Int(floor(t/t_inc));
    index2 = Int(ceil(t/t_inc));
    if index2 > len
        println("WARNING: Requested noise outside time domain.");
        index2 = len;
    end
    delta = (t - index1*t_inc)/t_inc;
    return w_xmax*(wx_data[index1+1,:]*(1-delta) + wx_data[index2+1,:]*delta);
end

function set_w_umax(x); global w_umax = x; end
function set_w_xmax(x); global w_xmax = x; end
function set_t_inc(x); global t_inc = x; end
function set_d_u(x); global d_u = x; end
function set_d_x(x); global d_x = x; end

end

nothing