

module Disturbances

wx_data = nothing;
wu_data = nothing;

w_xmax = nothing;
w_umax = nothing;

len = 100001;
t_inc = 10;

function write_data()
    f = open("WX.dat", "w");
    write(f, rand(len*6).*sign.(randn(len*6)));
    close(f);

    f = open("WU.dat", "w");
    write(f, rand(len*3).*sign.(randn(len*3)));
    close(f);
end

function read_data()
    f = open("WX.dat", "r");
    global wx_data = zeros(len, 6);
    for i=1:len
        for j=1:6
            global wx_data[i,j] = read(f, Float64);
        end
    end
    close(f);
    
    f = open("WU.dat", "r");
    global wu_data = zeros(len, 3);
    for i=1:len
        for j=1:3
            global wu_data[i,j] = read(f, Float64);
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

end

nothing