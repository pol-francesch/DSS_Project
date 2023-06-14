function M = meanAnomAtTime(t,t0,M0,n)
% Returns the mean anomaly given a current time, initial time, mean anomaly
% at initial time, and mean motion
% Inputs: t - current time (sec)
%         M0 - initial mean anomaly (rad)
% Outputs: M - mean anomaly at time t

    M = n*(t - t0) + M0;

end