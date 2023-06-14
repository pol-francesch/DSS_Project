%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Initial Conditions - GPS satellites
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% GPS Almanac Week 204
prns = [ 1,  2,  3,  4,  5,  6,  7,  8,  9, 10, 11, 12, 13, 14, 15, 16, ...
         7, 18, 19, 20, 21, 22, 23, 24, 25, 26, 27, 28, 29, 30, 31 ];

gps_sma = 26599.800;
gps_ecc = 0.000001;
gps_inc = 55.0;
gps_raan = [-143.881352, -149.446614,  -84.778512,  -22.818990,  -87.385619, ...
            -144.357305,   35.349176,  153.734801,  -26.028092,  -84.931462, ...
            -141.883192,   99.225855,  -16.924009,   97.083800,  -32.683446, ...
             100.277002,  158.759866, -143.755546,  161.314101,  -94.320245, ...
            -149.686511,  -86.467123,   30.236113,   94.562051,   91.458499, ...
             155.031466,   33.155708,  159.555151,   35.770562,   36.473751, ...
             -25.294003];
gps_aop = [  53.765545,  -75.631986,   57.743969, -173.043165,   64.681857, ...
             -46.573319, -127.231250,   13.148360,  112.699728, -139.444742, ...
            -154.121597,   76.981738,   53.130934, -176.496584,   67.669709, ...
              44.494200,  -81.162593, -178.499508,  128.866153, -166.197417, ...
             -45.969586, -177.202756,   50.789967,   59.425199,   26.153212, ...
              40.480607,   96.978207,  136.728158, -151.319311,   29.645619, ...
            -128.456075];
gps_ta = [ -86.960177, -144.505126, -152.714424,   13.494215,   61.374757, ...
             -74.233353, -120.546970,   39.999955,   59.759746,  152.572997, ...
              -2.914681, -167.377052,   14.453008, -163.841686,  -15.388842, ...
              87.047746,   37.617981,  -80.040336,  160.861859,  -37.675960, ...
              36.690388, -139.207914,  -70.966229, -173.403954,  142.599792, ...
              40.680892,  147.604237,   43.801439, -129.846039, -157.112496, ...
              78.628936];

% Convert to meters and rads
gps_sma = gps_sma * 1000;
gps_inc = deg2rad(gps_inc);
gps_raan = deg2rad(gps_raan);
gps_aop = deg2rad(gps_aop);
gps_ta = deg2rad(gps_ta);

% Initialize GPS satellite osculating Keplerian elements
gps01oe = [gps_sma, gps_ecc, gps_inc, gps_raan( 1), gps_aop( 1), gps_ta( 1)];
gps02oe = [gps_sma, gps_ecc, gps_inc, gps_raan( 2), gps_aop( 2), gps_ta( 2)];
gps03oe = [gps_sma, gps_ecc, gps_inc, gps_raan( 3), gps_aop( 3), gps_ta( 3)];
gps04oe = [gps_sma, gps_ecc, gps_inc, gps_raan( 4), gps_aop( 4), gps_ta( 4)];
gps05oe = [gps_sma, gps_ecc, gps_inc, gps_raan( 5), gps_aop( 5), gps_ta( 5)];
gps06oe = [gps_sma, gps_ecc, gps_inc, gps_raan( 6), gps_aop( 6), gps_ta( 6)];
gps07oe = [gps_sma, gps_ecc, gps_inc, gps_raan( 7), gps_aop( 7), gps_ta( 7)];
gps08oe = [gps_sma, gps_ecc, gps_inc, gps_raan( 8), gps_aop( 8), gps_ta( 8)];
gps09oe = [gps_sma, gps_ecc, gps_inc, gps_raan( 9), gps_aop( 9), gps_ta( 9)];
gps10oe = [gps_sma, gps_ecc, gps_inc, gps_raan(10), gps_aop(10), gps_ta(10)];
gps11oe = [gps_sma, gps_ecc, gps_inc, gps_raan(11), gps_aop(11), gps_ta(11)];
gps12oe = [gps_sma, gps_ecc, gps_inc, gps_raan(12), gps_aop(12), gps_ta(12)];
gps13oe = [gps_sma, gps_ecc, gps_inc, gps_raan(13), gps_aop(13), gps_ta(13)];
gps14oe = [gps_sma, gps_ecc, gps_inc, gps_raan(14), gps_aop(14), gps_ta(14)];
gps15oe = [gps_sma, gps_ecc, gps_inc, gps_raan(15), gps_aop(15), gps_ta(15)];
gps16oe = [gps_sma, gps_ecc, gps_inc, gps_raan(16), gps_aop(16), gps_ta(16)];
gps17oe = [gps_sma, gps_ecc, gps_inc, gps_raan(17), gps_aop(17), gps_ta(17)];
gps18oe = [gps_sma, gps_ecc, gps_inc, gps_raan(18), gps_aop(18), gps_ta(18)];
gps19oe = [gps_sma, gps_ecc, gps_inc, gps_raan(19), gps_aop(19), gps_ta(19)];
gps20oe = [gps_sma, gps_ecc, gps_inc, gps_raan(20), gps_aop(20), gps_ta(20)];
gps21oe = [gps_sma, gps_ecc, gps_inc, gps_raan(21), gps_aop(21), gps_ta(21)];
gps22oe = [gps_sma, gps_ecc, gps_inc, gps_raan(22), gps_aop(22), gps_ta(22)];
gps23oe = [gps_sma, gps_ecc, gps_inc, gps_raan(23), gps_aop(23), gps_ta(23)];
gps24oe = [gps_sma, gps_ecc, gps_inc, gps_raan(24), gps_aop(24), gps_ta(24)];
gps25oe = [gps_sma, gps_ecc, gps_inc, gps_raan(25), gps_aop(25), gps_ta(25)];
gps26oe = [gps_sma, gps_ecc, gps_inc, gps_raan(26), gps_aop(26), gps_ta(26)];
gps27oe = [gps_sma, gps_ecc, gps_inc, gps_raan(27), gps_aop(27), gps_ta(27)];
gps28oe = [gps_sma, gps_ecc, gps_inc, gps_raan(28), gps_aop(28), gps_ta(28)];
gps29oe = [gps_sma, gps_ecc, gps_inc, gps_raan(29), gps_aop(29), gps_ta(29)];
gps30oe = [gps_sma, gps_ecc, gps_inc, gps_raan(30), gps_aop(30), gps_ta(30)];
gps31oe = [gps_sma, gps_ecc, gps_inc, gps_raan(31), gps_aop(31), gps_ta(31)];

gps_constellation = [ gps01oe; gps02oe; gps03oe; gps04oe; gps05oe; gps06oe; gps07oe; gps08oe; ...
                      gps09oe; gps10oe; gps11oe; gps12oe; gps13oe; gps14oe; gps15oe; gps16oe; ...
                      gps17oe; gps18oe; gps19oe; gps21oe; gps22oe; gps23oe; gps24oe; gps25oe; ...
                      gps25oe; gps26oe; gps27oe; gps28oe; gps29oe; gps30oe; gps31oe ];
