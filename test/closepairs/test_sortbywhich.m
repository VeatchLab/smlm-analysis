function tests = test_sortbywhich
tests = functiontests(localfunctions);
end

function test_xyt(testcase)
xx = 0:1:1000;
yy = zeros(1,1000);
tt = 0:1:1000;

% correct answer is the one with smallest (range of data)/(range
% considered), where range considered is rmax for the spatial coordinates
% and (taumax - taumin) for time. 
rmax = 10;
taumin = 0;
taumax = 1;
% correct answer is 3, for time
verifyEqual(testcase,sortbywhich(xx,yy,rmax,tt,taumin,taumax),3);

rmax = 1;
taumin = 0;
taumax = 10;
% correct answer is 1, for x (first spatial coordinate)
verifyEqual(testcase,sortbywhich(xx,yy,rmax,tt,taumin,taumax),1);
% correct answer is 2, for y (second spatial coordinate)
verifyEqual(testcase,sortbywhich(yy,xx,rmax,tt,taumin,taumax),2);
end

function test_xy(testcase)
xx = 0:1:1000;
yy = zeros(1,1000);

rmax = 1;
% correct answer is the one that is more spread out (i.e. xx)
verifyEqual(testcase,sortbywhich(xx,yy,rmax),1);
verifyEqual(testcase,sortbywhich(yy,xx,rmax),2);
end