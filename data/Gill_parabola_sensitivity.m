%% This code performs sensitivity test for SP case (w=50 km, d=500 m, h=600 m)
% 1) D_\infty = 1000 m
% 2) different choices of D results in different dimensionless geometry
% parameters as well as different q
% 3) but we want to demonstrate that the critical Q1 shouldn't matter, at
% least not greatly mattered


% I then realize it indeed does not matter, the dimensional Gill's model
% relates h_c, \alpha, D_\infty, Q, and Q1. All except for Q1 can be
% determined for the Samaon Passage case and Q1 can be later computed. The
% dimensionless potential vorticity q=-D/D_\infty only have one dependence,
% which is D_infty. Changing D will cause the dimensionless h_c, r_c, and Q
% to change and the resulted dimensionless Q1 under critical condition will 
% change as well. But after scaling, the dimensional Q1 should be the same.