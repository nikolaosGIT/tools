function [Q, DVA] = DVA_method(t_data, U_data, I_data, I_smooth, U_smooth, dU_smooth, DVA_smooth)
% Function for differential voltage analysis (DVA) of measurement data independent of the current/voltage resolution
% 14.04.2021 - Nikolaos Wassiliadis and Philipp Rosner
%
% Details:
% This function calculates a differential voltage over capacity of a CC charge or discharge
% sequence, which allows for differential voltage analysis (DVA) similiar to the methods introduced
% by Bloom et al. (2005): 10.1016/j.jpowsour.2004.07.021. Measurements (t, I, U) are processed and
% can be smoothed (cubic splines) to account for low resolution measurements, e.g., at module and
% pack level.
%
% Input:
% t_data = Time in min
% U_data = Voltage in mV
% I_data = Current in A
% I_smooth = Spline smmothing factor for the current signal. Range: [0 1]. Increase for less smoothing, decrease for more.
% U_smooth = Spline smmothing factor for the voltage signal. Range: [0 1]. Increase for less smoothing, decrease for more.
% dU_smooth = Spline smmothing factor for the differential voltage signal. Range: [0 1]. Increase for less smoothing, decrease for more.
% DVA_smooth = Spline smmothing factor for the DVA signal. Range: [0 1]. Increase for less smoothing, decrease for more.

%% Global parameter definition

warning('off','SPLINES:CHCKXYWP:NaNs'); % Turns off warnings about using csaps on data containing NaN fields
SplineSmoothingFactor_I=I_smooth; % Spline smmothing factor for the current signal. Range: [0 1]. Increase for less smoothing, decrease for more.
SplineSmoothingFactor_U=U_smooth; % Spline smmothing factor for the voltage signal. Range: [0 1]. Increase for less smoothing, decrease for more.
SplineSmoothingFactor_dU=dU_smooth; % Spline smmothing factor for the differential voltage signal. Range: [0 1]. Increase for less smoothing, decrease for more.
SplineSmoothingFactor_DVA=DVA_smooth; % Spline smmothing factor for the DVA signal. Range: [0 1]. Increase for less smoothing, decrease for more.

%% Transfer data to main procedure

DVAData.Time  = t_data;
DVAData.U     = U_data;
DVAData.I     = I_data;

%% Differential voltage calculation

% Smoothing input measurement signals
DVAData.U_smooth=csaps(DVAData.Time,DVAData.U,SplineSmoothingFactor_U,DVAData.Time);
DVAData.I_smooth=csaps(DVAData.Time,DVAData.I,SplineSmoothingFactor_I,DVAData.Time);

% Differentiation of the voltage signal and smoothing the output 
DVAData.dU=[NaN;diff(DVAData.U)./diff(DVAData.Time)];
DVAData.dU_smooth=csaps(DVAData.Time,DVAData.dU,SplineSmoothingFactor_dU,DVAData.Time);

% Actual DVA calculation (dU/dQ=(dU/dt)/(dQ/dt)=(dU/dt)/I) and smoothing the output
DVAData.DVA=DVAData.dU_smooth./DVAData.I;
DVAData.DVA_smooth=csaps(DVAData.Time,DVAData.DVA,SplineSmoothingFactor_DVA,DVAData.Time);

%% Calculate charge

DVAData.Charge = (cumtrapz(DVAData.Time, abs(DVAData.I)));

%% Return

Q   = DVAData.Charge;
DVA = DVAData.DVA_smooth;
end