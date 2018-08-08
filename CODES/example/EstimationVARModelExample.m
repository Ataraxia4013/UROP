%% Estimate VAR(4) Model
% Fit a VAR(4) model to consumer price index (CPI) and unemployment
% rate data.
%%
% Load the |DATA_USEconModel| data set.
load Data_USEconModel
%%
% Plot the two series on separate plots.
figure;
plot(DataTable.Time,DataTable.CPIAUCSL);
ylabel('Consumer price index');
xlabel('Date');
figure;
plot(DataTable.Time,DataTable.UNRATE);
ylabel('Unemployment rate');
xlabel('Date');
%%
% Stabilize CPI by converting it to a series of growth rates. Synchronize
% the unemployment rate series with the CPI growth rate by removing its
% first observation.
rcpi = price2ret(DataTable.CPIAUCSL);
unrate = DataTable.UNRATE(2:end);
%%
% Create a default VAR(4) model using the shorthand syntax.
Mdl = varm(2,4)
%%
% |Mdl| is a |varm| model object. All properties containing |NaN| values
% correspond to parameters to be estimated given data.
%%
% Estimate the model using the entire data set.
EstMdl = estimate(Mdl,[rcpi unrate])
%%
% |EstMdl| is an estimated |varm| model object.  It is fully specified in
% the sense that all parameters have known values.  The description
% indicates that the autoregressive polynomial is stationary.
%%
% Display a summary statistics from the estimation.
summarize(EstMdl)