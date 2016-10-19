function FCSScalePlot(t, x, k, j);

% FCSScalePlot(t, x, k, j) scales and plots fcs-curves, where t is the time-axis,
% x are the FCS-curves

if nargin<3
    k=10;
end
if nargin<4
    j=10;
end

semilogx(t(k:end), FCSScale(x, k, j),'linewidth',2);
xlabel('time [s]');
ylabel('autocorrelation');
axis tight
