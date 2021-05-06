function [out, iout] = chunk(in)
% function [out, iout] = chunk(in)
%
% Return the largest chunk of data which doesn't contain any NaN's.
%
% Modified from chunk.m, by Felix Tubiana

% Deborah Le Bel
%    06.30.2011
%
%
% --------------------------------------------------------------------

%
%
if nargin < 1, help chunk, return, end
in = in(:);
bad = find(isnan(in));
if isempty(bad) | (length(bad) == length(in))
   out = in;
   iout = [1:length(in)];
   return
end

start = unique([1; bad+1]);
stop  = unique([bad-1; length(in)]);
indx = find( (stop-start) == max( (stop-start)) );
rows = size(indx,2);
cols = length(start(indx(1)):stop(indx(1)));
iout = NaN*ones(rows,cols);
out = NaN*ones(rows,cols);
for i = 1:length(indx)
   iout(i,:) = start(indx(i)):stop(indx(i));
   out(i,:) = in(iout(i,:));
end

