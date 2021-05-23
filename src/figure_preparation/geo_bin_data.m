function [Xbin,Ybin,Ylo,Yhi,CI] = geo_bin_data(X,Y,edges,minN)

Xbin = (edges(1:end-1)+edges(2:end) )./2;
[N,dummy,bin] = histcounts(X,edges);

for ii = 1:length(edges)-1
    
    if any( bin==ii)
        %Ybin(ii) = geomean( Y( (bin==ii) & ~isnan(Y) ) );
        %Ystd(ii) = nanstd(  Y( (bin==ii) ) );
        if sum(bin==ii)>20
            %Ybin(ii) = exp( trimmean(log(Y((bin == ii) & isreal(Y))),20) );
            Ybin(ii) = exp( nanmean(log(Y((bin == ii) & isreal(Y)))) );
        else
            Ybin(ii) = exp( nanmean(log(Y((bin == ii) & isreal(Y)))) );
        end
        Ylo(ii) = exp( log(Ybin(ii)) - nanstd(log(Y((bin == ii) & isreal(Y)))));
        Yhi(ii) = exp( log(Ybin(ii)) + nanstd(log(Y((bin == ii) & isreal(Y)))));
        CI(ii) = exp( 1.96 * nanstd(log(Y((bin == ii) & ~isnan(Y))))./ sqrt( sum((bin == ii) & isreal(Y))));
    else
        Ybin(ii) = NaN;
        Ylo(ii) = NaN;
        Yhi(ii) = NaN;
        CI(ii) = NaN;
    end
end

Ybin( N < minN ) = NaN;
Yhi( N < minN ) = NaN;
Ylo( N < minN ) = NaN;
CI(N<minN) = NaN;