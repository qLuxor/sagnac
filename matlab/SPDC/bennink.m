function out = bennink( lp,ls,li,wp,ws,wi,axp,axs,axi,L,Lambda )
%BENNINK Summary of this function goes here
%   Detailed explanation goes here

    function out = func(z,lp,ls,li,wp,ws,wi,axp,axs,axi,Lambda)


        function out = qp(z)
            out = wp.^2 + 1i.*z.*lp./pi./n(lp,axp);
        end

        function out = qs(z)
            out = ws.^2 + 1i.*z.*ls./pi./n(ls,axs);
        end

        function out = qi(z)
            out = wi.^2 + 1i.*z.*li./pi./n(li,axi);
        end

        out = wp.*ws.*wi.*exp(1i.* 2.*pi.*(n(lp,axp)./lp - n(ls,axs)./ls - n(li,axi)./li + 1./Lambda ).*z) ./ (conj(qs(z)).*conj(qi(z)) + qp(z).*conj(qi(z)) + qp(z).*conj(qs(z)));
    end

    f = @(z) func(z,lp,ls,li,wp,ws,wi,axp,axs,axi,Lambda);
    out = integral(f,-L/2,L/2);


end

