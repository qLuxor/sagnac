lp = 405e-9;
% ls = 810e-9; li = ls;
wp = 35e-6;
ws = 50e-6; wi = ws;
axp = 'y'; axi = axp;
axs = 'z';
L = 0.03;
Lambda = 10e-6;


ls = 809.9e-9:5e-12:810.1e-9;
li = 809.9e-9:5e-12:810.1e-9;

psi = zeros(size(ls,2),size(li,2));

for i = 1:length(ls)
    disp(i)
    for j = 1:length(li)
        psi(i,j) = bennink(1/(1/ls(i) + 1/li(j)),ls(i),li(j),wp,ws,wi,axp,axs,axi,L,Lambda);
    end
end

%%
imagesc(ls,li,abs(psi));
colormap(jet);
colorbar;