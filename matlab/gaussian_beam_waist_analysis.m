clear all; % clear all variables
close all; % close all plots

files=dir('d*.bmp');

%% Process images
for i=1:length(files)
    % get file name
    filecurr=files(i).name;

    % get distance from the file name
    L(i)=str2num(filecurr(2:4))/1000;

    % apply gaussian filter and normalize
    fig=imread(filecurr);
    fig1=double(fig);%rgb2gray(fig);
    PSF = fspecial('gaussian',15,5);
    fig2 = imfilter(fig1,PSF,'conv');
    fig2=fig2/max(max(fig2));

    % find coordinates of the max
    [c,ind]=max(fig2);
    x=find(c==1);
    y=ind(x);

    % get the subfig with width L_x and height L_y
    L_x = 500;
    L_y = 500;
    x_i = x - L_x/2;
    if x_i < 1
        x_i = 1;
    end
    x_f = x + L_x/2;
    if x_f > 1280
        x_f = 1280;
    end
    y_i = y - L_y/2;
    if y_i < 1
        y_i = 1;
    end
    y_f = y + L_y/2;
    if y_f > 1024
        y_f = 1024;
    end
    subfig2 = fig2(y_i:y_f,x_i:x_f); 

    % plot subfig
    h2=figure();
    imagesc(subfig2);
    title(['Image ' num2str(L(i))]);
    saveas(h2,[num2str(i) 'Image ' num2str(L(i))],'png')

    % define fit function on the subfig
    size_subfig = size(subfig2);
    [X,Y] = meshgrid(1:size_subfig(2),1:size_subfig(1));
    gaussian=@(par)sum(sum((exp(-2.*(X-par(1)).^2./(par(2).^2)-2.*(Y-par(3)).^2./(par(4).^2))-subfig2).^2));

    % starting parameters of the fit
    par0=[x-x_i,5,y-y_i,5];

    % fit
    [var,chi]=fminsearch(gaussian,par0);

    % plot final fit
    h1=figure();
    par=var;
    imagesc((exp(-2.*(X-par(1)).^2./(par(2).^2)-2.*(Y-par(3)).^2./(par(4).^2))));
    title(['Fit ' num2str(L(i))]);
    saveas(h1,[num2str(i) 'Image ' num2str(L(i)) '-fit'],'png')
    
    % get waist
    Wx(i)=par(2)*5.2e-6;
    Wy(i)=par(4)*5.2e-6;
    
    % close plots
    close(h1);
    close(h2);
end
save('data','L','Wx','Wy');

%% fit waists
close all;
figure()
lambda = 810e-9;
x=0:0.001:0.4;
wai=@(x,par)par(1).*sqrt(1+((x-par(2)).*lambda./(pi.*par(1).^2)).^2);
par0=[1e-5,-0.2];       
[varWX,chi]=fminsearch(@(par)sum((wai(L,par)-Wx).^2),par0);
plot(L,Wx,'+',x,wai(x,varWX));
text(0,0,['$W_{0X}$=' num2str(varWX(1))],'Interpreter','latex','FontSize',20);
hold on
[varWY,chi]=fminsearch(@(par)sum((wai(L,par)-Wy).^2),par0);
plot(L,Wy,'x',x,wai(x,varWY))
text(0,1e-4,['$W_{0Y}$=' num2str(varWY(1))],'Interpreter','latex','FontSize',20);
ylabel('Waist','FontSize',20);
xlabel('L(m)','FontSize',20);
grid on

disp(['W_{0X}=' num2str(varWX(1))])
disp(['W_{0Y}=' num2str(varWY(1))])
%end

