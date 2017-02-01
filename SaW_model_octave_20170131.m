%Authors: Igor G. Olaizola, Angel Martin, Mikel Zorrilla, Iñigo Tamayo
%Main 2 is based on main. It includes the case of GPUs that wasn't modelled
%in "main.m"
%close all;

close all;
clc; clear; 
main_fig = figure();
%# figure size displayed on screen (50% scaled, but same aspect ratio)
set(main_fig, 'Units','centimeters', 'Position',[0 0 19 23])
%# figure size printed on paper
set(main_fig, 'PaperUnits','centimeters')
%set(main_fig, 'PaperSize',[480 600])
set(main_fig, 'PaperPosition',[0 0 19 23])
%set(main_fig, 'PaperOrientation','portrait')

disp('INIT')
W=5e17; %Work Load
gm = 1e13; gt = gm; gp = gm; gs = gm; %Data communication cost (Server side)
mg = 1000; mgs = mg/10; %management cost on the server side
fpm = 0.15; fpt = 0.15; fpp = .15; fps = 1; %CPU Processing utilization factor
fpgm = 0.30; fpgt = 0.30; fpgp = 0.5; %GPU Processing utilization factor
fbm = 0.15; fbt = fbm; fbp = .15; fbs = 1; %Bandwidth utilization factor
ffpm = 0.5; ffpt = ffpm; ffpp = 0.8; %  factor applied to CPU when GPU available
fbgm = 0.15; fbgt = fbm; fbgp = .15;  %Bandwidth utilization factor

Tpm = 2:200:20000; Tpt = 5:50:5000; Tpp = 1:10:1000; 

ps = 1:1:100; %number of processors
pm_gpu_frac = 0.1; %fraction of mobile device with gpu inside
pt_gpu_frac = 0.3; %fraction of tablet device with gpu inside
pp_gpu_frac = 0.5; %fraction of PC with gpu inside

pgm =  pm_gpu_frac * Tpm; 
pm = (1-pm_gpu_frac) * Tpm;

pgt =  pt_gpu_frac * Tpt; 
pt = (1-pt_gpu_frac) * Tpt; %fraction of tablets with gpu inside

pgp =  pp_gpu_frac * Tpp; 
pp = (1-pp_gpu_frac) * Tpp; %fraction of PCs with gpu inside



Fm = 0.05e9; Ft = 0.08e9; Fp = 2.5e8; Fs = 85e9; %FLOPS
Fgm =  4*Fm; Fgt = 4* Ft ; Fgp = 4*Fp; %FLOPS


ggm = gm; ggt = gm; ggp = gm;  %Data communication cost (Server side)
bm = 3e6; bt = 8e6; bp = 20e6; bs = 6e9; %bandwidth
bgm = bm; bgt = bt; bgp = bp; %we asume that the internal bandiwth is 
%higher than the one reached by the modem

Dm = 1/68;Dt = 1/8;Dp = 1/20;Dgm = 1/25;Dgt = 1/2;Dgp = (1-(Dm+Dt+Dp+Dgm+Dgt));


disp('CALC')

pos = 80;
c = 1;
C = zeros(1,6);
delta = 0.01.*ones(1,5);
s = zeros(1,5); %signs of delta to decide to decrease it
D=1/6.*ones(1,6);
thr = 1; %threshold in %
itermax = 500;
counter = 1;
while (c)
    fprintf('iter:%d',counter)
    Cm =  D(1) .* (W./(fpm .*   pm .* Fm) + gm./(fbm .* bm*pm)+mg.*pm);
    Ct =  D(2) * (W./(fpt .*   pt .* Ft) + gt./(fbt .* bt*pt)+mg.*pt);
    Cp =  D(3) * (W./(fpp .*   pp .* Fp) + gp./(fbp .* bp*pp)+mg.*pp);
    Cgm =  D(4) * (W./(fpgm .*   pgm .* Fgm + ffpm* fpm .*  pgm .* Fm) + ggm./(fbgm .* bgm*pgm)+mg.*pgm);
    Cgt =  D(5) * (W./(fpgt .*   pgt .* Fgt + ffpt*fpt .*   pgt .* Ft) + ggt./(fbgt .* bgt*pgt)+mg.*pgt);
    Cgp =  D(6) * (W./(fpgp .*   pgp .* Fgp + ffpp*fpp .*   pgp .* Fp) + ggp./(fbgp .* bgp*pgp)+mg.*pgp);
        
    C(1) = Cm(pos);  C(2) = Ct(pos); C(3) = Cp(pos); 
    C(4) = Cgm(pos); C(5) = Cgt(pos);C(6) = Cgp(pos);
    
    for i1 = 1:5
        Dev(i1) = abs(C(i1)-mean(C))/mean(C);
        if( (Dev(i1)>(thr/100)) )
            if C(i1)>mean(C)
                D(i1) = D(i1)-delta(i1);
                if s(i1) == 1
                    delta = delta/2;
                    fprintf(' delta/2  ');
                end
                s(i1) = -1;
            elseif C(i1)<mean(C)
                D(i1) = D(i1)+delta(i1);
                if s(i1) == -1
                    delta = delta/2;
                    fprintf(' delta/2  ');
                end
                s(i1) = 1;
            end
            
        end
    end
    Dev(6) = abs(C(6)-mean(C))/mean(C);
    D(6) = (1-sum(D(1:5)));
   %improvements under threshold
   fprintf('\tmax Dev:%f ',max(Dev));
   if(max(Dev)<(thr/100))
       c = 0;
   end
    
   %safety break 
   if counter >= itermax
       c = 0;
   end
   counter = counter + 1;
   fprintf('\n');
end



Cs =  W./(fps .*   ps .* Fs) +  gs./(fbs .* bs./ps) +mgs.*ps ;
CT = (Cm.^-1+Ct.^-1+Cp.^-1+Cgm.^-1+Cgt.^-1+Cgp.^-1).^-1;
disp('PLOT')
%m = plot(1:numel(Cm),log(Cm));
% hold on;
% t = plot(log(Ct),'r')
% p = plot(log(Cp),'g')
% s = plot(log(Cs),'m')
% T = plot(log(CT),'k')
% hold off;
%aux_fig = figure();
%plot(log(CT))
%hold on;
%plot(log(Cs),'m')

m = semilogy(...
    1:numel(Cm),(Cm),':','Color', [0.3 0.3 0.8],'MarkerSize',6,'LineWidth',1,...
    1:numel(Ct),(Ct),'--','Color', [0.3 0.3 0.8],'MarkerSize',6,'LineWidth',1,...
    1:numel(Cp),(Cp),'-.','Color', [0.3 0.3 0.8],'MarkerSize',6,'LineWidth',1,...
    1:numel(Cgm),(Cgm),':','Color', [0.3 0.8 0.3],'MarkerSize',6,'LineWidth',2,...
    1:numel(Cgt),(Cgt),'--','Color', [0.3 0.8 0.3],'MarkerSize',6,'LineWidth',2,...
    1:numel(Cgp),(Cgp),'-.','Color', [0.3 0.8 0.3],'MarkerSize',6,'LineWidth',2,...
    1:numel(Cs),(Cs),'-','Color', [0.0 0.0 0.0],'MarkerSize',6,'LineWidth',3,...
    1:numel(CT),(CT),'-','Color', [0.8 0.3 0.3],'MarkerSize',6,'LineWidth',3
    );
axis([1 100 8e4 1e10])
%legend(sprintf("Mobile Phone from %i to %i units", min(pm), max(pm)), sprintf("Tablet from %i to %i units", min(pt), max(pt)), sprintf("PC from %i to %i units", min(pp), max(pp)), sprintf("Mobile with GPU from %i to %i units", min(pgm), max(pgm)), sprintf("Tablet with GPU from %i to %i units", min(pgt), max(pgt)), sprintf("PC with GPU from %i to %i units", min(pgp), max(pgp)), sprintf("Local Server from %i to %i processors", min(ps), max(ps)), "Total Distributed among a custom zoo of different devices:");
leg = legend(sprintf("Mobile Phone from %i to %i units", 1, 18000), sprintf("Tablet from %i to %i units", 1, 3500), sprintf("PC from %i to %i units", 1, 500), sprintf("Mobile with GPU from %i to %i units", 1, 2000), sprintf("Tablet with GPU from %i to %i units", 1, 1500), sprintf("PC with GPU from %i to %i units", 1, 500), sprintf("Local Server from %i to %i processors", 1, 100), "Total Distributed among a custom zoo of different devices:");
set(leg,'box','off');
%grid on
%set(gca,'XTickLabel',ps)
xlabel('Percentage of computing nodes employed for each device type distribution');
ylabel('Computational Cost to perform W (Log Scale)');
title(sprintf("W=%.1e, g=%.1e mg=%.2f ms=%.2f",W,gm,mg,mgs));

set (main_fig, "visible", "off");
print(main_fig, sprintf("performance_%i_%i_%i.pdf",W/1e9,gm/1e6,round(mg*100)),"-dpdf");
set (main_fig, "visible", "on");

[M,I] = min(Cs)

fprintf('Cm :%i\nCt :%i\nCp :%i\nCgm:%i\nCgt:%i\nCgp:%i\nCT :%d\nCs :%i\n',round(Cm(pos)),round(Ct(pos)),round(Cp(pos)),round(Cgm(pos)),...
    round(Cgt(pos)),round(Cgp(pos)),round(CT(pos)),round(Cs(pos)))
%dist_fig = figure();
sprintf('names = {''Mobile CPU: %.2f%%'',''Tablet CPU: %.2f%%'',''PC CPU: %.2f%%'',''Mobile GPU: %.2f%%'',''Tablet GPU: %.2f%%'',''PC GPU: %.2f%%''}',100*D(1),100*D(2),100*D(3),100*D(4),100*D(5),100*D(6))
%pie(D,names);