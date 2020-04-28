function [ output_args ] = chernobyl_fallout( ~ ) %#ok<STOUT>
% model of chernobyl fission product decay
% compare to figures found at https://en.wikipedia.org/wiki/Chernobyl_disaster

% decay chains (g is gamma radiation)
%  I-131 -> Xe-131 + g
% Xe-140 -> Cs-140 -> Ba-140 -> La-140 -> Ce-140 + g
% Xe-137 -> Cs-137 + g -> Ba-137 + g
% Cs-134 -> Ba-134 + g
% Sn-132 -> Sb-132 + g -> Te-132 ->  I-132 + g -> Xe-132 + g
% Zr-95  -> Nb-95  + g -> Mo-95
% Ru-103 -> Rh-103 + g
% Ru-106 -> Rh-106 + g
% Xe-133 -> Cs-133

% to match relative contribution units of data, two
% "other" contributions are adjusted empirically

% other1 -> fast decay
% other2 -> slow decay

%Rate constants (/day)
log2 = log(2);
k_I131 = log2 / 8.0197;         % I-131 decay
k_Xe140 = log2 * 86400 / 14;       %Xe-140 decay chain
k_Cs140 = log2 * 86400 / 64;       %Cs-140
k_Ba140 = log2 / 13;             %Ba-140
k_La140 = log2 * 24 / 40;          %La-140
k_Xe137 = log2 * 1440 / 3.818;     %Xe-137 decay chain
k_Cs137 = log2 / 30.17 / 365.25;   %Cs-137
k_Cs134 = log2 / 2.0652 / 365.25;  %Cs-134 decay
k_Sn132 = log2 * 86400 / 39.7;     %Sn-132 decay chain
k_Sb132 = log2 * 1440 / 2.79;     %Sb-132
k_Te132 = log2 / 3.204;           %Te-132
k_I132 = log2 * 24 / 2.295;      % I-132
k_Zr95 = log2 / 64.032;          %Zr-95  decay chain
k_Nb95 = log2 / 34.991;          %Nb-95
k_Ru103 = log2 / 39.26;           %Ru-103 decay
k_Ru106 = log2 / 373.59;          %Ru-106 decay
k_Xe133 = log2 / 5.2475;          %Xe-133 decay
k_other_fast = log2 / 3;
k_other_slow = log2 / 100;


% parameters for bateman equation
time = linspace(0, 10000, 10000)';

k_Xe140_chain = [k_Xe140; k_Cs140; k_Ba140; k_La140];
k_Xe137_chain = [k_Xe137; k_Cs137];
k_Sn132_chain = [k_Sn132; k_Sb132; k_Te132; k_I132];
k_Zr95_chain = [k_Zr95; k_Nb95];
    
% calculate time-dependent isotope amounts
nt_I131 = bateman(22.4, k_I131, time);
nt_Xe140 = bateman(15.5, k_Xe140_chain, time);
nt_Xe137 = bateman(1.7, k_Xe137_chain, time);
nt_Cs134 = bateman(1.18, k_Cs134, time);
nt_Sn132 = bateman(34.46, k_Sn132_chain, time);
nt_Zr95 = bateman(11.6, k_Zr95_chain, time);
nt_Ru103 = bateman(3.05, k_Ru103, time);
nt_Ru106 = bateman(0.41, k_Ru106, time);
nt_Xe133 = bateman(4.85, k_Xe133, time);
nt_other_fast = bateman(4.6, k_other_fast, time);
nt_other_slow = bateman(0.25, k_other_slow, time);

% convert to relative percent contribution for every time point
per_tot = (nt_I131 + sum(nt_Xe140, 2) + sum(nt_Xe137, 2) + nt_Cs134 ...
    + sum(nt_Sn132, 2) + sum(nt_Zr95, 2) + nt_Ru103 + nt_Ru106 ...
    + nt_Xe133 + nt_other_fast + nt_other_slow)/100;
rp_I131 = nt_I131 ./ per_tot;
rp_Xe140 = nt_Xe140 ./ repmat(per_tot, 1, 4);
rp_Xe137 = nt_Xe137 ./ repmat(per_tot, 1, 2);
rp_Cs134 = nt_Cs134 ./ per_tot;
rp_Sn132 = nt_Sn132 ./ repmat(per_tot, 1, 4);
rp_Zr95 = nt_Zr95 ./ repmat(per_tot, 1, 2);
rp_Ru103 = nt_Ru103 ./ per_tot;
rp_Ru106 = nt_Ru106 ./ per_tot;
rp_Xe133 = nt_Xe133 ./ per_tot;
rp_other_fast = nt_other_fast ./ per_tot;
rp_other_slow = nt_other_slow ./ per_tot;

% rate gamma radiation emission
gamma = 100 * k_I131 * rp_I131 + 100 * k_La140 * rp_Xe140(:,4) ...
    + 2000 * k_Cs137 * rp_Xe137(:,2) + 1000 * k_Cs134 * rp_Cs134 ...
    + k_Sn132 * rp_Sn132(:,1) + 100 * k_Te132 * rp_Sn132(:,3) ...
    + 100 * k_I132 * rp_Sn132(:,4) + 100 * k_Zr95 * rp_Zr95(:,1) ...
    + 800 * k_Ru103 * rp_Ru103 + 800 * k_Ru106 * rp_Ru106 ...
    + 100 * k_other_fast * rp_other_fast + 100 * k_other_slow * rp_other_slow;
disp([size(time), size(gamma)])
% figure 1 (gamma emission)
scrsz = get(0, 'ScreenSize');
figure('Position', [5, scrsz(4) / 2 - 70, scrsz(3) / 2 - 68, scrsz(4) / 2]);
loglog(time(2:end), gamma(2:end));
title('figure 1: chernobyl gamma radiation');
xlabel('time (days)');
ylabel('relative gamma dose rate');
xlim([1 10000]);
ylim([10 10000]);

% figure 2 (isotope detection)
figure('Position', [scrsz(3) / 2 - 55, scrsz(4) / 2 - 70, ...
    scrsz(3) / 2 + 55, scrsz(4) / 2]);
semilogx(time, [rp_Xe133, rp_I131, rp_Sn132(:,3) + rp_Sn132(:,4), ...
    rp_Xe140(:,3) + rp_Xe140(:,4), rp_Zr95(:,1) + rp_Zr95(:,2), ...
    rp_other_fast + rp_other_slow, rp_Ru103 + rp_Ru106, rp_Cs134, ...
    rp_Xe137(:,2)]);
title('figure 2: chernobyl isotope fallout');
legend('Xe', 'I-131', 'Te-132/I-132', 'Ba-140/La-140', 'Zr-95/Nb-95', ...
    'Others', 'Ru', 'Cs-134', 'Cs-137', 'Location', 'NorthWest');
xlabel('time (days)');
ylabel('% contribution to detecable signal');
xlim([1 10000]);
ylim([0 100]);
end