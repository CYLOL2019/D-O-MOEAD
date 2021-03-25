SetPATH;
problems = 20;
exetime = zeros(problems,1);
for pn = problems : problems
    %% GAMS input.csv
    pfile = sprintf('../../problem/portfolio problem/port%d.txt',pn);
%     pfile = fopen(pfile);
    input = textread(pfile);
%     input = textscan(pfile,'%f');
    [NoA u Covariance] = DataInput(input);
    Original_r = u;
    Original_Q = Covariance;
    mu = u;
    cov = Covariance;
    pfile = sprintf('../../problem/portfolio problem/portef%d.txt',pn);
%     pfile = fopen(pfile);
    input = textread(pfile);
%     input = textscan(pfile,'%f');
    nor = zeros(4,1);
    nor(3) = input(1,1);
    nor(4) = input(end,1);
    nor(1) = input(1,2);
    nor(2) = input(end,2);
    iwgdx('MtoG','nor','mu','cov');
    gdxWhos('MtoG');
    t1 = clock;
    disp(datestr(now));
    system 'gams AUGMECON2_MIQP2I_Portfolio --gdxin=MtoG logOption=0';
%     system 'gams AUGMECON2_MIQP2I_Portfolio --gdxin=MtoG';
    t2 = clock;
    exetime(pn) = etime(t2,t1);
    exetime(pn)
    timefile = sprintf('extime%d.mat',pn);
    save(timefile,'exetime');
end
% end
