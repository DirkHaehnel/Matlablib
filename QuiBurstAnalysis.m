dirname = 'm:\MTusers\Qui\110601_WW_HTH_ctrl\';


fname = '110601_WW_6MGdnHCl+DTT_25uW_dilute.ht3';

[sync, tcspc, chan, special, num, head] = ht3v2read([dirname fname],[0 inf]);

res = MultiFocus2FCS([dirname fname],4,[],[],0.8e-9);
[y,t]=FCSCrossRead(res);

binwidth = 1e3;
bb = tttr2bin(sync, binwidth);