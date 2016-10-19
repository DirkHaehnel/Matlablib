emin = 0.02;
emax = 0.04;

load viscosity
t = (10:0.1:90)+273.15;

c = polyfit(T,log(visc(1,:)),3);
y = t./exp(polyval(c,t));
z1 = mean([y(1:end-1);y(2:end)])./diff(y,1)*mean(diff(t));

rhoda = [10	2120
20	1650
30	1280
40	980
50	755
24	1550
35	1180
37	1120
45	870
51	740
60	585
71	450
74	400
79	350
81	335
85	280
94	240];
[ind ind] = sort(rhoda(:,1));
rhoda = rhoda(ind,:);

c = polyfit(rhoda(:,1)+273.15,log(rhoda(:,2)),3);
y = exp(polyval(c,t));
%z2 = abs(mean([y(1:end-1);y(2:end)])./diff(y,1)*mean(diff(t)));
z2 = 50/emin*abs(mean(diff(t))./diff(y,1));

gel = [9.945342381     158.4190476
10.94120333     157.8738095
11.96389095     156.8914286
12.97198143     155.6271429
13.97837476     154.3828571
14.98514        154.0609524
15.99437857     153.1690476
17.00112952     152.107619
18.00396762     150.9719048
19.00751333     149.8542857
20.00593762     149.4319048
21.0027519      148.2495238
21.99825524     145.9161905
22.99573762     145.1461905
23.99695952     144.0914286
24.99801952     142.1280952
26.00732143     140.227619
27.00567476     138.8528571
28.00636        136.4319048
29.00804476     133.2980952
30.0344119      128.4042857
31.04304619     122.9333333
32.0423319      103.6009524
33.03634048     66.927
34.03615667     60.5172381
35.02775143     58.20333333
36.03708048     56.9192381
37.02666571     56.1007619
38.02305905     55.59471429
39.02521667     55.20152381
40.01416429     54.82047619
41.00779381     54.57261905
42.00812238     54.4637619
42.99332571     54.43719048
43.98613952     54.18395238
44.98656714     54.11742857
45.96985952     54.1242381
46.96571476     53.9687619
47.96214619     53.94419048
48.94484143     53.76852381
49.94795048     53.73171429
];

t3 = gel(:,1)'+273.15;
close; hp = Simplex('HillPoly',[33 1],[],[],[],[],gel(:,1),gel(:,2));
hp = Simplex('HillPoly',hp,[],[],[],[],gel(:,1),gel(:,2),2);
[z z z] = Hill(hp,gel(:,1),gel(:,2),t-273.15);
c = polyfit(T,log(visc(1,:)),3);
y = t./exp(polyval(c,t))./z';
z3 = abs(mean([y(1:end-1);y(2:end)])./diff(y)*mean(diff(t)));

close
plot(mean([t(1:end-1);t(2:end)]),emin*z2,'Color','b')
line(mean([t(1:end-1);t(2:end)]),emin*z1,'Color','r')
line(mean([t(1:end-1);t(2:end)]),emin*z3,'Color','g')
patch([t(1:end-1);t(2:end);t(2:end);t(1:end-1)],[emin*z1; emin*z1; emax*z1; emax*z1],'r','EdgeColor','none','FaceAlpha',0.3)
line(mean([t(1:end-1);t(2:end)]),emax*z1,'Color','r')
patch([t(1:end-1);t(2:end);t(2:end);t(1:end-1)],[emin*z2; emin*z2; emax*z2; emax*z2],'b','EdgeColor','none','FaceAlpha',0.3)
line(mean([t(1:end-1);t(2:end)]),emax*z2,'Color','b')
patch([t(1:end-1);t(2:end);t(2:end);t(1:end-1)],[emin*z3; emin*z3; emax*z3; emax*z3],'g','EdgeColor','none','FaceAlpha',0.3)
line(mean([t(1:end-1);t(2:end)]),emax*z3,'Color','g')

axis([10+273.15 90+273.15 0 4])
xlabel('temperature (K)');
ylabel('accuracy (K)')
legend({'Lifetime (Rhodamine B)','2fFCS (arbitrary dye)','2fFCS (\mugel)'},4)

