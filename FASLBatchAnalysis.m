%DATA 1 :
    cellNameFile(1,1) = {'D:\FASL\SOFI\20110819\Data\'};
    cellDirectoryOut(1,1) = {'D:\FASL\SOFI\20110819\Tests\'};
    DataTest(1,:) = {'01_Q625_GABAR1' '01_Q625_GABAR1_X1' '01_Q625_GABAR1_X2' '02_Q625_GABAR1' '03_Q625_GABAR1' '03_Q625_GABAR1_X1' '03_Q625_GABAR1_X2' '04_Q625_GABAR1' '04_Q625_GABAR1_X1' '04_Q625_GABAR1_X2' '06_Q625_GABAR1_X2' '06_Q625_GABAR1_X2_X1' '06_Q625_GABAR1_X2_X2' '07_Q625_GABAR1' '07_Q625_GABAR1_X1' '07_Q625_GABAR1_X2'}; 

%DATA 2 :
    cellNameFile(1,2) = {'D:\FASL\SOFI\20110901\Data\'};
    cellDirectoryOut(1,2) = {'D:\FASL\SOFI\20110901\Tests\'};
    DataTest(2,:)= { '01_KDEL_R1_525nm' '01_KDEL_R1_625nm' '02_KDEL_R1_625nm' '03_KDEL_R1_625nm' '04_KDEL_R1_525nm' '05_KDEL_R1_525nm' '06_KDEL_R1_625nm' '07_KDEL_R1_625nm' '08_KDEL_R1_525nm' '09_KDEL_R1_625nm' '10_KDEL_R1_525nm' '11_KDEL_R1_625nm' '12_KDEL_R1_525nm' }; 

%DATA 3 :
    cellNameFile(1,3) = {'D:\FASL\SOFI\20110902\Data\'};
    cellDirectoryOut(1,3) = {'D:\FASL\SOFI\20110902\Tests\'};
    DataTest(3,:) = {'01_KDEL_R1_625nm' '02_KDEL_R1_525nm' '03_KDEL_R1_625nm' '04_KDEL_R1_625nm' '05_KDEL_R1_525nm' '06_KDEL_R1_525nm' '07_KDEL_R1_625nm' '08_KDEL_R1_625nm' '09_KDEL_R1_525nm' '10_KDEL_R1_525nm' '11_KDEL_R1_625nm' '12_KDEL_R1_525nm' '13_KDEL_R1_625nm' '14_KDEL_R1_625nm' '15_KDEL_R1_525nm' '16_KDEL_R1_525nm' '17_KDEL_R1_625nm'}; 

    %LAST
%DATA 4 :
    cellNameFile(1,4) = {'D:\FASL\SOFI\20110818\Data\'};
    cellDirectoryOut(1,4) = {'D:\FASL\SOFI\20110818\Tests\'};
    DataTest(4,:) = {'test1' 'test2' 'test3' 'test4' 'test5' 'test6' 'test7'};
%DATA 5 :
    cellNameFile(1,5) = {'D:\FASL\SOFI\20110826\Data\'};
    cellDirectoryOut(1,5) = {'D:\FASL\SOFI\20110826\Tests\'};
    DataTest(5,:) = {'01_KDEL_ERGIC_525nm' '02_KDEL_ERGIC_525nm' '03_KDEL_ERGIC_625nm'};
%DATA 6 :
    cellNameFile(1,6) = {'D:\FASL\SOFI\20110829\Data\'};
    cellDirectoryOut(1,6) = {'D:\FASL\SOFI\20110829\Tests\'};
    DataTest(6,:) = {'01_KDEL_R1_625nm' '02_KDEL_R1_525nm' '03_KDEL_R1_525nm' '04_KDEL_R1_625nm' '05_KDEL_R1_625nm' '06_KDEL_R1_525nm' '07_KDEL_R1_525nm' '08_KDEL_R1_625nm'};
%DATA 7 :
    cellNameFile(1,7) = {'D:\FASL\SOFI\20110830A\Data\'};
    cellDirectoryOut(1,7) = {'D:\FASL\SOFI\20110830A\Tests\'};
    DataTest(7,:) = {'01_KDEL_R1_625nm' '02_KDEL_R1_525nm' '03_KDEL_R1_525nm' '04_KDEL_R1_625nm' '05_KDEL_R1_625nm' '06_KDEL_R1_525nm'};
%DATA 8 :
    cellNameFile(1,8) = {'D:\FASL\SOFI\20110830B\Data\'};
    cellDirectoryOut(1,8) = {'D:\FASL\SOFI\20110830B\Tests\'};
    DataTest(8,:) = {'01_KDEL_ERGIC_525nm' '02_KDEL_ERGIC_625nm' '03_KDEL_ERGIC_625nm' '04_KDEL_ERGIC_525nm' '05_KDEL_ERGIC_525nm' '06_KDEL_ERGIC_625nm'};
%DATA 9 :
    cellNameFile(1,9) = {'D:\FASL\SOFI\20110831A\Data\'};
    cellDirectoryOut(1,9) = {'D:\FASL\SOFI\20110831A\Tests\'};
    DataTest(9,:) = {'01_Golgi_ERGIC_625nm' '02_Golgi_ERGIC_525nm' '03_Golgi_ERGIC_625nm' '04_Golgi_ERGIC_525nm' '05_Golgi_ERGIC_525nm' '06_Golgi_ERGIC_625nm' '07_Golgi_ERGIC_625nm' '08_Golgi_ERGIC_625nm' '09_Golgi_ERGIC_625nm' '10_Golgi_ERGIC_525nm'};
%DATA 10 :
    cellNameFile(1,10) = {'D:\FASL\SOFI\20110831B\Data\'};
    cellDirectoryOut(1,10) = {'D:\FASL\SOFI\20110831B\Tests\'};
    DataTest(10,:) = {'01_Golgi_R1_525nm' '02_Golgi_R1_625nm' '03_Golgi_R1_625nm' '04_Golgi_R1_525nm'};

for c = 4:size(cellNameFile,2)
    externNameFile = char(cellNameFile(1,c));
    externDirectoryOut = char(cellDirectoryOut(1,c));
    for f=1:size(DataTest,2)
        externActualFile = char(DataTest(c,f));
        FASLAnalysis;
    end
end