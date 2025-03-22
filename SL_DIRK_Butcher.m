function [A,b,c] = SL_DIRK_Butcher(s,p,q,scheme_no)
%=========================================================================%
% SL_DIRK_Butcher: Returns the Butcher tableau for specified DIRK methods.
%
% Inputs:
%   s          - Number of stages
%   p          - Order of accuracy
%   q          - Weak stage order
%   scheme_no  - Scheme identifier (determines the specific DIRK method)
%
% Outputs:
%   A          - Coefficient matrix (Butcher tableau)
%   b          - Weights vector
%   c          - Nodes vector
%=========================================================================%
%-------------------------------------------------------------------------%
% 4th order A stable stiffly accurate DIRK for comparison: DIRK(5,4,1)
%-------------------------------------------------------------------------%
    if s==5 && p ==4 && q==1 && scheme_no == 1
       A = [.25,0,0,0,0;
            .50,.25,0,0,0;
            .34,-.04,.25,0,0;
            371/1360,-137/2720,15/544,.25,0;
            25/24,-49/48,125/16,-85/12,.25];
        b = [25/24,-49/48,125/16,-85/12,.25];
        c = [0.2500;0.7500;0.5500;0.5000;1.0000];
%-------------------------------------------------------------------------%
% 5th order stiffly accurate DIRK for comparison: DIRK(5,5,1)
%-------------------------------------------------------------------------%
    elseif s==5 && p==5 && q==1 && scheme_no == 2
        A = [ 4024571134387/14474071345096  0                             0                             0                             0;
              9365021263232/12572342979331  4024571134387/14474071345096  0                             0                             0;
              2144716224527/9320917548702  -397905335951/4008788611757    4024571134387/14474071345096  0                             0;
             -291541413000/6267936762551    226761949132/4473940808273   -1282248297070/9697416712681   4024571134387/14474071345096  0; 
             -2481679516057/4626464057815  -197112422687/6604378783090    3952887910906/9713059315593   4906835613583/8134926921134   4024571134387/14474071345096];
        b = [-2522702558582/12162329469185, 1018267903655/12907234417901, 4542392826351/13702606430957, 5001116467727/12224457745473, 1509636094297/3891594770934];
        c = [4024571134387/14474071345096; 5555633399575/5431021154178; 5255299487392/12852514622453; 3/20; 10449500210709/14474071345096];
%--------------------------------------------------------------------------%
% SL order 3 for ODEs and SL order 4 for PDEs. 
% A stable stiffly accurate: ESDIRK(8,4,3)
%-------------------------------------------------------------------------%
    elseif s ==8 && p ==4 && q ==3 && scheme_no == 3
        % Old
        A = [0,0,0,0,0,0,0,0;0.24800000000000000E0,0.24800000000000000E0,0,0,0,0,0,0;0.29270060483870968E0,0.33429939516129032E0,0.24800000000000000E0,0,0,0,0,0;(-0.51362481734263786E-1),(-0.51362481734263786E-1),0,0.24800000000000000E0,0,0,0,0;0.35369973435638459E0,0.35132507700050000E0,(-0.14157022829091814E-1),(-0.92142825059265208E-1),0.24800000000000000E0,0,0,0;0.35192957487563523E0,0.28906373575611147E0,(-0.37478801620774382E0),(-0.13213263862149485E0),0.33333333333333333E0,0.24800000000000000E0,0,0;0.82414187342077589E0,0.82787555937212080E0,0.22259159671547227E-1,(-0.10397878035287728E1),(-0.84615384615384615E0),0.96366505721817505E0,0.24800000000000000E0,0;(-0.22913012140578189E-1),0.26472048593920004E-1,0.29441950035521513E0,0.40530401868189481E0,(-0.63930341710591595E0),0.78826133225267898E0,(-0.10024047063721478E0),0.24800000000000000E0];
        b = [(-0.22913012140578189E-1),0.26472048593920004E-1,0.29441950035521513E0,0.40530401868189481E0,(-0.63930341710591595E0),0.78826133225267898E0,(-0.10024047063721478E0),0.24800000000000000E0];
        c = [0,0.49600000000000000E0,0.87500000000000000E0,0.14527503653147243E0,0.84672496346852757E0,0.71540598913584137E0,0.10000000000000000E1,0.10000000000000000E1]';
        % Latest from Steven
        % A = [0,0,0,0,0,0,0,0;0.248000000000000000E0,0.248000000000000000E0,0,0,0,0,0,0;0.243935483870967742E0,0.508064516129032258E0,0.248000000000000000E0,0,0,0,0,0;(-0.513624817342637861E-1),(-0.513624817342637861E-1),0,0.248000000000000000E0,0,0,0,0;0.277174770246897791E0,0.281127459876945706E0,0.371131874399222943E-2,0.367114146006918458E-1,0.248000000000000000E0,0,0,0;0.656430932516777865E0,0.983116140854959422E0,0.306736184868580828E0,(-0.764591554818762461E0),(-0.714285714285714286E0),0.248000000000000000E0,0,0;0.776303285539662910E0,0.871092929875633608E0,0.890013172514585640E-1,(-0.985099626921967170E0),(-0.962962962962962963E0),0.963665057218175051E0,0.248000000000000000E0,0;(-0.276299824073057127E-1),0.486470165237278563E-2,0.305103970506605148E-1,0.416600580608291813E0,(-0.360366558519483595E0),0.788261332252678975E0,(-0.100240470637214781E0),0.248000000000000000E0];
        % b = [(-0.276299824073057127E-1),0.486470165237278563E-2,0.305103970506605148E-1,0.416600580608291813E0,(-0.360366558519483595E0),0.788261332252678975E0,(-0.100240470637214781E0),0.248000000000000000E0];
        % c = [0,0.496000000000000000E0,0.100000000000000000E1,0.145275036531472428E0,0.846724963468527572E0,0.715405989135841368E0,0.100000000000000000E1,0.100000000000000000E1]';

%--------------------------------------------------------------------------%
% SL order 4 for ODEs and SL order 4 for PDEs. Stage order = 1.
% A stable stiffly accurate: ESDIRK(7,4,4).
%-------------------------------------------------------------------------%
    elseif s==7 && p==4 && q==4 && scheme_no == 4
        % A =  [ 0.000000050675404                   0                   0                   0                   0                   0                   0
        %        0.714129386047127   0.714129335371719                   0                   0                   0                   0                   0
        %        0.110476243329833  -0.029602006604617   0.301826571358646                   0                   0                   0                   0
        %       -5.974636867642467   1.829731241101338   7.564196029481168   0.821730121929775                   0                   0                   0
        %      -10.161124731847998  -3.676913695391693  14.046214582218884  -0.043611983848188   0.386230239553503                   0                   0
        %       14.662471220721153  -3.790771578863373 -10.706249908456313   0.038219014098556  -0.203999027826617   0.561450587039992                   0
        %        0.125099251818450  -0.064493475313523   0.571013965875813   0.000402281931904  -0.008685988240820  -0.000000000000000   0.376663963928177];
        % b = A(end,:);
        % c = sum(A,2);
        A = [0,0,0,0,0,0,0;0.635590872820565575E0,0.635590872820565575E0,0,0,0,0,0;0.983262502142798822E-1,(-0.263464393396968004E-1),0.268632311303142846E0,0,0,0,0;0.764458877470564230E1,0.179116098091253414E1,(-0.656175202937676871E1),0.548945472526618387E0,0,0,0;0.909651966561551601E1,0.219452704150501351E1,(-0.842175699107812196E1),0.181803934754497280E0,0.371849547971121284E0,0,0;(-0.771769834340603586E0),0.576964753493517400E1,(-0.112350363897065485E1),(-0.209709693509362848E0),0.204279033401788590E0,0.105010351556378430E1,0;0.988874617234384450E-1,(-0.103949735712193650E0),0.561554315511425850E0,(-0.882214420498338574E-1),0.859368520628129828E-1,0.738514422533704999E-3,0.445054034041816525E0];
        b = [0.988874617234384450E-1,(-0.103949735712193650E0),0.561554315511425850E0,(-0.882214420498338574E-1),0.859368520628129828E-1,0.738514422533704999E-3,0.445054034041816525E0];
        c = [0,0.127118174564113115E1,0.340612122177725928E0,0.342294319876802612E1,0.342294319876802612E1,0.491904691708012560E1,0.100000000000000000E1]';
%-------------------------------------------------------------------------%
% SL order 4 for ODEs and SL order 5 for PDEs. Stage order = 2.
% A stable stiffly accurate: ESDIRK(10,5,4).
%-------------------------------------------------------------------------%
    elseif s==10 && p ==5 && q==4 && scheme_no == 5
        A = [0,0,0,0,0,0,0,0,0,0;0.278053841136452325E0,0.278053841136452325E0,0,0,0,0,0,0,0,0;(-0.575868360343262785E-1),(-0.575868360343262785E-1),0.278053841136452325E0,0,0,0,0,0,0,0;(-0.180642456048994587E-1),(-0.180642456048994587E-1),0.508074650073346593E0,0.278053841136452325E0,0,0,0,0,0,0;0.265854318479767113E0,0.265854318479767113E0,0.106904188570680116E0,0,0.278053841136452325E0,0,0,0,0,0;(-0.325466509175089601E-1),(-0.325466509175089601E-1),(-0.701033964442915475E-1),0,0,0.278053841136452325E0,0,0,0,0;(-0.167904264379527886E0),(-0.167904264379527886E0),0.108254891523028329E1,(-0.651582215997644616E-1),0.332786922620662875E-1,(-0.421486126841410240E0),0.278053841136452325E0,0,0,0;(-0.109670417213868564E0),(-0.109670417213868564E0),0.595390642560852413E0,0.125116691580959791E0,(-0.416368027944150149E-1),(-0.386714019015184213E0),(-0.128647296818705952E0),0.278053841136452325E0,0,0;0.190877500952266521E0,0.190877500952266521E0,0.638300555626834943E0,0.198175542907826854E1,(-0.604991443212251072E0),(-0.160424853428720174E1),(-0.193306077833444137E1),0.186243592808780533E1,0.278053841136452325E0,0;(-0.124204710873885030E0),(-0.124204710873885030E0),0,(-0.122471035700209359E1),0.117062718337373206E1,0.121793319273673784E1,0.149047633502855261E1,(-0.102182754466399741E1),(-0.662143228861613777E0),0.278053841136452325E0];
        b = [(-0.124204710873885030E0),(-0.124204710873885030E0),0,(-0.122471035700209359E1),0.117062718337373206E1,0.121793319273673784E1,0.149047633502855261E1,(-0.102182754466399741E1),(-0.662143228861613777E0),0.278053841136452325E0];
        c = [0,0.556107682272904650E0,0.162880169067799768E0,0.750000000000000000E0,0.916666666666666667E0,0.142857142857142857E0,0.571428571428571429E0,0.222222222222222222E0,0.100000000000000000E1,0.100000000000000000E1]';
%-------------------------------------------------------------------------%
% SL order 4 for ODEs and SL order 5 for PDEs. Stage order = 2. This is NOT
% stiffly accurate, A stable: ESDIRK(19,5,4).
%-------------------------------------------------------------------------%
    elseif s==19 && p ==5 && q==4 && scheme_no == 6
        A = [0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0;0.635590872820565575E0,0.635590872820565575E0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0;0.983262502142798822E-1,(-0.263464393396968004E-1),0.268632311303142846E0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0;0.764458877470564230E1,0.179116098091253414E1,(-0.656175202937676871E1),0.548945472526618387E0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0;0.909651966561551601E1,0.219452704150501351E1,(-0.842175699107812196E1),0.181803934754497280E0,0.371849547971121284E0,0,0,0,0,0,0,0,0,0,0,0,0,0,0;(-0.771769834340603586E0),0.576964753493517400E1,(-0.112350363897065485E1),(-0.209709693509362848E0),0.204279033401788590E0,0.105010351556378430E1,0,0,0,0,0,0,0,0,0,0,0,0,0;0.988874617234384450E-1,(-0.103949735712193650E0),0.561554315511425850E0,(-0.882214420498338574E-1),0.859368520628129828E-1,0.738514422533704999E-3,0.445054034041816525E0,0,0,0,0,0,0,0,0,0,0,0,0;0.317795436410282787E0,0,0,0,0,0,0,0.317795436410282787E0,0,0,0,0,0,0,0,0,0,0,0;0.491631251071399411E-1,0,0,0,0,0,0,(-0.131732196698484002E-1),0.134316155651571423E0,0,0,0,0,0,0,0,0,0,0;0.382229438735282115E1,0,0,0,0,0,0,0.895580490456267072E0,(-0.328087601468838435E1),0.274472736263309193E0,0,0,0,0,0,0,0,0,0;0.454825983280775800E1,0,0,0,0,0,0,0.109726352075250675E1,(-0.421087849553906098E1),0.909019673772486398E-1,0.185924773985560642E0,0,0,0,0,0,0,0,0;(-0.385884917170301793E0),0,0,0,0,0,0,0.288482376746758700E1,(-0.561751819485327426E0),(-0.104854846754681424E0),0.102139516700894295E0,0.525051757781892150E0,0,0,0,0,0,0,0;0.494437308617192225E-1,0,0,0,0,0,0,(-0.519748678560968252E-1),0.280777157755712925E0,(-0.441107210249169287E-1),0.429684260314064914E-1,0.369257211266852500E-3,0.222527017020908262E0,0,0,0,0,0,0;0.494437308617192225E-1,0,0,0,0,0,0,(-0.519748678560968252E-1),0.280777157755712925E0,(-0.441107210249169287E-1),0.429684260314064914E-1,0.369257211266852500E-3,0.540322453431191050E0,0.317795436410282787E0,0,0,0,0,0;0.494437308617192225E-1,0,0,0,0,0,0,(-0.519748678560968252E-1),0.280777157755712925E0,(-0.441107210249169287E-1),0.429684260314064914E-1,0.369257211266852500E-3,0.271690142128048203E0,(-0.131732196698484002E-1),0.134316155651571423E0,0,0,0,0;0.494437308617192225E-1,0,0,0,0,0,0,(-0.519748678560968252E-1),0.280777157755712925E0,(-0.441107210249169287E-1),0.429684260314064914E-1,0.369257211266852500E-3,0.404482140437372941E1,0.895580490456267072E0,(-0.328087601468838435E1),0.274472736263309193E0,0,0,0;0.494437308617192225E-1,0,0,0,0,0,0,(-0.519748678560968252E-1),0.280777157755712925E0,(-0.441107210249169287E-1),0.429684260314064914E-1,0.369257211266852500E-3,0.477078684982866627E1,0.109726352075250675E1,(-0.421087849553906098E1),0.909019673772486398E-1,0.185924773985560642E0,0,0;0.494437308617192225E-1,0,0,0,0,0,0,(-0.519748678560968252E-1),0.280777157755712925E0,(-0.441107210249169287E-1),0.429684260314064914E-1,0.369257211266852500E-3,(-0.163357900149393531E0),0.288482376746758700E1,(-0.561751819485327426E0),(-0.104854846754681424E0),0.102139516700894295E0,0.525051757781892150E0,0;0.494437308617192225E-1,0,0,0,0,0,0,(-0.519748678560968252E-1),0.280777157755712925E0,(-0.441107210249169287E-1),0.429684260314064914E-1,0.369257211266852500E-3,0.271970747882627485E0,(-0.519748678560968252E-1),0.280777157755712925E0,(-0.441107210249169287E-1),0.429684260314064914E-1,0.369257211266852500E-3,0.222527017020908262E0];
        b = [0.461474821376046077E-1,0.692998238081291003E-2,(-0.374369543674283900E-1),0.588142946998892382E-2,(-0.572912347085419885E-2),(-0.492342948355803333E-4),(-0.296702689361211016E-1),(-0.554398590465032802E-1),0.299495634939427120E0,(-0.470514357599113906E-1),0.458329877668335908E-1,0.393874358684642666E-3,0.290102131074802651E0,(-0.554398590465032802E-1),0.299495634939427120E0,(-0.470514357599113906E-1),0.458329877668335908E-1,0.393874358684642666E-3,0.237362151488968813E0];
        c = [0,0.127118174564113115E1,0.340612122177725928E0,0.342294319876802612E1,0.342294319876802612E1,0.491904691708012560E1,0.100000000000000000E1,0.635590872820565575E0,0.170306061088862964E0,0.171147159938401306E1,0.171147159938401306E1,0.245952345854006280E1,0.500000000000000000E0,0.113559087282056557E1,0.670306061088862964E0,0.221147159938401306E1,0.221147159938401306E1,0.295952345854006280E1,0.100000000000000000E1]';
    end 
end
%=========================================================================%




