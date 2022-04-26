 function genIRNSS=IRNSS_PRN(PRN_ID, signal_type)
L=2^10-1;
%  signal_type = 'L5';
%  PRN_ID = 1;
%% generate
%%
    if signal_type == 'L5'
        if PRN_ID == 1
            loadReg2 = [0 0 1 0 0 1 1 0 0 0];
        elseif PRN_ID == 2
            loadReg2 = [0 0 0 0 1 0 0 1 1 0];
        elseif PRN_ID == 3
            loadReg2 = [1 0 0 0 1 1 0 1 0 0];
        elseif PRN_ID == 4
            loadReg2 = [0 1 0 1 1 1 0 0 1 0];
        elseif PRN_ID == 5
            loadReg2 = [1 1 1 0 1 1 0 0 0 0];
        elseif PRN_ID == 6
            loadReg2 = [0 0 0 1 1 0 1 0 1 1];
        elseif PRN_ID == 7
            loadReg2 = [0 0 0 0 0 1 0 1 0 0];
        else
            disp('PRN_ID is incorrect. Must be from 1 to 7');
        end;
    elseif signal_type == 'S'    
        if PRN_ID == 1
            loadReg2 = [0 0 1 1 1 0 1 1 1 1];
        elseif PRN_ID == 2
            loadReg2 = [0 1 0 1 1 1 1 1 0 1];
        elseif PRN_ID == 3
            loadReg2 = [1 0 0 0 1 1 0 0 0 1];
        elseif PRN_ID == 4
            loadReg2 = [0 0 1 0 1 0 1 0 1 1];
        elseif PRN_ID == 5
            loadReg2 = [1 0 1 0 0 1 0 0 0 1];
        elseif PRN_ID == 6
            loadReg2 = [0 1 0 0 1 0 1 1 0 0];
        elseif PRN_ID == 7
            loadReg2 = [0 0 1 0 0 0 1 1 1 0];
        else
            disp('PRN_ID is incorrect. Must be from 1 to 7');
        end;
    else 
        disp('signal_type is incorrect. Must be L5 or S');
    end;
loadReg2 = loadReg2(end:-1:1);
    
goldseq = comm.GoldSequence('FirstPolynomial',[10 7 0],...
    'SecondPolynomial',[10, 8, 7, 4, 2, 1, 0],...
    'FirstInitialConditions', [1 1 1 1 1 1 1 1 1 1],...
    'SecondInitialConditions',loadReg2,...
    'Index',0,'SamplesPerFrame',L);
genIRNSS = (step(goldseq))*2-1;
