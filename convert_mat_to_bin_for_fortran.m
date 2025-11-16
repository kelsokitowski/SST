function convert_mat_to_bin_for_fortran(bfValues)
%---------------------------------------------------------------------
%  Converts existing .mat files in results/param_<bf>/ into .bin
%  files readable by data_loader_all.f90.
%
%  Required files in each folder:
%     checkpoint.mat   (contains v_1..v_6,t,tStar,counter)
%     weightStuff.mat  (contains weight,triadFlag,outsideCutCell,
%                       insideCutCell,CxVals,CyVals,Q11)
%     ETDRKcoeffs.mat  (contains EX1_*,EX2_*,Q_*,f1_*,f2_*,f3_*)
%
%  Output:
%     checkpoint.bin, weightStuff.bin, ETDRKcoeffs.bin
%
%  Each .bin begins with int32(kLength) header, then arrays in
%  column-major (Fortran) order, as expected by data_loader_all.
%---------------------------------------------------------------------

for ib = 1:numel(bfValues)
    bf = bfValues(ib);
    folder = sprintf('results/param_%g', bf);
    fprintf('\n=== Processing bf = %.3f ===\n', bf);
    if ~isfolder(folder)
        warning('Folder %s not found, skipping.', folder);
        continue;
    end

    % ---------------- checkpoint ----------------
    chkfile = fullfile(folder,'checkpoint.mat');
    if isfile(chkfile)
        S = load(chkfile);
        kLength = numel(S.v_1);
        fid = fopen(fullfile(folder,'checkpoint.bin'),'w');
        fwrite(fid,int32(kLength),'int32');
        fwrite(fid,[S.v_1 S.v_2 S.v_3 S.v_4 S.v_5 S.v_6 ...
                    S.t S.tStar double(S.counter)],'double');
        fclose(fid);
        fprintf('  ✓ checkpoint.bin written (kLength=%d)\n',kLength);
    else
        warning('  checkpoint.mat missing.');
    end

    % ---------------- weightStuff ----------------
    wfile = fullfile(folder,'weightStuff.mat');
    if isfile(wfile)
        W = load(wfile);
        kLength = size(W.weight,1);
        fid = fopen(fullfile(folder,'weightStuff.bin'),'w');
        fwrite(fid,int32(kLength),'int32');
        fwrite(fid,W.weight,'double');
        fwrite(fid,int32(W.triadFlag),'int32');
        fwrite(fid,int32(W.outsideCutCell),'int32');
        fwrite(fid,int32(W.insideCutCell),'int32');
        fwrite(fid,W.CxVals,'double');
        fwrite(fid,W.CyVals,'double');
        fwrite(fid,int32(W.Q11),'int32');
        fclose(fid);
        fprintf('  ✓ weightStuff.bin written\n');
    else
        warning('  weightStuff.mat missing.');
    end

    % ---------------- ETDRKcoeffs ----------------
    efile = fullfile(folder,'ETDRKcoeffs.mat');
    if isfile(efile)
        E = load(efile);
        names = {'EX1_1','EX2_1','Q_1','f1_1','f2_1','f3_1', ...
                 'EX1_2','EX2_2','Q_2','f1_2','f2_2','f3_2', ...
                 'EX1_3','EX2_3','Q_3','f1_3','f2_3','f3_3', ...
                 'EX1_4','EX2_4','Q_4','f1_4','f2_4','f3_4', ...
                 'EX1_5','EX2_5','Q_5','f1_5','f2_5','f3_5', ...
                 'EX1_6','EX2_6','Q_6','f1_6','f2_6','f3_6'};
        kLength = numel(E.EX1_1);
        fid = fopen(fullfile(folder,'ETDRKcoeffs.bin'),'w');
        fwrite(fid,int32(kLength),'int32');
        for i = 1:numel(names)
            fwrite(fid,E.(names{i}),'double');
        end
        fclose(fid);
        fprintf('  ✓ ETDRKcoeffs.bin written\n');
    else
        warning('  ETDRKcoeffs.mat missing.');
    end
end
fprintf('\nAll conversions complete.\n');
end
