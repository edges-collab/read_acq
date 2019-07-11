function acq2mat(filename, outpath)        
  % Input file name
  % Output directory
    
%matlabpool open 8;

  [inpath, inname, inext] = fileparts(filename);
  files = dir(filename);
    
  if (length(files) == 0)
    % No input file found
    disp('Error: No input file found');
    return;
  end
        
  for n=1:length(files)

    thisfile = [inpath '/' files(n).name];
    [dummy, inname] = fileparts(thisfile);
    outfile_anc = [params.outpath params.inname '_anc.mat'];
    outfile_data = [params.outpath params.inname '_data.mat'];

    disp(sprintf('Importing: %s', thisfile));

    params.tcal = 400;
    params.tload = 300;
    params.nchannels = 16384*4;
    params.max_freq = 200;
    fstep = max_freq/params.nchannels;
    params.freqs = 0:fstep:(params.max_freq .* (1 - params.nchannels.^(-1)));
    nta = 0; 
    params.nanc = 13;
    line_ancillary = zeros(1,nanc);
    nchunk = 2000;
    anc = zeros(params.nanc, nchunk);
    p0 = zeros(params.nchannels, nchunk);
    p1 = zeros(params.nchannels, nchunk);
    p2 = zeros(params.nchannels, nchunk);
    
    bSpec = 0;
    nLine = 0;
    fid = fopen(thisfile, 'r');
    while 1
        % Read lines
        if (bSpec)
          if (nLine > 0)
            % fgets is too slow on big lines, so we have a hack using fread
            line = fread(fid, nLine, '*char');
          else
            % This is the first time we are reading a line containing a spectrum
            % and need to get its size
            line = fgets(fid);
            nLine = length(line);            
          end
          bSpec = 0;
        else
          line = fgets(fid);
          bSpec = 1;
        end
        
        if (~ischar(line) || length(line)==0)
            break;
        end
        
        if ((line(1) ~= '*') & (line(1) ~= '\n') & (line(1) ~= '#'))

            % Read the ancillary info
            try 
              [line_ancillary(1:10), count, errmsg, index] = sscanf(line, '%g:%g:%g:%g:%g %g %g %g %g %g spectrum ', 10);  
            catch
              
              % Uh oh, we probably have an incomplete file.
              disp(sprintf('Error: Failed to parse ancillary data at spectrum# %d', nta));              
              break;
             
            end             
              
            % Read the spectrum
            spec = decode(line((index+params.nskip*4):end));
            total_power = sum(spec);
            
            % Check the switch state and act accordingly
            swpos = line_ancillary(6);
            switch(swpos)
                case 0
                    p0 = spec;
                    a0 = line_ancillary;
                case 1
                    p1 = spec;
                    ancillary(12) = total_power;
                case 2
                    p2 = spec;
                    
                    if ((length(p0) == length(p1)) && (length(p0) == length(p2)) && (length(p1) == length(p2)))         
                      
                      ancillary(13) = total_power;
                      ta = (p0-p1) ./ (p2-p1) .* params.tcal + params.tload;
                      
                      nta = nta + 1;
                      if (mod(nta, nchunk)==0)
                                                
                        newspec = zeros(size(p0) + [0 nchunk]);
                        newspec(:,1:nta) = p0;
                        p0 = newspec;
                        
                        newspec = zeros(size(p1) + [0 nchunk]);
                        newspec(:,1:nta) = p1;
                        p1 = newspec;
                        
                        newspec = zeros(size(p2) + [0 nchunk]);
                        newspec(:,1:nta) = p2;
                        p2 = newspec;     
                                             
                      end
                      
                      waterfall(:,nta) = ta;
                      waterfall_0(:,nta) = p0;
                      waterfall_1(:,nta) = p1;
                      waterfall_2(:,nta) = p2;
                      
                      if (mod(nta,100)==0)
                        disp(sprintf('%04d:%03d:%02d:%02d:%02d - line=%d', ancillary(1:5), nta));
                      end
                    else
                      % Uh oh, we probably have an incomplete file.
                      disp(sprintf('Error: Failed to calibrate at spectrum# %d', nta));
                      break;
                    end
            end
        end
    end
    
    fclose(fid);
    
    % Do some additional metrics and write to the disk the current waterfall (if it exists)
    if (bWriteEach == 1 && nta > 0)
      do_plots(params, waterfall(:,1:nta), waterfall_0(:,1:nta), waterfall_1(:,1:nta), waterfall_2(:,1:nta));
    end
  
  end % loop through all files  
     
  % Write to the disk the current waterfall (if it exists)
  if (bWriteEach == 0 && nta > 0)
    do_plots(params, waterfall(:,1:nta), waterfall_0(:,1:nta), waterfall_1(:,1:nta), waterfall_2(:,1:nta));
  end

  disp(sprintf('Done. %d lines processed', nta));
    
  matlabpool close;

end


        








