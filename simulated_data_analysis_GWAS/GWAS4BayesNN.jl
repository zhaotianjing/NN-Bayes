
function GWAS4BayesNN(mme,map_file,marker_effects_file1,marker_effects_file2,weights_file,ebv_file;
              #mapfile
              #window
              window_size = "1 Mb",sliding_window = false,
              #GWAS
              GWAS = true, threshold = 0.001,
              #genetic correlation
              genetic_correlation = false,
              #misc
              header = true, output_winVarProps = false)

    if typeof(window_size) == String
        if split(window_size)[2] != "Mb"
            error("The format for window_size is \"1 Mb\".")
        end
    end

    if map_file == false && typeof(window_size) <: Integer
        println("The map file is not provided. A fake map file is generated with $window_size markers in each 1 Mb window.")
        nmarkers=length(readdlm(marker_effects_file1,',',header=true)[2])
        nWindow=Int(floor(nmarkers/window_size))
        println("In total $nWindow windows, $nmarkers markers.")
        #Part1: window1 to windiw[end-1]
        nmarkers_part1 = window_size*(nWindow-1)
        start_part1    = 1
        end_part1      = (nWindow-1)*10^6
        skip_part1     = end_part1/nmarkers_part1
        position_part1=Int.(floor.(collect(start_part1:skip_part1:end_part1)))
        #Part2: window[end]
        nmarkers_part2 = nmarkers - nmarkers_part1
        start_part2    = end_part1 + 1
        end_part2      = nWindow*10^6
        skip_part2     = (end_part2-start_part2+1)/nmarkers_part2
        position_part2=Int.(floor.(collect(start_part2:skip_part2:end_part2)))
        #Combine part1&part2
        mapfile = DataFrame(markerID=1:nmarkers,
                            chromosome=fill(1,nmarkers),
                            position=[position_part1;position_part2])
        CSV.write("mapfile.temp",mapfile)
        window_size_bp=1_000_000
    else
        mapfile = (header == true ? readdlm(map_file,',',header=true)[1] : readdlm(map_file,','))
        window_size_bp = map(Int64,parse(Float64,split(window_size)[1])*1_000_000)
    end


    chr     = map(string,mapfile[:,2])
    pos     = map(Int64,mapfile[:,3])

    window_size_nSNPs   = Array{Int64,1}()  #save number of markers in ith window for all windows
    window_chr          = Array{String,1}() #1
    window_pos_start    = Array{Int64,1}()  #1_000_000
    window_pos_end      = Array{Int64,1}()  #2_000_000
    window_snp_start    = Array{Int64,1}()  #1_314_314
    window_snp_end      = Array{Int64,1}()  #1_999_003
    window_column_start = Array{Int64,1}()  #101
    window_column_end   = Array{Int64,1}()  #200

    index_start = 1
    for i in unique(chr)
      pos_on_chri     = pos[chr.== i] #assume chr and pos are sorted
      if sliding_window == false
          nwindow_on_chri = ceil(Int64,pos_on_chri[end]/window_size_bp)
      else
          nwindow_on_chri = findfirst(x -> x >= pos_on_chri[end] - window_size_bp, pos_on_chri)
      end

      for j in 1: nwindow_on_chri
        if sliding_window == false
            thisstart = window_size_bp*(j-1)
        else
            thisstart = pos_on_chri[j]
        end
        thisend  = thisstart + window_size_bp
        snps_window = thisstart .<= pos_on_chri .< thisend
        snps_window_sizej = sum(snps_window)
        #empty windows exist in non-sliding window; no empty window in sliding windows
        #empty windows were deleted
        if snps_window_sizej!=0
            push!(window_snp_start,pos_on_chri[findfirst(snps_window)])
            push!(window_snp_end,pos_on_chri[findlast(snps_window)])
            push!(window_column_start,index_start)
            push!(window_column_end,index_start+snps_window_sizej-1)
            push!(window_chr,i)
            push!(window_pos_start,thisstart)
            push!(window_pos_end,thisend)
            push!(window_size_nSNPs,snps_window_sizej)
        end
        if sliding_window == false
            index_start += snps_window_sizej
        else
            index_start += 1
        end
      end
    end

    out=[]
    if GWAS == true
        println("Compute the posterior probability of association of the genomic window that explains more than ",threshold," of the total genetic variance.")

            #using marker effect files
            output_hn1        = readdlm(marker_effects_file1,',',header=true)[1]
            output_hn2        = readdlm(marker_effects_file2,',',header=true)[1]
            weights           = readdlm(weights_file,',')
            ebv_sample        = readdlm(ebv_file,header=true,',')[1]

            nsamples,nMarkers = size(output_hn1)
            nWindows          = length(window_size_nSNPs)
            winVarProps       = zeros(nsamples,nWindows)
            winVar            = zeros(nsamples,nWindows)
            #window_mrk_start ID and window_mrk_end ID are not provided now
            X = (typeof(mme) <: Array ? mme : mme.M[1].output_genotypes)
            nInd=size(X,1)
            @showprogress "running GWAS..." for i=1:nsamples
                α_hn1  = output_hn1[i,:]                           #i-th MCMC iteration
                α_hn2  = output_hn2[i,:]
                w      = weights[i,:]

                g_all = ebv_sample[i,:]
                genVar = var(g_all)
                for winj = 1:length(window_column_start)
                  wStart = window_column_start[winj]
                  wEnd   = window_column_end[winj]
                  BV_winj= tanh.(X[:,wStart:wEnd]*[α_hn1[wStart:wEnd] α_hn2[wStart:wEnd]])*w[2:end]
                  var_winj = var(BV_winj)
                  winVar[i,winj]      = var_winj
                  winVarProps[i,winj] = var_winj/genVar
                end
            end
            winVarProps[isnan.(winVarProps)] .= 0.0 #replace NaN caused by situations no markers are included in the model
            WPPA, prop_genvar = vec(mean(winVarProps .> threshold,dims=1)), vec(mean(winVarProps,dims=1))
            prop_genvar = round.(prop_genvar*100,digits=2)
            winVarmean = vec(mean(winVar,dims=1))
            winVarstd  = vec(std(winVar,dims=1))

            srtIndx = sortperm(WPPA,rev=true)

            outi = DataFrame(trait  = fill(1,length(WPPA))[srtIndx],
                            window = (1:length(WPPA))[srtIndx],
                            chr    = window_chr[srtIndx],
                            wStart = window_pos_start[srtIndx],
                            wEnd   = window_pos_end[srtIndx],
                            start_SNP = window_snp_start[srtIndx],
                            end_SNP   = window_snp_end[srtIndx],
                            numSNP  = window_size_nSNPs[srtIndx],
                            estimateGenVar  = winVarmean[srtIndx],
                            stdGenVar     = winVarstd[srtIndx],
                            prGenVar = prop_genvar[srtIndx],
                            WPPA     = WPPA[srtIndx],
                            PPA_t  = cumsum(WPPA[srtIndx]) ./ (1:length(WPPA)))
             push!(out,outi)
             outfile="GWAS_WPPA4Bayesnn.txt"
             CSV.write(outfile, outi)
    end

    return output_winVarProps ? (Tuple(out),winVarProps) : Tuple(out)
end
