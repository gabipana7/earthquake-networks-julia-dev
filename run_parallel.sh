echo `date`
julia --project motifs_analysis_tsallis.jl Italy totalenergy &
julia --project motifs_analysis_tsallis.jl Italy meanenergy &
julia --project motifs_analysis_tsallis.jl California totalenergy &
julia --project motifs_analysis_tsallis.jl California meanenergy &
julia --project motifs_analysis_tsallis.jl Japan totalenergy &
julia --project motifs_analysis_tsallis.jl Japan meanenergy &
wait 
echo `date`