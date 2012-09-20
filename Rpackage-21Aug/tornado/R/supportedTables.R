supportedTables <-
function(genome){
	message("getting supported tables. may take several minutes.")
	require(GenomicFeatures)
	genome.tracks = supportedUCSCFeatureDbTracks(genome)
	intersect(genome.tracks,rownames(supportedUCSCtables()))
	}
