version 1.0

struct IndexData {
	File data
	File data_index
}

struct Sample {
	String sample_id
	Array[IndexData] movie_bams
}
