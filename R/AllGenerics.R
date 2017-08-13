

setGeneric("unpack_transcripts", 
           function(object, ...) standardGeneric("unpack_transcripts"))

setGeneric("subset_transcripts",
           function(window, object, ...) standardGeneric("subset_transcripts"))

setGeneric("add_stepping",
           function(object,  ...) standardGeneric("add_stepping"))

setGeneric("make_annotation_track",
           function(window, object, ...) standardGeneric("make_annotation_track"))

setGeneric("make_signal_track",
           function(window, object, ...) standardGeneric("make_signal_track"))

setGeneric("make_locus_summary",
           function(object, ...) standardGeneric("make_locus_summary"))

setGeneric("make_locus_summaries",
           function(object, ...) standardGeneric("make_locus_summaries"))

setGeneric("make_summary_plotter",
           function(object, ...) standardGeneric("make_summary_plotter"))

setGeneric("make_track_plotter",
           function(object, ...) standardGeneric("make_track_plotter"))

setGeneric("single_locus_view",
           function(window, object, ...) standardGeneric("single_locus_view"))


setGeneric("multi_locus_view",
           function(windows, object, ...) standardGeneric("multi_locus_view"))


setGeneric("relative_position",
           function(view, ...) standardGeneric("relative_position"))

setGeneric("get_layout",
           function(object, ...) standardGeneric("get_layout"))



setGeneric("normalize_coverage_matrix",
           function(mats, ...) standardGeneric("normalize_coverage_matrix"))

