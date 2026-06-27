if (getRversion() >= "2.15.1") {
    utils::globalVariables(c(
        "Cell",
        "Coefficient",
        "Color_Val",
        "Detector",
        "Expression",
        "File",
        "Fluorophore",
        "FSC-A",
        "Intensity",
        "Marker",
        "Marker1",
        "Marker2",
        "Max_NPS",
        "MedianResidual",
        "metric_value",
        "NPS",
        "Receiving_Marker",
        "Residual",
        "RMS",
        "Sample",
        "Similarity",
        "SSC-A",
        "Signature",
        "Spilling_Marker",
        "Spread",
        "scatter_x",
        "scatter_y",
        "bin_idx",
        "ch_idx",
        "color",
        "count",
        "fill",
        "x",
        "x_bin",
        "x_high",
        "x_low",
        "y",
        "y_bin",
        "y_high",
        "y_low",
        "y_orig"
    ))
}

# Package-level cache for static data files (spectral libraries, dictionaries)
.spectreasy_cache <- new.env(parent = emptyenv())

